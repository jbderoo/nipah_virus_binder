import glob
import sys
import re
import argparse
import gemmi
import pandas             as pd
import numpy              as np
from   typing             import List, Tuple, Optional
from   collections        import defaultdict
from   Bio.PDB.MMCIF2Dict import MMCIF2Dict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Parse CIF predictions, extract sequences, and rank by ipSAE."
    )
    parser.add_argument(
        "--structure_preds", required=True,
        help="Path to the structure predictions directory (contains /predictions/*)"
    )
    parser.add_argument(
        "--out_prefix", required=True,
        help="Prefix for output CSV files"
    )
    return parser.parse_args()


def cif_chain_sequences(cif_path, prefer_auth_ids=False):
    """Return per-chain protein sequences from an mmCIF/ModelCIF file."""
    d = MMCIF2Dict(cif_path)

    eids = d.get("_entity_poly_seq.entity_id", [])
    nums = d.get("_entity_poly_seq.num", [])
    mons = d.get("_entity_poly_seq.mon_id", [])
    if isinstance(eids, str): eids = [eids]
    if isinstance(nums, str): nums = [nums]
    if isinstance(mons, str): mons = [mons]

    aa3to1 = {
        'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
        'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
        'THR':'T','TRP':'W','TYR':'Y','VAL':'V','SEC':'U','PYL':'O'
    }

    entity_res = defaultdict(list)
    for eid, num, mon in zip(eids, nums, mons):
        try:
            n = int(str(num))
        except Exception:
            n = float("inf")
        entity_res[str(eid)].append((n, aa3to1.get(mon.upper(), 'X')))

    entity_to_seq = {}
    for eid, items in entity_res.items():
        items.sort(key=lambda t: t[0])
        entity_to_seq[eid] = "".join(one for _, one in items)

    asym_ids = d.get("_struct_asym.id", [])
    asym_eids = d.get("_struct_asym.entity_id", [])
    auth_ids = d.get("_struct_asym.pdbx_auth_asym_id", [])
    if isinstance(asym_ids, str): asym_ids = [asym_ids]
    if isinstance(asym_eids, str): asym_eids = [asym_eids]
    if isinstance(auth_ids, str) or auth_ids is None:
        auth_ids = [] if auth_ids in (None, "") else [auth_ids]

    chain_to_seq = {}
    for i, (asym, eid) in enumerate(zip(asym_ids, asym_eids)):
        chain_id = (
            auth_ids[i] if (prefer_auth_ids and auth_ids and i < len(auth_ids) and auth_ids[i])
            else asym
        )
        chain_to_seq[str(chain_id).strip()] = entity_to_seq.get(str(eid), "")
    return chain_to_seq

### begin huge occlusion check ###

def get_chain(structure: gemmi.Structure, chain_id: str) -> gemmi.Chain:
    """
    Extract a specific chain from a structure.

    Args:
        structure: gemmi Structure object
        chain_id: Single-letter chain identifier

    Returns:
        The requested chain

    Raises:
        ValueError: If chain not found
    """
    for model in structure:
        for chain in model:
            if chain.name == chain_id:
                return chain
    raise ValueError(f"Chain '{chain_id}' not found in structure")


def get_ca_atoms(chain: gemmi.Chain, skip_first_n: int = 0) -> List[gemmi.Atom]:
    """
    Extract all CA (alpha carbon) atoms from a chain.

    Args:
        chain: gemmi Chain object
        skip_first_n: Number of N-terminal residues to skip (default: 0)

    Returns:
        List of CA atoms
    """
    ca_atoms = []
    residue_count = 0

    for residue in chain:
        residue_count += 1
        if residue_count <= skip_first_n:
            continue

        for atom in residue:
            if atom.name == "CA":
                ca_atoms.append(atom)
    return ca_atoms


def kabsch_algorithm(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute optimal rotation matrix and translation vector using Kabsch algorithm.

    This aligns mobile structure P onto fixed structure Q.

    Args:
        P: Nx3 array of mobile structure coordinates
        Q: Nx3 array of fixed (reference) structure coordinates

    Returns:
        (R, t) where R is 3x3 rotation matrix and t is translation vector
        Apply as: aligned_P = (R @ P.T).T + t
    """
    # Center both coordinate sets
    centroid_P = np.mean(P, axis=0)
    centroid_Q = np.mean(Q, axis=0)

    P_centered = P - centroid_P
    Q_centered = Q - centroid_Q

    # Compute covariance matrix
    H = P_centered.T @ Q_centered

    # SVD
    U, S, Vt = np.linalg.svd(H)

    # Compute rotation matrix
    R = Vt.T @ U.T

    # Handle reflection case (ensure proper rotation)
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T

    # Compute translation
    t = centroid_Q - (R @ centroid_P)

    return R, t


def find_dummy_atom(structure: gemmi.Structure) -> Optional[gemmi.Atom]:
    """
    Find the dummy atom marking the pore entrance.
    Looks for atoms with element 'X' or names like 'DUM', 'DU', or 'X'.

    Args:
        structure: gemmi Structure to search

    Returns:
        The dummy atom, or None if not found
    """
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    element = atom.element.name.upper()
                    atom_name = atom.name.strip().upper()

                    if element == "X" or atom_name in ("DUM", "DU", "X"):
                        return atom
    return None


def check_occlusion(
    reference_pdb: str,
    target_cif: str,
    alignment_chain: str = "A",
    contact_cutoff: float = 6.0,
    min_contacts: int = 3,
    skip_n_terminal: int = 175,
    include_alignment_chain: bool = False,
    save_aligned: str = None,
    verbose: bool = False
) -> int:
    """
    Check if a binder protein occludes a pore entrance.

    Workflow:
    1. Load reference pore structure (with dummy atom at pore entrance)
    2. Load target binder structure
    3. Align structures based on chain A CA atoms using Kabsch algorithm
       (skipping first N residues to avoid flexible regions)
    4. Apply rotation + translation to entire target structure
    5. Count CA atoms in target within cutoff distance of dummy atom
    6. Return PASS (0) or FAIL (1) based on contact threshold

    Args:
        reference_pdb: Path to reference structure with dummy atom
        target_cif: Path to target binder structure
        alignment_chain: Chain ID to use for alignment (default: "A")
        contact_cutoff: Distance cutoff for CA-dummy contacts in Angstroms (default: 6.0)
        min_contacts: Minimum contacts required to pass (default: 3)
        skip_n_terminal: Number of N-terminal residues to skip in alignment (default: 175)
        include_alignment_chain: Include alignment chain in contact counting (default: False)
        save_aligned: Path to save aligned target structure (default: None)
        verbose: Print detailed output (default: False)

    Returns:
        0 = PASS (sufficient occlusion)
        1 = FAIL (insufficient occlusion)

    Raises:
        RuntimeError: If CA atoms or dummy atom not found
        ValueError: If alignment chain not found
    """

    # ========================================================================
    # 1. LOAD STRUCTURES
    # ========================================================================
    if verbose:
        print("=" * 70)
        print("PORE OCCLUSION CHECK")
        print("=" * 70)
        print(f"\nLoading structures:")
        print(f"  Reference (pore):  {reference_pdb}")
        print(f"  Target (binder):   {target_cif}")

    ref_structure = gemmi.read_structure(reference_pdb)
    target_structure = gemmi.read_structure(target_cif)

    # ========================================================================
    # 2. EXTRACT ALIGNMENT CHAINS
    # ========================================================================
    if verbose:
        print(f"\nExtracting alignment chain '{alignment_chain}' from both structures...")

    ref_chain = get_chain(ref_structure, alignment_chain)
    target_chain = get_chain(target_structure, alignment_chain)

    ref_cas = get_ca_atoms(ref_chain, skip_first_n=skip_n_terminal)
    target_cas = get_ca_atoms(target_chain, skip_first_n=skip_n_terminal)

    if not ref_cas or not target_cas:
        raise RuntimeError(
            f"No CA atoms found on chain '{alignment_chain}' in reference or target"
        )

    if verbose:
        print(f"  Reference chain {alignment_chain}: {len(ref_cas)} CA atoms (skipped first {skip_n_terminal} residues)")
        print(f"  Target chain {alignment_chain}:    {len(target_cas)} CA atoms (skipped first {skip_n_terminal} residues)")

    # ========================================================================
    # 3. PREPARE COORDINATES FOR KABSCH ALIGNMENT
    # ========================================================================
    if verbose:
        print("\nPreparing coordinates for Kabsch alignment...")
        print("  (Target is MOBILE, Reference is FIXED)")

    # We need to match CA atoms by residue number for proper alignment
    # Build dictionaries keyed by sequence position (after skipping)
    def build_ca_dict(ca_list, skip_n):
        ca_dict = {}
        for ca in ca_list:
            # Get the residue this CA belongs to
            # The residue number after skipping
            for model in [ref_structure, target_structure]:
                for chain in model:
                    for i, residue in enumerate(chain):
                        if i < skip_n:
                            continue
                        for atom in residue:
                            if atom is ca:
                                ca_dict[i - skip_n] = ca
                                break
        return ca_dict

    # Simpler approach: just use the CA lists in order
    # Assume they correspond by position after skipping
    min_len = min(len(ref_cas), len(target_cas))
    ref_cas = ref_cas[:min_len]
    target_cas = target_cas[:min_len]

    if verbose:
        print(f"  Using {min_len} corresponding CA pairs for alignment")

    # Convert to numpy arrays
    ref_coords = np.array([[ca.pos.x, ca.pos.y, ca.pos.z] for ca in ref_cas])
    target_coords = np.array([[ca.pos.x, ca.pos.y, ca.pos.z] for ca in target_cas])

    if verbose:
        print(f"  Reference CA center: ({ref_coords.mean(axis=0)[0]:.2f}, "
              f"{ref_coords.mean(axis=0)[1]:.2f}, {ref_coords.mean(axis=0)[2]:.2f})")
        print(f"  Target CA center:    ({target_coords.mean(axis=0)[0]:.2f}, "
              f"{target_coords.mean(axis=0)[1]:.2f}, {target_coords.mean(axis=0)[2]:.2f})")

    # ========================================================================
    # 4. COMPUTE KABSCH ALIGNMENT
    # ========================================================================
    if verbose:
        print("\nComputing Kabsch alignment (rotation + translation)...")

    R, t = kabsch_algorithm(target_coords, ref_coords)

    # Compute RMSD before and after
    initial_rmsd = np.sqrt(np.mean(np.sum((target_coords - ref_coords)**2, axis=1)))
    aligned_coords = (R @ target_coords.T).T + t
    final_rmsd = np.sqrt(np.mean(np.sum((aligned_coords - ref_coords)**2, axis=1)))

    if verbose:
        print(f"  Initial RMSD: {initial_rmsd:.3f} Å")
        print(f"  Final RMSD:   {final_rmsd:.3f} Å")
        print(f"  Rotation matrix determinant: {np.linalg.det(R):.6f} (should be ~1.0)")
        print(f"  Translation vector: ({t[0]:.2f}, {t[1]:.2f}, {t[2]:.2f})")

    # ========================================================================
    # 5. APPLY TRANSFORMATION TO ENTIRE TARGET STRUCTURE
    # ========================================================================
    if verbose:
        print("\nApplying rotation and translation to entire target structure...")

    for model in target_structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    # Get current position as numpy array
                    pos = np.array([atom.pos.x, atom.pos.y, atom.pos.z])

                    # Apply rotation and translation
                    new_pos = R @ pos + t

                    # Update atom position
                    atom.pos = gemmi.Position(new_pos[0], new_pos[1], new_pos[2])

    # ========================================================================
    # 6. SAVE ALIGNED STRUCTURE (if requested)
    # ========================================================================
    if save_aligned:
        if verbose:
            print(f"\nSaving aligned target structure to: {save_aligned}")
        target_structure.write_pdb(save_aligned)
        if verbose:
            print("  ✓ Saved! Load in PyMOL with:")
            print(f"    load {reference_pdb}")
            print(f"    load {save_aligned}")
            print("  (Target structure is now aligned to reference frame)")

    # ========================================================================
    # 7. FIND DUMMY ATOM IN REFERENCE
    # ========================================================================
    if verbose:
        print("\nLocating dummy atom in reference structure...")

    dummy_atom = find_dummy_atom(ref_structure)
    if dummy_atom is None:
        raise RuntimeError("Dummy atom not found in reference structure")

    dummy_pos = dummy_atom.pos

    if verbose:
        print(f"  Dummy atom found at: ({dummy_pos.x:.2f}, {dummy_pos.y:.2f}, {dummy_pos.z:.2f})")

    # ========================================================================
    # 8. COLLECT CA ATOMS FROM TARGET (AFTER TRANSFORMATION)
    # ========================================================================
    ca_atoms_for_contacts = []

    for model in target_structure:
        for chain in model:
            # Skip alignment chain if requested
            if chain.name == alignment_chain and not include_alignment_chain:
                continue

            for residue in chain:
                for atom in residue:
                    if atom.name == "CA":
                        ca_atoms_for_contacts.append((atom, residue, chain.name))

    if verbose:
        print(f"\nCA atoms available for contact analysis:")
        print(f"  Total: {len(ca_atoms_for_contacts)}")
        if not include_alignment_chain:
            print(f"  (Excluding alignment chain '{alignment_chain}')")

    # ========================================================================
    # 9. COMPUTE DISTANCES TO DUMMY ATOM
    # ========================================================================
    contacts = []

    for atom, residue, chain_name in ca_atoms_for_contacts:
        distance = (atom.pos - dummy_pos).length()
        if distance <= contact_cutoff:
            contacts.append((atom, residue, chain_name, distance))

    if verbose:
        print(f"\n{'─' * 70}")
        print(f"CONTACT ANALYSIS (cutoff = {contact_cutoff} Å)")
        print(f"{'─' * 70}")
        print(f"\nContacts found: {len(contacts)}")

        if contacts:
            print(f"\nDetailed contact list:")
            print(f"  {'Chain':<6} {'Residue':<8} {'Number':<8} {'Distance (Å)':<12}")
            print(f"  {'-'*6} {'-'*8} {'-'*8} {'-'*12}")
            for atom, residue, chain_name, dist in sorted(contacts, key=lambda x: x[3]):
                print(f"  {chain_name:<6} {residue.name:<8} {residue.seqid.num:<8} {dist:>6.2f}")

    # ========================================================================
    # 10. DETERMINE PASS/FAIL
    # ========================================================================
    if verbose:
        print(f"\n{'=' * 70}")
        print(f"RESULT")
        print(f"{'=' * 70}")
        print(f"Contacts found:    {len(contacts)}")
        print(f"Minimum required:  {min_contacts}")
        print()

    if len(contacts) >= min_contacts:
        if verbose:
            print("✔ PASS - Binder sufficiently occludes the pore entrance")
        return 0
    else:
        if verbose:
            print("✘ FAIL - Insufficient occlusion (off-target binding likely)")
        return 1


def main():
    """Command-line interface."""
    parser = argparse.ArgumentParser(
        description="Check if a binder protein occludes a pore entrance",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s occlude_me.pdb model_rank2.cif
  %(prog)s reference.pdb binder.cif --verbose --cutoff 8.0
  %(prog)s pore.pdb binder.cif --save-aligned aligned_binder.pdb -v
  %(prog)s pore.pdb binder.cif --min-contacts 5 --include-alignment-chain
        """
    )

    parser.add_argument(
        "reference_pdb",
        help="Reference structure (pore) with dummy atom at entrance"
    )
    parser.add_argument(
        "target_cif",
        help="Target structure (binder protein)"
    )
    parser.add_argument(
        "--chain", "-c",
        default="A",
        help="Chain ID for alignment (default: A)"
    )
    parser.add_argument(
        "--skip-n-terminal",
        type=int,
        default=175,
        help="Number of N-terminal residues to skip in alignment (default: 175)"
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=6.0,
        help="Distance cutoff for contacts in Angstroms (default: 6.0)"
    )
    parser.add_argument(
        "--min-contacts",
        type=int,
        default=3,
        help="Minimum contacts required to pass (default: 3)"
    )
    parser.add_argument(
        "--include-alignment-chain",
        action="store_true",
        help="Include alignment chain in contact counting"
    )
    parser.add_argument(
        "--save-aligned",
        type=str,
        default=None,
        help="Save aligned target structure to this file (for PyMOL visualization)"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Print detailed output"
    )

    args = parser.parse_args()

    try:
        result = check_occlusion(
            reference_pdb=args.reference_pdb,
            target_cif=args.target_cif,
            alignment_chain=args.chain,
            contact_cutoff=args.cutoff,
            min_contacts=args.min_contacts,
            skip_n_terminal=args.skip_n_terminal,
            include_alignment_chain=args.include_alignment_chain,
            save_aligned=args.save_aligned,
            verbose=args.verbose
        )

        exit(result)

    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        exit(2)

### end huge occlusion check ###

def main():
    args = parse_args()
    structure_preds = args.structure_preds
    out_prefix = args.out_prefix

    dirs = glob.glob(f'{structure_preds}/predictions/*')
    rows = []

    for d in dirs:
        txts = glob.glob(f'{d}/*txt')
        for txt in txts:
            if "byres" not in txt:
                df = pd.read_csv(txt, sep=r"\s+", engine='python', skip_blank_lines=True)
                vals = df.loc[df["Type"].eq("max"), "ipSAE"].squeeze()

        # cif = glob.glob(f'{d}/*model_0*.cif')[0] # old, hopefully not bad?
        cifs = glob.glob(f'{d}/*model_0.cif')
        if len(cifs) != 1:
            raise RuntimeError(f"Ambiguous or missing model_0.cif in {d}: {cifs}")
        cif = cifs[0]

        seqs = cif_chain_sequences(cif)
        kept_seq = ""
        for chain, seq in seqs.items():
            if len(seq) < 251:
                kept_seq = seq

        in_pore = check_occlusion("full_nipah_occlusion.pdb", cif, verbose=False)
        name    = d.split('/')[-1]
        rows.append({"name": name, "ipSAE": vals, "in_pore": in_pore, "sequence": kept_seq})
        if "0420" in name: # specific candidate that is failing:
            print(f'my input file is: {cif}')
            in_pore = check_occlusion("full_nipah_occlusion.pdb", cif, verbose=True)
            print(f'quick result: {in_pore}')

    df = pd.DataFrame(rows)
    df["ipSAE"] = pd.to_numeric(df["ipSAE"], errors="coerce")
    df_sorted = df.sort_values("ipSAE", ascending=False).reset_index(drop=True)
    df_sorted.to_csv(f"{out_prefix}_full_rank_sort.csv", index=False)

    sub = (
        df_sorted.loc[df_sorted["ipSAE"] > 0, ["name", "sequence"]]
        .head(10)
        .reset_index(drop=True)
    )
    sub.insert(0, "index", range(1, len(sub) + 1))
    sub.to_csv(f"{out_prefix}_submission.csv", index=False)


if __name__ == "__main__":
    main()

