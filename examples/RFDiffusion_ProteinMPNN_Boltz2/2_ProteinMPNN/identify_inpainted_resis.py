#!/usr/bin/env python3

import sys
import os
import pickle

def extract_pdb_sequences(pdb_path):
    """
    Extracts the amino acid sequences of all chains in a PDB file.
    Returns a dictionary: {chain_id: sequence}
    """

    aa_map = {
        'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q',
        'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
        'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V','UNK':'X'
    }

    sequences = {}
    last_res = {}

    with open(pdb_path, 'r') as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue

            resname = line[17:20].strip()
            chain_id = line[21].strip()
            resnum = line[22:26].strip()

            if resname not in aa_map:
                continue

            if chain_id not in sequences:
                sequences[chain_id] = ""
                last_res[chain_id] = None

            if last_res[chain_id] != resnum:
                sequences[chain_id] += aa_map[resname]
                last_res[chain_id] = resnum

    return sequences


def get_false_positions(pdb_path, binder_chain, design_bool):
    """
    Returns a list of chain+resnum (e.g., 'A12') for residues where design_bool is False
    in the binder chain only. It indexes design_bool only over the binder residues.
    """
    positions = []
    binder_index = 0
    seen = set()

    with open(pdb_path, 'r') as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue

            chain = line[21].strip()
            resnum = line[22:26].strip()

            if chain == binder_chain and resnum not in seen:
                if binder_index < len(design_bool):
                    if not design_bool[binder_index]:
                        positions.append(f"{chain}{resnum}")
                seen.add(resnum)
                binder_index += 1

    return positions


def main():
    if len(sys.argv) != 2:
        print("\nUsage:\n  python find_fixed_resis.py <path/to/structure.pdb>\n")
        sys.exit(1)

    pdb_path = sys.argv[1]

    if not os.path.exists(pdb_path):
        print(f"Error: File not found: {pdb_path}")
        sys.exit(1)

    # Build matching TRB path
    trb_path = os.path.splitext(pdb_path)[0] + ".trb"

    if not os.path.exists(trb_path):
        print(f"Error: Matching TRB not found: {trb_path}")
        sys.exit(1)

    # Extract sequences to identify binder
    seqs = extract_pdb_sequences(pdb_path)

    # Known receptor sequence
    nipah_seq = ('PKLISYTLPVVGQSGTCITDPLLAMDEGYFAYSHLERIGSCSRGVSKQRIIGVGEVLDRGDEVPSLFMTNVWTPPNPNTVYHCSAVYNNEFYYVLCAVSTVGDPILNSTYWSGSLMMTRLAVKPKSNGGGYNQHQLALRSIEKGRYDKVMPYGPSGIKQGDTLYFPAVGFLVRTEFKYNDSNCPITKCQYSKPENCRLSMGIRPNSHYILRSGLLKYNLSDGENPKVVFIEISDQRLSIGSPSKIYDSLGQPVFYQASFSWDTMIKFGDVLTVNPLVVNWRNNTVISRPGQSQCPRFNTCPEICWEGVYNDAFLIDRINWISAGVFLDSNQTAENPVFTVFKDNEILYRAQLASEDTNAQKTITNCFLLKNKIWCISLVEIYDTGDNVIRPKLFAVKIPEQC')

    for chain, seq in seqs.items():
        if seq == nipah_seq:
            receptor = chain
        else:
            binder = chain

    # Load design booleans
    with open(trb_path, "rb") as f:
        data = pickle.load(f)

    design_bool = data['inpaint_seq']

    # Get residues that remain fixed
    false_positions = get_false_positions(pdb_path, binder, design_bool)

    print(" ".join(false_positions))


if __name__ == "__main__":
    main()

