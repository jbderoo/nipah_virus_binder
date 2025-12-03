import glob
import sys
import re
import argparse
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import pandas as pd
from collections import defaultdict


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

        cif = glob.glob(f'{d}/*model_0*.cif')[0]
        seqs = cif_chain_sequences(cif)
        kept_seq = ""
        for chain, seq in seqs.items():
            if len(seq) < 251:
                kept_seq = seq

        name = d.split('/')[-1]
        rows.append({"name": name, "ipSAE": vals, "sequence": kept_seq})

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

