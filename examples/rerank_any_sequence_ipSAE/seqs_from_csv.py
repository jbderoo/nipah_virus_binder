#!/usr/bin/env python3
import pandas as pd
import argparse
import os
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert a CSV with 'id' and 'designed_sequence' columns to a FASTA file."
    )
    parser.add_argument(
        "csv",
        help="Input CSV file containing 'id' and 'designed_sequence' columns"
    )
    parser.add_argument(
        "fasta",
        help="Output FASTA filename"
    )
    return parser.parse_args()


def csv_to_fasta(csv_path, fasta_path):
    # Check CSV exists
    if not os.path.exists(csv_path):
        sys.exit(f"Error: Input CSV not found: {csv_path}")

    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        sys.exit(f"Error reading CSV: {e}")

    # Validate required columns
    required_cols = {"id", "designed_sequence"}
    missing = required_cols - set(df.columns)
    if missing:
        sys.exit(f"Error: Missing required columns: {', '.join(missing)}")

    seqs = df["designed_sequence"]
    ids = df["id"]

    if len(seqs) == 0:
        sys.exit("⚠️No sequences found in the input CSV.")

    with open(fasta_path, "w") as f:
        for seq, identifier in zip(seqs, ids):
            if pd.isna(seq) or pd.isna(identifier):
                continue
            f.write(f">{identifier}\n{seq}\n\n")

    print(f"Wrote {len(seqs)} sequences to {fasta_path}")


def main():
    args = parse_args()
    csv_to_fasta(args.csv, args.fasta)


if __name__ == "__main__":
    main()

