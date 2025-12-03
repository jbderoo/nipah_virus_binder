#!/usr/bin/env python3

import argparse
import os
import glob


def parse_header(header):
    """
    Convert header such as:
       >nipah_motif_bb2_9, id=4, T=0.1, ...
    to:
       >nipah_motif_bb2_9_4
    """
    header = header.lstrip(">")
    parts = header.split(",")

    name = parts[0].strip()  # nipah_motif_bb2_9

    # find id=#
    id_val = None
    for p in parts:
        p = p.strip()
        if p.startswith("id="):
            id_val = p.split("=")[1]
            break

    if id_val is None:
        raise ValueError(f"Could not find id= in header: {header}")

    return f">{name}_{id_val}"


def strip_sequence(seq, remove_seq):
    """
    Strip sequence N (1-based) from colon-separated sequences.
    Example:
        seq = "AAA:BBB"
        remove_seq = 2  →  returns "AAA"
    """
    if remove_seq is None:
        return seq

    parts = seq.split(":")

    if len(parts) == 1:
        return seq  # nothing to remove

    idx = remove_seq - 1  # convert to 0-based

    if idx < 0 or idx >= len(parts):
        return seq  # invalid index → leave unchanged

    # keep all except the removed part
    remaining = parts[:idx] + parts[idx+1:]
    return ":".join(remaining)


def process_fasta(filepath, seen_sequences, out_handle, remove_seq):
    """
    Skip first entry; rewrite others; deduplicate by processed sequence.
    """
    with open(filepath, "r") as f:
        lines = f.read().splitlines()

    entries = []
    current_header = None
    current_seq = []

    # Group into (header, seq)
    for line in lines:
        if line.startswith(">"):
            if current_header:
                entries.append((current_header, "".join(current_seq)))
            current_header = line
            current_seq = []
        else:
            current_seq.append(line.strip())

    if current_header:
        entries.append((current_header, "".join(current_seq)))

    # skip the first entry
    for header, seq in entries[1:]:
        seq = seq.strip()

        # strip out one of the paired sequences if requested
        seq_processed = strip_sequence(seq, remove_seq)

        # dedupe post-processing
        if seq_processed in seen_sequences:
            continue

        new_header = parse_header(header)

        out_handle.write(new_header + "\n")
        out_handle.write(seq_processed + "\n")

        seen_sequences.add(seq_processed)


def main():
    parser = argparse.ArgumentParser(description="Condense FASTA files.")
    parser.add_argument("input_dir", help="Directory with .fa files")
    parser.add_argument("output_fasta", help="Output merged FASTA")
    parser.add_argument(
        "--remove_seq",
        type=int,
        default=None,
        help="1-based index of colon-separated sequence to remove (optional)"
    )

    args = parser.parse_args()

    fa_files = glob.glob(os.path.join(args.input_dir, "*.fa"))
    if not fa_files:
        raise RuntimeError("No .fa files found.")

    seen_sequences = set()

    with open(args.output_fasta, "w") as out:
        for f in sorted(fa_files):
            process_fasta(f, seen_sequences, out, args.remove_seq)

    print(f"Done. Wrote {len(seen_sequences)} unique sequences to {args.output_fasta}")


if __name__ == "__main__":
    main()

