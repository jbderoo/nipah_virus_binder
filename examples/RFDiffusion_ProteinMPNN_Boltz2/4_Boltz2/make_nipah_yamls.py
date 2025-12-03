#!/usr/bin/env python3
import yaml
import glob
import sys
import os
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="Match binder sequences to MSAs and write paired YAML input files."
    )
    parser.add_argument(
        "input_fasta",
        help="Input FASTA file containing designed binder sequences."
    )
    parser.add_argument(
        "input_msa_dir",
        help="Directory containing MSA (.a3m) files to match sequences against."
    )
    parser.add_argument(
        "yaml_out_dir",
        help="Output directory for generated YAML files."
    )
    return parser.parse_args()


nipah_seq = (
    "MPAENKKVRFENTTSDKGKIPSKVIKSYYGTMDIKKINEGLLDSKILSAFNTVIALLGSIVIIVMNIMIIQNYTRSTDNQAVIKDALQGIQQQIKGLADKIGTEIGPKVSLIDTSSTITIPANIGLLGSKISQSTASINENVNEKCKFTLPPLKIHECNISCPNPLPFREYRPQTEGVSNLVGLPNNICLQKTSNQILKPKLISYTLPVVGQSGTCITDPLLAMDEGYFAYSHLERIGSCSRGVSKQRIIGVGEVLDRGDEVPSLFMTNVWTPPNPNTVYHCSAVYNNEFYYVLCAVSTVGDPILNSTYWSGSLMMTRLAVKPKSNGGGYNQHQLALRSIEKGRYDKVMPYGPSGIKQGDTLYFPAVGFLVRTEFKYNDSNCPITKCQYSKPENCRLSMGIRPNSHYILRSGLLKYNLSDGENPKVVFIEISDQRLSIGSPSKIYDSLGQPVFYQASFSWDTMIKFGDVLTVNPLVVNWRNNTVISRPGQSQCPRFNTCPEICWEGVYNDAFLIDRINWISAGVFLDSNQTAENPVFTVFKDNEILYRAQLASEDTNAQKTITNCFLLKNKIWCISLVEIYDTGDNVIRPKLFAVKIPEQCT"
)


def get_sequences(fasta):
    my_dict = {}
    with open(fasta) as f:
        for line in f:
            if line[0] == '>':
                name = line.strip().replace('>', '')
                my_dict[name] = None
            elif line == '\n':
                continue
            else:
                my_dict[name] = line.strip()
    return my_dict


def seq_matches_msa(seq, msa):
    with open(msa) as f:
        f.readline()
        seq_de_msa = f.readline().rstrip("\n")
    return seq_de_msa == seq.strip()


def main():
    args = parse_args()
    input_fasta = args.input_fasta
    input_msa_dir = args.input_msa_dir
    yaml_out_dir = args.yaml_out_dir

    seq_dict = get_sequences(input_fasta)
    msas = glob.glob(f"{input_msa_dir}/*a3m")

    if not msas:
        print(f"No MSA files found in {input_msa_dir}")
        sys.exit(1)

    match_msa = {}
    for name, seq in seq_dict.items():
        for msa in msas:
            if seq_matches_msa(seq, msa):
                match_msa[name] = msa
                break

    for name, seq in seq_dict.items():
        yaml_dict = {
            "version": 1,
            "sequences": [
                {
                    "protein": {
                        "id": "A",
                        "sequence": nipah_seq,
                        "msa": "full_nipah.a3m"
                    }
                },
                {
                    "protein": {
                        "id": "B",
                        "sequence": seq,
                        "msa": match_msa.get(name, "NO_MATCH_FOUND")
                    }
                }
            ]
        }

        if "NO_MATCH_FOUND" in yaml_dict["sequences"][1]["protein"]["msa"]:
            print(f"error! Couldn't find MSA match for {name}")
            sys.exit(1)

        os.makedirs(yaml_out_dir, exist_ok=True)
        out_name = f"{yaml_out_dir}/{name}.yaml"
        with open(out_name, "w") as f:
            yaml.dump(yaml_dict, f, sort_keys=False, width=120, indent=4)

        # Uncomment for verbose output:
        # print(f"Wrote {out_name}")


if __name__ == "__main__":
    main()

