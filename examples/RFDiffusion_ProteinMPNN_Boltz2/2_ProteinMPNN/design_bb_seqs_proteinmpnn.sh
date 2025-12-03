#!/bin/bash

for pdb in ../1_RFDiffusion/outputs/*.pdb; do
    design_resis=$(python identify_inpainted_resis.py "$pdb")

    echo "Running: $pdb  fixed_residues=$fixed_resis"

    python /path/to/LigandMPNN/run.py \
        --seed 111 \
        --pdb_path "$pdb" \
        --out_folder "./nipah_outputs/fix_residues" \
        --redesigned_residues "$design_resis" \
        --omit_AA "C" \
        --model_type "protein_mpnn" \
    	--verbose 0 \
    	--batch_size 4
done
