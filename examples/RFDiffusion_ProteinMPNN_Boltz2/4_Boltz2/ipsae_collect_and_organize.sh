boltz2_outdir=boltz_results_nipah_yamls
out_prefix=nipah_binders

for dir in "$boltz2_outdir"/predictions/*; do
    base=$(basename "$dir")
    python ipsae.py \
        "$dir/pae_${base}_model_0.npz" \
        "$dir/${base}_model_0.cif" \
        10 10
done

# After ipsae calc, organize it all into a pretty csv
python parse_ipsae_data_and_nipah_occlusion.py --structure_preds $boltz2_outdir --out_prefix $out_prefix
