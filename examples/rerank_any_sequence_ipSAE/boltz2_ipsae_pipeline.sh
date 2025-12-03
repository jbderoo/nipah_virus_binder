#!/usr/bin bash
# @jderoo

design_csv=boltzgen_final_designs_4.csv               # input csv; this is the final design csv from boltzgen
fasta_file=boltzgen_example.fasta                     # intermediate file name; a name for the fasta file with all seqs extracted from ^
msa_dir=boltzgen_example_msas                         # intermediate dir name; what rime_mmseqs2 will produce
yaml_dir=boltzgen_example_yamls                       # intermediate dir name; what boltz2 will predict from
out_prefix=boltzgen_example                           # out prefix for csv files; submission (top 10) and full rank 
organizing_dir_name=boltzgen_example                  # dir name where all files end up after run



#### begin pipeline
boltz2_outdir="boltz_results_$yaml_dir" # don't change! boltz2 output dir name

# make sure your conda env is boltz221!
python seqs_from_csv.py $design_csv $fasta_file  # extract sequences from input csv to fasta file for rime interpretability
/path/to/your/colabfold_search $fasta_file $msa_dir # make msa's via MMseqs2 
python make_nipah_yamls.py $fasta_file $msa_dir $yaml_dir # boltz2 needs yamls for prediction. Create the yamls of the binder seq, binder a3m, nipah virus seq, nipah virus a3m
boltz predict $yaml_dir --write_full_pae # do boltz prediction. Make sure you're in correct conda env! All my scripts have super low reqs; json, yaml, numpy, etc. Should be able to do in default boltz env with no changes

# process the ipsae calcs for each structure pred. THIS IS NOT MY SCRIPT, this comes from the dunbrack lab
for dir in "$boltz2_outdir"/predictions/*; do
    base=$(basename "$dir")
    python ipsae.py \
        "$dir/pae_${base}_model_0.npz" \
        "$dir/${base}_model_0.cif" \
        10 10
done

# After ipsae calc, organize it all into a pretty csv
#python parse_ipsae_data.py --structure_preds $boltz2_outdir --out_prefix $out_prefix
#python parse_ipsae_data_simple.py --structure_preds $boltz2_outdir --out_prefix $out_prefix ## depending on which version, this line and the above line are the exact same. _simple.py is more up to date and more generic; the next line includes a pore occlusion calculation to ensure we're binding to the correct region of the nipah G protein and not off to the side on the "beta barrel".
python parse_ipsae_data_and_nipah_occlusion.py --structure_preds $boltz2_outdir --out_prefix $out_prefix

# organize
mkdir -p $organizing_dir_name
mv $design_csv       $organizing_dir_name/
mv $fasta_file       $organizing_dir_name/
mv $yaml_dir         $organizing_dir_name/
mv $msa_dir          $organizing_dir_name/
mv $boltz2_outdir    $organizing_dir_name/
mv ${out_prefix}*csv $organizing_dir_name/

