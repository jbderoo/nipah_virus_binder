python make_nipah_yamls.py \
    ../3_MMseqs2/all_seqs.fasta \
    ../3_MMseqs2/all_msas \
    nipah_yamls 

# this script is "hard coded" to take the nipah virus sequence and a3m, then pairs it with every designed binder sequence and it's corresponding MSA. These are both inputs; the giant fasta file we produced previously that was an input into MMseqs2 to make the .a3m's, and then that produced directory of .a3m's. Thirdly we provide an output directory for the input yamls to go to for folding with boltz2.
