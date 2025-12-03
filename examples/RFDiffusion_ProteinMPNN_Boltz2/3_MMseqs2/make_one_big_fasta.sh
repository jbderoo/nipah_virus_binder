#!/bin/bash

python merge_and_clean_proteinmpnn_fastas.py \
    ../2_ProteinMPNN/simplified_LigandMPNN_output \
    all_seqs.fasta \
    --remove_seq 2 # custom python script to read all the .fa's produced by ProteinMPNN, simplify them, and condense them into a single .fasta file to pass to MMseqs2 to generate many MSA's very effectively. When ProteinMPNN makes a fasta, it does it in the pattern "$binder_seq:$nipah_virus_seq". This "--remove_seq 2" flag specifies to remove the $nipah_virus_seq because boltz2 wants unpaired MSAs (i.e. separate; one MSA for the nipah virus and one MSA for the binder). 

 
