#!/bin/bash

# Activate Alpine HMMer environment
module load miniforge
mamba activate /pl/active/jravilab/HMMer

# Set input/output paths
INPUT_DIR="hmmer_inputs"
OUTPUT_DIR="hmmer_results"
mkdir -p $OUTPUT_DIR

# Placeholder HMM databases
COG_DB="/pl/active/jravilab/HMMer/databases/COG_database2024.hmm"
ARG_DB="/pl/active/jravilab/HMMer/databases/NCBI_ResFinder.hmm"

# Run HMMer on gene sequences
hmmscan --cpu 16 --tblout $OUTPUT_DIR/genes_ARG.tbl $ARG_DB $INPUT_DIR/Aba_genes_for_hmmer.fasta
hmmscan --cpu 16 --tblout $OUTPUT_DIR/genes_COG.tbl $COG_DB $INPUT_DIR/Aba_genes_for_hmmer.fasta

# Run HMMer on protein sequences
hmmscan --cpu 16 --tblout $OUTPUT_DIR/proteins_ARG.tbl $ARG_DB $INPUT_DIR/Aba_proteins_for_hmmer.fasta
hmmscan --cpu 16 --tblout $OUTPUT_DIR/proteins_COG.tbl $COG_DB $INPUT_DIR/Aba_proteins_for_hmmer.fasta

echo "HMMer time complete. Results saved in $OUTPUT_DIR"
