#!/bin/bash

QUERY_FILE= '/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/data/pcs_to_seq.faa'
DB_FILE='/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/VFDB_setB_pro.fas'
DB_NAME="virulence_factors"
OUTPUT_FILE="/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/blast_results_PCs.txt"

# Create the BLAST database if it doesn't already exist
if [ ! -f "${DB_NAME}.pin" ]; then
    echo "Creating BLAST database..."
    makeblastdb -in $DB_FILE -dbtype prot -out $DB_NAME
fi

# Run BLASTP
echo "Running BLASTP..."
blastp -query $QUERY_FILE -db $DB_NAME -out $OUTPUT_FILE -outfmt 6 -evalue 1e-5 -num_threads 8

echo "BLASTP search completed. Results saved to $OUTPUT_FILE."