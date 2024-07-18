#!/bin/bash
#$ -cwd
#$ -j y
#$ -o /u/scratch/m/madelain/blast.log
#$ -l h_rt=12:00:00
#$ -l h_vmem=128G  # Request 128B of memory
#$ -m bea

#Set paths
BLAST_DIR="/u/home/m/madelain/ncbi-blast-2.16.0+/bin"
DB_DIR="/u/scratch/m/madelain/VFDB"
QUERY="/u/scratch/m/madelain/pcs_to_seq.faa"
DB="${DB_DIR}/VFDB_setB_pro"
OUT="/u/scratch/m/madelain/blast_results.txt"

# Create BLAST database
$BLAST_DIR/makeblastdb -in /u/scratch/m/madelain/VFDB_setB_pro.fas -dbtype prot -out $DB

# Run BLASTP
$BLAST_DIR/blastp -query $QUERY -db $DB -out $OUT -evalue 1e-5 -outfmt 6 -num_threads 8
