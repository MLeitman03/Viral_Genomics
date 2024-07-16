### run_prokka_vf.sh ###

#!/bin/bash

#$ -cwd
# error = Merged with joblog
#$ -o /u/scratch/m/madelain/prokka_log.$JOB_ID

#$ -j y
#$ -l h_rt=24:00:00,h_data=8G
#$ -pe shared 8
#$ -M $USER
#$ -m bea

# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

# Define variables
PROTEINS_FILE="u/scratch/m/madelain/pcs_to_seq.faa"
VIRULENCE_DB="u/scratch/m/madelain/VFDB_setB_pro.fas"
OUTPUT_DIR="/u/home/m/madelain/prokka_output/"
PREFIX="annotated"

# Load the job environment
. /u/local/Modules/default/init/modules.sh
module load anaconda3

# Activate the Prokka environment
conda activate my_env_name

# Run Prokka with the custom database
prokka --proteins $VIRULENCE_DB --outdir $OUTPUT_DIR --prefix $PREFIX --cpus 8 $PROTEINS_FILE

# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "

### run_prokka_vf.sh END ###
