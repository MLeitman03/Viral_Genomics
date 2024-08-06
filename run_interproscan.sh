### run_interproscan.sh START ###
#! /bin/bash

#$ -cwd
# error = Merged with joblog
#$ -o /u/scratch/m/madelain/interproscan_logs/interproscan_log.$JOB_ID

#$ -j y

#$ -l h_rt=24:00:00,h_data=16G # request 16G of memory

#$ -pe shared 8
#$ -M $USER
#$ -m bea

# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

conda activate /u/home/m/madelain/miniconda3/envs/openjdk_11.0.1
cd $SCRATCH/interproscan-5.68-100.0
./interproscan.sh -appl AntiFam-7.0,CDD-3.20,Coils-2.2.1,FunFam-4.3.0,Gene3D-4.3.0,Hamap-2023_05,MobiDBLite-2.0,NCBIfam-14.0,PANTHER-18.0,Pfam-37.0,PIRSF-3.10,PIRSR-2023_05,PRINTS-42.0,SFLD-4,SMART-9.0,SUPERFAMILY-1.75, -dp -f tsv -goterms -i /u/home/m/madelain/pcs_to_seq_noasterick.faa -o pcs_interproscan_alldb.tsv
# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "

### run_interproscan.sh END ###
