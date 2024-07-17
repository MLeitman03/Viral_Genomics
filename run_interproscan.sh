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

# Load Java 11 module
module load java/jdk-11.0.14

# Check Java version to ensure correct version is loaded
java -version

# Run InterProScan with adjusted memory settings
export JAVA_OPTS="-Xms2g -Xmx4g"
./interproscan-5.68-100.0/interproscan.sh -i predicted_proteins.faa -f tsv -o interproscan_output.tsv

# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "

### run_interproscan.sh END ###
