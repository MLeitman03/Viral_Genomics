### run_vcontact2.sh START ###
#! /bin/bash

#$ -cwd

# error = Merged with joblog
#$ -o /u/scratch/m/madelain/vcontact2_logs/vc_log.$JOB_ID

#$ -j y

#$ -l h_rt=24:00:00,h_data=8G # runtime and memory

#$ -pe shared 8
#$ -M $USER
#$ -m bea

# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

# read in arguments by flag and assign them to variables
while getopts "p:" opt; do
  case $opt in
    p) proteins_file="$OPTARG";;
    *) echo "Unknown option: -$opt" >&2; exit 1;;
  esac
done

# load the job environment:
. /u/local/Modules/default/init/modules.sh
module load anaconda3
# To see which versions of anaconda are available use: module av anaconda

conda activate vcontact2

vDir="/u/scratch/m/madelain/vcontact2/vcontact2_output_${JOB_ID}"
mkdir -p $vDir

# Run vContact2
vcontact2_gene2genome -p ${proteins_file} -o ${vDir}/protein_map.csv -s Prodigal-FAA
vcontact2 --raw-proteins ${proteins_file} --rel-mode 'Diamond' --proteins-fp ${vDir}/protein_map.csv \
--db 'ProkaryoticViralRefSeq201-Merged' --output-dir ${vDir} --pcs-mode MCL \
--vcs-mode ClusterONE --c1-bin ~/miniconda3/envs/vcontact2_env/bin/cluster_one-1.0.jar

# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "

### run_vcontact2.sh END ###
