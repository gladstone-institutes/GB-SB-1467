#!/bin/bash
#$ -cwd    
#$ -o /gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/tmp/
#$ -e /gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/tmp/
#$ -pe smp 1
#$ -l mem_free=30G
#$ -l scratch=50G
#$ -l h_rt=02:00:00
#$ -j yes

# define variables
data_dir=/gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025
script_dir=$data_dir/scripts/GB-SB-1467
container_dir=$data_dir/assets
export APPTAINER_BINDPATH="$data_dir"

# run doublet finder for each sample
apptainer exec $container_dir/seurat-v5-2-1_soupx-1-6-2_doubletfinder.sif Rscript \
$script_dir/02_doublet_removal/02b_summarize_doublet_removal_mtcutoff15.R 

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"