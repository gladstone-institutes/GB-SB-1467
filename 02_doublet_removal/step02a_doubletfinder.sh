#!/bin/bash
#$ -cwd    
#$ -o /gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/tmp/
#$ -e /gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/tmp/
#$ -pe smp 1
#$ -l mem_free=40G
#$ -l scratch=50G
#$ -l h_rt=08:00:00
#$ -t 1-16 # 16 samples
#$ -j yes

# define variables
data_dir=/gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025
script_dir=$data_dir/scripts/GB-SB-1467
container_dir=$data_dir/assets
export APPTAINER_BINDPATH="$data_dir"

# get the approriate sample details
input_file="$script_dir/02_doublet_removal/input.csv"
row_number=$((SGE_TASK_ID + 1))
row_value=$(awk -F, "NR==$row_number {print \$1}" "$input_file")

task_sample=$(basename "$row_value")
echo -e "Working on sample $task_sample"

# create the output directory
mkdir -p $data_dir/results/02_doublet_removal/

# run doublet finder for each sample
apptainer exec $container_dir/seurat-v5-2-1_soupx-1-6-2_doubletfinder.sif Rscript \
$script_dir/02_doublet_removal/02a_doubletfinder.R \
--cellranger_h5 $row_value/count/sample_filtered_feature_bc_matrix.h5 \
--doublet_rate_estimate $script_dir/02_doublet_removal/10X_doublet_rate_estimation.csv \
--output $data_dir/results/02_doublet_removal/$task_sample


## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"