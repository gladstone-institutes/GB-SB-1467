#!/usr/bin/env bash
#$ -cwd
#$ -o /gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/tmp/
#$ -e /gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/tmp/
#$ -pe smp 1
#$ -l mem_free=10G
#$ -l scratch=10G
#$ -l h_rt=00:30:00
#$ -t 1-3
#$ -j yes

data_dir=/gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025
result_dir=$data_dir/results/05_cluster_resolution_optimization
script_dir=$data_dir/scripts/GB-SB-1467
container_dir=$data_dir/assets
export APPTAINER_BINDPATH="$data_dir"

# Dynamically get the list of folders, excluding one
folder_list=($(find $result_dir -mindepth 1 -maxdepth 1 -type d | sort))

# Get the target folder for this task
folder="${folder_list[$((SGE_TASK_ID - 1))]}"
folder_basename=$(basename "$folder")

# Echo for debugging
echo "Processing folder: $folder"

# Run the script
singularity exec $container_dir/seurat-v5-2-1_soupx-1-6-2_doubletfinder.sif \
Rscript $script_dir/05_cluster_resolution_optimization/05c_cluster_resolution_optimization_summarize.R \
    --input $folder/per_sample_results \
    --output $folder \
    --output_prefix "gb_sb_1467_${folder_basename}"



## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"

######################## END ########################