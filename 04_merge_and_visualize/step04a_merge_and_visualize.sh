#!/bin/bash
#$ -cwd    
#$ -o /gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/tmp/
#$ -e /gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/tmp/
#$ -pe smp 2
#$ -l mem_free=200G
#$ -l scratch=100G
#$ -l h_rt=05:00:00
#$ -j yes
#$ -t 1-3

# Define paths
data_dir=/gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025
script_dir=$data_dir/scripts/GB-SB-1467
container_dir=$data_dir/assets
export APPTAINER_BINDPATH="$data_dir"

# Define the list of PCs
pcs_list=(15 20 30)

# Map the task ID (1-based) to the correct value in pcs_list (0-based indexing in bash arrays)
npcs=${pcs_list[$((SGE_TASK_ID - 1))]}
output_dir="$data_dir/results/04_merge_and_visualize/${npcs}PCs"
output_prefix="gb_sb_1467_${npcs}pcs"

# Run the R script
singularity exec $container_dir/seurat-v5-2-1_soupx-1-6-2_doubletfinder.sif Rscript \
$script_dir/04_merge_and_visualize/04a_merge_and_visualize.R \
--input_dir $data_dir/results/02_doublet_removal \
--metadata $data_dir/data/V25-03_Samples.csv \
--output_dir "$output_dir" \
--output_prefix "$output_prefix" \
--project "gb_sb_1467" \
--npcs "$npcs"


## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"