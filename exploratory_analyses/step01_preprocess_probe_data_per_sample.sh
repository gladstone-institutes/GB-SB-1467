#!/bin/bash
#$ -cwd    
#$ -o /gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/tmp/
#$ -e /gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/tmp/
#$ -pe smp 2
#$ -l mem_free=100G
#$ -l scratch=100G
#$ -l h_rt=15:00:00
#$ -j yes

# Define paths
data_dir=/gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025
script_dir=$data_dir/scripts/GB-SB-1467
container_dir=$data_dir/assets
export APPTAINER_BINDPATH="$data_dir"

# Run the R script
singularity exec $container_dir/seurat-v5-2-1_soupx-1-6-2_doubletfinder.sif Rscript \
$script_dir/exploratory_analyses/01_preprocess_probe_data_per_sample.R \
--samplesheet $script_dir/02_doublet_removal/input.csv \
--metadata $data_dir/data/V25-03_Samples.csv \
--output_dir $data_dir/results/exploratory_analyses/01_preprocess_probe_data_per_sample \
--output_prefix "gb_sb_1467_probes_30pcs" \
--subset_cells $data_dir/results/04_merge_and_visualize/30PCs/gb_sb_1467_30pcs_post_qc_cell_ids.csv


## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
