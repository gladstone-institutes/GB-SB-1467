#!/usr/bin/env bash
#$ -cwd
#$ -o /gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/tmp/
#$ -e /gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/tmp/
#$ -pe smp 1
#$ -l mem_free=200G
#$ -l scratch=100G
#$ -l h_rt=10:00:00
#$ -t 1-21
#$ -j yes


# define paths
data_dir=/gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025
script_dir=$data_dir/scripts/GB-SB-1467
container_dir=$data_dir/assets
export APPTAINER_BINDPATH="$data_dir"

# Define array of resolution-ndim combos
param_list=(
  "30 0.5"
  "15 0.08"
  "15 0.1"
  "15 0.2"
  "15 0.3"
  "15 0.4"
  "15 0.5"
  "15 0.6"
  "20 0.08"
  "20 0.1"
  "20 0.2"
  "20 0.3"
  "20 0.4"
  "20 0.5"
  "20 0.6"
  "30 0.08"
  "30 0.1"
  "30 0.2"
  "30 0.3"
  "30 0.4"
  "30 0.6"
)

# Get combo for this task
params="${param_list[$SGE_TASK_ID-1]}"
ndim=$(echo $params | cut -d' ' -f1)
res=$(echo $params | cut -d' ' -f2)

echo "Running for PCs=$ndim and resolution=$res"

# Run the R script via Singularity
singularity exec $container_dir/seurat-v5-2-1_soupx-1-6-2_doubletfinder-presto.sif \
Rscript $script_dir/06_clustering_cell_type/06b_plots_findallmarkers.R \
    --input $data_dir/results/06_clustering_cell_type/${ndim}PC_${res}res/gb_sb_1467_${ndim}PC_${res}res_clustered_and_cell_typed.rds \
    --output $data_dir/results/06_clustering_cell_type/${ndim}PC_${res}res \
    --output_prefix "gb_sb_1467_${ndim}PC_${res}res"


## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"

######################## END ########################