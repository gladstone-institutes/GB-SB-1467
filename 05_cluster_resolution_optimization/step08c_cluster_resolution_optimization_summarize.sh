#!/usr/bin/env bash
#$ -cwd
#$ -o /gladstone/bioinformatics/adnetworksppg/Project_1_Huang/ZL02/tmp/
#$ -e /gladstone/bioinformatics/adnetworksppg/Project_1_Huang/ZL02/tmp/
#$ -pe smp 1
#$ -l mem_free=20G
#$ -l scratch=50G
#$ -l h_rt=02:00:00
#$ -j yes
#$ -P neuroppg


data_dir=/gladstone/bioinformatics/adnetworksppg/Project_1_Huang/ZL02
script_dir=$data_dir/scripts/YH_ZL02
container_dir=$data_dir/assets/containers
export APPTAINER_BINDPATH="$data_dir"

singularity exec $container_dir/seurat-v5_soupx-1-6-1_doubletfinder.sif \
Rscript $script_dir/08_cluster_resolution_optimization/08c_cluster_resolution_optimization_summarize.R \
    --input $data_dir/results/08_cluster_resolution_optimization_no_cellbender/per_sample_results \
    --output $data_dir/results/08_cluster_resolution_optimization_no_cellbender \
    --output_prefix "YH_ZL02_no_cellbender"


## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"

######################## END ########################