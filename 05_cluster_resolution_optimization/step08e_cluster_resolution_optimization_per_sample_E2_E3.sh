#!/usr/bin/env bash
#$ -cwd
#$ -o /gladstone/bioinformatics/adnetworksppg/Project_1_Huang/ZL02/tmp/
#$ -e /gladstone/bioinformatics/adnetworksppg/Project_1_Huang/ZL02/tmp/
#$ -pe smp 4
#$ -l mem_free=100G
#$ -l scratch=150G
#$ -l h_rt=40:00:00
#$ -t 1-65
#$ -j yes

data_dir=/gladstone/bioinformatics/adnetworksppg/Project_1_Huang/ZL02
script_dir=$data_dir/scripts/YH_ZL02
container_dir=$data_dir/assets/containers
export APPTAINER_BINDPATH="$data_dir"

singularity exec $container_dir/seurat-v5_soupx-1-6-1_doubletfinder.sif \
Rscript $script_dir/08_cluster_resolution_optimization/08b_cluster_resolution_optimization_per_sample.R \
    --input $data_dir/results/08_cluster_resolution_optimization_no_cellbender/E2_E3/YH_ZL02_no_cellbender_E2_E3_sketched_data.rds \
    --output $data_dir/results/08_cluster_resolution_optimization_no_cellbender/E2_E3 \
    --output_prefix "YH_ZL02_no_cellbender_E2_E3" \
    --mode "scRNA" \
    --cores ${NSLOTS} \
    --ndim 20 \
    --metadata "MouseID" \
    --sample ${SGE_TASK_ID}


## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"

######################## END ########################