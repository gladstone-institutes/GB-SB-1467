#!/usr/bin/env bash
#$ -cwd
#$ -o /gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/tmp/
#$ -e /gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/tmp/
#$ -pe smp 4
#$ -l mem_free=100G
#$ -l scratch=150G
#$ -l h_rt=40:00:00
#$ -t 1-16
#$ -j yes


# define paths
data_dir=/gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025
script_dir=$data_dir/scripts/GB-SB-1467
container_dir=$data_dir/assets
export APPTAINER_BINDPATH="$data_dir"

# Capture the arguments
pcs_to_use=$1

# run script
singularity exec $container_dir/seurat-v5-2-1_soupx-1-6-2_doubletfinder.sif \
Rscript $script_dir/05_cluster_resolution_optimization/05b_cluster_resolution_optimization_per_sample.R \
    --input $data_dir/results/05_cluster_resolution_optimization/gb_sb_1467_sketched_data.rds \
    --output $data_dir/results/05_cluster_resolution_optimization/${pcs_to_use}PCs \
    --output_prefix "gb_sb_1467_${pcs_to_use}pcs" \
    --mode "scRNA" \
    --cores ${NSLOTS} \
    --ndim $pcs_to_use \
    --metadata "SampleID" \
    --sample ${SGE_TASK_ID}


## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"

######################## END ########################