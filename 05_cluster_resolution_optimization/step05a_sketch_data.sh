#!/usr/bin/env bash
#$ -cwd
#$ -o /gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/tmp/
#$ -e /gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/tmp/
#$ -pe smp 2
#$ -l mem_free=50G
#$ -l scratch=100G
#$ -l h_rt=02:00:00
#$ -q gpu.q
#$ -j yes


# define paths
data_dir=/gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025
script_dir=$data_dir/scripts/GB-SB-1467
container_dir=$data_dir/assets
export APPTAINER_BINDPATH="$data_dir"


# run sketching 
singularity exec $container_dir/seurat-v5-2-1_soupx-1-6-2_doubletfinder.sif \
Rscript $script_dir/05_cluster_resolution_optimization/05a_sketch_data.R \
    --input $data_dir/results/04_merge_and_visualize/30PCs/gb_sb_1467_30pcs_merged_data_sct_pca_umap.rds \
    --output $data_dir/results/05_cluster_resolution_optimization \
    --output_prefix "gb_sb_1467" \
    --mode "scRNA" \
    --cores ${NSLOTS} \
    --sketch_percent 25


## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"

######################## END ########################