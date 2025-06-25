#!/bin/bash
#$ -cwd    
#$ -o /gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/tmp/
#$ -e /gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/tmp/
#$ -l mem_free=200G
#$ -l scratch=50G
#$ -l h_rt=05:00:00
#$ -j yes


# Define paths
data_dir=/gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025
script_dir=$data_dir/scripts/GB-SB-1467
container_dir=$data_dir/assets
export APPTAINER_BINDPATH="$data_dir"



# Run the R script
singularity exec $container_dir/seurat-v5-2-1_soupx-1-6-2_doubletfinder.sif Rscript \
$script_dir/08_between_cluster_association_genotype/08_between_cluster_association_genotype_ref_WT_wo_B3.R
 




## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"