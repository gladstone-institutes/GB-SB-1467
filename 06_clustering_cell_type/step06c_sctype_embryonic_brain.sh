#!/usr/bin/env bash
#$ -cwd
#$ -o /gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/tmp/
#$ -e /gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/tmp/
#$ -pe smp 2
#$ -l mem_free=200G
#$ -l scratch=100G
#$ -l h_rt=04:00:00
#$ -j yes

# define paths
data_dir=/gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025
script_dir=$data_dir/scripts/GB-SB-1467
container_dir=$data_dir/assets
export APPTAINER_BINDPATH="$data_dir"

# set variables
ndim=30
res=0.08

echo "Running for PCs=$ndim and resolution=$res"

# Run the R script via Singularity
singularity exec $container_dir/seurat-v5-2-1_soupx-1-6-2_doubletfinder.sif \
Rscript $script_dir/06_clustering_cell_type/06c_sctype_embryonic_brain.R \
    --input $data_dir/results/06_clustering_cell_type/30PC_0.08res/gb_sb_1467_30PC_0.08res_clustered_and_cell_typed.rds \
    --output $data_dir/results/06_clustering_cell_type \
    --output_prefix "gb_sb_1467" \
    --ndim $ndim \
    --resolution $res \
    --tissue "Embryonic brain" \
    --custom_db $container_dir/Cell_marker_Mouse_20Apr2025.xlsx


## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"

######################## END ########################