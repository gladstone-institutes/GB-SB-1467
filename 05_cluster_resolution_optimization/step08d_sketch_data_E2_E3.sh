#!/usr/bin/env bash
#$ -cwd
#$ -o /gladstone/bioinformatics/adnetworksppg/Project_1_Huang/ZL02/tmp/
#$ -e /gladstone/bioinformatics/adnetworksppg/Project_1_Huang/ZL02/tmp/
#$ -pe smp 1
#$ -l mem_free=100G
#$ -l scratch=100G
#$ -l h_rt=02:00:00
#$ -q gpu.q
#$ -j yes
#$ -P neuroppg

data_dir=/gladstone/bioinformatics/adnetworksppg/Project_1_Huang/ZL02
script_dir=$data_dir/scripts/YH_ZL02
container_dir=$data_dir/assets/containers
export APPTAINER_BINDPATH="$data_dir"


singularity exec $container_dir/seurat-v5_soupx-1-6-1_doubletfinder.sif \
Rscript $script_dir/08_cluster_resolution_optimization/08a_sketch_data.R \
    --input $data_dir/results/07_merge_and_visualize/no_cellbender/YH_ZL02_merged_data_sct_pca_umap.rds \
    --output $data_dir/results/08_cluster_resolution_optimization_no_cellbender/E2_E3 \
    --output_prefix "YH_ZL02_no_cellbender_E2_E3" \
    --metadata $data_dir/data/Final_YH_ZL02_snRNAseq_MAPT_fApoE_mouse_list_added_corrected_240816_AA.csv \
    --genotypes_to_keep "MAPT-fE2,MAPT-fE3" \
    --mode "scRNA" \
    --cores ${NSLOTS} \
    --sketch_size 202998 #use 27.44% of 739,866 cells for sketching


## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"

######################## END ########################