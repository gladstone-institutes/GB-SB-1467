#!/usr/bin/env bash

###############################################################################
## Project ID: GB-SB-1467
## Authors: Ayushi Agrawal
##
## Script Goal: Wrapper script to perform pseudobulk differential expression
## run locally on Ayushi's laptop
###############################################################################

# create the output folder
mkdir -p /Volumes/Jain-Boinformatics-Collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/results/07_pseudobulk_de

cd /Volumes/Jain-Boinformatics-Collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/scripts/GB-SB-1467/07_pseudobulk_de

# run analysis for subcluster 2 ---------------------------------
Rscript 07_pseudobulk_de.R \
--input '~/Downloads/Isha_proj/gb_sb_1467_30PC_0.08res_clustered_and_cell_typed.rds' \
--output '/Volumes/Jain-Boinformatics-Collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/results/07_pseudobulk_de' \
--output_prefix 'gb_sb_1467_30PC_0.08res' \
--predictor "Condition" \
--sampleid "SampleID" \
--cell_annotation "seurat_clusters" \
--pathway_db '~/Downloads/Mm_20250415.RData' \
--p_val_cutoff 0.01 \
--fold_change_cutoff 1 \
--species "mouse" > /Volumes/Jain-Boinformatics-Collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/results/07_pseudobulk_de/logs.txt 2>&1


## END 