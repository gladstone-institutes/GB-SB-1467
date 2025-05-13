#!/usr/bin/env Rscript

###############################################################################
## Project ID: GB-SB-1467
## Authors: Ayushi Agrawal
## Script Goal: Generate figures for paper
## Note: This script was run locally on Ayushi's laptop
###############################################################################

# Load libraries
library(Seurat)
library(dplyr)
library(ggplot2)

# set the working directory
setwd(file.path("/Volumes/Jain-Boinformatics-Collaboration",
                "sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/results"))

# load the Seurat object
data <- readRDS("~/Documents/downloads/Isha_proj/gb_sb_1467_30PC_0.08res_clustered_and_cell_typed_pathway_modules_added.rds")

# get module score columns from meta data
cols <- colnames(data@meta.data)[grep(".v2024.1.Mm1$", colnames(data@meta.data))]

pdf(file = "exploratory_analyses/evaluate_cell_death_pathways/gb_sb_1467_30PC_0.08res_module_score_featureplots.pdf",
    width = 10,
    height = 7)
for (c in cols){
  print(FeaturePlot(data, 
                    features = c, 
                    raster = FALSE, 
                    order = TRUE, 
                    label = TRUE, 
                    reduction = "umap"))
}
dev.off()

