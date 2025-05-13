#!/usr/bin/env Rscript

###############################################################################
## Project ID: GB-SB-1467
## Authors: Ayushi Agrawal
## Script Goal: Save the top 50 marker genes per cluster from Seurat's 
##              FindAllMarkers results to a CSV.
## Note: This script was run locally on Ayushi's laptop
###############################################################################


# load required packages
library(Seurat)
library(dplyr)

# set working directory
setwd(file.path("/Volumes/Jain-Boinformatics-Collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025",
                "results/06_clustering_cell_type/30PC_0.08res"))

# load the FindAllMarkers results 
markers_list <- read.csv("gb_sb_1467_30PC_0.08res_FindAllMarkers_marker_genes_per_cluster.csv")

# find top 50 markers genes per cluster
top50_per_cluster2 <- markers_list %>%
  filter(p_val_adj < 0.05) %>%              
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  arrange(desc(avg_log2FC), p_val) %>%
  slice_head(n = 50) %>%
  ungroup()

# save the marker list
top50_per_cluster %>%
  write.csv(., 
            "gb_sb_1467_30PC_0.08res_FindAllMarkers_top50_marker_genes_per_cluster.csv", 
            row.names = FALSE)


#END