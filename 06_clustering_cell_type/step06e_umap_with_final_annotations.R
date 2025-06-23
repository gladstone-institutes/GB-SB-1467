#!/usr/bin/env Rscript

###############################################################################
## Project ID: GB-SB-1467
## Authors: Ayushi Agrawal
## Script Goal: Generate UMAP with confirmed annotations
## Note: This script was run locally on Ayushi's laptop
###############################################################################


# load required packages
library(Seurat)
library(dplyr)

# set working directory
setwd(file.path("/Volumes/Jain-Boinformatics-Collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025",
                "results/06_clustering_cell_type/30PC_0.08res"))

# create output folder
outdir <- file.path("/Volumes/Jain-Boinformatics-Collaboration",
                    "sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/results/paper/figures")
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)


# load the Seurat object 
data <- readRDS("gb_sb_1467_30PC_0.08res_clustered_and_cell_typed.rds")

# Full cell type annotations (per cluster)
celltype_annotations <- c(
  "Upper-Layer Excitatory Neurons",      # Cluster 0
  "Hypothalamic/Thalamic Neurons",       # Cluster 1
  "Astrocytes",                          # Cluster 2
  "Striatum Medium Spiny Neurons",       # Cluster 3
  "Progenitors",                         # Cluster 4
  "Oligodendrocyte Progenitor Cells",    # Cluster 5
  "MGE-interneuron",                     # Cluster 6
  "Deeper-Layer Excitatory Neurons",     # Cluster 7
  "Mural Cells",                         # Cluster 8
  "Microglia",                           # Cluster 9
  "LGE/CGE Interneurons",                # Cluster 10
  "Ependymal Cells",                     # Cluster 11
  "Endothelial Cells",                   # Cluster 12
  "Hippocampal Neurons",                 # Cluster 13
  "Unknown",                             # Cluster 14
  "Purkinje Cells"                       # Cluster 15
)

# Abbreviated cluster labels
celltype_abbr <- c(
  "0_UL-ExN",      # Cluster 0
  "1_HypoTh_Neur", # Cluster 1
  "2_Astro",       # Cluster 2
  "3_SMSN",        # Cluster 3
  "4_Prog",        # Cluster 4
  "5_OPC",         # Cluster 5
  "6_MGE-Int",     # Cluster 6
  "7_DL-ExN",      # Cluster 7
  "8_Mural",       # Cluster 8
  "9_Micro",       # Cluster 9
  "10_LGE-Int",    # Cluster 10
  "11_Epend",      # Cluster 11
  "12_Endo",       # Cluster 12
  "13_HPC-N",      # Cluster 13
  "14_Unknown",    # Cluster 14
  "15_Purk"        # Cluster 15
)

# Write legend to CSV
legend_df <- data.frame(
  Cluster = 0:15,
  Abbreviation = celltype_abbr,
  Full_Name = celltype_annotations,
  stringsAsFactors = FALSE
)
write.csv(legend_df, "celltype_abbreviation_legend.csv", row.names = FALSE)

# Assign full annotations
names(celltype_annotations) <- as.character(0:15)
data$cell_type <- unname(celltype_annotations[as.character(Idents(data))])
data$seurat_clusters_cell_type <- paste0(data$seurat_clusters, "_", data$cell_type)

# Factor with sorted cluster order
data$seurat_clusters_cell_type <- factor(
  data$seurat_clusters_cell_type,
  levels = unique(data$seurat_clusters_cell_type[order(as.numeric(data$seurat_clusters))])
)

# Assign abbreviated annotations
names(celltype_abbr) <- as.character(0:15)
data$cell_type_abbr <- unname(celltype_abbr[as.character(Idents(data))])
data$seurat_clusters_abbr <- factor(data$cell_type_abbr, levels = celltype_abbr)


# ---- UMAP Plots with Full Annotations ----
pdf("gb_sb_1467_30PC_0.08res_finalized_annotations_labeled_UMAP.pdf", 
    width = 14, 
    height = 10)
print(DimPlot(data, 
              group.by = "seurat_clusters_cell_type", 
              raster = FALSE, 
              label = TRUE, 
              reduction = "umap"
))
dev.off()

pdf("gb_sb_1467_30PC_0.08res_finalized_annotations_labeled_condition_split_UMAP.pdf", 
    width = 30, 
    height = 10)
print(DimPlot(data, 
              group.by = "seurat_clusters_cell_type", 
              split.by = "Condition", 
              raster = FALSE, 
              label = TRUE, 
              reduction = "umap"
))
dev.off()


# ---- UMAP Plots with Abbreviated Annotations ----
pdf("gb_sb_1467_30PC_0.08res_finalized_abbreviated_annotations_labeled_UMAP.pdf", 
    width = 12, 
    height = 10)
print(DimPlot(data, 
              group.by = "seurat_clusters_abbr", 
              raster = FALSE, 
              label = TRUE,
              reduction = "umap"
))
dev.off()

pdf("gb_sb_1467_30PC_0.08res_finalized_abbreviated_annotations_labeled_condition_split_UMAP.pdf", 
    width = 28, 
    height = 10)
print(DimPlot(data, 
              group.by = "seurat_clusters_abbr", 
              split.by = "Condition", 
              raster = FALSE, 
              label = TRUE, 
              reduction = "umap"
))
dev.off()

# Subfigure A (Abbreviated Annotations)
pdf(file.path(outdir,"subfig_A.pdf"), 
    width = 12, 
    height = 10)
print(DimPlot(data, 
              group.by = "seurat_clusters_abbr", 
              raster = FALSE, 
              label = TRUE, 
              reduction = "umap"
))
dev.off()

pdf(file.path(outdir,"subfig_supp_A.pdf"), 
    width = 28, 
    height = 10)
print(DimPlot(data, 
              group.by = "seurat_clusters_abbr", 
              split.by = "Condition", 
              raster = FALSE, 
              label = TRUE, 
              reduction = "umap"))
dev.off()


# Subfigure A (Abbreviated Annotations + lower point size)
pdf(file.path(outdir,"subfig_A_v2.pdf"), 
    width = 7, 
    height = 6)
plt <- DimPlot(data, 
               group.by = "seurat_clusters_abbr", 
               raster = FALSE, 
               label = TRUE, 
               label.size = 3,
               reduction = "umap") 
plt[[1]]$layers[[1]]$aes_params$alpha = 0.2
plt[[1]]$layers[[1]]$aes_params$size = 0.2
plt[[1]]$layers[[1]]$aes_params$shape = 16
print(plt)
dev.off()

pdf(file.path(outdir,"subfig_supp_A_v2.pdf"), 
    width = 18, 
    height = 6)
plt <- DimPlot(data, 
               group.by = "seurat_clusters_abbr", 
               split.by = "Condition", 
               raster = FALSE, 
               label = TRUE, 
               label.size = 3,
               reduction = "umap") 
plt[[1]]$layers[[1]]$aes_params$alpha = 0.2
plt[[1]]$layers[[1]]$aes_params$size = 0.2
plt[[1]]$layers[[1]]$aes_params$shape = 16
print(plt)
dev.off()


# END