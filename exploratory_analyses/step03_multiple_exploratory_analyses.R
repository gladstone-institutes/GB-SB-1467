#!/usr/bin/env Rscript

###############################################################################
## Project ID: GB-SB-1467
## Authors: Ayushi Agrawal
## Script Goal: Exploratory analyses
## Note: This script was run locally on Ayushi's laptop
###############################################################################


# load required packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(purrr)

# set working directory
setwd("/Volumes/Jain-Boinformatics-Collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025")

# load the input data 
data <- readRDS("results/06_clustering_cell_type/30PC_0.08res/gb_sb_1467_30PC_0.08res_clustered_and_cell_typed.rds")




# 1. Are the samples male or female? -------------------------------------------
# Reference: 1. https://doi.org/10.1038/s41598-019-50731-x
#            2. https://doi.org/10.2144/btn-2022-0077
# Kdm5c (X-chromosome-specific) and Kdm5d (Y-chromosome-specific)

# get y-chr genes from 10X genomics bed files
bed <- read.delim("assets/Chromium_Mouse_Transcriptome_Probe_Set_v1.1.1_GRCm39-2024-A.bed", 
                  header = FALSE, sep = "\t", stringsAsFactors = FALSE) %>%
  as.data.frame()
genes_chrY <- bed[(bed$V1 == "chrY"),"V4"]
genes_chrY <- unique(sub("^[^|]*\\|([^|]*)\\|.*$", "\\1", genes_chrY))

# generate dotplot
pdf(file = "results/exploratory_analyses/gb_sb_1467_30PC_0.08res_chrx_chry_genes_group_sampleid_dotplot.pdf",
    width = 10)
print(DotPlot(data, 
              features = c(genes_chrY, "Kdm5c", "Kdm6a"),
              group.by = "SampleID") + 
        RotatedAxis())   
dev.off()

genes_chrX <- bed[(bed$V1 == "chrX"),"V4"]
genes_chrX <- unique(sub("^[^|]*\\|([^|]*)\\|.*$", "\\1", genes_chrX))
# too many genes to plot
length(genes_chrX)
#[1] 747

# Conclusions
# 9 Male mice: samples A, C, D, F-K
# 7 Female mice: samples B, E, L-P




# 2. Identify catecholaminergic neurons  ---------------------------------------
catecholaminergic_neuron_markers <- c("Th", "Dbh", "Pnmt")

pdf(file = "results/exploratory_analyses/gb_sb_1467_30PC_0.08res_catecholaminergic_neuron_markers_dotplot.pdf",
    height = 7,
    width = 20)
print(DotPlot(data, split.by = "Condition",
              features = catecholaminergic_neuron_markers,
              cols = "RdBu") + 
        RotatedAxis() + 
        coord_flip() +
        theme(legend.position = "bottom"))   
dev.off()

pdf(file = "results/exploratory_analyses/gb_sb_1467_30PC_0.08res_catecholaminergic_neuron_markers_featureplot.pdf",
    width = 25,
    height = 14)
print(FeaturePlot(data, 
                  features= catecholaminergic_neuron_markers, 
                  split.by = "Condition",
                  raster = FALSE, 
                  order = TRUE, 
                  label = TRUE, 
                  reduction = "umap"))
print(FeaturePlot(data, 
                  features= catecholaminergic_neuron_markers, 
                  raster = FALSE, 
                  order = TRUE, 
                  label = TRUE, 
                  reduction = "umap"))
dev.off()




# 3. Add module scores for cell death ------------------------------------------
# Path to MSigDB gene sets
pathways <- list.files(path = "assets/MSigDB_gene_sets", pattern = "\\.tsv$", full.names = TRUE)

for (p in pathways){
  # Read gene set file
  meta <- read.delim(p, header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
  kv <- setNames(meta$V2, meta$V1)
  gene_symbols <- unique(trimws(unlist(strsplit(kv["GENE_SYMBOLS"], ","))))
  gene_symbols <- gene_symbols[nzchar(trimws(gene_symbols))]
  
  # Extract base name
  meta_colname <- basename(p) |> tools::file_path_sans_ext()
  
  # Add module score (creates e.g., HALLMARK_XYZ.v2024.1.Mm12345_1)
  data <- AddModuleScore(
    data,
    features = list(gene_symbols),
    ctrl = 100,
    name = meta_colname
  )
}

# Identify pathway score columns
pathway_score_cols <- paste0(basename(pathways) |> tools::file_path_sans_ext(), "1")
cleaned_names <- sub("\\.v2024\\.1\\.Mm[0-9]+$", "", pathway_score_cols)

# create output directory
if (!(dir.exists("results/exploratory_analyses/evaluate_cell_death_pathways"))) {
  dir.create("results/exploratory_analyses/evaluate_cell_death_pathways", 
             recursive = T, 
             showWarnings = F)
}

# save the metadata
write.csv(data@meta.data,
          file = paste0("results/exploratory_analyses/evaluate_cell_death_pathways/",
                        "gb_sb_1467_30PC_0.08res_clustered_and_cell_typed_pathway_modules_added_metadata.csv")
)

# create a dotplot f the pathway scores
pdf(file = paste0("results/exploratory_analyses/evaluate_cell_death_pathways/",
                  "gb_sb_1467_30PC_0.08res_cell_death_pathways_dotplot.pdf"),
    height = 10,
    width = 20)
print(DotPlot(data,
              pathway_score_cols) +
        RotatedAxis() +
        scale_x_discrete(labels = setNames(cleaned_names, pathway_score_cols))
)
dev.off()

# save the rds object
saveRDS(data,
        file = paste0("results/exploratory_analyses/evaluate_cell_death_pathways/",
                      "gb_sb_1467_30PC_0.08res_clustered_and_cell_typed_pathway_modules_added.rds")
)




# Write session info
writeLines(
  capture.output(sessionInfo()),
  file.path("results/exploratory_analyses", "sessionInfo.txt")
)
print("***** Script Complete *****")

# END