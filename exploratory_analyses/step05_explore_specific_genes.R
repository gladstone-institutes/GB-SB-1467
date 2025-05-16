#!/usr/bin/env Rscript

###############################################################################
## Project ID: GB-SB-1467
## Authors: Reuben Thomas, Ayushi Agrawal
## Script Goal: Gene expression of specific genes and overall DEG comparisons
## Usage: Rscript step05_explore_specific_genes.R
## Note: This script was run locally on Ayushi's laptop
###############################################################################


# Setup ------------------------------------------------------------------------
# set working directory
setwd(file.path("/Volumes/Jain-Boinformatics-Collaboration",
                "sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025",
                "results/07_pseudobulk_de"))

# load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(tidyverse)
library(pheatmap)
library(readxl)
library(stringr)
library(muscat)
library(SummarizedExperiment)
library(SingleCellExperiment)


# load gene sets of known pathways ---------------------------------------------
parse_msigdb <- function(file) {
  df <- read.delim(file, header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
  kv <- setNames(df$V2, df$V1)
  genes <- unique(trimws(unlist(strsplit(kv["GENE_SYMBOLS"], ","))))
  genes[nzchar(genes)]
}

serine_bio_genes <- parse_msigdb("../../assets/MSigDB_gene_sets/REACTOME_SERINE_BIOSYNTHESIS.v2024.1.Mm.tsv")
apoptosis_genes <- parse_msigdb("../../assets/MSigDB_gene_sets/HALLMARK_APOPTOSIS.v2024.1.Mm.tsv")
atf4_genes <- str_to_title(parse_msigdb("../../assets/MSigDB_gene_sets/ATF4_Q2.v2024.1.Hs.tsv"))
atf3_genes <- str_to_title(parse_msigdb("../../assets/MSigDB_gene_sets/ATF3_Q6.v2024.1.Hs.tsv"))

serine_deprivation_induced_genes <- read_xlsx("../../assets/NIHMS1742915-supplement-2.xlsx", sheet = 1, skip = 3) %>%
  mutate(Symbol = str_to_title(Symbol)) %>%
  pull(Symbol)


# list of all gene lists
genelist <- list(mouse_reactome_serine_biosynthesis_genes = serine_bio_genes,
                 mouse_hallmark_apoptosis_genes = apoptosis_genes,
                 human_ATF4_target_genes = atf4_genes, 
                 human_ATF3_target_genes = atf3_genes, 
                 PMC8491098_supp2_serine_deprivation_induced_genes = serine_deprivation_induced_genes
)


# load DEG lists and Seurat object ---------------------------------------------
# Find all DEG csv files 
deg_dir <- "gene_expression_associations"
deg_files <- list.files(deg_dir, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)

# Function to read and attach cluster and comparison info
read_deg_file <- function(file) {
  df <- read.csv(file)
  if (!all(c("p_val", "logFC") %in% colnames(df))) return(NULL)
  parts <- strsplit(file, "/")[[1]]
  df$comparison <- parts[length(parts) - 2]
  df
}

# Load all DEGs
all_deg_data <- map_dfr(deg_files, read_deg_file)

# read in the data
data <- readRDS("~/Downloads/Isha_proj/gb_sb_1467_30PC_0.08res_clustered_and_cell_typed.rds")


# Prepare the pseudobulk data --------------------------------------------------
predictors <- "Condition"
sampleid <- "SampleID"
clusters <- "seurat_clusters"

# Store the meta-data for each cell in the PhenoData object
PhenoData <- data@meta.data
PhenoData <- PhenoData[,!grepl("^DF\\.classifications_", colnames(PhenoData))]
PhenoData <- PhenoData[,!grepl("^pANN_", colnames(PhenoData))]

# Create SingleCellExperiment object
sce <- SummarizedExperiment(assays = list(counts = data[["RNA"]]$counts), 
                            colData = PhenoData)
sce <- as(sce, "SingleCellExperiment")

# Prep this object for subsequent aggregation analyses
sce <- prepSCE(sce,
               kid = clusters, # sub-population assignments
               gid = predictors, # group IDs
               sid = sampleid, # sample IDs
               drop = FALSE
)

# Aggregate counts across cells for each mouse (sample_id) within each cluster (cluster_id)
pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))

# Flatten pseudobulk array into matrix
cluster_ids <- names(assays(pb))

pb_mat <- do.call(cbind, lapply(cluster_ids, function(cl) {
  mat <- assays(pb)[[cl]]  # This is a matrix: genes × samples
  colnames(mat) <- paste0("cl", cl, "_", colnames(mat))  # Rename columns
  mat
}))

# Filter out columns (cluster-samples) with all zeros
lib_sizes <- colSums(pb_mat)
pb_mat_filtered <- pb_mat[, lib_sizes > 0]

# Compute logCPM
dge <- edgeR::DGEList(counts = pb_mat_filtered)
logCPM <- edgeR::cpm(dge, log = TRUE, prior.count = 1)


# Generate plots for each gene set ---------------------------------------------
# Main loop
for (idx in 1:length(genelist)) {
  output_prefix <- names(genelist)[idx]
  outdir <- paste0("../exploratory_analyses/evaluate_", output_prefix)
  if (!dir.exists(outdir)) dir.create(outdir)
  
  genes_to_plot <- genelist[[idx]]
  genes_found <- intersect(genes_to_plot, rownames(data))
  genes_found <- intersect(genes_found, rownames(logCPM))
  genes_missing <- setdiff(genes_to_plot, genes_found)
  
  if (length(genes_found) == 0) {
    warning(paste("No valid genes found in Seurat object for", output_prefix, "- skipping plot."))
    next
  }
  
  # If too many genes, select the 100 most variable ones based on logCPM
  if (length(genes_found) > 100) {
    message(paste("Gene list too long for", output_prefix, "- selecting top 100 most variable genes based on logCPM"))
    # Select top 100 genes by variance
    sub_logCPM <- logCPM[genes_found, , drop = FALSE]
    gene_variances <- rowVars(as.matrix(sub_logCPM))
    top_variable_genes <- names(sort(gene_variances, decreasing = TRUE))[1:100]
    genes_found <- top_variable_genes
  }
  
  if (length(genes_missing) > 0) {
    warning(paste("Some genes not found in Seurat object for", output_prefix, ":", paste(genes_missing, collapse = ", ")))
  }
  
  
  # ---- Dotplot for pseudobulk DE results
  filtered_deg <- all_deg_data %>%
    filter(gene %in% genes_to_plot) %>%
    filter(p_adj.loc < 0.01 & abs(logFC) > 1) %>%
    mutate(Significance = -log10(p_adj.loc))
  
  gheight <- max(6, length(unique(filtered_deg$gene)) * 0.3)
  
  pdf(file.path(outdir, paste0(output_prefix, "_pseudobulkDE_results_dotplot.pdf")),
      width = 10, height = gheight)
  print(
    ggplot(filtered_deg, aes(y = factor(gene, levels = unique(gene)), x = as.factor(cluster_id),
                             size = Significance, color = logFC)) +
      geom_point(stroke = 0.4, alpha = 0.95) +
      scale_color_gradientn(colors = c("#2c7bb6", "#abd9e9", "gray", "#fdae61", "#d7191c")) +
      facet_grid(. ~ comparison, scales = "free_y", space = "free_y") +
      labs(
        title = paste0("Differential Expression of ", str_replace_all(output_prefix, "_", " ")),
        y = str_replace_all(output_prefix, "_", " "),
        x = "Comparison",
        color = "logFC",
        size = "-log10(FDR)"
      ) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            strip.background = element_rect(fill = "gray90", color = "black", linewidth = 1),
            panel.spacing = unit(1, "lines"))
  )
  dev.off()
  
  
  # ---- Heatmap of pseudobulk logCPM
  annotation_df <- data.frame(cluster_sample = colnames(logCPM)) %>%
    separate(cluster_sample, into = c("cluster", "sample"), sep = "_", extra = "merge") %>%
    left_join(as.data.frame(colData(sce)) %>% rownames_to_column("cell") %>% distinct(sample_id, .keep_all = TRUE), 
              by = c("sample" = "sample_id")) %>%
    distinct(sample, cluster, group_id) 
  
  # Fix rownames to match cluster_sample
  annotation_df$cluster_sample <- paste0(annotation_df$cluster, "_", annotation_df$sample)
  rownames(annotation_df) <- annotation_df$cluster_sample
  
  heatmap_mat <- logCPM[genes_found, , drop = FALSE]
  annotation_df <- annotation_df[colnames(heatmap_mat), , drop = FALSE]
  
  # Define annotation color mappings
  group_levels <- unique(annotation_df$group_id)
  cluster_levels <- unique(annotation_df$cluster)
  
  ann_colors <- list(
    group_id = setNames(colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(group_levels)), group_levels),
    cluster = setNames(colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(cluster_levels)), cluster_levels)
  )
  
  # Heatmap with color-blind friendly colors
  heatmap_height <- max(5, length(genes_found) * 0.12)
  pdf(file.path(outdir, paste0(output_prefix, "_heatmap_pseudobulk_logCPM.pdf")),
      width = 24,
      height = heatmap_height)
  print(pheatmap(heatmap_mat,
                 annotation_col = annotation_df[, c("group_id", "cluster")],
                 annotation_colors = ann_colors,
                 scale = "row", 
                 cluster_rows = TRUE,
                 cluster_cols = FALSE,
                 show_colnames = FALSE,
                 fontsize_row = 6,
                 clustering_distance_rows = "euclidean", 
                 #clustering_distance_cols = "euclidean",
                 clustering_method = "complete",
                 color = colorRampPalette(c("#2c7bb6", "#ffffbf", "#d7191c"))(40),
                 angle_col = "45",
                 main = paste0("LogCPM Heatmap of ", str_replace_all(output_prefix, "_", " ")))
  )
  dev.off()
  
  
  # ---- DotPlot
  dot_height <- ifelse(length(unique(data$seurat_clusters)) < 10, 7,
                       7 + ceiling(length(unique(data$seurat_clusters)) / 10))
  dot_width <- max(8, length(genes_found) * 0.3)
  pdf(file.path(outdir, paste0(output_prefix, "_dotplot.pdf")),
      width = dot_width, height = dot_height)
  print(DotPlot(data, features = genes_found) + RotatedAxis())
  dev.off()
  
  
  # ---- Violin plots
  n_features <- length(genes_found)
  n_per_plot <- 8
  n_panels <- ceiling(n_features / n_per_plot)
  vln_height <- max(6, n_panels * 3)
  
  pdf(file.path(outdir, paste0(output_prefix, "_violinplot.pdf")), width = 25, height = vln_height)
  for (j in seq(1, n_features, by = n_per_plot)) {
    print(VlnPlot(data,
                  features = genes_found[j:min(j + n_per_plot - 1, n_features)],
                  raster = FALSE, ncol = 4))
  }
  dev.off()
}


# Additional plots for reactome serine biosynthesis pathway genes --------------
output_prefix <- "mouse_reactome_serine_biosynthesis_genes"
outdir <- "../exploratory_analyses/evaluate_mouse_reactome_serine_biosynthesis_genes"

# ---- DotPlot
dot_width <- ifelse(length(unique(data$seurat_clusters)) < 10, 7,
                     16 + ceiling(length(unique(data$seurat_clusters)) / 10))
dot_height <- max(8, length(genelist$mouse_reactome_serine_biosynthesis_genes) * 0.5)
pdf(file.path(outdir, paste0(output_prefix, "_splitby_condition_cluster_dotplot.pdf")),
    width = dot_width, height = dot_height)
print(DotPlot(data, 
              features = genelist$mouse_reactome_serine_biosynthesis_genes, 
              split.by = "Condition",
              cols = "RdBu"
              ) + 
        RotatedAxis() +
        coord_flip())
dev.off()


# ---- Violin plots
pdf(file = file.path(outdir, paste0(output_prefix,"_splitby_condition_cluster_violinplot.pdf")),
    width = 25,
    height = 7)
for(i in 1:length(genelist$mouse_reactome_serine_biosynthesis_genes)){
  print(VlnPlot(data, 
                features= genelist$mouse_reactome_serine_biosynthesis_genes[i], 
                raster = FALSE,
                split.by = "Condition"))
}
dev.off()



print("***** Script Complete! *****")

# END