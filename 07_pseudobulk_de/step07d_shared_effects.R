#!/usr/bin/env Rscript

###############################################################################
## Project ID: GB-SB-1467
## Authors: Ayushi Agrawal
## Script Goal: Summarize pseudobulk DEGs and visualize shared effects
## Note: This script was run locally on Ayushi's laptop
###############################################################################

# Load required libraries
library(pheatmap)
library(purrr)
library(tidyverse)
library(ggplot2)
library(viridis)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

# Set working directory
setwd(file.path("/Volumes/Jain-Boinformatics-Collaboration",
                "sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025",
                "results/07_pseudobulk_de_v2"))

# annotations for clusters
celltype_abbr <- read.csv("../06_clustering_cell_type/30PC_0.08res/celltype_abbreviation_legend.csv")

# Define comparison of interest
comparison_of_interest <- "KO-VitaminB3_vs_WT-VitaminB3"

# Define path to DE result files
deg_dir <- "gene_expression_associations"

# Find all CSV files recursively
deg_files <- list.files(deg_dir, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)

# Function to read and annotate DEG files
read_deg_file <- function(file) {
  df <- read.csv(file)
  
  if (!all(c("p_adj.loc", "logFC", "gene") %in% colnames(df))) return(NULL)
  
  parts <- strsplit(file, "/")[[1]]
  comp <- parts[length(parts) - 2]
  
  cluster <- celltype_abbr$Abbreviation[celltype_abbr$Cluster == unique(df$cluster_id)]
  
  df %>%
    mutate(
      Comparison = comp,
      Cluster = cluster,
      Gene = gene
    )
}

# Load all DEG data
all_deg_data <- deg_files %>%
  map(read_deg_file) %>%
  compact() %>%
  bind_rows()

# Summarize significant DEGs per comparison × cluster
deg_summary <- all_deg_data %>%
  filter(p_adj.loc < 0.01, abs(logFC) > 1) %>%
  group_by(Comparison, Cluster) %>%
  summarise(
    DEG_Count = n(),
    Genes = I(list(Gene)),
    .groups = "drop"
  ) %>%
  filter(DEG_Count > 0)

# Order Cluster factor numerically for consistent plotting
extract_numeric_cluster <- function(cluster_label) {
  as.numeric(gsub("[^0-9]", "", cluster_label))
}
deg_summary <- deg_summary %>%
  mutate(Cluster_num = extract_numeric_cluster(Cluster)) %>%
  arrange(Cluster_num) %>%
  mutate(Cluster = factor(Cluster, levels = unique(Cluster))) %>%
  select(-Cluster_num)

# Flatten summary to one gene list per row
deg_summary_flat <- deg_summary %>%
  mutate(Gene_List = sapply(Genes, function(x) paste(x, collapse = ";"))) %>%
  select(-Genes)

# Filter for a specific comparison of interest
deg_subset <- deg_summary_flat %>%
  filter(Comparison == comparison_of_interest)

# Split gene list into long format
gene_cluster_df <- deg_subset %>%
  separate_rows(Gene_List, sep = ";") %>%
  rename(gene = Gene_List)

# Count number of clusters each gene appears in
gene_counts <- gene_cluster_df %>%
  distinct(gene, Cluster) %>%
  count(gene, name = "n_clusters")

# Bin genes by how many clusters they are shared across
summary_df <- gene_counts %>%
  count(n_clusters, name = "n_genes") %>%
  mutate(n_clusters_bin = ifelse(n_clusters > 10, ">10", as.character(n_clusters))) %>%
  group_by(n_clusters_bin) %>%
  summarise(n_genes = sum(n_genes)) %>%
  ungroup() %>%
  mutate(n_clusters_bin = factor(n_clusters_bin, 
                                 levels = c(as.character(1:10), ">10")))

# Plot 1: Bar plot of gene sharing across clusters
p1 <- ggplot(summary_df, aes(x = n_clusters_bin, y = n_genes)) +
  geom_bar(stat = "identity", fill = "gray70", width = 0.7) +
  geom_text(aes(label = n_genes), vjust = -0.5, size = 3) +
  theme_classic(base_size = 14) +
  labs(
    title = paste("DEG overlap across clusters"),
    subtitle = paste("Comparison:", comparison_of_interest),
    x = "Number of Cell Types Sharing Effect",
    y = "Number of Genes"
  ) +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold")
  )

# Save bar plot
ggsave(paste0("deg_overlap_across_clusters_", comparison_of_interest, ".pdf"), 
       p1, 
       width = 7, 
       height = 5)


# Identify commonly shared genes
common_genes <- gene_counts %>%
  filter(n_clusters > 2) %>%
  pull(gene)

# Prepare heatmap input: filter logFC values and cap extremes
heatmap_df <- all_deg_data %>%
  filter(gene %in% common_genes & Comparison == comparison_of_interest) %>%
  mutate(
    gene = factor(gene, levels = common_genes),
    sig_marker = ifelse(p_adj.loc < 0.01 & abs(logFC) > 1, "+", ""),
    logFC_capped = pmax(pmin(logFC, 2), -2),
    Cluster_num = as.numeric(sub("_.*", "", Cluster))
  ) %>%
  arrange(Cluster_num) %>%
  mutate(Cluster = factor(Cluster, levels = unique(Cluster)))  # set consistent order here

# Create heatmap matrix
heatmap_matrix <- heatmap_df %>%
  select(gene, Cluster, logFC_capped) %>%
  pivot_wider(names_from = gene, values_from = logFC_capped) %>%
  column_to_rownames("Cluster") %>%
  as.matrix()

# Create significance marker matrix
sig_marker_matrix <- heatmap_df %>%
  select(gene, Cluster, sig_marker) %>%
  pivot_wider(names_from = gene, values_from = sig_marker) %>%
  column_to_rownames("Cluster") %>%
  as.matrix()


# Define custom color scale: blue-white-red
custom_colors <- colorRampPalette(c("#053061", "#f7f7f7", "#B2182B"))(100)

# Plot heatmap 
# no column clustering
p2 <- pheatmap(
  mat = heatmap_matrix,
  cluster_rows = FALSE,        # keep clusters in specified order
  cluster_cols = TRUE,         # cluster genes
  color = custom_colors,
  display_numbers = sig_marker_matrix,
  number_color = "black",
  main = paste("Common DEGs across >2 clusters (colors = log2FC)\n",
               "+ indicates FDR < 0.01 & abs(log2FC) > 1",
               "\nComparison:", comparison_of_interest),
  legend = TRUE,
  fontsize_row = 6,
  fontsize_col = 4,
  fontsize_number = 3,
  border_color = "grey90",
  angle_col = 45,
  labels_col = colnames(heatmap_matrix),
  labels_row = rownames(heatmap_matrix)
)

# with column clustering
p3 <- pheatmap(
  mat = heatmap_matrix,
  cluster_rows = TRUE,        # keep clusters in specified order
  cluster_cols = TRUE,         # cluster genes
  color = custom_colors,
  display_numbers = sig_marker_matrix,
  number_color = "black",
  main = paste("Common DEGs across >2 clusters (colors = log2FC)\n",
               "+ indicates FDR < 0.01 & abs(log2FC) > 1",
               "\nComparison:", comparison_of_interest),
  legend = TRUE,
  fontsize_row = 6,
  fontsize_col = 4,
  fontsize_number = 3,
  border_color = "grey90",
  angle_col = 45,
  labels_col = colnames(heatmap_matrix),
  labels_row = rownames(heatmap_matrix)
)

# Save heatmap as PDF
pdf(paste0("heatmap_common_genes_", comparison_of_interest, ".pdf"), 
    width = 6, 
    height = 4)
print(p2)
dev.off()

pdf(paste0("heatmap_common_genes_", comparison_of_interest, "_cluster_cols.pdf"), 
    width = 6, 
    height = 4)
print(p3)
dev.off()

# Save data used for heatmap
write.csv(heatmap_df, paste0("heatmap_common_genes_", comparison_of_interest, ".csv"), 
          row.names = FALSE)




# Save session info for reproducibility
writeLines(capture.output(sessionInfo()), 
           "sessionInfo_deg_summary.txt")

# Print message
print("***** Script Complete *****")

# END