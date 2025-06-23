#!/usr/bin/env Rscript

###############################################################################
## Project ID: GB-SB-1467
## Authors: Ayushi Agrawal
## Script Goal: Summarize the number of DEGs from pseudobulk DE analysis
## Note: This script was run locally on Ayushi's laptop
###############################################################################


# Load required libraries
library(dplyr)
library(ggplot2)
library(purrr)
library(tidytext)
library(tidyr)
library(UpSetR)
library(tibble)
library(readr)
library(forcats)

# set working directory
setwd(file.path("/Volumes/Jain-Boinformatics-Collaboration",
                "sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025",
                "results/07_pseudobulk_de_v2"))

# annotations for clusters
celltype_abbr <- read.csv("../06_clustering_cell_type/30PC_0.08res/celltype_abbreviation_legend.csv")

# comparisons to summarize
comparisons_of_interest <- c("KO-VitaminB3_vs_WT-VitaminB3",
                             "KO+VitaminB3_vs_WT-VitaminB3")

# Define path to DE result files
deg_dir <- "gene_expression_associations"

# Find all CSV files recursively
deg_files <- list.files(deg_dir, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)

# Filter function
count_deg <- function(file) {
  df <- read.csv(file)
  if (!all(c("p_adj.loc", "logFC", "gene") %in% colnames(df))) return(NULL)
  
  sig_degs <- df %>%
    filter(p_adj.loc < 0.01, abs(logFC) > 1)
  
  parts <- strsplit(file, "/")[[1]]
  comp <- parts[length(parts) - 2]
  
  cluster <- celltype_abbr$Abbreviation[celltype_abbr$Cluster == unique(df$cluster_id)]
  
  data.frame(
    Comparison = comp,
    Cluster = cluster,
    DEG_Count = nrow(sig_degs),
    Genes = I(list(sig_degs$gene))  # Store as list-column
  )
}

# Apply to all files
deg_summary <- map_dfr(deg_files, count_deg)

# Filter to keep only comparisons of interest
deg_summary <- deg_summary %>% filter(Comparison %in% comparisons_of_interest)

# Order Cluster factor numerically
extract_numeric_cluster <- function(cluster_label) {
  as.numeric(gsub("[^0-9]", "", cluster_label))
}
deg_summary <- deg_summary %>%
  mutate(Cluster_num = extract_numeric_cluster(Cluster)) %>%
  arrange(Cluster_num) %>%
  mutate(Cluster = factor(Cluster, levels = unique(Cluster))) %>%
  select(-Cluster_num)

# Save summary if needed
deg_summary_flat <- deg_summary %>%
  mutate(Gene_List = sapply(Genes, function(x) paste(x, collapse = ";"))) %>%
  select(-Genes)
write.csv(deg_summary_flat, "DEG_counts_by_comparison_cluster_summary.csv", row.names = FALSE)

# Plot barplot with text labels only for DEG_Count > 0
ggplot(deg_summary, aes(x = Cluster, y = DEG_Count, fill = Comparison)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = DEG_Count),
            position = position_dodge(width = 0.9),
            vjust = -0.3, size = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Number of DEGs per Cluster and Comparison",
       y = "# DEGs (padj < 0.01 & abs(log2FC) > 1)", x = "Cluster")

# save plot
ggsave("DEG_counts_barplot.pdf", width = 14, height = 6)



# Choose the reference comparison to order by
reference_comparison <- "KO-VitaminB3_vs_WT-VitaminB3"

# Get the cluster ordering based on the reference comparison
cluster_order <- deg_summary %>%
  filter(Comparison == reference_comparison) %>%
  arrange(DEG_Count) %>%
  pull(Cluster)

# Convert Cluster to a factor with the desired order
deg_summary$Cluster <- factor(deg_summary$Cluster, levels = cluster_order)


# Plot: horizontal, faceted, ordered by DEG count
ggplot(deg_summary, aes(x = DEG_Count, y = Cluster, fill = Comparison)) +
  geom_col() +
  geom_text(aes(label = ifelse(DEG_Count > 0, DEG_Count, "")),
            hjust = -0.1, size = 3) +
  scale_y_reordered() +
  facet_wrap(~ Comparison) +
  theme_bw() +
  theme(strip.text = element_text(face = "bold"),
        axis.text.y = element_text(size = 10)) +
  labs(title = "DEG Counts per Cluster (Faceted by Comparison)",
       x = "# DEGs (padj < 0.01 & abs(log2FC) > 1)", y = "Cluster")
ggsave("DEG_counts_barplot_facet_by_comparison.pdf", width = 15, height = 5)



# upset plot
# Create a binary matrix of gene presence across each comparison-cluster
deg_long <- deg_summary %>%
  unnest(Genes) %>%
  mutate(set_id = paste(Comparison, Cluster, sep = "_"))

# Reformat for UpSetR: binary presence matrix
deg_binary_matrix <- deg_long %>%
  distinct(Genes, set_id) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = set_id, values_from = value, values_fill = 0) %>%
  column_to_rownames("Genes")

# Plot
png("DEG_upset_plot.png", width = 3400, height = 3000, res = 200)
(upset(deg_binary_matrix, 
       sets = colnames(deg_binary_matrix),  # Ensure all comparisons are included
       nsets = length(deg_binary_matrix),  # Show all sets
       order.by = "freq", sets.bar.color = "#56B4E9",
       mainbar.y.label = "Number of overlapping DE genes",
       sets.x.label = "DE comparisons", text.scale = 1.2))
dev.off()



# END