#!/usr/bin/env Rscript

###############################################################################
## Project ID: GB-SB-1467
## Authors: Ayushi Agrawal
## Script Goal: Correlation analyses for the pseudobulk DE results
## Note: This script was run locally on Ayushi's laptop
###############################################################################


# Load required libraries
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(reshape2)
library(pheatmap)

# set working directory
setwd(file.path("/Volumes/Jain-Boinformatics-Collaboration",
                "sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025",
                "results/07_pseudobulk_de"))

# create output directory
if (!(dir.exists("correlation_bw_comparisons"))) {
  dir.create("correlation_bw_comparisons", recursive = T, showWarnings = F)
}

# Define path to DE result files
deg_dir <- "gene_expression_associations"

# Find all CSV files recursively
deg_files <- list.files(deg_dir, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)

# Function to read and attach cluster and comparison info
read_deg_file <- function(file) {
  df <- read.csv(file)
  if (!all(c("p_val", "logFC") %in% colnames(df))) return(NULL)
  parts <- strsplit(file, "/")[[1]]
  cluster <- parts[length(parts) - 1]
  comp <- parts[length(parts) - 2]
  df$Comparison <- comp
  df$Cluster <- cluster
  df$Gene <- df$gene  # ensure gene column exists for merging
  df
}

# Load all DEGs
all_deg_data <- map_dfr(deg_files, read_deg_file)

# Filter for raw p < 0.05
filtered_deg <- all_deg_data %>% filter(p_val < 0.05)

# Get list of clusters and comparisons
clusters <- unique(filtered_deg$Cluster)
comparisons <- unique(filtered_deg$Comparison)

# Prepare empty list to collect correlation values per cluster
cor_matrix_list <- list()

# For annotation: store a mapping from long names to short labels
comparison_pairs <- combn(comparisons, 2, simplify = FALSE)
short_names <- paste0("P", seq_along(comparison_pairs))
pair_labels <- sapply(comparison_pairs, function(x) paste(x, collapse = " vs "))
names(short_names) <- pair_labels

# Iterate over each cluster
for (clus in clusters) {
  df_clus <- filtered_deg %>% filter(Cluster == clus)
  cor_row <- c()
  
  for (pair in comparison_pairs) {
    c1 <- pair[1]
    c2 <- pair[2]
    
    df_c1 <- df_clus %>% filter(Comparison == c1) %>% select(Gene, logFC)
    df_c2 <- df_clus %>% filter(Comparison == c2) %>% select(Gene, logFC)
    
    merged <- inner_join(df_c1, df_c2, by = "Gene", suffix = c("_C1", "_C2"))
    
    if (nrow(merged) > 2) {
      cor_val <- cor(merged$logFC_C1, merged$logFC_C2, method = "pearson")
      comp_label <- paste(c1, "vs", c2)
      cor_row[short_names[comp_label]] <- cor_val
      
      # Save scatterplot
      p <- ggplot(merged, aes(x = logFC_C1, y = logFC_C2)) +
        geom_point(alpha = 0.6) +
        geom_smooth(method = "lm", se = FALSE, color = "red") +
        theme_bw() +
        labs(title = paste("Cluster", clus, "-", c1, "vs", c2),
             x = paste("log2FC", c1),
             y = paste("log2FC", c2),
             subtitle = paste("Pearson r =", round(cor_val, 3)))
      
      print(p)
      ggsave(filename = paste0("correlation_bw_comparisons/scatter_cluster_", 
                               clus, "_", c1, "_vs_", c2, ".pdf"), plot = p, width = 6, height = 5)
    }
  }
  
  cor_matrix_list[[clus]] <- cor_row
}

# Combine all rows into a matrix and plot heatmap
cor_matrix <- bind_rows(cor_matrix_list) %>% as.data.frame()
rownames(cor_matrix) <- names(cor_matrix_list)

# Add annotation to explain short column labels
annotation_col <- data.frame(Comparison = names(short_names))
rownames(annotation_col) <- short_names

pdf("correlation_bw_comparisons/correlation_heatmap.pdf", width = 15, height = 6)
(pheatmap(cor_matrix,
          cluster_rows = TRUE, cluster_cols = TRUE,
          main = "Pearson Correlation of log2FC between Comparisons by Cluster",
          display_numbers = TRUE,
          number_color = "black",
          annotation_col = annotation_col,
          legend = TRUE,
          legend_position = "bottom")
)
dev.off()




# Write session info
writeLines(
  capture.output(sessionInfo()),
  "correlation_bw_comparisons/sessionInfo.txt"
)
print("***** Script Complete *****")


# END