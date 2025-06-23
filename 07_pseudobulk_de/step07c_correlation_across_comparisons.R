#!/usr/bin/env Rscript

###############################################################################
## Project ID: GB-SB-1467
## Authors: Ayushi Agrawal
## Script Goal: Correlation analyses for the pseudobulk DE results
## Note: This script was run locally on Ayushi's laptop
###############################################################################


# Load required libraries
library(tidyr)
library(stringr)
library(dplyr)
library(ggplot2)
library(purrr)
library(reshape2)
library(pheatmap)
library(ggrepel)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")


# set working directory
setwd(file.path("/Volumes/Jain-Boinformatics-Collaboration",
                "sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025",
                "results/07_pseudobulk_de_v2"))

# create output directory
if (!(dir.exists("correlation_bw_comparisons"))) {
  dir.create("correlation_bw_comparisons", recursive = T, showWarnings = F)
}

# annotations for clusters
celltype_abbr <- read.csv("../06_clustering_cell_type/30PC_0.08res/celltype_abbreviation_legend.csv")

# Define path to DE result files
deg_dir <- "gene_expression_associations"

# Find all CSV files recursively
deg_files <- list.files(deg_dir, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)

# Function to read and attach cluster and comparison info
read_deg_file <- function(file) {
  df <- read.csv(file)
  if (!all(c("p_val", "logFC") %in% colnames(df))) return(NULL)
  
  parts <- strsplit(file, "/")[[1]]
  comp <- parts[length(parts) - 2]
  
  cluster <- celltype_abbr$Abbreviation[celltype_abbr$Cluster == unique(df$cluster_id)]
  
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
  print(paste("Processing ", clus, "......"))
  df_clus <- filtered_deg %>% filter(Cluster == clus)
  cor_row <- c()
  
  for (pair in comparison_pairs) {
    print(paste("compairson pair:", pair))
    c1 <- pair[1]
    c2 <- pair[2]
    
    df_c1 <- df_clus %>% dplyr::filter(Comparison == c1) %>% dplyr::select(Gene, logFC)
    df_c2 <- df_clus %>% dplyr::filter(Comparison == c2) %>% dplyr::select(Gene, logFC)
    
    merged <- inner_join(df_c1, df_c2, by = "Gene", suffix = c("_C1", "_C2"))
    print(paste("dimensions of merged data frame:", dim(merged)))
    
    if (nrow(merged) > 2) {
      cor_val <- cor(merged$logFC_C1, merged$logFC_C2, method = "pearson")
      comp_label <- paste(c1, "vs", c2)
      cor_row[short_names[comp_label]] <- cor_val
      
      # Define threshold for labeling (e.g., difference > 1)
      label_threshold <- 2
      
      # Add label and label status (with explicit, interpretable values)
      merged$label <- ifelse(abs(merged$logFC_C1 - merged$logFC_C2) > label_threshold,
                             merged$Gene, NA)
      merged$Label_Status <- factor(
        ifelse(!is.na(merged$label), "|Δlog2FC| > 2", "|Δlog2FC| ≤ 2"),
        levels = c("|Δlog2FC| > 2", "|Δlog2FC| ≤ 2")
      )
      
      # Custom axis label logic for the first comparison pair
      if (setequal(pair, c("KO-VitaminB3_vs_WT-VitaminB3", "KO+VitaminB3_vs_KO-VitaminB3"))) {
        x_label <- paste("Knockout effect\n(log2FC", c1, ")")
        y_label <- paste("Rescue effect\n(log2FC", c2, ")")
      } else {
        x_label <- paste("log2FC", c1)
        y_label <- paste("log2FC", c2)
      }
      
      # Create the plot
      p <- ggplot(merged, aes(x = logFC_C1, y = logFC_C2)) +
        geom_smooth(method = "lm", se = FALSE, color = "#DC3220") +
        geom_point(aes(color = Label_Status), alpha = 0.6) +
        {
          if (any(!is.na(merged$label))) {
            geom_text_repel(aes(label = label), size = 2, max.overlaps = 20,
                            box.padding = 0.5,
                            point.padding = 0.3,
                            min.segment.length = 0,
                            force = 2,
                            force_pull = 0.1,
                            nudge_y = 0.1,
                            nudge_x = 0.1,
                            segment.color = "grey50",
                            segment.alpha = 0.5)
          }
        } +
        theme_bw() +
        labs(
          title = paste("Cluster", clus, "-", c1, "vs", c2),
          x = x_label,
          y = y_label,
          subtitle = paste("Pearson r =", round(cor_val, 3)),
          caption = "Red line: linear regression fit"
        ) +
        scale_color_manual(
          values = c("|Δlog2FC| > 2" = "#0C7BDC", "|Δlog2FC| ≤ 2" = "black"),
          labels = c(
            "|Δlog2FC| > 2" = expression("|"~Delta~log[2]*"FC| > 2"),
            "|Δlog2FC| ≤ 2" = expression("|"~Delta~log[2]*"FC| ≤ 2")
          ),
          name = expression("|"~Delta~log[2]*"FC| threshold")
        )
      
      print(p)
      ggsave(filename = paste0("correlation_bw_comparisons/scatter_cluster_", 
                               clus, "_", c1, "_vs_", c2, ".pdf"), 
             plot = p, width = 8, height = 6)
    }
  }
  
  cor_matrix_list[[clus]] <- cor_row
}

# Combine all rows into a matrix and plot heatmap
cor_matrix <- bind_rows(cor_matrix_list) %>% as.data.frame()
rownames(cor_matrix) <- names(cor_matrix_list)

# Add annotation to explain short column labels
annotation_col <- data.frame(Comparison = paste(names(short_names),
                                                paste0("(",short_names,")")))
rownames(annotation_col) <- short_names

pdf("correlation_bw_comparisons/correlation_heatmap_all.pdf", width = 15, height = 6)
(pheatmap(cor_matrix,
          cluster_rows = TRUE, cluster_cols = TRUE,
          main = "Pearson Correlation of log2FC between Comparisons by Cluster",
          display_numbers = TRUE,
          number_color = "black",
          annotation_col = annotation_col,
          legend = TRUE,
          angle_col = "0")
)
dev.off()


# Rename the columns of cor_matrix
long_names <- setNames(names(short_names), short_names)
colnames(cor_matrix) <- long_names[colnames(cor_matrix)]


# plot heatmap for comparison pairs of interest
comparisons_of_interest <- c( "KO-VitaminB3_vs_WT-VitaminB3",
                              "KO+VitaminB3_vs_WT-VitaminB3",
                              "KO+VitaminB3_vs_KO-VitaminB3"
                              )
col_names <- colnames(cor_matrix)

# Find those with at least 2 of the targets
matched_cols <- col_names[sapply(col_names, 
                                 function(x) sum(str_detect(x, fixed(comparisons_of_interest))) >= 2)]

# View result
cor_matrix_subset <- cor_matrix[,matched_cols]
cor_matrix_subset <- cor_matrix_subset[rowSums(is.na(cor_matrix_subset)) == 0, ]

# Add annotation to explain short column labels
annotation_col <- data.frame(Comparison = matched_cols)
rownames(annotation_col) <- matched_cols

pdf("correlation_bw_comparisons/correlation_heatmap_comparisons_of_interest.pdf", 
    width = 12, height = 6)
(pheatmap(cor_matrix_subset,
          cluster_rows = TRUE, 
          cluster_cols = TRUE,
          labels_col = rep("",3),
          main = "Pearson Correlation of log2FC between Comparisons by Cluster",
          display_numbers = TRUE,
          number_color = "black",
          annotation_col = annotation_col,
          legend = TRUE,
          angle_col = "0")
)
dev.off()


# clean condition labels
colnames(cor_matrix) <- gsub("-", "wo", gsub("\\+", "w", gsub(" ", "_", colnames(cor_matrix))))

# output the matrix
write.csv(cor_matrix, 
          file = "correlation_bw_comparisons/correlation_heatmap_data.csv",
          row.names = TRUE)





# Write session info
writeLines(
  capture.output(sessionInfo()),
  "correlation_bw_comparisons/sessionInfo.txt"
)
print("***** Script Complete *****")


# END