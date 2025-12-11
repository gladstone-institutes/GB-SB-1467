#!/usr/bin/env Rscript

###############################################################################
## Project ID: GB-SB-1467
## Authors: Ayushi Agrawal
## Script Goal: Generate supplementary figures for paper. Increase font sizes 
##              for all the figures (gene names, axis numbers and and labels) 
##              and make them all square shaped
## Note: This script was run locally on Ayushi's laptop
###############################################################################

# Setup ------------------------------------------------------------------------
# Load libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(clusterProfiler)
library(EnhancedVolcano)

# set the working directory
setwd("~/Dropbox (Gladstone)/GB-SB-1467/paper/revision_Dec2025/analysis_suppFig_v2/")
basedir <- file.path("/Volumes/Jain-Boinformatics-Collaboration",
                     "sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/results")


# Functions --------------------------------------------------------------------
# Volcano plots
plot_volcano <- function(de,
                         cell_type,
                         comparison,
                         height = 10,
                         width = 15,
                         p_value_cutoff = 0.05,
                         fold_change_cutoff = 0.25,
                         output) {
  number_of_de_genes <- de %>%
    filter(p_adj.loc < p_value_cutoff & abs(logFC) > fold_change_cutoff) %>%
    nrow()
  
  filename <- paste0(cell_type, "_", comparison)
  subtitle_expr <- paste0(comparison, 
                          "\np_adj.loc cutoff = ", p_value_cutoff, 
                          "; log2 fold change cutoff = ", fold_change_cutoff,
                          "\nNumber of DE genes = ", number_of_de_genes)
  
  p_value_cutoff_new <- min(de$p_val[de$p_adj.loc >= p_value_cutoff])
  
  #set the limits for x-axis
  x_pad  <- 0.2
  x_min  <- round(min(de$logFC, na.rm = TRUE),1) - x_pad
  x_max  <- round(max(de$logFC, na.rm = TRUE),1) + x_pad
  
  plot <- EnhancedVolcano(de,
                          title = cell_type,
                          subtitle = subtitle_expr,
                          lab = de$gene,
                          pCutoff = p_value_cutoff_new,
                          FCcutoff = fold_change_cutoff,
                          x = 'logFC',
                          y = 'p_val',
                          xlim = c(x_min, x_max),
                          ylab = bquote(~-Log[10]~italic(p_val)),
                          drawConnectors = TRUE,
                          widthConnectors = 0.3,
                          boxedLabels = TRUE,
                          labSize = 3.5,
                          colConnectors = "black",
                          legendLabSize = 10,
                          legendIconSize = 4.0,
                          legendLabels = c("NS", 
                                           bquote("Log"[2]~"FC > " ~ .(fold_change_cutoff)),
                                           bquote("p_adj.loc < " ~ .(p_value_cutoff)),
                                           bquote("p_adj.loc < " ~ .(p_value_cutoff) ~ " & Log"[2]~"FC > " ~ .(fold_change_cutoff)))
  )
  ggsave(
    plot = plot,
    filename = file.path(output, paste0(filename, "_volcano_plot.pdf")),
    width = width,
    height = height,
    units = "in"
  )
}


# Volcano plots ----------------------------------------------------------------
deDir <- "07_pseudobulk_de_v2/gene_expression_associations/KO-VitaminB3_vs_WT-VitaminB3"

# Define cluster numbers and corresponding figure labels
volcano_configs <- list(
  list(cluster = 2, fig_label = "A"),
  list(cluster = 8, fig_label = "B"),
  list(cluster = 12, fig_label = "C")
)

# Loop through clusters
for (config in volcano_configs) {
  cluster_num <- config$cluster
  
  dat <- read.csv(file.path(basedir, 
                            deDir,
                            paste0("cluster_", cluster_num),
                            paste0("cluster_", cluster_num, "_KO-VitaminB3_vs_WT-VitaminB3.csv")))
  
  # Generate the plot
  plot_volcano(dat,
               cell_type = paste("Cluster", cluster_num),
               comparison = "KO-VitaminB3_vs_WT-VitaminB3",
               height = 8,
               width = 8,
               p_value_cutoff = 0.01,
               fold_change_cutoff = 1,
               output = getwd()
  )
  
  # Rename the output file to use the desired naming convention
  old_filename <- file.path(getwd(), paste0("Cluster ", cluster_num, "_KO-VitaminB3_vs_WT-VitaminB3_volcano_plot.pdf"))
  new_filename <- file.path(getwd(), paste0("analysis_suppFig_", config$fig_label, ".pdf"))
  file.rename(old_filename, new_filename)
}
rm(dat)


# Dot plot for pathway enrichment ----------------------------------------------
gseaDir <- "07_pseudobulk_de_v2/GSEA/KO-VitaminB3_vs_WT-VitaminB3"

# Define cluster numbers and corresponding figure labels
cluster_configs <- list(
  list(cluster = 2, fig_label = "D", width = 7),
  list(cluster = 8, fig_label = "E", width = 7),
  list(cluster = 12, fig_label = "F", width = 7.6)
)

# Loop through clusters
for (config in cluster_configs) {
  cluster_num <- config$cluster
  
  df <- read.csv(file.path(basedir, gseaDir,
                           paste0("cluster_", cluster_num),
                           "wp_mm_20250415",
                           paste0("cluster_", cluster_num, "_KO-VitaminB3_vs_WT-VitaminB3_GSEA_WikiPathways_results_table.csv")), 
                 stringsAsFactors = FALSE, row.names = 1)
  er <- new("gseaResult")
  er@result <- df
  er@result <- dplyr::mutate(er@result, Description = str_trunc(Description, 80))
  
  pdf(paste0("analysis_suppFig_", config$fig_label, ".pdf"), 
      height = 7, 
      width = config$width)
  print( enrichplot::dotplot(er, showCategory = 20, label_format=50) + 
           theme_bw(base_size = 14) +
           labs(title = paste0("Cluster ", cluster_num, "\nKO-B3_vs_WT-B3"), 
                subtitle = "\nGSEA Terms for WikiPathways"))
  dev.off()
}

rm(df, er)


# Correlation plots ------------------------------------------------------------
# set working directory
corDir <- file.path(basedir, "07_pseudobulk_de_v2/gene_expression_associations")

# annotations for clusters
celltype_abbr <- read.csv(file.path(basedir, "/06_clustering_cell_type/30PC_0.08res/celltype_abbreviation_legend.csv"))

# Define cluster-to-figure-label mapping
cluster_fig_labels <- list(
  "2" = "G",
  "8" = "H",
  "12" = "I"
)

# All CSV files list
deg_files <- file.path(corDir, 
                       c("KO-VitaminB3_vs_WT-VitaminB3/cluster_2/cluster_2_KO-VitaminB3_vs_WT-VitaminB3.csv",
                         "KO+VitaminB3_vs_KO-VitaminB3/cluster_2/cluster_2_KO+VitaminB3_vs_KO-VitaminB3.csv",
                         "KO-VitaminB3_vs_WT-VitaminB3/cluster_8/cluster_8_KO-VitaminB3_vs_WT-VitaminB3.csv",
                         "KO+VitaminB3_vs_KO-VitaminB3/cluster_8/cluster_8_KO+VitaminB3_vs_KO-VitaminB3.csv",
                         "KO-VitaminB3_vs_WT-VitaminB3/cluster_12/cluster_12_KO-VitaminB3_vs_WT-VitaminB3.csv",
                         "KO+VitaminB3_vs_KO-VitaminB3/cluster_12/cluster_12_KO+VitaminB3_vs_KO-VitaminB3.csv"))

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

# For annotation: store a mapping from long names to short labels
comparison_pairs <- combn(comparisons, 2, simplify = FALSE)

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
      
      # Define threshold for labeling (e.g., difference > 1)
      label_threshold <- 2
      
      # Add label and label status (with explicit, interpretable values)
      merged$label <- ifelse(abs(merged$logFC_C1 - merged$logFC_C2) > label_threshold,
                             merged$Gene, NA)
      merged$Label_Status <- factor(
        ifelse(!is.na(merged$label), "|Δlog2FC| > 2", "|Δlog2FC| ≤ 2"),
        levels = c("|Δlog2FC| > 2", "|Δlog2FC| ≤ 2")
      )
      
      # Set axis labels based on comparison type
      # Since we only have 2 comparisons, there's only 1 pair to compare
      x_label <- if (c1 == "KO-VitaminB3_vs_WT-VitaminB3") {
        paste("Knockout effect\n(log2FC", c1, ")")
      } else {
        paste("Rescue effect\n(log2FC", c1, ")")
      }
      
      y_label <- if (c2 == "KO+VitaminB3_vs_KO-VitaminB3") {
        paste("Rescue effect\n(log2FC", c2, ")")
      } else {
        paste("Knockout effect\n(log2FC", c2, ")")
      }
      
      # Create the plot
      p <- ggplot(merged, aes(x = logFC_C1, y = logFC_C2)) +
        geom_smooth(method = "lm", se = FALSE, color = "#DC3220") +
        geom_point(aes(color = Label_Status), alpha = 0.6) +
        {
          if (any(!is.na(merged$label))) {
            geom_text_repel(aes(label = label), size = 3, max.overlaps = 20,
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
        theme_bw(base_size = 14) +
        labs(
          title = paste("Cluster", clus),
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
        ) +
        theme(
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.key.size = unit(0.4, "cm")
        )
      
      # Extract cluster number from the cluster abbreviation
      cluster_num <- celltype_abbr$Cluster[celltype_abbr$Abbreviation == clus]
      fig_label <- cluster_fig_labels[[as.character(cluster_num)]]
      
      # Save with the appropriate figure label
      ggsave(filename = paste0("analysis_suppFig_", fig_label, ".pdf"), 
             plot = p, width = 6.5, height = 5.5)
    }
  }
}




## END