#!/usr/bin/env Rscript

###############################################################################
## Project ID: GB-SB-1467
## Authors: Ayushi Agrawal
## Script Goal: Custom volcano plots 
## Note: This script was run locally on Ayushi's laptop
###############################################################################


# Load required libraries
library(dplyr)
library(EnhancedVolcano)

# set working directory
setwd(file.path("/Volumes/Jain-Boinformatics-Collaboration",
                "sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025",
                "results/07_pseudobulk_de_v2",
                "gene_expression_associations/KO-VitaminB3_vs_WT-VitaminB3"))


# define variables
comparison <- "KO-VitaminB3_vs_WT-VitaminB3"
height <- 10
width <- 15
p_value_cutoff <- 0.01
fold_change_cutoff <- 1
output <- file.path("/Volumes/Jain-Boinformatics-Collaboration",
                    "sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025",
                    "results/paper/figures/subfig_I")
clusters_of_interst <- c("cluster_2",
                         "cluster_4",
                         "cluster_8",
                         "cluster_12")
genes_of_interest <- c("Mt1", "Mt2", "Srr","Serinc3","Serinc1","Serinc5",
                       "Serinc2","Psat1","Psph","Serinc4","Phgdh")


# create output folder
if(!dir.exists(output)) dir.create(output)


# generate volcano plots
for (i in clusters_of_interst){
  deg <- read.csv(paste0(i, "/", i, "_KO-VitaminB3_vs_WT-VitaminB3.csv"))
  deg$label <- ifelse(deg$gene %in% genes_of_interest, deg$gene, NA)
  
  number_of_de_genes <- deg %>%
    filter(p_adj.loc < p_value_cutoff & abs(logFC) > fold_change_cutoff) %>%
    nrow()
  
  filename <- paste0(i, "_", comparison)
  subtitle_expr <- paste0(comparison, 
                          "\np_adj.loc cutoff = ", p_value_cutoff, 
                          "; log2 fold change cutoff = ", fold_change_cutoff,
                          "\nNumber of DE genes = ", number_of_de_genes)
  
  p_value_cutoff_new <- min(deg$p_val[deg$p_adj.loc >= p_value_cutoff])
  
  p1 <- EnhancedVolcano(deg,
                        title = i,
                        subtitle = subtitle_expr,
                        lab = deg$gene,
                        selectLab = genes_of_interest,
                        pCutoff = p_value_cutoff_new,
                        FCcutoff = fold_change_cutoff,
                        x = 'logFC',
                        y = 'p_val',
                        ylab = bquote(~-Log[10]~italic(p_val)),
                        drawConnectors = TRUE,
                        widthConnectors = 0.3,
                        boxedLabels = TRUE,
                        labSize = 3.0,
                        colConnectors = "black",
                        legendLabSize = 10,
                        legendIconSize = 3.0,
                        legendLabels = c("NS", 
                                         bquote("Log"[2]~"FC > " ~ .(fold_change_cutoff)),
                                         bquote("p_adj.loc < " ~ .(p_value_cutoff)),
                                         bquote("p_adj.loc < " ~ .(p_value_cutoff) ~ " & Log"[2]~"FC > " ~ .(fold_change_cutoff)))
  )
  ggsave(
    plot = p1,
    filename = file.path(output, paste0(filename, "_volcano_plot.pdf")),
    width = width,
    height = height,
    units = "in"
  )
}


# END