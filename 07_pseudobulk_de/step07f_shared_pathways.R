#!/usr/bin/env Rscript

###############################################################################
## Project ID: GB-SB-1467
## Authors: Ayushi Agrawal
## Script Goal: Summarize shared pathways from GSEA and ORA
## Note: This script was run locally on Ayushi's laptop
###############################################################################

# Load required libraries
library(purrr)
library(dplyr)
#library(ggplot2)
#library(viridis)
#library(conflicted)
#conflict_prefer("filter", "dplyr")
#conflict_prefer("select", "dplyr")

# Set working directory
setwd(file.path("/Volumes/Jain-Boinformatics-Collaboration",
                "sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025",
                "results/07_pseudobulk_de_v2"))
comparisons_of_interest <- c("KO-VitaminB3_vs_WT-VitaminB3")

# ================================
# Define helper function
# ================================

get_common_pathways <- function(comp, db, data, min_clusters = 2) {
  df_sub <- data %>%
    filter(Comparison == comp, Pathway_db == db)
  
  if (nrow(df_sub) == 0) {
    return(NULL)
  }
  
  method <- unique(df_sub$Method)
  if (length(method) > 1) method <- paste(method, collapse = ";")
  
  # Count how many clusters each pathway appears in
  pathway_cluster_counts <- df_sub %>%
    group_by(Description) %>%
    summarise(
      N_clusters = n_distinct(Cluster),
      Clusters = paste(sort(unique(Cluster)), collapse = ", "),
      .groups = "drop"
    ) %>%
    filter(N_clusters >= min_clusters) %>%
    mutate(
      Comparison = comp,
      Pathway_db = db,
      Method = method
    ) %>%
    select(Comparison, Pathway_db, Method, Description, N_clusters, Clusters)
  
  return(pathway_cluster_counts)
}


# ================================
# Step 1: Load GSEA results
# ================================
gsea_dir <- "GSEA"

gsea_files <- list.files(gsea_dir, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)

read_gsea_file <- function(file) {
  df <- read.csv(file)
  
  if (!all(c("p.adjust", "NES", "Description") %in% colnames(df))) return(NULL)
  
  parts <- strsplit(file, "/")[[1]]
  pathway_db <- parts[length(parts) - 1]
  cluster <- parts[length(parts) - 2]
  comp <- parts[length(parts) - 3]
  enrichment_method <- parts[length(parts) - 4]
  
  df %>%
    mutate(
      Comparison = comp,
      Cluster = cluster,
      Pathway_db = pathway_db,
      Method = enrichment_method
    )
}

gsea_data <- gsea_files %>%
  map(read_gsea_file) %>%
  compact() %>%
  bind_rows()

gsea_sig <- gsea_data %>% 
  filter(p.adjust < 0.05) %>%
  filter(Comparison %in% comparisons_of_interest)

gsea_combos <- gsea_sig %>%
  select(Comparison, Pathway_db) %>%
  distinct()

gsea_common_summary <- map2_dfr(
  gsea_combos$Comparison,
  gsea_combos$Pathway_db,
  ~ get_common_pathways(.x, .y, gsea_sig)
)

gsea_common_serine <- gsea_common_summary %>%
  filter(!is.na(Description) & grepl("serine", Description, ignore.case = TRUE))


# ================================
# Step 2: Load ORA results
# ================================
ora_dir <- "ORA"

ora_files <- list.files(ora_dir, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)

read_ora_file <- function(file) {
  df <- read.csv(file)
  
  if (!all(c("p.adjust", "FoldEnrichment", "Description") %in% colnames(df))) return(NULL)
  
  parts <- strsplit(file, "/")[[1]]
  pathway_db <- parts[length(parts) - 1]
  cluster <- parts[length(parts) - 2]
  comp <- parts[length(parts) - 3]
  enrichment_method <- parts[length(parts) - 4]
  
  df %>%
    mutate(
      Comparison = comp,
      Cluster = cluster,
      Pathway_db = pathway_db,
      Method = enrichment_method
    )
}

ora_data <- ora_files %>%
  map(read_ora_file) %>%
  compact() %>%
  bind_rows()

ora_sig <- ora_data %>%
  filter(p.adjust < 0.05) %>%
  filter(Comparison %in% comparisons_of_interest)

ora_combos <- ora_sig %>%
  select(Comparison, Pathway_db) %>%
  distinct()

ora_common_summary <- map2_dfr(
  ora_combos$Comparison,
  ora_combos$Pathway_db,
  ~ get_common_pathways(.x, .y, ora_sig)
)

ora_common_serine <- ora_common_summary %>%
  filter(!is.na(Description) & grepl("serine", Description, ignore.case = TRUE))


# ================================
# Save results
# ================================

# Save each file with a distinct name
if (nrow(gsea_common_summary) > 0) {
  write.csv(gsea_common_summary, 
            "common_pathways_gsea_summary_KO-VitaminB3_vs_WT-VitaminB3.csv", row.names = F)
}

if (nrow(gsea_common_serine) > 0) {
  write.csv(gsea_common_serine, 
            "common_pathways_gsea_serine_only_KO-VitaminB3_vs_WT-VitaminB3.csv", row.names = F)
}

if (nrow(ora_common_summary) > 0) {
  write.csv(ora_common_summary, 
            "common_pathways_ora_summary_KO-VitaminB3_vs_WT-VitaminB3.csv", row.names = F)
}

if (nrow(ora_common_serine) > 0) {
  write.csv(ora_common_serine, 
            "common_pathways_ora_serine_only_KO-VitaminB3_vs_WT-VitaminB3.csv", row.names = F)
}


# END