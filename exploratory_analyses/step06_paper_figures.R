#!/usr/bin/env Rscript

###############################################################################
## Project ID: GB-SB-1467
## Authors: Ayushi Agrawal
## Script Goal: Dotplots for ANOVA and Tukey HSD results
## Note: This script was run locally on Ayushi's laptop
###############################################################################

# Load libraries
library(ggplot2)
library(dplyr)

# set working directory
setwd(file.path("/Volumes/Jain-Boinformatics-Collaboration",
                "sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/results"))

# read in the data
celltype_abbr <- read.csv("06_clustering_cell_type/30PC_0.08res/celltype_abbreviation_legend.csv")

tukey_df <- read.csv(file.path("exploratory_analyses/evaluate_pathways",
                               "gb_sb_1467_30PC_0.08res_tukey_results_module_scores.csv"))

anova_df <- read.csv(file.path("exploratory_analyses/evaluate_pathways",
                               "gb_sb_1467_30PC_0.08res_anova_summary_module_scores.csv"))


# Summary dotplot for significant Tukey results
# Define significant module × cluster pairs from ANOVA
signif_anova <- anova_df %>%
  filter(p_adj < 0.05) %>%
  mutate(module_cluster = paste(module, cluster, sep = "_"))

# Filter Tukey results based on significant ANOVA hits
sig_tukey <- tukey_df %>%
  mutate(
    comparison = paste(group2, "vs", group1),
    module = sub("\\.v2024\\.1\\.Mm[0-9]+$", "", module)
  ) %>%
  filter(comparison %in% c("KO_wo_Vitamin_B3 vs WT_wo_Vitamin_B3",
                           "KO_w_Vitamin_B3 vs KO_wo_Vitamin_B3")) %>%
  left_join(celltype_abbr, by = c("cluster" = "Cluster")) %>%
  mutate(
    cell_type = Abbreviation,
    module_cluster = paste(module, cluster, sep = "_"),
    signif = ifelse(p.adj < 0.05 & module_cluster %in% signif_anova$module_cluster, "1", "0")
  )

sig_tukey$module_cluster <- factor(sig_tukey$module_cluster,
                                   levels = sig_tukey %>%
                                     arrange(module, cluster) %>%
                                     pull(module_cluster) %>% unique())

# Order Cluster factor numerically for consistent plotting
extract_numeric_cluster <- function(cluster_label) {
  as.numeric(gsub("[^0-9]", "", cluster_label))
}
sig_tukey <- sig_tukey %>%
  mutate(Cluster_num = extract_numeric_cluster(cell_type)) %>%
  arrange(Cluster_num) %>%
  mutate(cell_type = factor(cell_type, levels = unique(cell_type))) %>%
  select(-Cluster_num)

p1 <- ggplot(sig_tukey %>% filter(module == "REACTOME_SERINE_BIOSYNTHESIS"), 
             aes(x = comparison, y = cell_type)) +
  geom_point(aes(color = cohen_d, size = -log10(p.adj), shape = signif)) +
  scale_color_gradient2(
    low = "#053061",   # darker blue
    mid = "#f7f7f7",   # light gray
    high = "#B2182B",  # strong red-orange
    midpoint = 0,
    name = "Cohen's d"
  ) +
  scale_shape_manual(
    values = c("0" = 16, "1" = 17),
    name = "Significant (FDR < 0.05)",
    labels = c("No", "Yes")
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 8)) +
  labs(
    title = "Condition Comparisons (Tukey HSD)", 
    subtitle = "Module: REACTOME_SERINE_BIOSYNTHESIS",
    x = "Condition Comparison", 
    y = "Module × Cluster"
  )

ggsave(file.path("exploratory_analyses/evaluate_pathways",
                 "gb_sb_1467_adj_pval_significant_tukey_comparisons_serine_only_dotplot.pdf"), 
       plot = p1, width = 6, height = 7)


p2 <- ggplot(sig_tukey %>% filter(module == "HALLMARK_APOPTOSIS"), 
             aes(x = comparison, y = cell_type)) +
  geom_point(aes(color = cohen_d, size = -log10(p.adj), shape = signif)) +
  scale_color_gradient2(
    low = "#053061",   # darker blue
    mid = "#f7f7f7",   # light gray
    high = "#B2182B",  # strong red-orange
    midpoint = 0,
    name = "Cohen's d"
  ) +
  scale_shape_manual(
    values = c("0" = 16, "1" = 17),
    name = "Significant (FDR < 0.05)",
    labels = c("No", "Yes")
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 8)) +
  labs(
    title = "Condition Comparisons (Tukey HSD)", 
    subtitle = "Module: HALLMARK_APOPTOSIS",
    x = "Condition Comparison", 
    y = "Module × Cluster"
  )

ggsave(file.path("exploratory_analyses/evaluate_pathways",
                 "gb_sb_1467_adj_pval_significant_tukey_comparisons_apoptosis_only_dotplot.pdf"), 
       plot = p2, width = 6, height = 7)


p3 <- ggplot(sig_tukey %>% filter(module == "Rummagene_PMC9314667_Table_1_xlsx_Regulons_ATF4_ATF4"), 
             aes(x = comparison, y = cell_type)) +
  geom_point(aes(color = cohen_d, size = -log10(p.adj), shape = signif)) +
  scale_color_gradient2(
    low = "#053061",   # darker blue
    mid = "#f7f7f7",   # light gray
    high = "#B2182B",  # strong red-orange
    midpoint = 0,
    name = "Cohen's d"
  ) +
  scale_shape_manual(
    values = c("0" = 16, "1" = 17),
    name = "Significant (FDR < 0.05)",
    labels = c("No", "Yes")
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 8)) +
  labs(
    title = "Condition Comparisons (Tukey HSD)", 
    subtitle = "Module: Rummagene_PMC9314667_Table_1_xlsx_Regulons_ATF4_ATF4",
    x = "Condition Comparison", 
    y = "Module × Cluster"
  )

ggsave(file.path("exploratory_analyses/evaluate_pathways",
                 "gb_sb_1467_adj_pval_significant_tukey_comparisons_rummagene_atf4_only_dotplot.pdf"), 
       plot = p3, width = 6, height = 7)


# END