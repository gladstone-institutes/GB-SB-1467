###############################################################################
## Project ID: GB-SB-1467
## Authors: Ayushi Agrawal
## Script Goal: Run ANOVA and Tukey HSD on per-sample module scores across 
##              conditions within each cluster, save results, and plot significance.
## Note: This script was run locally on Ayushi's laptop
###############################################################################

# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(rstatix)
library(ggpubr)
library(ggsignif)

# set working directory
setwd("/Volumes/Jain-Boinformatics-Collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/results")

# load metadata
meta <- read.csv(file.path("exploratory_analyses/evaluate_cell_death_pathways",
                           "gb_sb_1467_30PC_0.08res_clustered_and_cell_typed_pathway_modules_added_metadata.csv"))

# clean condition labels
meta$Condition <- gsub("-", "wo", gsub("\\+", "w", gsub(" ", "_", meta$Condition)))

# identify module score columns
module_cols <- grep("v2024.1.Mm1$", colnames(meta), value = TRUE)

# Aggregate per sample per cluster
agg_scores <- meta %>%
  select(all_of(module_cols), cluster = seurat_clusters, SampleID, Condition) %>%
  group_by(cluster, SampleID, Condition) %>%
  summarise(across(all_of(module_cols), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

# Output directory
plot_dir <- "exploratory_analyses/evaluate_cell_death_pathways/per_cluster_plots"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Initialize results lists
anova_results <- list()
tukey_results <- list()

# Run ANOVA and Tukey HSD per module × cluster
for (module in module_cols) {
  module_clean <- sub("\\.v2024\\.1\\.Mm[0-9]+$", "", module)
  
  for (cl in unique(agg_scores$cluster)) {
    df <- agg_scores %>% filter(cluster == cl)
    df$Condition <- factor(df$Condition)
    
    aov_res <- anova_test(df, formula = as.formula(paste(module, "~ Condition")))
    tukey_res <- tukey_hsd(df, as.formula(paste(module, "~ Condition")))
    
    aov_res$module <- module_clean
    aov_res$cluster <- cl
    tukey_res$module <- module_clean
    tukey_res$cluster <- cl
    
    anova_results[[paste0(module, "_cl", cl)]] <- aov_res
    tukey_results[[paste0(module, "_cl", cl)]] <- tukey_res
    
    tukey_res_sig <- tukey_res %>%
      add_xy_position(x = "Condition") %>%
      filter(p.adj < 0.05) %>%
      mutate(label = paste0(p.adj.signif, " (Effect size: ", round(estimate, 3), ")"))
    
    p <- ggboxplot(df, x = "Condition", y = module, add = "jitter") +
      stat_pvalue_manual(tukey_res_sig, label = "label", hide.ns = TRUE) +
      labs(
        title = paste("Cluster", cl, "|", module_clean),
        subtitle = get_test_label(aov_res, detailed = TRUE),
        caption = "Tukey HSD: ns (>0.05), * (<=0.05), ** (<=0.01), *** (<=0.001), **** (<=0.0001)\nValues in parentheses are effect sizes (mean difference)",
        y = "Module Score", x = "Condition"
      ) +
      theme_bw()
    
    ggsave(file.path(plot_dir, paste0("gb_sb_1467_30PC_0.08res_cluster", cl, "_", module_clean, "_anova_res.pdf")),
           plot = p, width = 8, height = 8)
  }
}

# Combine ANOVA and Tukey results
anova_df <- map_df(anova_results, ~ as_tibble(as.data.frame(.x)))
tukey_df <- map_df(tukey_results, ~ as_tibble(as.data.frame(.x)))

# Adjust ANOVA p-values (Benjamini-Hochberg)
anova_df$p_adj <- p.adjust(anova_df$p, method = "BH")

# Save results
write.csv(anova_df, 
          "exploratory_analyses/evaluate_cell_death_pathways/gb_sb_1467_30PC_0.08res_anova_summary_module_scores.csv", 
          row.names = FALSE)
write.csv(tukey_df, 
          "exploratory_analyses/evaluate_cell_death_pathways/gb_sb_1467_30PC_0.08res_tukey_results_module_scores.csv", 
          row.names = FALSE)



# Summary plots: faceted boxplots for raw and adjusted p-values
agg_long <- agg_scores %>%
  pivot_longer(cols = all_of(module_cols), names_to = "module", values_to = "score") %>%
  mutate(
    module = sub("\\.v2024\\.1\\.Mm[0-9]+$", "", module),
    module_cluster = paste(module, cluster, sep = "_")
  )

agg_long$cluster <- factor(agg_long$cluster, levels = sort(unique(agg_long$cluster)), 
                           labels = paste0("Cluster ", sort(unique(agg_long$cluster))))

signif_anova_raw <- anova_df %>% filter(p < 0.05) %>% mutate(module_cluster = paste(module, cluster, sep = "_"))
signif_anova_adj <- anova_df %>% filter(p_adj < 0.05) %>% mutate(module_cluster = paste(module, cluster, sep = "_"))

agg_long_raw <- agg_long %>% filter(module_cluster %in% signif_anova_raw$module_cluster)
agg_long_adj <- agg_long %>% filter(module_cluster %in% signif_anova_adj$module_cluster)

p_raw <- ggplot(agg_long_raw, aes(x = Condition, y = score, fill = Condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
  facet_grid(cluster ~ module, scales = "free_y", switch = "y") +
  labs(title = "Module Scores: Raw ANOVA p < 0.05", x = "Condition", y = "Module Score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        strip.text = element_text(size = 8), 
        panel.spacing = unit(1, "lines"))

ggsave("exploratory_analyses/evaluate_cell_death_pathways/gb_sb_1467_raw_pval_summary_module_scores.pdf", 
       plot = p_raw, width = 20, height = 16)

p_adj <- ggplot(agg_long_adj, aes(x = Condition, y = score, fill = Condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
  facet_grid(cluster ~ module, scales = "free_y", switch = "y") +
  labs(title = "Module Scores: FDR-Adjusted ANOVA p < 0.05", x = "Condition", y = "Module Score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        strip.text = element_text(size = 8), 
        panel.spacing = unit(1, "lines"))

ggsave("exploratory_analyses/evaluate_cell_death_pathways/gb_sb_1467_adj_pval_summary_module_scores.pdf", 
       plot = p_adj, width = 16, height = 12)


# Summary dotplot for significant Tukey results
# Define significant module × cluster pairs from ANOVA
signif_anova <- anova_df %>%
  filter(p_adj < 0.05) %>%
  mutate(module_cluster = paste(module, cluster, sep = "_"))

# Filter Tukey results based on significant ANOVA hits
sig_tukey <- tukey_df %>%
  mutate(
    comparison = paste(group1, "vs", group2),
    module = sub("\\.v2024\\.1\\.Mm[0-9]+$", "", module),
    module_cluster = paste(module, cluster, sep = "_")
  ) %>%
  filter(p.adj < 0.05 & module_cluster %in% signif_anova$module_cluster)

sig_tukey$module_cluster <- factor(sig_tukey$module_cluster,
                                   levels = sig_tukey %>%
                                     arrange(module, cluster) %>%
                                     pull(module_cluster) %>% unique())

p3 <- ggplot(sig_tukey, aes(x = comparison, y = module_cluster)) +
  geom_point(aes(size = abs(estimate), color = -log10(p.adj))) +
  scale_color_viridis_c(option = "D", name = "-log10(p.adj)") +
  scale_size_continuous(name = "|Mean Difference|") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size = 8)) +
  labs(title = "Significant Condition Comparisons (Tukey HSD): p.adj < 0.05", 
       x = "Condition Comparison", 
       y = "Module × Cluster")

ggsave("exploratory_analyses/evaluate_cell_death_pathways/gb_sb_1467_adj_pval_significant_tukey_comparisons_dotplot.pdf", 
       plot = p3, width = 8, height = 8)




# Write session info
writeLines(
  capture.output(sessionInfo()),
  file.path("exploratory_analyses/evaluate_cell_death_pathways/", "sessionInfo.txt")
)
print("***** Script Complete *****")

# END