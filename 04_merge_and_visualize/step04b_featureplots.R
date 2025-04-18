# this script was run locally on Ayushi's laptop

# load libraries
library(Seurat)
library(ggplot2)

# set working directory
setwd(file.path("/Volumes/Jain-Boinformatics-Collaboration",
                "sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025",
                "results/04_merge_and_visualize/30PCs"
                )
      )

# load data
dat <- readRDS("gb_sb_1467_30pcs_merged_data_sct_pca_umap.rds")

# generate feature plots
FeaturePlot(dat, 
            raster = FALSE, 
            order = TRUE, 
            label = FALSE, 
            features = "Naxd", 
            reduction = "umap") %>%
  ggsave(file = "gb_sb_1467_30pcs_Naxd_featureplot.pdf",
         plot = .,
         width = 12,
         height = 7)

((FeaturePlot(dat, 
            raster = FALSE, 
            order = TRUE, 
            label = FALSE, 
            features = "Naxd",
            split.by = "Condition",
            reduction = "umap") & theme(legend.position = "right")) +
    patchwork::plot_layout(ncol = 3, nrow = 2)) %>%
  ggsave(file = "gb_sb_1467_30pcs_Naxd_featureplot_splitby_condition.pdf",
         plot = .,
         width = 24,
         height = 12)

(FeaturePlot(dat, 
             raster = FALSE, 
             order = TRUE, 
             label = FALSE, 
             features = "Naxd",
             split.by = "Genotype",
             reduction = "umap") & theme(legend.position = "right")) %>%
  ggsave(file = "gb_sb_1467_30pcs_Naxd_featureplot_splitby_genotype.pdf",
         plot = .,
         width = 25,
         height = 7)


# END