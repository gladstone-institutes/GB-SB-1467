#!/usr/bin/env Rscript

###############################################################################
## Project ID: GB-SB-1467
## Authors: Ayushi Agrawal, Natalie Elphick, Michela Traglia
## Script Goal: Generate plots of QC metrics for a series of scRNA-seq samples
## Usage example: Rscript 02b_summarize_doublet_removal.R
###############################################################################

setwd("/gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/results/02_doublet_removal/")

library(tidyverse)

results <- data.frame()
sample_col <- character()

samples <- list.dirs(".", full.names = FALSE, recursive = F)

for(i in 1:length(samples)){
  tmp <- read_csv(paste0(samples[i],"/",samples[i],"_QC_metrics_and_doublets_removed.csv"))
  results <- as.data.frame(rbind(results,tmp))
  sample_col <- c(sample_col, samples[i])
}

results <- as.data.frame(cbind(Sample=sample_col, results))
write_csv(x = results,file = "doublet_removal_summary.csv")

results <- results %>%
  mutate(Percent_doublets_removed = round((Doublet/ (Doublet + Singlet))*100,
                                          digits = 3),
         Percent_filtered = round(((prefilter_ncells-postfilter_ncells)/prefilter_ncells)*100,
                                  digits = 2))

p1 <- results %>%
  select(Sample,prefilter_ncells,postfilter_ncells) %>%
  pivot_longer(!Sample, values_to = "Cells") %>%
  mutate(name = ifelse(name == "prefilter_ncells","Pre-filtering","Post-filtering")) %>%
  ggplot(aes(x = Sample, y = Cells, 
             fill = factor(name, levels=c("Pre-filtering","Post-filtering")))) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank()
  ) +
  xlab("Sample") +
  ylab("Number of cells") +
  ggtitle("Number of cells per sample before and after filtering")

write_csv(x = results,file = "pre_post_filtering_metrics.csv")
ggsave(filename = "pre_post_percent_mt_nFeature_RNA_filtering.png",
       dpi = 300,
       units = "in",
       width = 20,
       height = 7)


#record logs
print("***** Script completed! *****")


########################## END ########################## 