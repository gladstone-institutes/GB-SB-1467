#!/usr/bin/env Rscript

###############################################################################
# Project ID: GB-SB-1467
# Authors: Ayushi Agrawal
# Project goal: Analysis of proteomics data
# Script goal: Generate PCA plot of samples using cell proteome data 
###############################################################################

#set working directory
working_dir <- "/Volumes/Jain-Boinformatics-Collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/"
setwd(working_dir)

#load libraries
library(dplyr)
library(ggplot2)
library(readxl)
library(reshape2)

# create output folder
if(!dir.exists("results/09_proteomics_data_analysis")){
  dir.create("results/09_proteomics_data_analysis")
}

# load data
# read the proteomics data
input_file <- "data/Proteomics_NAXD_15samples_Forcore.xlsx"
cell_data <- input_file %>%
  read_xlsx(.) %>% as.data.frame()

#rownames(cell_data) <- cell_data$`Master accession number`
#keep only samples as columns
cell_data <- cell_data[,(colnames(cell_data) %in% LETTERS[1:16])]
colnames(cell_data) <- paste0("Sample", colnames(cell_data))
cell_data <- as.data.frame(lapply(cell_data, function(x) as.numeric(as.character(x))))

# Calculate log intensities
log_intensities <- log2(cell_data + 1)  # Adding 1 to avoid log(0)

# Create box plots to visualize the data
pdf("results/09_proteomics_data_analysis/proteomics_data_boxplot.pdf", width = 14)
print(ggplot(data = melt(log_intensities), aes(x=variable, y=value)) + 
        geom_boxplot(fill = "lightblue", color = "black") +
        labs(title = "Box Plots of Log Intensities per Sample (Proteomics Data)",
             x = "Sample", y = "Log Intensity") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
              plot.title = element_text(face = "bold", hjust = 0.5)))
dev.off()


# END