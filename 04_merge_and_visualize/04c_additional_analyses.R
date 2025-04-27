#!/usr/bin/env Rscript

###############################################################################
## Project ID: GB-SB-1467
## Authors: Ayushi Agrawal
##
## Script Goal: Additional analyses
##
## Usage example:
## Rscript 04c_additional_analyses.R
##  --input 'input_seurat_object.RDS' \  # Input Seurat object 
##  --output '/output_directory' \       # Location for output files
##  --output_prefix "outputs_scRNA_"     # Prefix for output files
##
## Run "Rscript 04c_additional_analyses.R --help" for more information
###############################################################################

# Get input arguments -----------------------------------------------------
library(optparse)
option_list <- list(
  make_option(c("-i", "--input"),
              action = "store", default = NA, type = "character",
              help = "Input Seurat object in RDS format (required)"
  ),
  make_option(c("-o", "--output"),
              action = "store", default = NA, type = "character",
              help = "Output directory, will create if it doesn't exist (required)"
  ),
  make_option(c("--output_prefix"),
              action = "store",
              default = "clustering_cell_type",
              type = "character",
              help = "Prefix for output files, [default %default]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# Check if required args are provided
if (is.na(opt$input) | is.na(opt$output) | is.na(opt$output_prefix)) {
  stop("Missing one or more required arguments")
}


# Setup ------------------------------------------------------------------------
# Load required packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)

#set seed
set.seed(42)

#increase RAM usage so the futurized Seurat functions... 
#can access larger global variables
options(future.globals.maxSize = 60 * 1024^3)


# Main analysis ----------------------------------------------------------------
# load data
data <- readRDS(opt$input)

# Fix metadata
data$Condition_all <- data$Condition
data$Condition <- gsub("^Het", "WT", data$Condition)
data$Genotype_all <- data$Genotype
data$Genotype <- gsub("^Het", "WT", data$Genotype)


# generate feature plots
VlnPlot(data, 
        features= c("Naxd", "Naxe"), 
        group.by = "SampleID",
        raster = FALSE,
        ncol = 1) %>%
  ggsave(file = file.path(opt$output, paste0(opt$output_prefix,"_Naxd_Naxe_per_sample_violinplot.pdf")),
         plot = .,
         width = 25,
         height = 15)

VlnPlot(data, 
        features= c("Naxd", "Naxe"), 
        group.by = "Condition",
        raster = FALSE,
        ncol = 2) %>%
  ggsave(file = file.path(opt$output, paste0(opt$output_prefix,"_Naxd_Naxe_per_condition_violinplot.pdf")),
         plot = .,
         width = 15,
         height = 7)

# extract metadata columns 
meta_cols <- names(data@meta.data)[!sapply(data@meta.data, is.numeric)]
# Exclude "orig.ident"
meta_cols <- meta_cols[meta_cols != "orig.ident"]
# Exclude "DF.classifications_" columns
meta_cols <- c(meta_cols[!grepl("^DF\\.classifications_", meta_cols)], 
               "GEM.Well",
               "cell_id",
               "cell_barcode")

# get cell ids
input <- data@meta.data |>
  rownames_to_column("cell_id")

# ArchR format is sample-name#cell-id
input <- input |>
  mutate(cell_barcode = sub(".*_", "", cell_id)) 

input |>
  select(all_of(meta_cols)) |>
  write_csv(file.path(opt$output,
                      paste0(opt$output_prefix,"_post_qc_cell_ids.csv")))



# Write sessionInfo
writeLines(capture.output(sessionInfo()),
           file.path(opt$output, "sessionInfo.txt"))


## END