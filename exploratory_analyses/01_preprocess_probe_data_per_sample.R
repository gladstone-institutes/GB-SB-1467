#!/usr/bin/env Rscript
#####################################################################################################
## Project ID: GB-SB-1467
## Authors: Ayushi Agrawal
##
## Script Goal: Pre-process probe level quantification for each sample
##
## Usage example:
## Rscript 01_preprocess_probe_data_per_sample.R \
## --samplesheet input.csv \              # Csv file containing the cellranger multi output location
## --metadata sample_metadata.csv \       # Sample metadata sheet
## --output_dir 07_merge_and_visualize \  # Output directory
## --output_prefix "GB-SB-1467" \         # Prefix for all output data and files
## --subset_cells cell_ids.csv            # cell ids of filtered cells per sample
## 
## Run "Rscript 01_preprocess_probe_data_per_sample.R --help" for more information
#####################################################################################################

# Get input arguments -------------------------------------------------------------------------------
library(optparse)

option_list <- list(
  make_option("--samplesheet",
              action = "store", default = NA, type = "character",
              help = "Path to csv file containing the cellranger multi output location (required)"
  ),
  make_option("--metadata",
              action = "store", default = NA, type = "character",
              help = "A CSV sample metadata sheet (required)"
  ),
  make_option("--output_dir",
              action = "store", default = NA, type = "character",
              help = "Output directory (required)"
  ),
  make_option("--output_prefix",
              action = "store", default = NA, type = "character",
              help = "Prefix for all output data and files (required)"
  ),
  make_option("--subset_cells",
              action = "store", default = NA, type = "character",
              help = "Path to subset file (CSV format) containing cell ids of filtered cells per sample (required)."
  )
)

# Read in the arguments
opt <- parse_args(OptionParser(option_list = option_list))

# Check that all required arguments are present
if (is.na(opt$samplesheet) | is.na(opt$metadata) | is.na(opt$output_dir) | is.na(opt$output_prefix) | is.na(opt$subset_cells)) {
  stop("***** ERROR: Missing required arguments! *****")
}


# Load required packages
library(Seurat)
library(dplyr)

# Create the output directory if it doesn't exist
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = T)
}

#increase RAM usage so the futurized Seurat functions... 
#can access larger global variables
options(future.globals.maxSize = 100 * 1024^3)

# Read in the metadata
meta_data <- read.csv(opt$metadata)
meta_data$SampleID <- paste0("sample", meta_data$SampleID)

# Read in the samplesheet
samples_list <- read.csv(opt$samplesheet)
samples_list <- samples_list$input

# Read in the cell ids to keep
keep_cell_ids <- read.csv(opt$subset_cells)




# Process the data and save Seurat object --------------------------------------------------------------
for (i in 1:length(samples_list)) {
  this_sample_name <- basename(samples_list[i])
  print(paste0("Processing sample ", this_sample_name))
  
  # Read in the Seurat object
  pre <- Read10X_h5(file.path(samples_list[i], "count/sample_raw_probe_bc_matrix.h5"))
  this_obj <- CreateSeuratObject(pre, min.cells = 4)
  rm(pre)
  
  # Remove any previous normalization
  DefaultAssay(this_obj) <- "RNA"
  this_obj <- DietSeurat(this_obj, assay = "RNA")
  
  # rename cells
  this_obj <- RenameCells(object = this_obj, add.cell.id = this_sample_name)
  
  # Calculate the percentage of mitochondrial genes
  this_obj[["percent.mt"]] <- PercentageFeatureSet(this_obj, pattern = "^mt-")
  
  # filter for QC'ed cells from main analysis
  keep <- keep_cell_ids[keep_cell_ids$SampleID == this_sample_name, "cell_id"]
  this_obj <- subset(this_obj, cells = keep)
  
  # Add the sample metadata to the Seurat object
  this_meta <- as.data.frame(meta_data[meta_data$SampleID == this_sample_name,])
  this_meta <- this_meta[rep(1, ncol(this_obj)), ]
  rownames(this_meta) <- colnames(this_obj)
  this_obj <- AddMetaData(object = this_obj, metadata =this_meta)
  
  # Fix metadata
  this_obj$Condition_all <- this_obj$Condition
  this_obj$Condition <- gsub("^Het", "WT", this_obj$Condition)
  this_obj$Genotype_all <- this_obj$Genotype
  this_obj$Genotype <- gsub("^Het", "WT", this_obj$Genotype)
  
  # save object
  saveRDS(this_obj, 
          file = file.path(opt$output_dir, 
                           paste0(this_sample_name, "_preprocessed.rds")))
  
  rm(this_obj)
  gc()
}

print("***** All samples processed individually and saved! *****")



# Write session info
writeLines(
  capture.output(sessionInfo()),
  file.path(opt$output_dir, "sessionInfo.txt")
)

#record logs
print("***** Script completed! *****")

##END 