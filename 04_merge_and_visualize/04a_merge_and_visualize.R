#!/usr/bin/env Rscript
#####################################################################################################
## Project ID: GB-SB-1467
## Authors: Ayushi Agrawal
##
## Script Goal: Merge and normalize all Seurat objects; and visualize metadata for batch effects
##
## Usage example:
## Rscript 04_merge_and_visualize.R \
## --input_dir 05_doublet_removal/ \  # Directory containing the QC filtered Seurat objects
## --metadata sample_metadata.csv \ # Sample metadata sheet
## --output_dir 07_merge_and_visualize \  # Output directory
## --output_prefix "GB-SB-1467" \  # Prefix for all output data and files
## --project "GB-SB-1467" \ # Project ID
## --exclude_samples sample1,sample2,sample3 \ # Comma separated list of samples to exclude (optional)
## --npcs 15 \ # Number of PCs to use for UMAP (optional)
## 
## Run "Rscript 04_merge_and_visualize.R --help" for more information
#####################################################################################################

# Get input arguments -------------------------------------------------------------------------------
library(optparse)
option_list <- list(
  make_option("--input_dir",
              action = "store", default = NA, type = "character",
              help = "Input directory (required)"
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
  make_option("--project",
              action = "store", default = NA, type = "character",
              help = "Project ID for the Seurat objects (required)"
  ),
  make_option("--exclude_samples",
              action = "store", default = NA, type = "character",
              help = "Comma separated list of samples to exclude, the sample names should match the RDS file names (optional)"
  ),
  make_option("--npcs",
              action = "store", default = 20, type = "integer",
              help = "Number of PCs to use for UMAP, default is 20 (optional)"
  )
)
# Read in the arguments
opt <- parse_args(OptionParser(option_list = option_list))

# Check that all required arguments are present
if (is.na(opt$input_dir) | is.na(opt$metadata) | is.na(opt$output_dir) | is.na(opt$output_prefix) | is.na(opt$project)) {
  stop("***** ERROR: Missing required arguments! *****")
}


# Load required packages
library(Seurat)
library(dplyr)
library(ggplot2)

# Set the working directory
setwd(opt$input_dir)

# Create the output directory if it doesn't exist
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = T)
}

#increase RAM usage so the futurized Seurat functions... 
#can access larger global variables
options(future.globals.maxSize = 50 * 1024^3)

# Read in the sample sheet CSV
meta_data <- read.csv(opt$metadata)
meta_data$SampleID <- paste0("sample", meta_data$SampleID)




# Merge the data into a singel Seurat object --------------------------------------------------------------
# Read in all the sample datasets
samples_list <- list.files(".", pattern = "\\.rds$", full.names = TRUE, recursive=TRUE)

for(i in 1:length(samples_list)){
  this_sample_name <- basename(dirname(samples_list[i]))
  print(paste0("Processing sample ", this_sample_name))
  
  # Read in the Seurat object
  this_obj <- readRDS(samples_list[i])
  
  # Remove any previous normalization
  DefaultAssay(this_obj) <- "RNA"
  this_obj <- DietSeurat(this_obj, assay = "RNA")
  
  # rename cells
  this_obj <- RenameCells(object = this_obj, add.cell.id = this_sample_name)
  
  # Add the sample metadata to the Seurat object
  this_meta <- as.data.frame(meta_data[meta_data$SampleID == this_sample_name,])
  this_meta <- this_meta[rep(1, ncol(this_obj)), ]
  rownames(this_meta) <- colnames(this_obj)
  this_obj <- AddMetaData(object = this_obj, metadata =this_meta)
  
  # merge into a single Seurat object
  if(i == 1){
    data <- this_obj
  } else {
    data <- merge(data, 
                  this_obj,
                  project = opt$project, 
                  merge.data=FALSE)
    data[["RNA"]] <- JoinLayers(data[["RNA"]])
  }
  rm(this_obj)
}

print("***** Dataset merge completed! *****")




# SCTransform normalization, PCA and UMAP --------------------------------------------------------------
data <- SCTransform(data, vst.flavor = "v2", verbose = TRUE) %>%
  RunPCA(assay = "SCT", npcs = 50) %>%
  RunUMAP(reduction = "pca", dims =  1:opt$npcs, verbose = TRUE)

saveRDS(data, 
        file = file.path(opt$output_dir, 
                         paste0(opt$output_prefix, 
                                "_merged_data_sct_pca_umap.rds")))
print("***** SCT normalization, PCA and UMAP completed! *****")




# Cell cycle scoring ------------------------------------------------------------------------------------
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
data <- CellCycleScoring(data, 
                         s.features = s.genes, 
                         g2m.features = g2m.genes, 
                         set.ident = TRUE)
data$Phase <- ifelse(data$Phase == "S", "G1S", data$Phase)
data$Phase <- ifelse(data$Phase == "G1", "non-cycling", data$Phase)




# Visualize the data ------------------------------------------------------------------------------------
# PCA plots
pdf(file = file.path(opt$output_dir, paste0(opt$output_prefix, "_PCA_plots.pdf")),
    width = 12, height = 7)
print(ElbowPlot(data, ndims = 50))
print(ElbowPlot(data, ndims = 20))
print(VizDimLoadings(data, dims = 1:2, reduction = "pca"))
print(DimPlot(data, reduction = "pca", raster = FALSE))
print(DimHeatmap(data, dims = 1:2, cells = 500, balanced = TRUE))
print(DimHeatmap(data, dims = 3:4, cells = 500, balanced = TRUE))
dev.off()

# Function to generate visualizations for the metadata variables
plot_metadata <- function(sc_dat, metadata_variable){
  DimPlot(sc_dat, 
          raster = FALSE, 
          order = TRUE, 
          label = FALSE, 
          group.by=metadata_variable, 
          reduction = "umap") %>%
    ggsave(file = file.path(opt$output_dir, 
                            paste0(opt$output_prefix, "_", metadata_variable, "_umap.pdf")),
           plot = .,
           width = 12,
           height = 7)
  pdf_width <- (min(length(unique(sc_dat@meta.data[[metadata_variable]])),4)) * 6  
  pdf_height <- ifelse(length(unique(sc_dat@meta.data[[metadata_variable]])) <= 4,
                       7,
                       5 * ceiling(length(unique(sc_dat@meta.data[[metadata_variable]]))/4))
  (DimPlot(sc_dat, 
           raster = FALSE, 
           order = TRUE, 
           label = FALSE, 
           group.by=metadata_variable,
           split.by=metadata_variable, 
           reduction = "umap",
           ncol = 4) + NoLegend()) %>%
    ggsave(file = file.path(opt$output_dir, 
                            paste0(opt$output_prefix, "_", metadata_variable, "_split_umap.pdf")),
           plot = .,
           width = pdf_width,
           height = pdf_height,
           limitsize = FALSE)
}

# UMAPs for the metadata variables
meta_cols <- c(colnames(meta_data), "Phase")
data@meta.data[meta_cols] <- lapply(data@meta.data[meta_cols], as.factor)
lapply(meta_cols, function(col_name) {
  plot_metadata(data, col_name)
})

# Feature plots for numeric data
for(meta_feature in c("nCount_RNA","nFeature_RNA","percent.mt","nCount_SCT","nFeature_SCT")){
  FeaturePlot(data, 
              raster = FALSE, 
              order = TRUE, 
              label = FALSE, 
              features = meta_feature, 
              reduction = "umap") %>%
    ggsave(file = file.path(opt$output_dir, 
                            paste0(opt$output_prefix, "_", meta_feature, "_featureplot.pdf")),
           plot = .,
           width = 12,
           height = 7)
}

print("***** Visualization completed! *****")




# Write session info
writeLines(
  capture.output(sessionInfo()),
  file.path(opt$output_dir, "sessionInfo.txt")
)

#record logs
print("***** Script completed! *****")

########################## END ########################## 