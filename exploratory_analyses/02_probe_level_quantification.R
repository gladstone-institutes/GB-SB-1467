#!/usr/bin/env Rscript
#####################################################################################################
## Project ID: GB-SB-1467
## Authors: Ayushi Agrawal
##
## Script Goal: Check probe level quantification of the reads for Naxd for each sample
##
## Usage example:
## Rscript 02_probe_level_quantification.R \
## --input_dir probe_data_per_sample \    # Directory containing the QC filtered Seurat objects
## --output_dir 02_probe_quantification \ # Output directory
## --output_prefix "GB-SB-1467" \         # Prefix for all output data and files
## --project "GB-SB-1467" \               # Project ID
## --exclude_samples sample1,sample2 \    # Comma separated list of samples to exclude (optional)
## --npcs 15                              # Number of PCs to use for UMAP (optional)
## 
## Run "Rscript 02_probe_level_quantification.R --help" for more information
#####################################################################################################

# Get input arguments -------------------------------------------------------------------------------
library(optparse)
option_list <- list(
  make_option("--input_dir",
              action = "store", default = NA, type = "character",
              help = "Directory containing the QC filtered Seurat objects (required)"
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
if (is.na(opt$input_dir) | is.na(opt$output_dir) | is.na(opt$output_prefix) | is.na(opt$project)) {
  stop("***** ERROR: Missing required arguments! *****")
}


# Load required packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# Create the output directory if it doesn't exist
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = T)
}

#increase RAM usage so the futurized Seurat functions... 
#can access larger global variables
options(future.globals.maxSize = 60 * 1024^3)

# Read in the samplesheet
samples_list <- list.files(opt$input_dir, pattern = "\\.rds$", full.names = TRUE, recursive=FALSE)





# Merge the data into a single Seurat object --------------------------------------------------------------
for(i in 1:length(samples_list)){
  this_sample_name <- basename(samples_list[i])
  print(paste0("Processing ", this_sample_name))
  
  # Read in the Seurat object
  this_obj <- readRDS(samples_list[i])
  
  # Remove any previous normalization
  DefaultAssay(this_obj) <- "RNA"
  this_obj <- DietSeurat(this_obj, assay = "RNA")
  
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

print("Merged Seurat object:")
print(data)
print(table(data$SampleID))
print(table(data$Condition))

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
# get the metadata columns that not numeric
meta_cols <- c(names(data@meta.data)[!sapply(data@meta.data, is.numeric)],
               "Weight.mg",
               "GEM.Well")
# Exclude "orig.ident"
meta_cols <- meta_cols[meta_cols != "orig.ident"]
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


# generate feature plots and violin plots
VlnPlot(data, 
        features = grep("Naxd", rownames(data), value = T), 
        group.by = "SampleID",
        raster = FALSE,
        ncol = 1) %>%
  ggsave(file = file.path(opt$output_dir, paste0(opt$output_prefix,"_Naxd_per_sample_violinplot.pdf")),
         plot = .,
         width = 25,
         height = 20)

VlnPlot(data, 
        features= grep("Naxd", rownames(data), value = T), 
        group.by = "Condition",
        raster = FALSE,
        ncol = 3) %>%
  ggsave(file = file.path(opt$output_dir, paste0(opt$output_prefix,"_Naxd_per_condition_violinplot.pdf")),
         plot = .,
         width = 22,
         height = 7)

VlnPlot(data, 
        features = grep("Naxe", rownames(data), value = T), 
        group.by = "SampleID",
        raster = FALSE,
        ncol = 1) %>%
  ggsave(file = file.path(opt$output_dir, paste0(opt$output_prefix,"_Naxe_per_sample_violinplot.pdf")),
         plot = .,
         width = 25,
         height = 20)

VlnPlot(data, 
        features= grep("Naxe", rownames(data), value = T), 
        group.by = "Condition",
        raster = FALSE,
        ncol = 3) %>%
  ggsave(file = file.path(opt$output_dir, paste0(opt$output_prefix,"_Naxe_per_condition_violinplot.pdf")),
         plot = .,
         width = 22,
         height = 7)

genes_to_plot <- grep("Naxd|Naxe", rownames(data), value = T)
pdf(file = file.path(opt$output_dir, paste0(opt$output_prefix, 
                                        "_Naxd_Naxe_featureplots_splitby_condition.pdf")),
    width = 24, 
    height = 15)
# Loop through genes two at a time
for (i in seq(1, length(genes_to_plot), by = 2)) {
  genes_subset <- genes_to_plot[i:min(i+1, length(genes_to_plot))]  # 2 genes at a time
  plots <- lapply(genes_subset, function(gene) {
    FeaturePlot(dat, 
                raster = FALSE, 
                order = TRUE, 
                label = FALSE, 
                features = gene,
                split.by = "Condition",
                reduction = "umap") & 
      theme(legend.position = "right")
  })
  
  # Combine 2 plots vertically
  combined_plot <- wrap_plots(plots, ncol = 1, nrow = 2)
  
  # Print the combined plot to the PDF
  print(combined_plot)
}
dev.off()


print("***** Visualization completed! *****")



# Write session info
writeLines(
  capture.output(sessionInfo()),
  file.path(opt$output_dir, "sessionInfo.txt")
)

#record logs
print("***** Script completed! *****")

##END 