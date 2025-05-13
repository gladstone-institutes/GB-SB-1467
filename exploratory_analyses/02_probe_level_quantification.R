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
## --output_dir 02_probe_quantification   # Output directory
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
  )
)

# Read in the arguments
opt <- parse_args(OptionParser(option_list = option_list))

# Check that all required arguments are present
if (is.na(opt$input_dir) | is.na(opt$output_dir)) {
  stop("***** ERROR: Missing required arguments! *****")
}


# Load required packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpubr)

# Create the output directory if it doesn't exist
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = T)
}

#increase RAM usage so the futurized Seurat functions... 
#can access larger global variables
options(future.globals.maxSize = 80 * 1024^3)

# Read in the samplesheet
samples_list <- list.files(opt$input_dir, 
                           pattern = "\\.rds$", 
                           full.names = TRUE, 
                           recursive=FALSE)

# create plot lists
naxd_plot_list_sct <- NULL
naxd_plot_list_rna <- NULL
naxe_plot_list_sct <- NULL
naxe_plot_list_rna <- NULL




# Process each single Seurat object --------------------------------------------------------------
for(i in 1:length(samples_list)){
  this_sample_name <- basename(samples_list[i])
  this_sample_name <- sub("_.*$", "", this_sample_name)
  print(paste0("Processing ", this_sample_name))
  
  # Read in the Seurat object
  this_obj <- readRDS(samples_list[i])
  
  # Remove any previous normalization
  DefaultAssay(this_obj) <- "RNA"
  this_obj <- DietSeurat(this_obj, assay = "RNA")
  
  # Perform SCT normalization
  this_obj <- SCTransform(this_obj, vst.flavor = "v2", verbose = TRUE)
  
  # Generate plots for Naxd and Naxe probes (SCT and RNA assay)
  pl <- VlnPlot(this_obj, 
                features = grep("Naxd", rownames(this_obj), value = T), 
                group.by = "SampleID",
                raster = FALSE,
                ncol = 3) &
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          title = element_text(size = 9),
          plot.margin = margin(t = 10, b = 10, r = 3, l = 3))
  naxd_plot_list_sct[[this_sample_name]] <- pl
  
  pl <- VlnPlot(this_obj, 
                assay = "RNA",
                layer = "counts",
                features = grep("Naxd", rownames(this_obj), value = T), 
                group.by = "SampleID",
                raster = FALSE,
                ncol = 3) &
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          title = element_text(size = 9),
          plot.margin = margin(t = 10, b = 10, r = 3, l = 3))
  naxd_plot_list_rna[[this_sample_name]] <- pl
  
  
  pl <- VlnPlot(this_obj, 
                features = grep("Naxe", rownames(this_obj), value = T), 
                group.by = "SampleID",
                raster = FALSE,
                ncol = 3) &
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          title = element_text(size = 9),
          plot.margin = margin(t = 10, b = 10, r = 3, l = 3))
  naxe_plot_list_sct[[this_sample_name]] <- pl
  
  pl <- VlnPlot(this_obj, 
                assay = "RNA",
                layer = "counts",
                features = grep("Naxe", rownames(this_obj), value = T), 
                group.by = "SampleID",
                raster = FALSE,
                ncol = 3) &
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          title = element_text(size = 9),
          plot.margin = margin(t = 10, b = 10, r = 3, l = 3))
  naxe_plot_list_rna[[this_sample_name]] <- pl
  
  # save the normalized object
  saveRDS(this_obj,
          file = file.path(opt$output_dir,
                           paste0(this_sample_name,
                                  "_post_SCT.rds")))
  
  rm(this_obj)
}


naxd_rna_plots <- ggarrange(
  plotlist = naxd_plot_list_rna,
  ncol = 2,
  nrow = ceiling(length(naxd_plot_list_rna) / 2),
  labels = names(naxd_plot_list_rna),
  vjust = 1.3,
  label.x = .4,
  legend = "none"
)
ggsave(filename = file.path(opt$output,"naxd_probes_violinplot_raw_counts_expr.pdf"),
       plot = naxd_rna_plots,
       height = ceiling(length(naxd_plot_list_rna) / 2) * 3,
       width = 8, limitsize = FALSE)


naxd_sct_plots <- ggarrange(
  plotlist = naxd_plot_list_sct,
  ncol = 2,
  nrow = ceiling(length(naxd_plot_list_sct) / 2),
  labels = names(naxd_plot_list_sct),
  vjust = 1.3,
  label.x = .4,
  legend = "none"
)
ggsave(filename = file.path(opt$output,"naxd_probes_violinplot_sct_normalized_expr.pdf"),
       plot = naxd_sct_plots,
       height = ceiling(length(naxd_plot_list_sct) / 2) * 3,
       width = 8, limitsize = FALSE)


naxe_rna_plots <- ggarrange(
  plotlist = naxe_plot_list_rna,
  ncol = 2,
  nrow = ceiling(length(naxe_plot_list_rna) / 2),
  labels = names(naxe_plot_list_rna),
  vjust = 1.3,
  label.x = .4,
  legend = "none"
)
ggsave(filename = file.path(opt$output,"naxe_probes_violinplot_raw_counts_expr.pdf"),
       plot = naxe_rna_plots,
       height = ceiling(length(naxe_plot_list_rna) / 2) * 3,
       width = 8, limitsize = FALSE)


naxe_sct_plots <- ggarrange(
  plotlist = naxe_plot_list_sct,
  ncol = 2,
  nrow = ceiling(length(naxe_plot_list_sct) / 2),
  labels = names(naxe_plot_list_sct),
  vjust = 1.3,
  label.x = .4,
  legend = "none"
)
ggsave(filename = file.path(opt$output,"naxe_probes_violinplot_sct_normalized_expr.pdf"),
       plot = naxe_sct_plots,
       height = ceiling(length(naxe_plot_list_sct) / 2) * 3,
       width = 8, limitsize = FALSE)



print("***** Visualization completed! *****")



# Write session info
writeLines(
  capture.output(sessionInfo()),
  file.path(opt$output_dir, "sessionInfo.txt")
)

#record logs
print("***** Script completed! *****")

##END 
