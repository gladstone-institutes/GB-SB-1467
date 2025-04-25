#!/usr/bin/env Rscript

###############################################################################
## Project ID: GB-SB-1467
## Authors: Reuben Thomas, Natalie Elphick, Ayushi Agrawal
##
## Script Goal: Generate marker genes plots and run Seurat::FindAllMarkers().
##
## Usage example:
## Rscript 06b_plots_findallmarkers.R
##  --input 'input_seurat_object.RDS' \  # Input Seurat object 
##  --output '/output_directory' \       # Location for output files
##  --output_prefix "outputs_scRNA_"     # Prefix for output files
##
## Run "Rscript 06b_plots_findallmarkers.R --help" for more information
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
if (is.na(opt$input) | is.na(opt$output)) {
  stop("Missing one or more required arguments")
}

# Load required packages
library(dplyr)
library(Seurat)
library(ggplot2)

#set seed
set.seed(42)

#increase RAM usage so the futurized Seurat functions... 
#can access larger global variables
options(future.globals.maxSize = 60 * 1024^3)




# Gene expression plots --------------------------------------------------------
# Read in the Seurat object
data <- readRDS(opt$input)

# genes from proteomics results
proteomics_genes <- c("Pcna", "Cdk1", "Pbk", "Kpna2",
                      "Apoa4", "Dhcr24", "Dgat1", "Msmo1", "Gykl1",
                      "Ctsd", "Cst3", "Tpp1",
                      "Mt1", "Mt2", "Naxe")
# mitochondrial genes
mt_genes <- grep("^mt-", rownames(data), value = TRUE)
features_to_plot <- c("Naxd",
                      proteomics_genes[proteomics_genes %in% rownames(data)],
                      mt_genes)

#Generate featureplots
pdf(file = file.path(opt$output, paste0(opt$output_prefix, "_featureplot.pdf")),
    width = 25,
    height = 15)
for(i in 1:ceiling(ifelse(length(features_to_plot)%%2,
                          length(features_to_plot)+1,
                          length(features_to_plot))/8)){
  if((8*i) > length(features_to_plot)){
    print(FeaturePlot(data, 
                      features= features_to_plot[((8*(i-1))+1):(length(features_to_plot))], 
                      raster = FALSE, 
                      order = TRUE, 
                      label = TRUE, 
                      reduction = "umap", 
                      ncol = 4))
  }else{
    print(FeaturePlot(data, 
                      features= features_to_plot[((8*(i-1))+1):(8*i)], 
                      raster = FALSE, 
                      order = TRUE, 
                      label = TRUE, 
                      reduction = "umap", 
                      ncol = 4))
  }
}
dev.off()

#Generate featureplots split by sample and condition for all the marker genes of interest
pdf(file = file.path(opt$output, paste0(opt$output_prefix, "_split_condition_featureplot.pdf")),
    width = 25,
    height = 7)
for(i in 1:length(features_to_plot)){
  print(FeaturePlot(data, 
                    features= features_to_plot[i], 
                    split.by = "Condition",
                    raster = FALSE, 
                    order = TRUE, 
                    label = TRUE, 
                    reduction = "umap"))
}
dev.off()

#Generate violin plots for all the marker genes of interest
pdf(file = file.path(opt$output, paste0(opt$output_prefix, "_violinplot.pdf")),
    width = 25,
    height = 15)
for(i in 1:ceiling(ifelse(length(features_to_plot)%%2,
                          length(features_to_plot)+1,
                          length(features_to_plot))/8)){
  if((8*i) > length(features_to_plot)){
    print(VlnPlot(data, 
                  features= features_to_plot[((8*(i-1))+1):(length(features_to_plot))],
                  raster = FALSE,
                  ncol = 4))
  }else{
    print(VlnPlot(data, 
                  features= features_to_plot[((8*(i-1))+1):(8*i)], 
                  raster = FALSE,
                  ncol = 4))
  }
}
dev.off()

#Generate violin plots split by condition for all the marker genes of interest
pdf(file = file.path(opt$output, paste0(opt$output_prefix,"_split_condition_violinplot.pdf")),
    width = 25,
    height = 7)
for(i in 1:length(features_to_plot)){
  print(VlnPlot(data, 
                features= features_to_plot[i], 
                raster = FALSE,
                split.by = "Condition"))
  print(VlnPlot(data, 
                features= features_to_plot[i],
                raster = FALSE,
                group.by = "Condition"))
}
dev.off()

#Generate violin plots split by sample for all the marker genes of interest
pdf(file = file.path(opt$output, paste0(opt$output_prefix,"_split_sampleid_violinplot.pdf")),
    width = 35,
    height = 7)
for(i in 1:2){
  print(VlnPlot(data, 
                features= features_to_plot[i], 
                raster = FALSE,
                split.by = "SampleID"))
  print(VlnPlot(data, 
                features= features_to_plot[i],
                raster = FALSE,
                group.by = "SampleID"))
}
dev.off()

#Generate a dot plot for all marker genes of interest
h <- ifelse(length(unique(data$seurat_clusters)) < 10,
            7,
            7 + ceiling(length(unique(data$seurat_clusters))/10))
pdf(file = file.path(opt$output, paste0(opt$output_prefix, "_group_seurat_clusters_dotplot.pdf")),
    width = 25,
    height = h)
print(DotPlot(data, 
              features = features_to_plot) + 
        RotatedAxis())   
dev.off()

#Generate a dot plot for all marker genes of interest group by condition
pdf(file = file.path(opt$output, paste0(opt$output_prefix, "_group_condition_dotplot.pdf")),
    width = 20)
print(DotPlot(data, 
              features = features_to_plot,
              group.by = "Condition") + 
        RotatedAxis())   
dev.off()

#Generate a dot plot for all marker genes of interest group by sampleid
pdf(file = file.path(opt$output, paste0(opt$output_prefix, "_group_sampleid_dotplot.pdf")),
    width = 20)
print(DotPlot(data, 
              features = features_to_plot,
              group.by = "SampleID") + 
        RotatedAxis())   
dev.off()


#Generate a dot plot for all marker genes of interest group by seurat clusters and split by sampleid
#pdf_height <- (length(unique(data$seurat_clusters)) * length(unique(data$SampleID))) / 5
custom_colors <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
  "#666666", "#AD7700", "#1C91D4", "#007756", "#D5C711", "#005685", "#A04700",
  "#B14380", "#4D4D4D", "#FFBE2D", "#80C7EF", "#00F6B3", "#F4EB71", "#06A5FF",
  "#FF8320", "#D99BBD", "#8C8C8C", "#FFCB57", "#9AD2F2", "#2CFFC6", "#F6EF8E",
  "#38B7FF", "#FF9B4D", "#E0AFCA", "#A3A3A3", "#8A5F00", "#1674A9", "#005F45",
  "#AA9F0D", "#00446B", "#803800", "#8D3666", "#3D3D3D"
) 
pdf(file = file.path(opt$output, paste0(opt$output_prefix, "_group_seurat_clusters_split_sampleid_dotplot.pdf")),
    height = 45,
    width = 20)
print(DotPlot(data, features = features_to_plot, split.by = "SampleID", 
        cols = custom_colors) + 
  RotatedAxis())
dev.off()




# Find Markers --------------------------------------------------------------
print("**** FindAllMarkers start ****")
data <- PrepSCTFindMarkers(data)
# Find differentially expressed genes for identifying the clusters
marker_res <- FindAllMarkers(data, 
                             assay = "SCT", 
                             slot = "data", 
                             test.use = "wilcox",
                             min.pct = 0.25)
marker_res %>%
  write.csv(., paste0(opt$output, 
                      "/", 
                      opt$output_prefix, 
                      "_FindAllMarkers_marker_genes_per_cluster.csv"), 
            row.names = FALSE)
print("**** FindAllMarkers end ****")




# Write session info
writeLines(
  capture.output(sessionInfo()),
  file.path(opt$output, "sessionInfo_step06b.txt")
)
print("***** Script Complete *****")




# END