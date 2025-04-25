#!/usr/bin/env Rscript

###############################################################################
## Project ID: GB-SB-1467
## Authors: Reuben Thomas, Natalie Elphick, Ayushi Agrawal
##
## Script Goal: Normalize using SCTransform, use input PCA dim and resolution
## to cluster cells and use ScType to assign cell types to the clusters.
##
## Usage example:
## Rscript 06a_clustering_cell_type.R
##  --input 'input_seurat_object.RDS' \  # Input Seurat object 
##  --output '/output_directory' \       # Location for output files
##  --output_prefix "outputs_scRNA_" \   # Prefix for output files
##  --tissue "Brain" \                   # Tissue type to supply to ScType
##  --ndim 20 \                          # Number of PCs to use
##  --resolution 0.02 \                  # Resolution to use for Seurat clustering
##  --batch_var \                        # Batch variable for Harmony correction
##  --custom_db "custom_db.xlsx"         # Path to custom database for ScType
##
## Run "Rscript 06a_clustering_cell_type.R --help" for more information
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
  ),
  make_option(c("--tissue"),
              action = "store", default = NA, type = "character",
              help = "Tissue type to supply to ScType (required)"
  ),
  make_option(c("-n", "--ndim"),
              action = "store", default = NA, type = "numeric",
              help = "Number of principal components to use (required)"
  ),
  make_option(c("-r", "--resolution"),
              action = "store", default = NA, type = "numeric",
              help = "Resolution parameter for FindClusters (required)"
  ),
  make_option("--batch_var",
              action = "store", default = NA, type = "character",
              help = "If provided, the batch variable will be used as the grouping variable for harmony batch correction, sample_label must also be provided (optional)"
  ),
  make_option("--custom_db",
              action = "store", default = NA, type = "character",
              help = "Path to custom database for ScType, must be a xlsx in the same format as the original (optional)"
  ),
  make_option("--exclude_sample",
              action = "store", default = NA, type = "character",
              help = "Sample to remove (optional)"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# Check if required args are provided
if (is.na(opt$resolution) | is.na(opt$input) | is.na(opt$output) | is.na(opt$ndim) | is.na(opt$tissue)) {
  stop("Missing one or more required arguments")
}

# Load required packages
library(dplyr)
library(Seurat)
library(HGNChelper)
library(openxlsx)
library(magrittr)
library(ggplot2)
library(readr)
library(tidyverse)
library(readxl)

# update the name of the results folders
opt$output <- file.path(opt$output, paste0(opt$ndim, "PC_", opt$resolution, "res"))

# update the output prefix for the input resolution
opt$output_prefix <- paste0(opt$output_prefix, "_", opt$ndim, "PC_", opt$resolution, "res")

#set seed
set.seed(42)

#increase RAM usage so the futurized Seurat functions... 
#can access larger global variables
options(future.globals.maxSize = 60 * 1024^3)


# Functions --------------------------------------------------------------------
generate_custom_sctype_db <- function(path_to_db_file){
  # Authors: Natalie Elphick, Ayushi Agrawal
  # Goal: Create a custom ScType database using the CellMarker 2.0 mouse marker gene database
  
  # Read in CellMarker db
  cellmark <- read_xlsx(path_to_db_file) |>
    filter(cell_type != "Cancer cell") |>
    rename(cellName = cell_name,
           tissueType = tissue_type,
           gene = Symbol) |>
    select(cellName, tissueType, gene) |>
    drop_na() 
  
  # Collapse gene lists
  cellmark <- cellmark |>
    distinct() |>
    group_by(cellName, tissueType) |>
    summarise(geneSymbolmore1 = toString(gene), .groups = 'drop') |>
    ungroup() |>
    mutate(geneSymbolmore1 = gsub(" ", "", geneSymbolmore1),
           cellName = gsub("_",":",cellName),
           shortName = cellName,
           geneSymbolmore2 = "") |>
    select(tissueType, cellName, geneSymbolmore1, geneSymbolmore2, shortName)
  
  # Extract the file path and modify the file name
  file_path <- dirname(path_to_db_file) # Get the directory of the file
  file_name <- basename(path_to_db_file)
  
  # Construct the new file name
  new_file_name <- paste0("sctype_formatted_", file_name)
  
  write.xlsx(cellmark,
             file = file.path(file_path, new_file_name))
  
  return(file.path(file_path, new_file_name))
}


gene_sets_prepare_modified <- function(path_to_db_file, cell_type, 
                                       species = "human", gene_format = NA){
  # GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
  # Written by Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021
  #
  # gene_sets_prepare: prepare gene sets and calculate marker sensitivity from input Cell Type excel file
  #
  # @params: path_to_db_file - DB file with cell types
  # @cell_type - cell type (e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain)
  # @species - human or mouse
  # @gene_format - HGNC or NA
  # 
  # Note: This function was modified per issue: 
  # https://github.com/IanevskiAleksandr/sc-type/issues/4
  
  cell_markers = openxlsx::read.xlsx(path_to_db_file)
  cell_markers = cell_markers[cell_markers$tissueType == cell_type,] 
  cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1); cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  # correct gene symbols from the given DB (up-genes)
  cell_markers$geneSymbolmore1 = sapply(1:nrow(cell_markers), function(i){
    
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore1[i],",")))
    if(species == "human"){
      markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""]) # only for humans
    }
    markers_all = sort(markers_all)
    
    if(length(markers_all) > 0){
      if(is.na(gene_format)){
        # this is used if the symbols are not in the format of HGNC
        suppressMessages({markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))})
      }
      paste0(markers_all, collapse=",")
    } else {
      ""
    }
  })
  
  # correct gene symbols from the given DB (down-genes)
  cell_markers$geneSymbolmore2 = sapply(1:nrow(cell_markers), function(i){
    
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore2[i],",")))
    if(species == "human"){
      markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""]) # only for humans
    }
    markers_all = sort(markers_all)
    
    if(length(markers_all) > 0){
      if(is.na(gene_format)){
        # this is used if the symbols are not in the format of HGNC
        suppressMessages({markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))})
      }
      paste0(markers_all, collapse=",")
    } else {
      ""
    }
  })
  
  cell_markers$geneSymbolmore1 = gsub("///",",",cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1)
  
  cell_markers$geneSymbolmore2 = gsub("///",",",cell_markers$geneSymbolmore2)
  cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  gs = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]),","))))
  names(gs) = cell_markers$cellName
  
  gs2 = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]),","))))
  names(gs2) = cell_markers$cellName
  
  list(gs_positive = gs, gs_negative = gs2)
}


run_sc_type <- function(custom_db=NA, tissue, db_prefix){
  # This is the function to run sc-type annotation
  # If a custom database is provided, use that, otherwise use the default
  if (!is.na(custom_db)) {
    gs_list <- gene_sets_prepare_modified(custom_db, tissue, "mouse", "HGNC")
  } else {
    gs_list <- gene_sets_prepare("/opt/sc-type/ScTypeDB_full.xlsx", tissue)
  }
  
  # get cell-type by cell matrix
  es.max <- sctype_score(
    scRNAseqData = as.matrix(data[["SCT"]]@scale.data), scaled = TRUE,
    gs = gs_list$gs_positive, gs2 = gs_list$gs_negative,
    gene_names_to_uppercase = 0 # for mouse 
  )
  
  rm(gs_list)
  print("---ScType Score completed---")
  
  # merge by cluster
  cL_resutls <- do.call(
    "rbind",
    lapply(
      unique(data@meta.data$seurat_clusters),
      function(cl) {
        es.max.cl <- sort(rowSums(es.max[, rownames(data@meta.data[data@meta.data$seurat_clusters == cl, ])]), decreasing = !0)
        head(
          data.frame(
            cluster = cl,
            type = names(es.max.cl),
            scores = es.max.cl,
            ncells = sum(data@meta.data$seurat_clusters == cl)
          ),
          10
        )
      }
    )
  )
  # clear memory
  rm(es.max)
  
  sctype_scores <- cL_resutls %>%
    group_by(cluster) %>%
    top_n(n = 1, wt = scores)
  
  # set low-confident (low ScType score) clusters to "unknown"
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells / 4] <- "Unknown"
  print(sctype_scores[, 1:3])
  
  write_csv(sctype_scores[order(sctype_scores$cluster),],
            file = paste0(
              opt$output,
              "/",
              opt$output_prefix,
              "_sctype_",
              db_prefix,
              "_celltype_annotation_table.csv"
            )
  )
  
  return(sctype_scores)
}


# Cluster using optimized resolution -------------------------------------------

# Read in the Seurat object
data <- readRDS(opt$input)

# Fix metadata
data$Condition_all <- data$Condition
data$Condition <- gsub("^Het", "WT", data$Condition)
data$Genotype_all <- data$Genotype
data$Genotype <- gsub("^Het", "WT", data$Genotype)

# Check for sample exclusion
if (!is.na(opt$exclude_sample)) {
  print(table(data$Sample_name))
  print(paste0("Removing sample ", opt$exclude_sample))
  data <- subset(data, subset = Sample_name == opt$exclude_sample, invert = TRUE)
  print(table(data$Sample_name))
  
  # Drop previous normalization
  DefaultAssay(data) <- "RNA"
  data <- DietSeurat(data, assay = "RNA")
  
  # Update output names
  opt$output <- paste0(opt$output, "_exclude_", opt$exclude_sample)
  opt$output_prefix <- paste0(opt$output_prefix, "_exclude_", opt$exclude_sample)
}

# If SCT assay is missing, rerun normalization
if (!"SCT" %in% Assays(data)) {
  print("***** Running SCTransform normalization *****")
  data <- SCTransform(data, vst.flavor = "v2", verbose = FALSE)
}

# Run integration and clustering
if (!is.na(opt$batch_var)) {
  library(harmony)
  print(paste0("***** Correcting for batch effects using ", opt$batch_var, " *****"))
  data <- RunPCA(data, npcs = opt$ndim)
  data <- RunHarmony(data,
                     assay.use = "SCT",
                     group.by.vars = opt$batch_var,
                     kmeans_init_nstart = 20,
                     kmeans_init_iter_max = 100)
  data <- RunUMAP(data, reduction = "harmony", dims = 1:opt$ndim)
  data <- FindNeighbors(data, dims = 1:opt$ndim, reduction = "harmony")
} else {
  print("***** No batch variable provided, skipping batch correction *****")
  data <- RunPCA(data, npcs = 50, assay = "SCT")
  data <- RunUMAP(data, dims = 1:opt$ndim)
  data <- FindNeighbors(data, dims = 1:opt$ndim)
}

data <- FindClusters(data, resolution = opt$resolution)
print("***** Clustering completed! *****")


# create the output directory
if (!(dir.exists(opt$output))) {
  dir.create(opt$output,recursive = TRUE)
}


# Run ScType --------------------------------------------------------------
source("/opt/sc-type/R/sctype_score_.R")

# modify the custom to sctype format
opt$custom_db <- generate_custom_sctype_db(opt$custom_db)

# run sctype 
# https://sctype.app/database.php
sctype_scores <- run_sc_type(opt$custom_db, 
                             opt$tissue, 
                             tolower(gsub(" ", "_", opt$tissue)))
data@meta.data$cell_type_CellMarker2.0 <- ""
for (j in unique(sctype_scores$cluster)) {
  cl_type <- sctype_scores[sctype_scores$cluster == j, ]
  data@meta.data$cell_type_CellMarker2.0[data@meta.data$seurat_clusters == j] <- as.character(cl_type$type[1])
}
# adding additional metadata of interest
data@meta.data$seurat_clusters_celltype_CellMarker2.0 <- paste0(data@meta.data$seurat_clusters, 
                                                                "_", 
                                                                data@meta.data$cell_type_CellMarker2.0)
# Extract the numeric values from the strings
numeric_values <- as.numeric(gsub("[^0-9.]", "", unique(data@meta.data$seurat_clusters_celltype_CellMarker2.0)))
# Sort the numeric values
sorted_vec <- unique(data@meta.data$seurat_clusters_celltype_CellMarker2.0)[order(numeric_values)]
# Convert the metadata column to a factor with levels set accordingly
data@meta.data$seurat_clusters_celltype_CellMarker2.0 <- factor(data@meta.data$seurat_clusters_celltype_CellMarker2.0, 
                                                                levels = sorted_vec)



# save the rds object
saveRDS(data,
        file = paste0(
          opt$output,
          "/",
          opt$output_prefix,
          "_clustered_and_cell_typed.rds"
        )
)
print("***** sctype complete and saved Seurat object! *****")


# UMAPs for all metadata variables ---------------------------------------------
meta_cols <- names(data@meta.data)[!sapply(data@meta.data, is.numeric)]
# Exclude "orig.ident"
meta_cols <- meta_cols[meta_cols != "orig.ident"]
# Exclude "DF.classifications_" columns
meta_cols <- meta_cols[!grepl("^DF\\.classifications_", meta_cols)]
for (meta in meta_cols) {
  # Skip if there is only one value
  if (length(unique(data@meta.data[[meta]])) == 1){
    next
  }
  print(paste0("***** Plotting ", meta, " *****"))
  pdf(file.path(opt$output, paste0(opt$output_prefix, "_", meta, "_UMAP.pdf")),
      width = 10,
      height = 7)
  print(DimPlot(data,
                raster = FALSE,
                order = TRUE,
                label = FALSE,
                group.by = meta,
                reduction = "umap"
  ))
  dev.off()
  w <- (min(length(unique(data@meta.data[[meta]])),4)) * 6
  h <- ifelse(length(unique(data@meta.data[[meta]])) <= 4,
              7,
              5 * ceiling(length(unique(data@meta.data[[meta]]))/4))
  pdf(file.path(opt$output, paste0(opt$output_prefix, "_", meta, "_split_UMAP.pdf")),
      width = w,
      height = h)
  print(DimPlot(data,
                raster = FALSE,
                order = TRUE,
                label = FALSE,
                group.by = meta,
                split.by = meta,
                reduction = "umap",
                ncol = 4
  ) + NoLegend())
  dev.off()
  
}



# custom UMAPs with cell type labels -------------------------------------------
pdf(file.path(opt$output, paste0(opt$output_prefix, "_seurat_clusters_labeled_UMAP.pdf")),
    width = 10,
    height = 7
)
print(DimPlot(data,
              group.by = "seurat_clusters",
              raster = FALSE,
              order = TRUE,
              label = TRUE,
              reduction = "umap"
))
dev.off()

pdf(file.path(opt$output, paste0(opt$output_prefix, "_ScType_celltype_CellMarker2.0_labeled_UMAP.pdf")),
    width = 10,
    height = 7
)
print(DimPlot(data,
              group.by = "seurat_clusters_celltype_CellMarker2.0",
              raster = FALSE,
              order = TRUE,
              label = TRUE,
              reduction = "umap"
))
dev.off()

pdf(file.path(opt$output, paste0(opt$output_prefix, 
                                 "_seurat_clusters_split_sample_UMAP.pdf")),
    width = 26,
    height = 20)
print(DimPlot(data,
              raster = FALSE,
              order = TRUE,
              label = FALSE,
              group.by = "seurat_clusters",
              split.by = "SampleID",
              reduction = "umap",
              ncol = 4
) + theme(legend.position = "right"))
dev.off()

pdf(file.path(opt$output, paste0(opt$output_prefix, 
                                 "_seurat_clusters_split_condition_UMAP.pdf")),
    width = 26,
    height = 7)
print(DimPlot(data,
              raster = FALSE,
              order = TRUE,
              label = TRUE,
              group.by = "seurat_clusters",
              split.by = "Condition",
              reduction = "umap",
              ncol = 4
) + theme(legend.position = "right"))
dev.off()


# Check sample bias ------------------------------------------------------------
# counts per sample in each sub cluster
column_names <- c("SampleID", "Condition", "seurat_clusters")
sample_cell_count_data <- data@meta.data %>%
  select(all_of(column_names)) %>%
  group_by(across(all_of(column_names)), .drop = FALSE) %>%
  summarise(number_of_cells_per_sample_in_cluster = n(), .groups = "drop") %>%
  complete(nesting(SampleID, Condition), seurat_clusters,
           fill = list(number_of_cells_per_sample_in_cluster = 0))
colnames(sample_cell_count_data)[1:3] <- c("sample_id", "animal_model", "cluster_id")
sample_cell_count_data <- sample_cell_count_data[order(sample_cell_count_data$cluster_id), ]

# use the sample total counts from the entire data and not just clustered data
sample_total_cells <- data@meta.data %>%
  group_by(SampleID) %>%
  summarise(total_cells_per_sample = n())

sample_cell_count_data <- sample_cell_count_data %>%
  merge(., sample_total_cells, by.x = "sample_id", by.y = "SampleID", all.x = TRUE)

write.csv(sample_cell_count_data, 
          file = file.path(opt$output,
                           paste0(opt$output_prefix,
                                  "_counts_per_sample_per_cluster.csv")),
          row.names = FALSE)
# Compute Proportions per Sample per Cluster
smp_cluster_counts_avg <- sample_cell_count_data %>%
  group_by(cluster_id) %>%
  mutate(total_numbers_of_cells_per_cluster = sum(number_of_cells_per_sample_in_cluster)) %>%
  group_by(sample_id, .add = TRUE) %>%
  mutate(fraction_of_cluster_cells_from_sample = number_of_cells_per_sample_in_cluster / total_numbers_of_cells_per_cluster) %>%
  ungroup() %>%
  mutate(percent_category = ifelse(fraction_of_cluster_cells_from_sample == 0, "zero", "nonzero"))


write.csv(smp_cluster_counts_avg, 
          file = file.path(opt$output,
                           paste0(opt$output_prefix,
                                  "_sample_cluster_proportion_tileplot_data.csv")),
          row.names = FALSE)

# Dynamically adjust plot width based on number of clusters
n_clusters <- length(unique(smp_cluster_counts_avg$cluster_id))
plot_width <- max(12, n_clusters * 0.6) 

# Plot: Tile Heatmap
pdf(file.path(opt$output, paste0(opt$output_prefix, "_sample_cluster_proportion_tileplot.pdf")), 
    width = plot_width, 
    height = 12)
print(
  ggplot(smp_cluster_counts_avg, aes(x = factor(cluster_id), y = sample_id)) +
    # Grey tiles for zero proportions
    geom_tile(data = subset(smp_cluster_counts_avg, percent_category == "zero"),
              fill = "grey80", color = "black", linewidth = 0.2) +
    # Gradient tiles for non-zero proportions
    geom_tile(data = subset(smp_cluster_counts_avg, percent_category == "nonzero"),
              aes(fill = fraction_of_cluster_cells_from_sample),
              color = "black", linewidth = 0.2) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0.5, limits = c(0, 1),
      name = "Proportion",
      guide = guide_colorbar(title.position = "top", barwidth = 10)
    ) +
    labs(
      title = "Sample Representation Across Clusters",
      subtitle = paste0(
        "Proportion = (# cells from sample in cluster) ÷ (total cells in cluster).\n",
        "Color interpretation: blue = low contribution, white = ~50%, red = dominant sample, ",
        "grey = 0 cells in cluster"
      ),
      x = "Seurat Cluster ID",
      y = "Sample ID"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(margin = margin(b = 10)),
      legend.position = "bottom"
    )
)
dev.off()




# Write session info
writeLines(
  capture.output(sessionInfo()),
  file.path(opt$output, "sessionInfo.txt")
)
print("***** Script Complete *****")




# END