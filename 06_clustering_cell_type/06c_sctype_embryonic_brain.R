#!/usr/bin/env Rscript

###############################################################################
## Project ID: GB-SB-1467
## Authors: Reuben Thomas, Natalie Elphick, Ayushi Agrawal
##
## Script Goal: Normalize using SCTransform, use input PCA dim and resolution
## to cluster cells and use ScType to assign cell types to the clusters.
##
## Usage example:
## Rscript 06c_sctype_embryonic_brain.R
##  --input 'input_seurat_object.RDS' \  # Input Seurat object 
##  --output '/output_directory' \       # Location for output files
##  --output_prefix "outputs_scRNA_" \   # Prefix for output files
##  --tissue "Brain" \                   # Tissue type to supply to ScType
##  --ndim 20 \                          # Number of PCs to use
##  --resolution 0.02 \                  # Resolution to use for Seurat clustering
##  --custom_db "custom_db.xlsx"         # Path to custom database for ScType
##
## Run "Rscript 06c_sctype_embryonic_brain.R --help" for more information
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
  make_option("--custom_db",
              action = "store", default = NA, type = "character",
              help = "Path to custom database for ScType, must be a xlsx in the same format as the original (optional)"
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


# Run ScType --------------------------------------------------------------
# Read in the Seurat object
data <- readRDS(opt$input)

source("/opt/sc-type/R/sctype_score_.R")

# modify the custom to sctype format
opt$custom_db <- generate_custom_sctype_db(opt$custom_db)

# run sctype 
# https://sctype.app/database.php
sctype_scores <- run_sc_type(opt$custom_db, 
                             opt$tissue, 
                             tolower(gsub(" ", "_", opt$tissue)))
data@meta.data$cell_type_CellMarker2.0_embryonic_brain <- ""
for (j in unique(sctype_scores$cluster)) {
  cl_type <- sctype_scores[sctype_scores$cluster == j, ]
  data@meta.data$cell_type_CellMarker2.0_embryonic_brain[data@meta.data$seurat_clusters == j] <- as.character(cl_type$type[1])
}
# adding additional metadata of interest
data@meta.data$seurat_clusters_celltype_CellMarker2.0_embryonic_brain <- paste0(data@meta.data$seurat_clusters, 
                                                                                "_", 
                                                                                data@meta.data$cell_type_CellMarker2.0_embryonic_brain)
# Extract the numeric values from the strings
numeric_values <- as.numeric(gsub("[^0-9.]", "", unique(data@meta.data$seurat_clusters_celltype_CellMarker2.0_embryonic_brain)))
# Sort the numeric values
sorted_vec <- unique(data@meta.data$seurat_clusters_celltype_CellMarker2.0_embryonic_brain)[order(numeric_values)]
# Convert the metadata column to a factor with levels set accordingly
data@meta.data$seurat_clusters_celltype_CellMarker2.0_embryonic_brain <- factor(data@meta.data$seurat_clusters_celltype_CellMarker2.0_embryonic_brain, 
                                                                                levels = sorted_vec)



# save the rds object
saveRDS(data,
        file = paste0(
          opt$output,
          "/",
          opt$output_prefix,
          "_clustered_and_cell_typed_CellMarker2.0_embryonic_brain.rds"
        )
)
print("***** sctype complete and saved Seurat object! *****")


# custom UMAPs with cell type labels -------------------------------------------
pdf(file.path(opt$output, paste0(opt$output_prefix, "_ScType_celltype_CellMarker2.0_embryonic_brain_labeled_UMAP.pdf")),
    width = 10,
    height = 7
)
print(DimPlot(data,
              group.by = "seurat_clusters_celltype_CellMarker2.0_embryonic_brain",
              raster = FALSE,
              order = TRUE,
              label = TRUE,
              reduction = "umap"
))
dev.off()





# Write session info
writeLines(
  capture.output(sessionInfo()),
  file.path(opt$output, "sessionInfo.txt")
)
print("***** Script Complete *****")




# END