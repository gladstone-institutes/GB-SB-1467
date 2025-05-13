#!/usr/bin/env Rscript

###############################################################################
## Project ID: GB-SB-1467
## Authors: Ayushi Agrawal
##
## Script Goal: Perform pseudo-bulked differential gene expression analysis
##
## Usage example:
## Rscript 07_pseudobulk_de.R
##  --input 'input_seurat_object.RDS' \   # Input Seurat object
##  --output '/output_directory' \        # Location for output files
##  --output_prefix 'Interneurons' \      # Prefix for output files
##  --predictor c("Group", "Condition") \ # Vector of predictor variables - for each one of them pairwise differential analyses are run
##  --cell_annotation "seurat_clusters" \ # cell annotation over which pseudo-bulk analyses are performed
##  --sampleid "PID" \                    # Meta data column to identify biological replicates
##  --confounder "BATCH"                  # Condfounding variable to include in the model matrix
##
## Run "07_pseudobulk_de.R --help" for more information
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
              help = "Output directory where results will be saved. 
              Will be created if it doesn't exist (required)"
  ),
  make_option("--output_prefix",
              action = "store", default = "pseudobulk_de", type = "character",
              help = "File name prefix for ouptut files (optional)"
  ),
  make_option(c("-p", "--predictor"),
              action = "store", default = NA, type = "character",
              help = "Comma-separated predictor variables for the design formula (required)"
  ),
  make_option("--predictor_ref_level",
              action = "store", default = NA, type = "character",
              help = "Reference level for the predictor (optional)"
  ),
  make_option(c("-c", "--cell_annotation"),
              action = "store", default = "seurat_clusters", type = "character",
              help = "Cell-type or cluster annotation column in metadata used for 
              pseudo-bulk aggregation (default = 'seurat_clusters')"
  ),
  make_option(c("-s", "--sampleid"),
              action = "store", default = NA, type = "character",
              help = "Metadata column name that identifies biological replicates"
  ),
  make_option("--confounder",
              action = "store", default = NA, type = "character",
              help = "Condfounding variable to include in the model matrix (optional)"
  ),
  make_option("--minCells",
              action = "store", default = 10, type = "numeric",
              help = "Keep clusters with a median of at least 'minCells' cells 
              across all individual replicates (default = 10)"
  ),
  make_option("--minGSSize",
              action = "store", default = 5, type = "numeric",
              help = "Min set size for pathway enrichment (optional)"
  ),
  make_option("--maxGSSize",
              action = "store", default = 500, type = "numeric",
              help = "Max set size for pathway enrichment (optional)"
  ),
  make_option("--pathway_db",
              action = "store", default = NA, type = "character",
              help = "Path to RData file that contains pathway DBs in the same 
              format used by Interactive Enrichment Analysis (required)"
  ),
  make_option("--p_val_cutoff",
              action = "store", default = 0.05, type = "numeric",
              help = "DE adjusted p value cutoff for inclusion in ORA (optional)"
  ),
  make_option("--fold_change_cutoff",
              action = "store", default = 0.25, type = "numeric",
              help = "DE log fold change cutoff for inclusion in ORA (optional)"
  ),
  make_option("--species",
              action = "store", default = "human", type = "character",
              help = "Species - supported options are 'human' or 'mouse',
              default is 'human'"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# Check if required args are provided
if (is.na(opt$input) | is.na(opt$output) | is.na(opt$predictor) | is.na(opt$sampleid) | 
    is.na(opt$pathway_db) | is.na(opt$species)) {
  stop("Missing one or more required arguments")
}


# Setup ------------------------------------------------------------------------
# Load required packages
library(clusterProfiler)
library(Seurat)
library(HGNChelper)
library(muscat)
library(SummarizedExperiment)
library(emmeans)
library(gridExtra)
library(EnhancedVolcano)
library(edgeR)
library(ggpubr)
library(grid)
library(openxlsx)
library(tidyverse)
library(magrittr)

# set seed
set.seed(1234)

# create output directory
if (!(dir.exists(opt$output))) {
  dir.create(opt$output, recursive = T, showWarnings = F)
}

# Set default colors to Wong B (https://doi.org/10.1038/nmeth.1618)
# colorblindness friendly colors
custom_colors <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
  "#666666", "#AD7700", "#1C91D4", "#007756", "#D5C711", "#005685", "#A04700",
  "#B14380", "#4D4D4D", "#FFBE2D", "#80C7EF", "#00F6B3", "#F4EB71", "#06A5FF",
  "#FF8320", "#D99BBD", "#8C8C8C", "#FFCB57", "#9AD2F2", "#2CFFC6", "#F6EF8E",
  "#38B7FF", "#FF9B4D", "#E0AFCA", "#A3A3A3", "#8A5F00", "#1674A9", "#005F45",
  "#AA9F0D", "#00446B", "#803800", "#8D3666", "#3D3D3D"
)
options(
  ggplot2.continuous.colour = list(custom_colors),
  ggplot2.continuous.fill = list(custom_colors),
  ggplot2.discrete.colour = list(custom_colors),
  ggplot2.discrete.fill = list(custom_colors)
)


# Functions --------------------------------------------------------------------
# refine the metadata levels in a Seurat object
refine_metadata_levels <- function(seurat_data) {
  for (i in base::colnames(seurat_data@meta.data)) {
    if (base::is.factor(seurat_data@meta.data[[i]])) {
      base::print(base::paste("Re-evaluating levels for a factor column", i))
      base::print(
        base::paste(
          "before:", base::paste(base::levels(seurat_data@meta.data[[i]]), collapse = ", ")
        )
      )
      seurat_data@meta.data[[i]] <- base::droplevels(seurat_data@meta.data[[i]]) # need to drop levels of the removed values
      base::print(
        base::paste(
          "after:", base::paste(base::levels(seurat_data@meta.data[[i]]), collapse = ", ")
        )
      )
    }
  }
  return(seurat_data)
}

# plots on multiple pages of pdf
multipage_plot <- function(plot_list, per_page, filename, ncol = 1) {
  # Create a PDF file to save the plots
  pdf(file = filename, height = 12, width = 10)
  
  # Split the plots into groups of per_page per page
  num_plots <- length(plot_list)
  num_pages <- ceiling(num_plots / per_page)
  plot_indices <- split(
    seq_len(num_plots),
    rep(seq_len(num_pages), each = per_page)
  ) |>
    suppressWarnings()
  
  # Plot the plots on each page
  for (i in seq_len(num_pages)) {
    if (length(plot_indices[[i]]) < per_page) {
      # Fill in the rest of the last page with blank plots
      this_index <- max(plot_indices[[i]])
      for (j in seq_len(per_page - length(plot_indices[[i]]))) {
        this_index <- this_index + 1
        plot_indices[[this_index]] <- ggplot() +
          geom_blank()
        plot_indices[[i]] <- c(plot_indices[[i]], this_index)
      }
    }
    plots <- ggarrange(
      plotlist = plot_list[plot_indices[[i]]],
      ncol = ncol, nrow = ceiling(length(plot_indices[[i]]) / ncol)
    )
    plots_grob <- ggplotGrob(plots)
    grid.newpage()
    grid.draw(plots_grob)
  }
  # Close the PDF file
  dev.off()
}

# PCA plot
pca_analysis <- function(clusterID, pb_obj, nfeatures_for_pca){
  
  # check if sample has all 0s
  temp_pb_obj <- pb_obj@assays@data[clusterID]
  sample_totals <- colSums(temp_pb_obj@listData[[clusterID]])
  if(!(length(names(sample_totals[sample_totals==0])) == 0)){
    toremove <- names(sample_totals[sample_totals==0])
    print(paste0(
      "PCA analyses on sample ",
      paste(toremove, collapse = ","),
      " not performed. The number of counts across all genes is 0 for this sample in cluster",
      clusterID
    ))
    pb_obj <- subset(pb_obj, ,!(colnames(pb_obj) %in% toremove))
  }
  
  # Calculate log-CPM values with log transformation and prior count stabilization
  log_cpm_counts <- data.frame(cpm(as.matrix(pb_obj@assays@data[[clusterID]]), 
                                   log = TRUE, prior.count = 0.01))
  
  # identify the top variable features to use for PCA
  variance_per_gene <- rowVars(as.matrix(log_cpm_counts))
  use_features_for_pca <- variance_per_gene %>% sort(.,decreasing = TRUE) %>% names(.)
  use_features_for_pca <- use_features_for_pca[1:nfeatures_for_pca]
  
  # Run PCA on the log-CPM data
  log_cpm_counts_t <- t(log_cpm_counts[use_features_for_pca,])
  res.pca.norm<- prcomp(log_cpm_counts_t, 
                        scale = TRUE, 
                        rank. = min(10, ncol(log_cpm_counts_t)))
  
  # Extract the proportion of variance explained
  variance_explained <- res.pca.norm$sdev^2 / sum(res.pca.norm$sdev^2)
  percent_variance_explained <- variance_explained * 100
  
  
  # Create a data frame for visualization
  log_cpm_counts_PCs=data.frame(Sample=rownames(log_cpm_counts_t),res.pca.norm$x )
  log_cpm_counts_PCs[,opt$predictor] <- "NA"
  match_indices <- match(log_cpm_counts_PCs$Sample, rownames(colData(pb_obj)))
  meta.data <- as.data.frame(pb_obj@colData[match_indices,])
  log_cpm_counts_PCs[,opt$predictor] <- meta.data$group_id
  
  #plotting PCs
  if (!is.na(opt$confounder)) {
    log_cpm_counts_PCs[,opt$confounder] <- "NA"
    log_cpm_counts_PCs[,opt$confounder] <- meta.data[,opt$confounder]
    
    return(ggplot(log_cpm_counts_PCs, aes(x = PC1, y = PC2, 
                                          shape = eval(parse(text = opt$confounder)), 
                                          color = eval(parse(text = opt$predictor)))) +
             geom_point() +
             geom_text_repel(aes(label = Sample), vjust = -0.6, hjust = 0.6, size = 3, 
                             show.legend = FALSE) +  # Add text labels
             guides(color = guide_legend(paste(opt$predictor, sep = ""))) +
             scale_shape(name = opt$confounder) +
             xlab(paste0("PC1 (", round(percent_variance_explained[1], 2), "% variance)")) +
             ylab(paste0("PC2 (", round(percent_variance_explained[2], 2), "% variance)")) +
             ggtitle(paste0(opt$output_prefix, " Cluster ", clusterID)) +
             theme_bw())
  }else{
    return(ggplot(log_cpm_counts_PCs,aes(x=PC1,y=PC2, color=eval(parse(text=opt$predictor))))+
             geom_point() +
             geom_text_repel(aes(label = Sample), vjust = -0.6, hjust = 0.6, size = 3, 
                             show.legend = FALSE) +  # Add text labels
             guides(color=guide_legend(paste(opt$predictor, sep=""))) +
             xlab(paste0("PC1 (", round(percent_variance_explained[1], 2), "% variance)")) +
             ylab(paste0("PC2 (", round(percent_variance_explained[2], 2), "% variance)")) +
             ggtitle(paste0(opt$output_prefix, " Cluster ", clusterID)) +
             theme_bw())
  }
}

clean_char <- function(input_string) {
  # Replace all special characters with nothing
  cleaned_string <- gsub("[^a-zA-Z0-9 ]", "", input_string)
  
  # Replace spaces with underscores
  cleaned_string <- gsub(" ", "_", cleaned_string)
  
  return(cleaned_string)
}

# Volcano plots
plot_volcano <- function(de,
                         cell_type,
                         comparison,
                         height = 10,
                         width = 15,
                         p_value_cutoff = 0.05,
                         fold_change_cutoff = 0.25,
                         output) {
  number_of_de_genes <- de %>%
    filter(p_adj.loc < p_value_cutoff & abs(logFC) > fold_change_cutoff) %>%
    nrow()
  
  filename <- paste0(clean_char(cell_type), "_", comparison)
  subtitle_expr <- paste0(comparison, 
                          "\np_adj.loc cutoff = ", p_value_cutoff, 
                          "; log2 fold change cutoff = ", fold_change_cutoff,
                          "\nNumber of DE genes = ", number_of_de_genes)
  
  p_value_cutoff_new <- min(de$p_val[de$p_adj.loc >= p_value_cutoff])
  
  plot <- EnhancedVolcano(de,
                          title = cell_type,
                          subtitle = subtitle_expr,
                          lab = de$gene,
                          pCutoff = p_value_cutoff_new,
                          FCcutoff = fold_change_cutoff,
                          x = 'logFC',
                          y = 'p_val',
                          ylab = bquote(~-Log[10]~italic(p_val)),
                          drawConnectors = TRUE,
                          widthConnectors = 0.3,
                          boxedLabels = TRUE,
                          labSize = 3.0,
                          colConnectors = "black",
                          legendLabSize = 10,
                          legendIconSize = 3.0,
                          legendLabels = c("NS", 
                                           bquote("Log"[2]~"FC > " ~ .(fold_change_cutoff)),
                                           bquote("p_adj.loc < " ~ .(p_value_cutoff)),
                                           bquote("p_adj.loc < " ~ .(p_value_cutoff) ~ " & Log"[2]~"FC > " ~ .(fold_change_cutoff)))
  )
  ggsave(
    plot = plot,
    filename = file.path(output, paste0(filename, "_volcano_plot.pdf")),
    width = width,
    height = height,
    units = "in"
  )
}

# Add missing design info to colData
add_design_with_colData <- function(design_matrix, col_data) {
  # Ensure that the row names of design matrices and colData match
  if (!all(rownames(design_matrix) == rownames(col_data))) {
    stop("Row names of design matrices and colData do not match")
  }
  
  print("This function only works for pairwise comparison designs (2 columns)")
  
  not_ref <- design_matrix |>
    as.data.frame() |>
    pull(2) |>
    as.logical()
  
  col_data <- col_data |>
    mutate(
      group = !not_ref,
      group = ifelse(group,
                     colnames(design_matrix)[1],
                     colnames(design_matrix)[2]
      )
    )
  
  return(col_data)
}

# Box plot for genes with lowest p-values
plot_boxplots_lowest_50_p_val <- function(result, iter, output, filename) {
  if (!(dir.exists(output))) {
    dir.create(output, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Extract top 50 lowest p-value genes
  lowest_50_pval <- result$table[[1]][[iter]] |>
    as.data.frame() |>
    arrange(p_val) |>
    head(50) |>
    pull(gene)
  
  # Add group info
  group_info <- add_design_with_colData(
    design_matrix = result$data[[iter]]$design,
    col_data = result$data[[iter]]$samples
  ) |>
    dplyr::select(group) |>
    tibble::rownames_to_column(var = "Subject")
  
  # Prepare CPM data
  plot_df <- result$data[[iter]]$counts |>
    cpm(log = TRUE) |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "Subject") |>
    dplyr::select(Subject, all_of(lowest_50_pval)) |>
    full_join(group_info)
  
  # Add consistent jitter
  set.seed(123)
  jitter_width <- 0.25
  plot_df <- plot_df |>
    mutate(
      x_jitter = as.numeric(as.factor(group)) + runif(n(), -jitter_width, jitter_width)
    )
  
  # Plot top 50 genes one by one
  plots <- lowest_50_pval |>
    purrr::map(.f = function(gene) {
      y_expr <- bquote("Log"[2] ~ "CPM")
      ggplot(plot_df, aes(x = group, y = .data[[gene]], fill = group)) +
        geom_violin(alpha = 0.5, color = "black") +
        geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
        geom_point(aes(x = x_jitter), shape = 21, size = 4, color = "black") +
        geom_text_repel(
          aes(x = x_jitter, label = Subject),
          size = 3,
          segment.color = "gray30",
          segment.size = 0.3,
          min.segment.length = 0.1,
          box.padding = 0.5,
          point.padding = 0.25,
          max.overlaps = Inf,
          show.legend = FALSE
        ) +
        ggtitle(gene) +
        ylab(label = y_expr) +
        xlab(label = opt$predictor) +
        theme_bw()
    })
  
  # Save plots in multipage PDF
  multipage_plot(
    plot_list = plots,
    per_page = 3,
    filename = file.path(output, filename)
  )
}

# Convert ids
map_ids <- function(input, db, fromType) {
  toType.list <- c("ENTREZID", "SYMBOL")
  if (fromType == "ENTREZID") {
    toType.list <- c("SYMBOL")
  }
  
  # perform ID mapping
  output <- bitr(input$gene,
                 fromType = fromType,
                 toType = toType.list,
                 OrgDb = db
  )
  
  join_cols <- c("gene")
  names(join_cols) <- fromType
  output <- dplyr::left_join(output, input, by = join_cols)
  return(output)
}

# Function to extract the desired strings based on patterns
extract_category <- function(x) {
  if (grepl("^go_", x, ignore.case = TRUE)) {
    return("GO")
  } else if (grepl("^pfocr_", x, ignore.case = TRUE)) {
    return("PFOCR")
  } else if (grepl("^wp_", x, ignore.case = TRUE)) {
    return("WikiPathways")
  } else if (grepl("^KEGG", x, ignore.case = TRUE)) {
    return("KEGG")
  }
}

# Run ORA
run_ora <- function(gene_list,
                    output,
                    output_prefix,
                    fromType
) {
  
  # Identifier mapping
  ranked_genes_entrez <- map_ids(gene_list, org.db.name, fromType)
  
  gene_list_entrez <- ranked_genes_entrez %>%
    filter(p.value < opt$p_val_cutoff & abs(fold.change) > opt$fold_change_cutoff) %>%
    pull(ENTREZID)
  
  if ((length(gene_list) == 0) | (length(gene_list_entrez) == 0)) {
    return()
  }
  
  universe <- map_ids(data.frame(gene = background_genes), org.db.name, fromType)
  
  # Backup existing ggplot2 options
  old_options <- list(
    ggplot2.continuous.colour = getOption("ggplot2.continuous.colour"),
    ggplot2.continuous.fill = getOption("ggplot2.continuous.fill"),
    ggplot2.discrete.colour = getOption("ggplot2.discrete.colour"),
    ggplot2.discrete.fill = getOption("ggplot2.discrete.fill")
  )
  
  # Temporarily unset ggplot2 options for dotplot
  options(
    ggplot2.continuous.colour = NULL,
    ggplot2.continuous.fill = NULL,
    ggplot2.discrete.colour = NULL,
    ggplot2.discrete.fill = NULL
  )
  
  for (this_db in db) {
    # Object from string
    database <- eval(parse(text = this_db))
    db_name <- paste0("ORA_",extract_category(this_db))
    
    enrichResult <- enricher(
      gene_list_entrez,
      universe = universe$ENTREZID,
      TERM2GENE = database[, c("term", "gene")],
      TERM2NAME = database[, c("term", "name")],
      minGSSize = opt$minGSSize,
      maxGSSize = opt$maxGSSize,
      # pAdjustMethod="holm", #default is "BH"
      pvalueCutoff = 1 # to limit results
    )
    
    if (!is.null(enrichResult)) {
      enrichResult <- setReadable(enrichResult, org.db.name, keyType = "ENTREZID")
      
      if (nrow(enrichResult) > 0) {
        dir.create(file.path(output, this_db), showWarnings = F, recursive = T)
        
        write_csv(as.data.frame(enrichResult@result), 
                  file.path(output, this_db, 
                            paste0(output_prefix,
                                   "_",
                                   db_name, 
                                   ".csv")))
        
        #trim descriptions to 80 characters
        enrichResult@result <- dplyr::mutate(enrichResult@result, 
                                             Description = str_trunc(Description, 80))
        
        p <- barplot(enrichResult, 
                     showCategory = 20,
                     label_format=50) +
          labs(
            title = output_prefix,
            subtitle = paste0("\nORA Terms for ", this_db)
          )
        ggplot2::ggsave(p, file = file.path(output, 
                                            this_db,
                                            paste0(output_prefix,
                                                   "_",
                                                   db_name, 
                                                   "_enrichment_barplot",
                                                   ".pdf")), 
                        width = 2400, height = 2600, units = "px", device='pdf')
        
        p <- enrichplot::dotplot(enrichResult, 
                                 showCategory = 20,
                                 label_format=50) +
          labs(
            title = output_prefix,
            subtitle = paste0("\nORA Terms for ", this_db)
          )
        ggplot2::ggsave(p, file = file.path(output,  
                                            this_db,
                                            paste0(output_prefix,
                                                   "_",
                                                   db_name, 
                                                   "_enrichment_dotplot",
                                                   ".pdf")), 
                        width = 2400, height = 2600, units = "px", device='pdf')
        
      }
    }
  }
  
  #-----------------------------------
  #KEGG analysis
  #-----------------------------------
  print("Starting KEGG...")
  ekeggbp <- enrichKEGG(
    gene     = gene_list_entrez %>% subset(., !is.na(.)),
    universe = universe$ENTREZID %>% subset(., !is.na(.)),
    organism    = "mmu",
    minGSSize = 10,
    pvalueCutoff = 0.8,
    keyType = "ncbi-geneid"
  )
  
  #translating gene IDs to human readable symbols
  if (!is.null(ekeggbp)) {
    ekeggbp <- setReadable(ekeggbp, OrgDb = org.db.name, keyType="ENTREZID")
    if (nrow(ekeggbp) > 0) {
      #Visualize
      dir.create(file.path(output, "KEGG"), showWarnings = F, recursive = T)
      
      #save the list of enriched pathways
      write.csv(ekeggbp,
                file = file.path(output, "KEGG", 
                                 paste0( output_prefix,
                                         "_ORA_KEGG_results_table.csv")))
      
      ## save images
      #trim descriptions to 80 characters
      ekeggbp@result <- dplyr::mutate(ekeggbp@result, 
                                      Description = str_trunc(Description, 80))
      
      p <- barplot(ekeggbp, 
                   showCategory = 20,
                   label_format=50) +
        labs(
          title = output_prefix,
          subtitle = "\nORA Terms for KEGG"
        )
      ggplot2::ggsave(p, file = file.path(output, 
                                          "KEGG",
                                          paste0(output_prefix,
                                                 "_ORA_KEGG_enrichment_barplot.pdf")), 
                      width = 2400, height = 2600, units = "px", device='pdf')
      
      p <- enrichplot::dotplot(ekeggbp, 
                               showCategory = 20,
                               label_format=50) +
        labs(
          title = output_prefix,
          subtitle = "\nORA Terms for KEGG"
        )
      ggplot2::ggsave(p, 
                      file = file.path(output, 
                                       "KEGG",
                                       paste0(output_prefix,
                                              "_ORA_KEGG_enrichment_dotplot.pdf")), 
                      width = 2400, height = 2600, units = "px", device='pdf')
      
    }
  }
  
  # Restore the original ggplot2 options
  options(old_options)
  
}

# Run GSEA
run_gsea <- function(ranked_genes, 
                     output, 
                     fromType,
                     output_prefix, 
                     seed = 1) {
  set.seed(seed)
  minGSSize <- opt$minGSSize
  maxGSSize <- opt$maxGSSize
  
  # Identifier mapping
  ranked_genes_entrez <- map_ids(ranked_genes, org.db.name, fromType)
  
  # Record unmapped rows
  ranked_genes_unmapped <- ranked_genes %>%
    dplyr::filter(!gene %in% ranked_genes_entrez[[fromType]]) %>%
    dplyr::arrange(desc(rank))
  
  # Resolve duplicates (keep ENTREZID with largest abs(rank))
  ranked_genes_entrez_dedup <- ranked_genes_entrez %>%
    dplyr::mutate(absrank = abs(rank)) %>%
    dplyr::arrange(desc(absrank)) %>%
    dplyr::distinct(ENTREZID, .keep_all = T) %>%
    dplyr::select(-absrank) %>%
    dplyr::arrange(desc(rank))
  
  ranked_genes_entrez_l <- ranked_genes_entrez_dedup$rank
  names(ranked_genes_entrez_l) <- ranked_genes_entrez_dedup$ENTREZID
  ranked_genes_entrez_l <- na.omit(ranked_genes_entrez_l)
  
  geneList <- sort(ranked_genes_entrez_l, decreasing = T)
  
  # Backup existing ggplot2 options
  old_options <- list(
    ggplot2.continuous.colour = getOption("ggplot2.continuous.colour"),
    ggplot2.continuous.fill = getOption("ggplot2.continuous.fill"),
    ggplot2.discrete.colour = getOption("ggplot2.discrete.colour"),
    ggplot2.discrete.fill = getOption("ggplot2.discrete.fill")
  )
  
  # Temporarily unset ggplot2 options for dotplot
  options(
    ggplot2.continuous.colour = NULL,
    ggplot2.continuous.fill = NULL,
    ggplot2.discrete.colour = NULL,
    ggplot2.discrete.fill = NULL
  )
  
  for (this_db in db) {
    print(this_db)
    # Object from string
    database <- eval(parse(text = this_db))
    
    # Perform GSEA
    gseaResult <- clusterProfiler::GSEA(
      geneList,
      TERM2GENE = database[, c("term", "gene")],
      TERM2NAME = database[, c("term", "name")],
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      # pAdjustMethod="holm", #default is "BH"
      pvalueCutoff = 0.99, # to limit results
      verbose = F
    )
    
    
    if (!is.null(gseaResult) && nrow(gseaResult) > 0) {
      gseaResult <- setReadable(gseaResult,
                                OrgDb = org.db.name,
                                keyType = "ENTREZID"
      )
      
      dir.create(file.path(output, this_db), showWarnings = F, recursive = T)
      
      write_csv(as.data.frame(gseaResult@result), 
                file.path(output, 
                          this_db, 
                          paste0(output_prefix,
                                 "_GSEA_",
                                 extract_category(this_db), 
                                 "_results_table.csv")))
      
      #trim descriptions to 80 characters
      gseaResult@result <- dplyr::mutate(gseaResult@result, 
                                         Description = str_trunc(Description, 80))
      
      p <- enrichplot::dotplot(gseaResult, 
                               showCategory = 20,
                               label_format=50) +
        labs(
          title = output_prefix,
          subtitle = paste0("\nORA Terms for ", this_db)
        )
      ggplot2::ggsave(p, file = file.path(output,  
                                          this_db,
                                          paste0(output_prefix,
                                                 "_GSEA_",
                                                 extract_category(this_db), 
                                                 "_enrichment_dotplot.pdf")), 
                      width = 2400, height = 2600, units = "px", device='pdf')
    }
  }
  
  #-----------------------------------
  #KEGG analysis
  #------------------------------
  ekeggbp <- gseKEGG(geneList = geneList,
                     organism = 'mmu',
                     minGSSize = minGSSize,
                     maxGSSize = maxGSSize,
                     pvalueCutoff = 0.99,
                     verbose = FALSE)
  if (!is.null(ekeggbp)) {
    gseaResult <- setReadable(ekeggbp,
                              OrgDb = org.db.name,
                              keyType = "ENTREZID"
    )
  }
  
  if (nrow(gseaResult) > 0) {
    #Visualize
    dir.create(file.path(output, "KEGG"), showWarnings = F, recursive = T)
    
    write_csv(gseaResult@result, 
              file.path(output, 
                        "KEGG", 
                        paste0(output_prefix,
                               "_GSEA_KEGG_results_table.csv")))
    ## save images
    #trim descriptions to 80 characters
    gseaResult@result <- dplyr::mutate(gseaResult@result, 
                                       Description = str_trunc(Description, 80))
    
    p <- enrichplot::dotplot(gseaResult, 
                             showCategory = 20,
                             label_format=50) +
      labs(
        title = output_prefix,
        subtitle = "\nORA Terms for KEGG"
      )
    ggplot2::ggsave(p, 
                    file = file.path(output, 
                                     "KEGG",
                                     paste0(output_prefix,
                                            "_GSEA_KEGG_enrichment_dotplot.pdf")), 
                    width = 2400, height = 2600, units = "px", device='pdf')
    
  }
  
  # Restore the original ggplot2 options
  options(old_options)
  
}


# Process inputs ---------------------------------------------------------------
if(opt$species == "human"){
  org.db.name <- "org.Hs.eg.db"
} else if(opt$species == "mouse"){
  org.db.name <- "org.Mm.eg.db"
} else {
  stop("Species must be either 'human' or 'mouse'")
}

# load pathway dbs
db <- load(opt$pathway_db)

# Read in the Seurat object
data <- readRDS(opt$input)
data <- refine_metadata_levels(data) # Fix the levels

# Predictors
predictors <- opt$predictor
predictors <- strsplit(predictors, ",")[[1]] %>% trimws()

# Sample identifier in metadata
sampleid <- opt$sampleid %>% trimws()

# Cluster annotation column
clusters <- opt$cell_annotation

# Store the meta-data for each cell in the PhenoData object
PhenoData <- data@meta.data
PhenoData <- PhenoData[,!grepl("^DF\\.classifications_", colnames(PhenoData))]
PhenoData <- PhenoData[,!grepl("^pANN_", colnames(PhenoData))]

# Check if predictors are in colnames of PhenoData
if (sum(predictors %in% colnames(PhenoData)) == 0) {
  stop("Columns of the meta.data slot of the Seurat object should include names of the predictors")
}

# Check if sampleid is in colnames of PhenoData
if (!(sampleid %in% colnames(PhenoData))) {
  stop("Columns of the meta.data slot of the Seurat object should include sampleid")
}

# get the cluster levels
cluster_levels <- PhenoData[[clusters]] %>%
  unique() %>%
  sort()



# create output directories ----------------------------------------------------
# Define the output directories in a named list
dirs <- list(
  de_outputs = file.path(opt$output, "gene_expression_associations"),
  ora_outputs = file.path(opt$output, "ORA"),
  gsea_outputs = file.path(opt$output, "GSEA")
)

# Create the directories recursively
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# Assign the paths to variables in the global environment
invisible(list2env(dirs, envir = .GlobalEnv))


# Pseudo bulk for DE analysis --------------------------------------------------
# find the list of background genes for ORA analysis later
data <- PrepSCTFindMarkers(data)
all_data <- GetAssayData(data, assay = "SCT")
background_genes <- rownames(all_data[rowSums(all_data[])>0,])
write.csv(as.data.frame(background_genes),
          file = file.path(opt$output, 
                           "ORA",
                           paste0(opt$output_prefix, "_nonzero_background_gene_list.csv")),
          row.names = FALSE)

# Create SingleCellExperiment object
sce <- SummarizedExperiment(assays = list(counts = data[["RNA"]]$counts), 
                            colData = PhenoData)

sce <- as(sce, "SingleCellExperiment")
rm(data)

# Prep this object for subsequent aggregation analyses
sce <- prepSCE(sce,
               kid = clusters, # sub-population assignments
               gid = predictors[1], # group IDs
               sid = sampleid, # sample IDs
               drop = FALSE
)

# keep those clusters with a median of at least 10 cells across all individual replicates
toKeep <- table(sce$cluster_id, sce$sample_id) %>%
  t() %>%
  apply(., 2, function(x) median(x)) %>%
  subset(., is_greater_than(., opt$minCells)) %>%
  names()

if (length(toKeep) < length(cluster_levels)) {
  print("******************************************")
  print(paste0(
    "Pseudo-bulk differential expression analyses on clusters ",
    paste(setdiff(cluster_levels, toKeep), collapse = ","),
    " not performed. The median number of cells across all biological replicates is less than 10 for these clusters."
  ))
  print("******************************************")
  
  sce <- subset(sce, ,cluster_id %in% toKeep)
}

# Aggregate counts across cells for each mouse (sample_id) within each cluster (cluster_id)
pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))

# Visualize the results in a Multi-Dimensional Scaling (MDS) plot
pb_mds <- pbMDS(pb) 

pb_mds %>% ggsave(
  plot = .,
  filename = file.path(opt$output, 
                       paste0(opt$output_prefix,
                              "_pseudo_bulked_MDS_plots.svg")),
  width = 10,
  height = 10
)

# Perform DE analysis 
for (p in 1:length(predictors)) {
  if (p == 1) {
    predictor <- "group_id"
  } else {
    predictor <- predictors[p]
  }
  print(predictor)
  
  
  ## generate PCA plots
  pca_clusters <- names(pb@assays@data)
  pca_plot_list <- lapply( pca_clusters,
                           pca_analysis,
                           pb,
                           500)
  multipage_plot(plot_list = pca_plot_list,
                 per_page = 2,
                 ncol = 1,
                 filename = file.path(opt$output, 
                                      paste0(opt$output_prefix,
                                             "_PCA_within_clusters.pdf")))
  
  
  nlevels_predictor <- pb[[predictor]] %>%
    unique() %>%
    length()
  levels_predictor <- pb[[predictor]] %>%
    unique()
  n_comparisons_per_cluster <- choose(nlevels_predictor, 2)
  
  # list all pairwise comparisons
  comparisons <- matrix(
    c("WT - Vitamin B3", "WT + Vitamin B3", 
      "WT - Vitamin B3", "KO - Vitamin B3",
      "KO + Vitamin B3", "KO - Vitamin B3"
    ),
    nrow = 2
  )
  
  for (comp in 1:ncol(comparisons)) {
    this_comparison <- gsub(" ", "", paste0(comparisons[1, comp], "_vs_", comparisons[2, comp]))
    print(paste0("Processing comaprison ", this_comparison, "......"))
    temp_pb <- pb[, (pb[[predictor]] == comparisons[1, comp] |
                       pb[[predictor]] == comparisons[2, comp])]
    if (is.factor(temp_pb[[predictor]])) {
      temp_pb[[predictor]] <- droplevels(temp_pb[[predictor]])
    }
    
    ## Set up the design matrix.
    group_id <- colData(temp_pb)[[predictor]]
    group_id <- relevel(as.factor(group_id), ref = as.character(comparisons[2, comp]))
    print(levels(group_id))
    if (!is.na(opt$confounder)) {
      confounder <- colData(temp_pb)[[opt$confounder]] %>%
        as.factor()
      confounder <- relevel(confounder,
                            ref = as.character(sort(unique(confounder), decreasing = TRUE)[1]))
      mm <- model.matrix(~group_id+confounder)
      row.names(mm) <- row.names(colData(temp_pb))
      colnames(mm) <- c(levels(group_id),levels(confounder)[2:length(levels(confounder))])
      print("Using model matrix:")
      print(as.data.frame(mm))
      print("Running pseudobulk model")
      res <- pbDS(temp_pb,
                  coef = 2,
                  design = mm, method = "edgeR", verbose = TRUE, min_cells = 3, filter = "genes"
      )
    } else {
      mm <- model.matrix(~group_id)
      colnames(mm) <- levels(group_id)
      row.names(mm) <- row.names(colData(temp_pb))
      print("Using model matrix:")
      print(as.data.frame(mm))
      print("Running pseudobulk model")
      res <- pbDS(temp_pb,
                  coef = 2,
                  design = mm, method = "edgeR", verbose = TRUE, min_cells = 3, filter = "genes"
      )
    }
    
    n_results <- length(res$table[[1]])
    result_names <- names(res$table[[1]])
    
    for (r in 1:n_results) {
      print(paste0("Processing psuedobulk results for cluster ", result_names[r], "......"))
      this_output <- file.path(opt$output, 
                               "gene_expression_associations",
                               this_comparison, 
                               paste0("cluster_", result_names[r]))
      dir.create(this_output,
                 recursive = T)
      res$table[[1]][[r]] %>%
        as.data.frame() %>%
        arrange(p_val,-logFC) %>%
        write.csv(., file.path(this_output,
                               paste0("cluster_", 
                                      result_names[r], 
                                      "_",
                                      this_comparison, ".csv")), 
                  row.names = FALSE)
      
      plot_volcano(res$table[[1]][[r]],
                   cell_type = paste0("Cluster ", result_names[r]),
                   comparison = this_comparison,
                   p_value_cutoff = opt$p_val_cutoff,
                   fold_change_cutoff = opt$fold_change_cutoff,
                   output = this_output
      )
      
      # Plot lowest p-value box plots
      plot_boxplots_lowest_50_p_val(res, 
                                    iter = r, 
                                    output = this_output,
                                    filename = paste0("Cluster_", 
                                                      result_names[r], 
                                                      "_",
                                                      this_comparison, 
                                                      "_lowest_50_p_val_genes_boxplots.pdf")
      )
      
      # check if there are any DEGs
      number_of_de_genes <- res$table[[1]][[r]] %>%
        as.data.frame() %>%
        filter(p_adj.loc < opt$p_val_cutoff & abs(logFC) > opt$fold_change_cutoff) %>%
        nrow()
      if(number_of_de_genes > 0){
        ora_gene_list <- res$table[[1]][[r]] %>%
          as.data.frame() %>%
          dplyr::select(gene, logFC, p_adj.loc) %>%
          dplyr::rename(fold.change = logFC, p.value = p_adj.loc)
        
        # ORA for ranked genes
        run_ora(
          gene_list = ora_gene_list,
          output = file.path(opt$output, 
                             "ORA", 
                             this_comparison, 
                             paste0("cluster_", result_names[r])),
          output_prefix = paste0("cluster_", 
                                 result_names[r],
                                 "_",
                                 this_comparison),
          fromType = "SYMBOL"
        )
        
        gsea_gene_list <- res$table[[1]][[r]] %>%
          as.data.frame() %>%
          dplyr::select(gene, logFC, p_val) %>%
          mutate(rank = sign(logFC) * -log10(p_val)) %>%
          arrange(rank) %>%
          dplyr::rename(fold.change = logFC, p.value = p_val)
        
        run_gsea(
          ranked_genes = gsea_gene_list,
          output = file.path(opt$output, 
                             "GSEA", 
                             this_comparison, 
                             paste0("cluster_", result_names[r])),
          fromType = "SYMBOL",
          output_prefix = paste0("cluster_", 
                                 result_names[r],
                                 "_",
                                 this_comparison)
        )
      } # run ORA/GSEA if DE genes
      
    } # iterate over each result
    
  } # loop over all comparisons of interest
  
} # run DE analysis for each predictor




# save the session info
writeLines(
  capture.output(sessionInfo()),
  file.path(opt$output, "sessionInfo.txt")
)


# END
