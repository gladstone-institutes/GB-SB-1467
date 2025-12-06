#Quality control metrics per cluster
#supplementary figure

#run locally on Ayushi's laptop

# Setup ------------------------------------------------------------------------
#load required packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)

#set all directory paths
basedir <- "~/Dropbox (Gladstone)/GB-SB-1467/"
outdir <- file.path(basedir, "paper/figures/QC_suppFig/")
setwd(basedir)

# set seed
set.seed(1234)

# create output directory
if (!(dir.exists(outdir))) {
  dir.create(outdir, recursive = T, showWarnings = F)
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


# Load data and prepare metadata -----------------------------------------------
#load the Seurat object
dat <- readRDS("gb_sb_1467_30PC_0.08res_clustered_and_cell_typed.rds")

##the number of cells in each mouse in each cluster
all_metadata <- dat[[]]
all_metadata$MouseNumber_Condition <- paste0(all_metadata$SampleID, 
                                             "_", 
                                             gsub(" ", "", all_metadata$Condition))


# Calculate summary statistics -------------------------------------------------
# Calculate summary statistics by cluster
avg_qc_features_per_cluster <- all_metadata %>% group_by(seurat_clusters) %>% 
  summarise(total_cells = n(),
            avg_ngene = mean(nFeature_RNA),
            sd_ngene = sd(nFeature_RNA),
            se_ngene = sd(nFeature_RNA) / sqrt(length(nFeature_RNA)),
            avg_nUMI = mean(nCount_RNA),
            sd_nUMI = sd(nCount_RNA),
            se_nUMI = sd(nCount_RNA) / sqrt(length(nCount_RNA)),
            avg_mt =  mean(percent.mt),
            sd_mt =  sd(percent.mt),
            se_mt = sd(percent.mt) / sqrt(length(percent.mt))) %>% 
  as.data.frame()

# Calculate summary statistics by mouse
avg_qc_features_per_mouse <- all_metadata %>% group_by(SampleID, Condition) %>% 
  summarise(total_cells = n(),
            avg_ngene = mean(nFeature_RNA),
            sd_ngene = sd(nFeature_RNA),
            se_ngene = sd(nFeature_RNA) / sqrt(length(nFeature_RNA)),
            avg_nUMI = mean(nCount_RNA),
            sd_nUMI = sd(nCount_RNA),
            se_nUMI = sd(nCount_RNA) / sqrt(length(nCount_RNA)),
            avg_mt =  mean(percent.mt),
            sd_mt =  sd(percent.mt),
            se_mt = sd(percent.mt) / sqrt(length(percent.mt)),
            .groups = "drop") %>% 
  as.data.frame()
avg_qc_features_per_mouse$SampleID <- factor(avg_qc_features_per_mouse$SampleID, 
                                              levels = sort(unique(avg_qc_features_per_mouse$SampleID)))

# Export summary statistics as supplementary tables
write.csv(avg_qc_features_per_cluster, 
          file = file.path(outdir, "supp_table_qc_by_cluster.csv"),
          row.names = FALSE)
write.csv(avg_qc_features_per_mouse, 
          file = file.path(outdir, "supp_table_qc_by_mouse.csv"),
          row.names = FALSE)


# Create QC figures ------------------------------------------------------------
#QC fig a - Cell counts per cluster
pa_cluster_cells <- ggplot(avg_qc_features_per_cluster, 
                           aes(x=seurat_clusters, y=total_cells, fill=seurat_clusters)) +
  geom_bar(position=position_dodge(), stat="identity") + 
  scale_x_discrete(expand = c(0.04,0.04))+
  theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA, linewidth=2),
    text = element_text(size = 12),
    legend.position = "none") + 
  ggtitle("Number of Cells per Cluster") +
  xlab("Cell Cluster") + 
  ylab("Number of Cells")


#QC fig b - nFeature_RNA per cluster
pb_cluster_ngene <- ggplot(all_metadata, 
                           aes(x=seurat_clusters, y=nFeature_RNA, fill=seurat_clusters)) + 
  stat_boxplot(geom='errorbar') +
  geom_boxplot(outlier.color = NA) + 
  geom_point(alpha = 0.2, size = 0.2, position = position_jitter(width = 0.1, seed = 1234)) +
  scale_x_discrete(expand = c(0.04,0.04))+
  theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA, linewidth=2),
    text = element_text(size = 12),
    legend.position = "none") + 
  ggtitle("Number of Genes per Cell by Cluster") +
  xlab("Cell Cluster") + 
  ylab("Number of Genes per Cell") 


#QC fig c - nCount_RNA per cluster
pc_cluster_numi <- ggplot(all_metadata, 
                          aes(x=seurat_clusters, y=nCount_RNA, fill=seurat_clusters)) + 
  stat_boxplot(geom='errorbar') +
  geom_boxplot(outlier.color = NA) + 
  geom_point(alpha = 0.2, size = 0.2, position = position_jitter(width = 0.1, seed = 1234)) +
  scale_x_discrete(expand = c(0.04,0.04))+
  theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA, linewidth=2),
    text = element_text(size = 12),
    legend.position = "none") + 
  ggtitle("nUMI per Cell by Cluster") +
  xlab("Cell Cluster") + 
  ylab("nUMI per Cell") 


#QC fig d - percent.mt per cluster
pd_cluster_mt <- ggplot(all_metadata, 
                        aes(x=seurat_clusters, y=percent.mt, fill=seurat_clusters)) + 
  stat_boxplot(geom='errorbar') +
  geom_boxplot(outlier.color = NA) + 
  geom_point(alpha = 0.2, size = 0.2, position = position_jitter(width = 0.1, seed = 1234)) +
  scale_x_discrete(expand = c(0.04,0.04))+
  theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA, linewidth=2),
    text = element_text(size = 12),
    legend.position = "none") + 
  ggtitle("% Mitochondrial Genes per Cells by Cluster") +
  xlab("Cell Cluster") + 
  ylab("% Mito Genes per Cell") 


#QC fig e - nFeature_RNA per mouse
pe_mouse_ngene <- ggplot(all_metadata, 
                         aes(x=MouseNumber_Condition, y=nFeature_RNA, fill=MouseNumber_Condition)) + 
  stat_boxplot(geom='errorbar') +
  geom_boxplot(outlier.color = NA) + 
  geom_point(alpha = 0.2, size = 0.2, position = position_jitter(width = 0.1, seed = 1234)) +
  scale_x_discrete(expand = c(0.04,0.04))+
  theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA, linewidth=2),
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    plot.margin = margin(t = 5, r = 5, b = 40, l = 5, unit = "pt"),
    legend.position = "none") + 
  ggtitle("Number of Genes per Cell by individual animal") +
  xlab("Mouse number") + 
  ylab("Number of Genes per Cell") 


#QC fig f - nCount_RNA per mouse
pf_mouse_numi <- ggplot(all_metadata, 
                        aes(x=MouseNumber_Condition, y=nCount_RNA, fill=MouseNumber_Condition)) + 
  stat_boxplot(geom='errorbar') +
  geom_boxplot(outlier.color = NA) + 
  geom_point(alpha = 0.2, size = 0.2, position = position_jitter(width = 0.1, seed = 1234)) +
  scale_x_discrete(expand = c(0.04,0.04))+
  theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA, linewidth=2),
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    plot.margin = margin(t = 5, r = 5, b = 40, l = 5, unit = "pt"),
    legend.position = "none") + 
  ggtitle("nUMI per Cell by individual animal") +
  xlab("Mouse number") + 
  ylab("nUMI per Cell")




# Create 6-panel publication-ready QC figure
pdf(file.path(outdir, "supp_qc_figure_6panels.pdf"),
    width = 14, height = 18)
grid.arrange(pa_cluster_cells, pb_cluster_ngene, pc_cluster_numi, 
             pd_cluster_mt, pe_mouse_ngene, pf_mouse_numi, 
             ncol = 2) 
dev.off()

print("********** QC figures and tables exported successfully! **********")

## END