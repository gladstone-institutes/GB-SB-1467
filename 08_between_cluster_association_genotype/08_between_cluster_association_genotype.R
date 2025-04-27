#!/usr/bin/env Rscript

###############################################################################
## Project ID: SB-1467
## Authors: Michela Traglia
##
## Script Goal: Between -cluster comparison
##   For this analyses, we will use the **lme4** package in R to fit generalized 
##   linear mixed effects models. We are going the model the change in the chance 
##   (or more formally the odds) of cells from a given mouse belonging to a given 
##   cluster from MAPT-fE2 genotype to the MAPT-fE3 genotype. The random effects part
##   of these models captures the inherent correlation between the cells coming 
##   from the same mouse. Age is used as a confounder.
##
## Note: This script was run locally on Ayushi's laptop
##
## Usage example:
## Rscript step10_01_log_odds_calculation_between_cluster_comparison_AgeConfounder_E2E3.R
###############################################################################

# load required packages -------------------------------------------------------
library(lme4)
library(dplyr)
require(magrittr)
library(Seurat)
library(ggplot2)
require(lme4)
require(gdata)


# Define input variables -------------------------------------------------------
pca_dim <- c(15,20,30) 
cluster_res <- c(0.08, 
                 0.1, 
                 0.2,
                 0.3,
                 0.4,
                 0.5,
                 0.6
                 )

ref_genotype <- "WT - Vitamin B3"


# set directory paths ----------------------------------------------------------
basedir <- "/gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/results"
for(pc in pca_dim ){
	for(res in cluster_res){
	
indir <- file.path(basedir, "06_clustering_cell_type/",pc,"PC_",cluster_res,"res")
outdir <- file.path(basedir, 
                    "08_between_cluster_association_genotype/",pc,"PC_",res,"res")

#create the output directory if it doesn't exist
if(!(dir.exists(outdir))){
  dir.create(outdir, recursive = T, showWarnings = FALSE)
}
setwd(outdir)

# Load the Seurat object -------------------------------------------------------
dat <- readRDS(file.path(indir,
                         paste0("gb_sb_1467_",
                         		pc,
                                "PC_",
                                res,
                                "res_clustered_and_cell_typed.rds")))


# Main analysis ----------------------------------------------------------------
# read in the counts if they already exist
if(file.exists(paste0("gb_sb_1467_",pc,"PC_",res,"res_counts_per_sample_per_cluster.csv")){
  counts <- read.csv(paste0("gb_sb_1467_",pc,"PC_",res,"res_counts_per_sample_per_cluster.csv"), header = T)
} else{
  ##the number of cells in each mouse in each cluster
  all_metadata <- dat[[]]
  
  #make a table of the required data
  smp_cluster_counts <- unique(all_metadata %>%
                                 group_by(SampleID) %>%
                                 mutate(total_numbers_of_cells_per_sample = 
                                          n()) %>%
                                 group_by(seurat_clusters, .add=TRUE) %>%
                                 mutate(number_of_cells_per_sample_in_cluster = 
                                          n()) %>%
                                 select(SampleID,
                                        Condition,
                                        seurat_clusters,
                                        total_numbers_of_cells_per_sample,
                                        number_of_cells_per_sample_in_cluster))
  colnames(smp_cluster_counts)[1:3] <- c("sample_id","animal_model","cluster_id")
  smp_cluster_counts <- smp_cluster_counts[order(smp_cluster_counts$cluster_id),]
  write.csv(smp_cluster_counts, 
            file = paste0("gb_sb_1467_",pc,"PC_",res,"res_counts_per_sample_per_cluster.csv"),
            row.names = FALSE)
  counts <- read.csv(paste0("gb_sb_1467_",pc,"PC_",res,"res_counts_per_sample_per_cluster.csv"), header = T)
}

pheno <- counts %>%
  select(sample_id, animal_model, total_numbers_of_cells_per_sample) %>%
  unique()
rownames(pheno) <- seq(1,nrow(pheno))
Clusters <- unique(counts$cluster_id)

##function to estimate the change in the odds of cluster membership from the reference to the other genotypes
estimateCellStateChange <- function(k, counts, pheno) {
  print(paste("Cluster", k))
  cluster_counts <- counts %>% 
    filter(cluster_id == k)
  cluster_counts %<>% merge(., pheno, all.y=TRUE)
  cluster_counts[is.na(cluster_counts$number_of_cells_per_sample_in_cluster),
                 "number_of_cells_per_sample_in_cluster"] <- 0
  
  cluster_counts %<>% 
    arrange(animal_model) %>% 
    mutate(proportion=
             number_of_cells_per_sample_in_cluster/total_numbers_of_cells_per_sample)
  
   
  ##plot proportion of cells per genotype
  pdf(paste0("proportion_of_cells_per_genotype_cluster_",k,"_",pc,"PC_",res,"res.pdf"))
  print(ggplot(cluster_counts, 
               aes(x=animal_model, 
                   y=((number_of_cells_per_sample_in_cluster+0.02)/total_numbers_of_cells_per_sample),
                   colour = Age)) +
          geom_boxplot(outlier.color = NA) +
          geom_jitter(position=position_jitterdodge(0.1)) +
          scale_y_log10() +
          ylab("proportion of cells per genotype")+
          ggtitle(paste0("Cluster ", k)) +
          theme_bw())
  dev.off()
  
  cluster_counts %<>% mutate(animal_model = as.factor(animal_model))
  cluster_counts$animal_model <- relevel(cluster_counts$animal_model, ref=ref_genotype)
  
  formula1=as.formula(paste0("cbind(number_of_cells_per_sample_in_cluster, ",
                             "total_numbers_of_cells_per_sample - ",
                             "number_of_cells_per_sample_in_cluster) ~ ",
                             "(1 | sample_id) + animal_model"))
  glmerFit <- glmer(formula1, data = cluster_counts, family = binomial, nAGQ=10)
  
  sglmerFit1 <- summary(glmerFit)
  TempRes1 <- (sglmerFit1$coefficients[-1,])
  
  return(TempRes1[1,])
}

#run the log odds function for all clusters
ClusterRes <- sapply(Clusters, estimateCellStateChange, counts, pheno)
#reformat the results table
ClusterRes %<>% 
  as.data.frame() %>% 
  t() 
row.names(ClusterRes) <-  paste0("Cluster", Clusters)
ClusterRes <- data.frame(ClusterRes)
colnames(ClusterRes)[c(1:6, 10:12)] <- c("logOddsRatio_KO-VitaminB3_vs_WT-VitaminB3",
                                         "logOddsRatio_KO+VitaminB3_vs_WT-VitaminB3",
									     "logOddsRatio_WT+VitaminB3_vs_WT-VitaminB3",
									     "standardError_KO-VitaminB3_vs_WT-VitaminB3",
									     "standardError_KO+VitaminB3_vs_WT-VitaminB3",
									     "standardError_WT+VitaminB3_vs_WT-VitaminB3",
									     "pvalue_KO-VitaminB3_vs_WT-VitaminB3",
									     "pvalue_KO+VitaminB3_vs_WT-VitaminB3",
									     "pvalue_WT+VitaminB3_vs_WT-VitaminB3")
# make a vector of all p-values and p.adjust for all p-values together
p.adjust_all <- p.adjust(c(ClusterRes$`pvalue_KO-VitaminB3_vs_WT-VitaminB3`, 
                           ClusterRes$`pvalue_KO+VitaminB3_vs_WT-VitaminB3`,
                           ClusterRes$`pvalue_WT+VitaminB3_vs_WT-VitaminB3`), method = "BH")
##perform multiple-testing correction
ClusterRes[,"p.adjust-KO-VitaminB3_vs_WT-VitaminB3"] = p.adjust_all[1:nrow(ClusterRes)]
ClusterRes[,"p.adjust-KO+VitaminB3_vs_WT-VitaminB3"] = p.adjust_all[(nrow(ClusterRes) + 1):(nrow(ClusterRes)*2)]
ClusterRes[,"p.adjust-WT+VitaminB3_vs_WT-VitaminB3"] = p.adjust_all[(nrow(ClusterRes)*2 + 1):length(p.adjust_all)]

##output the results
ClusterRes <- ClusterRes[,!(colnames(ClusterRes) %in% c("X7","X8","X9"))]
print(ClusterRes)
write.csv(ClusterRes, file = paste0("log_odds_ratio_per_genotype_per_cluster", 
                                    "_reference_genotype_", "WT-VitaminB3",
                                    "_",pc,"PC",
                                    "_",res,"res.csv"))



#get data in long form
ClusterRes$cluster_id <- rownames(ClusterRes)
x1 <- ClusterRes[,c("cluster_id", 
                                 "logOddsRatio_KO-VitaminB3_vs_WT-VitaminB3", 
                                 "standardError_KO-VitaminB3_vs_WT-VitaminB3",
                                 "pvalue_KO-VitaminB3_vs_WT-VitaminB3",
                                 "p.adjust_KO-VitaminB3_vs_WT-VitaminB3")]
x1$genotype <- "KO-VitaminB3"
colnames(x1) <- c("cluster_id","logOddsRatio","standardError","pvalue","p.adjust","genotype")

x2 <- ClusterRes[,c("cluster_id", 
                    "logOddsRatio_KO+VitaminB3_vs_WT-VitaminB3", 
                    "standardError_KO+VitaminB3_vs_WT-VitaminB3",
                    "pvalue_KO+VitaminB3_vs_WT-VitaminB3",
                    "p.adjust-KO+VitaminB3_vs_WT-VitaminB3")]
x2$genotype <- "KO+VitaminB3"
colnames(x2) <- c("cluster_id","logOddsRatio","standardError","pvalue","p.adjust","genotype")


x3 <- ClusterRes[,c("cluster_id", 
                    "logOddsRatio_WT+VitaminB3_vs_WT-VitaminB3", 
                    "standardError_logOddsRatio_WT+VitaminB3_vs_WT-VitaminB3",
                    "pvalue_WT+VitaminB3_vs_WT-VitaminB3",
                    "p.adjust-WT+VitaminB3_vs_WT-VitaminB3")]

x3$genotype <- "-WT+VitaminB3"
colnames(x3) <- c("cluster_id","logOddsRatio","standardError","pvalue","p.adjust","genotype")


ClusterRes_plot <- rbind(x1,x2,x3)
rownames(ClusterRes_plot) <- seq(1,nrow(ClusterRes_plot))
ClusterRes_plot$cluster_id <- factor(ClusterRes_plot$cluster_id,
                                     levels = paste0("Cluster",seq(1,max(Clusters))))

#make the box plot for each cluster
#add label for p.adjusted < 0.05
label.df <- ClusterRes_plot[ClusterRes_plot$p.adjust < 0.05,]
if(nrow(label.df) >0 ){
  label.df$logOddsRatio <- 1.5
  print(paste0("Found ", nrow(label.df) ," associations with p.adjusted < 0.05."))
} else {
  print("No association has p.adjusted < 0.05.")
}

label.df2 <- ClusterRes_plot[ClusterRes_plot$pvalue < 0.05,]
if(nrow(label.df2) >0 ){
  label.df2$logOddsRatio <- 1.5
  print(paste0("Found ", nrow(label.df2) ," associations with pvalue < 0.05."))
} else {
  print("No association has pvalue < 0.05.")
}


pdf(paste0("log_odds_boxplot_with_padjusted_reference_genotype_", ref_genotype,
           "_",pc,"PC_",res,"res.pdf"),
    height = 30,
    width = 50
)
ggplot(ClusterRes_plot, 
       aes(x=genotype, y=logOddsRatio, fill=genotype)) +
  geom_hline(yintercept=0) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=logOddsRatio-standardError, ymax=logOddsRatio+standardError), width=.2,
                position=position_dodge(.9)) +
  facet_wrap(~ cluster_id) +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 40)) +
  geom_text(data = label.df, label = "*",
            vjust = 1.5, size = 30) +
  xlab("") +
  ylab(paste0("log odds ratio: animal model vs ", ref_genotype))
dev.off()


pdf(paste0("log_odds_boxplot_with_raw_pvalue_reference_genotype_", ref_genotype,
           "_",pc,"PC_",res,"res.pdf"),
    height = 30,
    width = 50
)
ggplot(ClusterRes_plot, 
       aes(x=genotype, y=logOddsRatio, fill=genotype)) +
  geom_hline(yintercept=0) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=logOddsRatio-standardError, ymax=logOddsRatio+standardError), width=.2,
                position=position_dodge(.9)) +
  facet_wrap(~ cluster_id) +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 40)) +
  geom_text(data = label.df2[label.df2$pvalue < 0.001,], label = "***",
            col = "darkgrey", vjust = 1.5,
            size = 30) +
  geom_text(data = label.df2[label.df2$pvalue > 0.001 & label.df2$pvalue < 0.01,], label = "**",
            col = "darkgrey", vjust = 1.5,
            size = 30) +
  geom_text(data = label.df2[label.df2$pvalue > 0.01 & label.df2$pvalue < 0.05,], label = "*",
            col = "darkgrey", vjust = 1.5,
            size = 30) +
  xlab("") +
  ylab(paste0("log odds ratio: animal model vs ", ref_genotype)) 
dev.off()
}
	}

#add session info
writeLines(capture.output(sessionInfo()), 
           "sessionInfo.txt")


#################### END ####################