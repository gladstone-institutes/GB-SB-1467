rm(list=ls())
library(lme4)
library(dplyr)
require(magrittr)
library(Seurat)
library(ggplot2)
require(lme4)
require(gdata)
library(cowplot)


basedir <- "~/Gladstone Dropbox/Michela Traglia/GB-SB-1467/"
setwd(basedir)
pca_dim <- c(#15,20,
             30) 
cluster_res <- c(0.08#, 
                 #0.1, 
                 #0.2,
                 #0.3,
                 #0.4,
                 #0.5,
                 #0.6
                 )
for(pc in 1:length(pca_dim)){
	for(res in 1:length(cluster_res)){

indir <- paste0(basedir, "06_clustering_cell_type/",pca_dim[pc],"PC_",cluster_res[res],"res/")
outdir <- paste0(basedir, 
                    "08_between_cluster_association_genotype/",pca_dim[pc],"PC_",cluster_res[res],"res/log_odds_wo_normalization/count+1/")


ClusterRes_WT <- read.csv( paste0(outdir,"log_odds_ratio_per_genotype_per_cluster", 
                                    "_reference_genotype_", "WT-VitaminB3",
                                    "_",pca_dim[pc],"PC",
                                    "_",cluster_res[res],"res.csv"))


cell_annot <- read.csv("~/Gladstone Dropbox/Michela Traglia/GB-SB-1467/06_clustering_cell_type_30PC_0.08res_noRDS_08May2025/gb_sb_1467_30PC_0.08res_sctype_embryonic_brain_celltype_annotation_table.csv")

cell_annot %<>% mutate(cluster_id = paste0("Cluster", cluster))
cell_annot %<>% filter(type != "Hippocampal granule cell")

full_results_WT <- merge(ClusterRes_WT, cell_annot, by.x="X", by.y="cluster_id")
full_results_WT <- full_results_WT[order(as.numeric(full_results_WT$cluster)),]


full_results_WT %<>% mutate(logOddsRatio_KO_wo_VitaminB3_vs_WT_wo_VitaminB3 = round(logOddsRatio_KO.VitaminB3_vs_WT.VitaminB3*100)/100,
                         logOddsRatio_KO_w_VitaminB3_vs_WT_wo_VitaminB3 = round(logOddsRatio_KO.VitaminB3_vs_WT.VitaminB3.1*100)/100,
                         logOddsRatio_WT_w_VitaminB3_vs_WT_wo_VitaminB3 = round(logOddsRatio_WT.VitaminB3_vs_WT.VitaminB3*100)/100,
                         pvalue_KO_wo_VitaminB3_vs_WT_wo_VitaminB3 = round(pvalue_KO.VitaminB3_vs_WT.VitaminB3*1000)/1000,
                         pvalue_KO_w_VitaminB3_vs_WT_wo_VitaminB3 = round(pvalue_KO.VitaminB3_vs_WT.VitaminB3.1*1000)/1000,
                         pvalue_WT_w_VitaminB3_vs_WT_wo_VitaminB3 = round(pvalue_WT.VitaminB3_vs_WT.VitaminB3*1000)/1000 ,
                         )

plotdir <- "~/Gladstone Dropbox/Michela Traglia/GB-SB-1467/08_between_cluster_association_genotype/30PC_0.08res/log_odds_wo_normalization/count+1/gb_sb_1467_proportion_between_condition_boxplots_scatterplots"

celltype_abbr <- c(
  "0_UL-ExN",    # Upper-Layer Excitatory Neurons
  "1_HypoTh-Neur",  # Hypothalamic/Thalamic Neurons HypoTh_Neur
  "2_Astro",     # Astrocytes
  "3_SMSN",       # Striatum Medium Spiny Neurons SMSN
  "4_Prog",      # Progenitors
  "5_OPC",       # Oligodendrocyte Progenitor Cells
  "6_MGE-Int",   # MGE-interneuron
  "7_DL-ExN",    # Deeper-Layer Excitatory Neurons
  "8_Mural",     # Mural Cells
  "9_Micro",     # Microglia
  "10_LGE-Int",  # LGE/CGE Interneurons
  "11_Epend",    # Ependymal Cells
  "12_Endo",     # Endothelial Cells
  "13_HPC-N",    # Hippocampal Neurons
  "14_Unknown",  # Unknown
  "15_Purk"      # Purkinje Cells
)

count_data <- read.csv("~/Gladstone Dropbox/Michela Traglia/GB-SB-1467/06_clustering_cell_type_30PC_0.08res_noRDS_08May2025/gb_sb_1467_30PC_0.08res_counts_per_sample_per_cluster.csv",
                       header = TRUE)

pheno_data <- count_data %>% dplyr::select(sample_id, animal_model,total_cells_per_sample) %>% unique(.)

##with sex info from Ayushi
male_mice <- paste0("sample",c("A", "C", "D", "F", "G", "H", "I", "J", "K"))
pheno_data[["sex"]] <- "female"
pheno_data$sex[pheno_data$sample_id %in% male_mice] <- "male"

orig_meta_data <- read.csv("~/Gladstone Dropbox/Michela Traglia/GB-SB-1467/08_between_cluster_association_genotype/30PC_0.08res/V25-03_Samples.csv",
                           header = T)

orig_meta_data %<>% mutate(sample_id = paste0("sample", SampleID))
pheno_data %<>% merge(., orig_meta_data)

counts <- count_data %>% dplyr::select(-animal_model, -total_cells_per_sample) %>% 
  tidyr::pivot_wider(.,names_from = sample_id, values_from = number_of_cells_per_sample_in_cluster) %>%
  as.data.frame(.) %>% arrange(cluster_id)

row.names(counts) <- paste0("cluster_", counts$cluster_id)
counts <- counts[,-1]

counts_res <- list(counts = counts)

pheno_data %<>% mutate(animal_model = as.factor(animal_model),
                       animal_model = relevel(animal_model, ref = "WT - Vitamin B3"))

pheno_data %<>% mutate(animal_model = as.factor(animal_model),
                       animal_model = relevel(animal_model, ref = "KO - Vitamin B3"))

full_count_data <- merge(count_data, pheno_data, by="sample_id")

#WT - B3 vs KO - B3
##WT - B3 vs WT + B3
##KO + B3 vs KO - B3
##KO + B3 vs WT - B3
cols <- c('orange','#66CCFF', '#009933', 'gold')
  
cond1="KO + Vitamin B3"
cond2="WT - Vitamin B3"
temp_data <- full_count_data %>% filter(animal_model.x %in% c(cond2,cond1))

old_name_celltype_abbr<- data.frame(cluster_id=unlist(strsplit(celltype_abbr, "_"))[2*(1:length(celltype_abbr))-1], renamed_cluster_id=celltype_abbr)
old_name_celltype_abbr[2,2] <- "1_HypoTh_Neur"
temp_data$cluster_id <- as.factor(temp_data$cluster_id)
temp_data_renamed <-  left_join(temp_data, old_name_celltype_abbr) 
temp_data_renamed$renamed_cluster_id[temp_data_renamed$renamed_cluster_id =="1_HypoTh_Neur" ] <- "1_HypoTh-Neur"
temp_data_renamed$renamed_cluster_id <- factor(temp_data_renamed$renamed_cluster_id, levels=celltype_abbr)
###logOdds plots 
pheno_data %<>% mutate(animal_model = as.factor(animal_model),
                       animal_model = relevel(animal_model, ref = "WT - Vitamin B3"))



old_name_comp_results_WT<- data.frame(cluster=full_results_WT$cluster, renamed_cluster_id=celltype_abbr)
old_name_comp_results_WT$renamed_cluster_id[old_name_comp_results_WT$renamed_cluster_id =="1_HypoTh_Neur" ] <- "1_HypoTh-Neur"

comp_results_renamed_WT <-  left_join(full_results_WT, old_name_comp_results_WT) 
comp_results_renamed_WT$renamed_cluster_id <- factor(comp_results_renamed_WT$renamed_cluster_id, levels=celltype_abbr)
comp_results_renamed_WT$dir <- "NA"
comp_results_renamed_WT$dir[comp_results_renamed_WT$logOddsRatio_KO_w_VitaminB3_vs_WT_wo_VitaminB3 < 0 ] <- "NEG"
comp_results_renamed_WT$dir[comp_results_renamed_WT$logOddsRatio_KO_w_VitaminB3_vs_WT_wo_VitaminB3 >= 0 ] <- "POS"


#plot
   p3<- (ggplot(temp_data_renamed, 
               aes(y=renamed_cluster_id, 
                   x=((number_of_cells_per_sample_in_cluster+0.02)/(total_cells_per_sample.x)),
                   fill=as.factor(animal_model.x))) +
          geom_boxplot() +
          xlab("proportion of cells per condition")+
          ylab("cell type") +
          ggtitle(paste0("")) +
          scale_fill_manual(values = cols[c(2,3)])+  
          labs(fill='Condition')+
          theme_bw()+ theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=8)))
 
 
  p4<- (ggplot(comp_results_renamed_WT, 
               aes(y=renamed_cluster_id, 
                   x=logOddsRatio_KO_w_VitaminB3_vs_WT_wo_VitaminB3,
                   color=as.numeric( pvalue_KO_w_VitaminB3_vs_WT_wo_VitaminB3)))+
           geom_point(aes(shape=dir, colour = cut( pvalue_KO_w_VitaminB3_vs_WT_wo_VitaminB3, c(0,0.1), include.lowest = T), size=abs(logOddsRatio_KO_w_VitaminB3_vs_WT_wo_VitaminB3), fill=ifelse(as.numeric(pvalue_KO_w_VitaminB3_vs_WT_wo_VitaminB3) <=0.1, "darkred","grey"))) +
           geom_vline(xintercept=0, lty="dashed", color="grey81")+
          xlab("logOR")+
          ylab("cell type") +
          ggtitle(paste0("")) +
          theme_bw()+ theme(legend.position = "none")
          )

  




  ### 

cond1="WT + Vitamin B3"
cond2="WT - Vitamin B3"
temp_data <- full_count_data %>% filter(animal_model.x %in% c(cond2,cond1))

old_name_celltype_abbr<- data.frame(cluster_id=unlist(strsplit(celltype_abbr, "_"))[2*(1:length(celltype_abbr))-1], renamed_cluster_id=celltype_abbr)
old_name_celltype_abbr[2,2] <- "1_HypoTh_Neur"
temp_data$cluster_id <- as.factor(temp_data$cluster_id)
temp_data_renamed <-  left_join(temp_data, old_name_celltype_abbr) 
temp_data_renamed$renamed_cluster_id[temp_data_renamed$renamed_cluster_id =="1_HypoTh_Neur" ] <- "1_HypoTh-Neur"
temp_data_renamed$renamed_cluster_id <- factor(temp_data_renamed$renamed_cluster_id, levels=celltype_abbr)
###logOdds plots 
pheno_data %<>% mutate(animal_model = as.factor(animal_model),
                       animal_model = relevel(animal_model, ref = "WT - Vitamin B3"))






old_name_comp_results_WT<- data.frame(cluster=full_results_WT$cluster, renamed_cluster_id=celltype_abbr)
old_name_comp_results_WT$renamed_cluster_id[old_name_comp_results_WT$renamed_cluster_id =="1_HypoTh_Neur" ] <- "1_HypoTh-Neur"

comp_results_renamed_WT <-  left_join(full_results_WT, old_name_comp_results_WT) 
comp_results_renamed_WT$renamed_cluster_id <- factor(comp_results_renamed_WT$renamed_cluster_id, levels=celltype_abbr)
comp_results_renamed_WT$dir <- "NA"
comp_results_renamed_WT$dir[comp_results_renamed_WT$logOddsRatio_WT_w_VitaminB3_vs_WT_wo_VitaminB3 < 0 ] <- "NEG"
comp_results_renamed_WT$dir[comp_results_renamed_WT$logOddsRatio_WT_w_VitaminB3_vs_WT_wo_VitaminB3 >= 0 ] <- "POS"


#plot
   p5<- (ggplot(temp_data_renamed, 
               aes(y=renamed_cluster_id, 
                   x=((number_of_cells_per_sample_in_cluster+0.02)/(total_cells_per_sample.x)),
                   fill=as.factor(animal_model.x))) +
          geom_boxplot() +
          xlab("proportion of cells per condition")+
          ylab("cell type") +
          ggtitle(paste0("")) +
          scale_fill_manual(values = cols[c(3,4)])+  
          labs(fill='Condition')+
          theme_bw()+ theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=8)))
 
 
  p6<- (ggplot(comp_results_renamed_WT, 
               aes(y=renamed_cluster_id, 
                   x=logOddsRatio_WT_w_VitaminB3_vs_WT_wo_VitaminB3,
                   color=as.numeric( pvalue_WT_w_VitaminB3_vs_WT_wo_VitaminB3)))+
           geom_point(aes(shape=dir, colour = cut( pvalue_WT_w_VitaminB3_vs_WT_wo_VitaminB3, c(0,0.1), include.lowest = T), size=abs(logOddsRatio_WT_w_VitaminB3_vs_WT_wo_VitaminB3), fill=ifelse(as.numeric(pvalue_WT_w_VitaminB3_vs_WT_wo_VitaminB3) <=0.1, "darkred","grey"))) +
           geom_vline(xintercept=0, lty="dashed", color="grey81")+
          xlab("logOR")+
          ylab("cell type") +
          ggtitle(paste0("")) +
          theme_bw()+ theme(legend.position = "none")
          )

  

 
  
cond1="KO - Vitamin B3"
cond2="WT - Vitamin B3"
temp_data <- full_count_data %>% filter(animal_model.x %in% c(cond2,cond1))

old_name_celltype_abbr<- data.frame(cluster_id=unlist(strsplit(celltype_abbr, "_"))[2*(1:length(celltype_abbr))-1], renamed_cluster_id=celltype_abbr)
old_name_celltype_abbr[2,2] <- "1_HypoTh_Neur"
temp_data$cluster_id <- as.factor(temp_data$cluster_id)
temp_data_renamed <-  left_join(temp_data, old_name_celltype_abbr) 
temp_data_renamed$renamed_cluster_id[temp_data_renamed$renamed_cluster_id =="1_HypoTh_Neur" ] <- "1_HypoTh-Neur"
temp_data_renamed$renamed_cluster_id <- factor(temp_data_renamed$renamed_cluster_id, levels=celltype_abbr)
###logOdds plots 
pheno_data %<>% mutate(animal_model = as.factor(animal_model),
                       animal_model = relevel(animal_model, ref = "WT - Vitamin B3"))



old_name_comp_results_WT<- data.frame(cluster=full_results_WT$cluster, renamed_cluster_id=celltype_abbr)
old_name_comp_results_WT$renamed_cluster_id[old_name_comp_results_WT$renamed_cluster_id =="1_HypoTh_Neur" ] <- "1_HypoTh-Neur"

comp_results_renamed_WT <-  left_join(full_results_WT, old_name_comp_results_WT) 
comp_results_renamed_WT$renamed_cluster_id <- factor(comp_results_renamed_WT$renamed_cluster_id, levels=celltype_abbr)

comp_results_renamed_WT$dir <- "NA"
comp_results_renamed_WT$dir[comp_results_renamed_WT$logOddsRatio_KO_wo_VitaminB3_vs_WT_wo_VitaminB3 < 0 ] <- "NEG"
comp_results_renamed_WT$dir[comp_results_renamed_WT$logOddsRatio_KO_wo_VitaminB3_vs_WT_wo_VitaminB3 >= 0 ] <- "POS"

#plot
   p7<- (ggplot(temp_data_renamed, 
               aes(y=renamed_cluster_id, 
                   x=((number_of_cells_per_sample_in_cluster+0.02)/(total_cells_per_sample.x)),
                   fill=as.factor(animal_model.x))) +
          geom_boxplot() +
          xlab("proportion of cells per condition")+
          ylab("cell type") +
          ggtitle(paste0("")) +
          scale_fill_manual(values = cols[c(1,3)])+  
          labs(fill='Condition')+
          theme_bw()+ theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=8)))
 
 
  p8<- (ggplot(comp_results_renamed_WT, 
               aes(y=renamed_cluster_id, 
                   x=logOddsRatio_KO_wo_VitaminB3_vs_WT_wo_VitaminB3,
                   color=as.numeric( pvalue_KO_wo_VitaminB3_vs_WT_wo_VitaminB3)))+
           geom_point(aes(shape=dir, colour = cut( pvalue_KO_wo_VitaminB3_vs_WT_wo_VitaminB3, c(0,0.1), include.lowest = T), size=abs(logOddsRatio_KO_wo_VitaminB3_vs_WT_wo_VitaminB3), fill=ifelse(as.numeric(pvalue_KO_wo_VitaminB3_vs_WT_wo_VitaminB3) <=0.1, "darkred","grey"))) +
           geom_vline(xintercept=0, lty="dashed", color="grey81")+
          xlab("logOR")+
          ylab("cell type") +
          ggtitle(paste0("")) +
          theme_bw()+ theme(legend.position = "none")
          )

  

###from here + add logodds with KO- as reference and add +1 to the logodds formula

ClusterRes_KO <- read.csv( paste0(outdir,"log_odds_ratio_per_genotype_per_cluster", 
                                    "_reference_genotype_", "KO-VitaminB3",
                                    "_",pca_dim[pc],"PC",
                                    "_",cluster_res[res],"res.csv"))

full_results_KO <- merge(ClusterRes_KO, cell_annot, by.x="X", by.y="cluster_id")
full_results_KO <- full_results_KO[order(as.numeric(full_results_KO$cluster)),]


full_results_KO %<>% mutate(logOddsRatio_KO_w_VitaminB3_vs_KO_wo_VitaminB3 = round(logOddsRatio_KO.VitaminB3_vs_KO.VitaminB3*100)/100,
                         logOddsRatio_WT_wo_VitaminB3_vs_KO_wo_VitaminB3 = round(logOddsRatio_WT.VitaminB3_vs_KO.VitaminB3*100)/100,
                         logOddsRatio_WT_w_VitaminB3_vs_KO_wo_VitaminB3 = round(logOddsRatio_WT.VitaminB3_vs_KO.VitaminB3.1*100)/100,
                         pvalue_KO_w_VitaminB3_vs_KO_wo_VitaminB3 = round(pvalue_KO.VitaminB3_vs_KO.VitaminB3*1000)/1000,
                         pvalue_WT_wo_VitaminB3_vs_KO_wo_VitaminB3 = round(pvalue_WT.VitaminB3_vs_KO.VitaminB3*1000)/1000,
                         pvalue_WT_w_VitaminB3_vs_KO_wo_VitaminB3 = round(pvalue_WT.VitaminB3_vs_KO.VitaminB3.1*1000)/1000 ,
                         )


cond1="KO + Vitamin B3"
cond2="KO - Vitamin B3"
temp_data <- full_count_data %>% filter(animal_model.x %in% c(cond2,cond1))

old_name_celltype_abbr<- data.frame(cluster_id=unlist(strsplit(celltype_abbr, "_"))[2*(1:length(celltype_abbr))-1], renamed_cluster_id=celltype_abbr)
old_name_celltype_abbr[2,2] <- "1_HypoTh_Neur"
temp_data$cluster_id <- as.factor(temp_data$cluster_id)
temp_data_renamed <-  left_join(temp_data, old_name_celltype_abbr) 
temp_data_renamed$renamed_cluster_id[temp_data_renamed$renamed_cluster_id =="1_HypoTh_Neur" ] <- "1_HypoTh-Neur"
temp_data_renamed$renamed_cluster_id <- factor(temp_data_renamed$renamed_cluster_id, levels=celltype_abbr)
###logOdds plots 
pheno_data %<>% mutate(animal_model = as.factor(animal_model),
                       animal_model = relevel(animal_model, ref = "KO - Vitamin B3"))


#rerun logODDS with KO- as reference


old_name_comp_results_KO<- data.frame(cluster=full_results_KO$cluster, renamed_cluster_id=celltype_abbr)
old_name_comp_results_KO$renamed_cluster_id[old_name_comp_results_KO$renamed_cluster_id =="1_HypoTh_Neur" ] <- "1_HypoTh-Neur"

comp_results_renamed_KO <-  left_join(full_results_KO, old_name_comp_results_KO) 
comp_results_renamed_KO$renamed_cluster_id <- factor(comp_results_renamed_KO$renamed_cluster_id, levels=celltype_abbr)

comp_results_renamed_KO$dir <- "NA"
comp_results_renamed_KO$dir[comp_results_renamed_KO$logOddsRatio_KO_w_VitaminB3_vs_KO_wo_VitaminB3 < 0 ] <- "NEG"
comp_results_renamed_KO$dir[comp_results_renamed_KO$logOddsRatio_KO_w_VitaminB3_vs_KO_wo_VitaminB3 >= 0 ] <- "POS"


#plot
   p1<- (ggplot(temp_data_renamed, 
               aes(y=renamed_cluster_id, 
                   x=((number_of_cells_per_sample_in_cluster+0.02)/(total_cells_per_sample.x)),
                   fill=as.factor(animal_model.x))) +
          geom_boxplot() +
          xlab("proportion of cells per condition")+
          ylab("cell type") +
          ggtitle(paste0("")) +
          scale_fill_manual(values = cols[c(1,2)])+  
          labs(fill='Condition')+
          theme_bw()+ theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=8)))
 
 
  p2<- (ggplot(comp_results_renamed_KO, 
               aes(y=renamed_cluster_id, 
                   x=logOddsRatio_KO_w_VitaminB3_vs_KO_wo_VitaminB3,
                   color=as.numeric( pvalue_KO_w_VitaminB3_vs_KO_wo_VitaminB3)))+
           geom_point(aes(shape=dir, colour = cut( pvalue_KO_w_VitaminB3_vs_KO_wo_VitaminB3, c(0,0.1), include.lowest = T), size=abs(logOddsRatio_KO_w_VitaminB3_vs_KO_wo_VitaminB3), fill=ifelse(as.numeric(pvalue_KO_w_VitaminB3_vs_KO_wo_VitaminB3) <=0.1, "darkred","grey"))) +
           geom_vline(xintercept=0, lty="dashed", color="grey81")+
          xlab("logOR")+
          ylab("cell type") +
          ggtitle(paste0("")) +
          theme_bw()+ theme(legend.position = "none")
          )


  ### 
pdf(paste0(plotdir, "/subfig_supp_B.pdf"), width=10, height=20)
plot_grid(p2, p1, p4, p3, p6, p5, p8, p7, align = "h", ncol = 2)
dev.off()