##support function 1
ineqfun_data <- function(pars, Total_Cells) {
  ineqres <- (pars*Total_Cells)
  return(ineqres)
}

##support function 2
estimate_variance_w_args <- function(pars, counts, Total_Cells, optim_clusters) {
  
  total_var.norm <- vector(mode = "numeric")
  nsamp <- ncol(counts)
  nclust <- nrow(counts)
  log_odds <- matrix(NA, nclust, nsamp)
  for(s in 1:ncol(counts)) {
    if(sum(pars[s]*Total_Cells[s] - counts[,s] < 0) > 0)
      print((pars[s]*Total_Cells[s] - counts[,s]))
    log_odds[,s] <- log((counts[,s]+1)/(pars[s]*Total_Cells[s] - counts[,s]))
  }
  for(c in optim_clusters) {
    temp_var_estimate <- (1/(nsamp - 1))*t(log_odds[c,])%*%(diag(nsamp) - (1/nsamp)*(rep(1, nsamp))%*%t(rep(1, nsamp)))%*%((log_odds[c,]))
    total_var.norm[c] <- temp_var_estimate[1,1]
  }
  return(sum(total_var.norm[optim_clusters]))
}



##New method for diff cell composition
logodds_optimized_norm_factors_likelihood_test_design <- function(counts_res, pheno_data) {
  
  found_stable_solution <- FALSE
  estimates <- vector(mode = "numeric")
  estimates_significance <- vector(mode = "numeric")
  
  optim_clusters <- 1:nrow(counts)
  old_outlier_clusters <- 1
  count_iter <- 1
  
  nsamp <- ncol(counts)
  while(!found_stable_solution & count_iter < 10) {
    print(count_iter)
    count_iter <- count_iter + 1
    
    
    counts <- counts_res$counts
    ##idea 1: sort by random effects
    #pheno_data <- data.frame(sample = colnames(counts), grp = c(rep(0, nsamp/2), rep(1, nsamp/2)))
    
    grp <- pheno_data$animal_model
    des <- model.matrix(~grp)
    row.names(des) <- colnames(counts)
    
    total_cells = data.frame(total_cells = colSums(counts),
                             sampleID = colnames(counts),
                             group = grp)
    
    
    
    counts %<>% as.data.frame() %>% tibble::rownames_to_column(var = "cluster")
    
    counts_long <- counts %>% tidyr::pivot_longer( cols = starts_with("s"), # Specify the columns to pivot
                                                   names_to = "sampleID",     # New column to hold variable names
                                                   values_to = "count"                # New column to hold values
    )
    
    
    counts_long_norm <- merge(counts_long, total_cells)
    counts_long_norm %<>% mutate(sampleID = as.factor(sampleID))
    
    Total_Cells <- total_cells$total_cells
    
    counts <- counts_res$counts
    
    
    
    estimate_variance <- function(pars) {
      estimate_variance_w_args(pars, counts, Total_Cells, optim_clusters)
    }
    
    ineqfun <- function(pars) {
      ineqres <- ineqfun_data(pars, Total_Cells)
      return(ineqres)
    }
    
    
    ui <- diag(nsamp)
    for(s in 1:nsamp) {
      ui[s,s] <- Total_Cells[s]
    }
    ci <- colMaxs(counts %>% as.matrix(.))
    
    temp_res <- solnp(pars = runif(nsamp, min = 0.7, max = 1.4), 
                      estimate_variance, 
                      LB = rep(0.1, nsamp),
                      UB = rep(5, nsamp),
                      ineqfun = ineqfun,
                      ineqLB = colMaxs(counts %>% as.matrix(.)),
                      ineqUB = rep(2*max(Total_Cells), nsamp),
                      control = list(tol = 1e-12,delta = 1e-9, trace = 0),
    )
    
    optim.pars <- temp_res$pars/exp(mean(log(temp_res$pars)))
    
    optim_factor <- data.frame(sampleID = colnames(counts),
                               optim.norm.factor = optim.pars)
    counts_long_norm %<>% merge(., optim_factor)
    
    total_var.norm <- vector(mode = "numeric")
    temp_count <- 0
    temp_count_cluster <- 1
    comparison <- vector(mode = "character")
    cluster_id <- vector(mode = "character")
    cluster_significance <- vector(mode = "numeric")
    for(clust in paste0("cluster_", 0:(nrow(counts) - 1))) {
      glmerFit <- glmer(cbind(count+1, (round(total_cells*optim.norm.factor) - count)) ~ (1|sampleID) + group,
                        data = counts_long_norm %>% filter(cluster == clust) ,
                        family = "binomial")
      
      sglmerFit <- summary(glmerFit)
      estimates[(temp_count+1):(temp_count+3)] <- sglmerFit$coefficients[2:4,1]
      
      glmerFit0 <- glmer(cbind(count+1, (round(total_cells*optim.norm.factor) - count)) ~ (1|sampleID),
                         data = counts_long_norm %>% filter(cluster == clust) ,
                         family = "binomial")
      
      anova_res <- anova(glmerFit, glmerFit0)
      
      cluster_significance[temp_count_cluster] <- anova_res$`Pr(>Chisq)`[2]
      estimates_significance[(temp_count+1):(temp_count+3)] <- sglmerFit$coefficients[2:4,4]
      comparison[(temp_count+1):(temp_count+3)] <- row.names(sglmerFit$coefficients)[-1]
      cluster_id[(temp_count+1):(temp_count+3)] <- rep(clust, 3)
      temp_count <- temp_count + 3
      temp_count_cluster <- temp_count_cluster + 1
    }
    
    results <- data.frame(cluster_id, comparison, estimates, estimates_significance)
    # names(estimates) <- paste0("logOR_estimates_", row.names(counts_res$counts))
    # names(estimates_significance) <- paste0("logOR_estimates_significance_", row.names(counts_res$counts))
    
    temp_res <- data.frame(cluster = 1:nrow(counts), cluster_significance) %>% arrange(cluster_significance)
    outlier_clusters <- temp_res %>% 
      filter(cluster_significance < 0.05) %>% 
      dplyr::slice(1:min(c(nrow(.), floor(nrow(counts)/2)))) %>%
      .$cluster %>%
      sort()
    
    if(identical(old_outlier_clusters, outlier_clusters)) {
      found_stable_solution <- TRUE
    }else{
      old_outlier_clusters <- outlier_clusters
      optim_clusters <- setdiff(1:nrow(counts), outlier_clusters)
    }
    
    #View(cbind(estimates, counts_res$observed_log_odds_ratios, change_mean, estimates_significance))
    
  }
  

  
  return(list(results=results, optim_factor = optim_factor))
  
}                                          

##variation of proposed method (can ignore)
 logodds_tmm_norm_factors_design <- function(counts_res, pheno_data) {
  
  counts <- counts_res$counts
  ##idea 1: sort by random effects
  #pheno_data <- data.frame(sample = colnames(counts), grp = c(rep(0, nsamp/2), rep(1, nsamp/2)))
  
  grp <- pheno_data$animal_model
  
  des <- model.matrix(~grp)
  row.names(des) <- colnames(counts)
  
  total_cells = data.frame(total_cells = colSums(counts),
                           sampleID = colnames(counts),
                           tmm.norm.factor = calcNormFactors(counts),
                           group = grp)
  
  
  total_cells %<>% mutate(total_cells.norm = round(total_cells*tmm.norm.factor))
  
  counts %<>% as.data.frame() %>% tibble::rownames_to_column(var = "cluster")
  
  counts_long <- counts %>% tidyr::pivot_longer( cols = starts_with("S"), # Specify the columns to pivot
                                                 names_to = "sampleID",     # New column to hold variable names
                                                 values_to = "count"                # New column to hold values
  )
  
  
  counts_long_norm <- merge(counts_long, total_cells)
  counts_long_norm %<>% mutate(sampleID = as.factor(sampleID))
  
  Total_Cells <- total_cells$total_cells
  
  counts <- counts_res$counts
  
  total_var.norm <- vector(mode = "numeric")
  temp_count <- 0
  temp_count_cluster <- 1
  comparison <- vector(mode = "character")
  cluster_id <- vector(mode = "character")
  cluster_significance <- vector(mode = "numeric")
  
  estimates <- vector(mode = "numeric")
  estimates_significance <- vector(mode = "numeric")
  for(clust in paste0("cluster_", 0:(nrow(counts) - 1))) {
    glmerFit <- glmer(cbind((count+1), (total_cells - count)) ~ (1|sampleID) + group,
                      data = counts_long_norm %>% filter(cluster == clust)  ,
                      family = "binomial",
                      nAGQ = 1,
                      control =glmerControl(optimizer = "bobyqa"))
    
    sglmerFit <- summary(glmerFit)
    cluster_significance[temp_count_cluster] <- anova_res$`Pr(>Chisq)`[2]
    estimates[(temp_count+1):(temp_count+3)] <- sglmerFit$coefficients[2:4,1]
    estimates_significance[(temp_count+1):(temp_count+3)] <- sglmerFit$coefficients[2:4,4]
    comparison[(temp_count+1):(temp_count+3)] <- row.names(sglmerFit$coefficients)[-1]
    cluster_id[(temp_count+1):(temp_count+3)] <- rep(clust, 3)
    temp_count <- temp_count + 3
    temp_count_cluster <- temp_count_cluster + 1
  }
  
  results <- data.frame(cluster_id, comparison, estimates, estimates_significance)
  
  # names(estimates) <- paste0("logOR_estimates_", paste0("c", 0:(nrow(counts) - 1)))
  # names(estimates_significance) <- paste0("logOR_estimates_significance_", paste0("c", 0:(nrow(counts) - 1)))
  # 
  # 
  # diff_res <- list(estimates = estimates,
  #                  estimates_significance = estimates_significance)
  
  return(results)
  
}                                          


require(edgeR)
require(lme4)
require(Rsolnp)
library(dplyr)
library(magrittr)
library(ggplot2)
library(matrixStats)




count_data <- read.csv("/Volumes/Jain-Boinformatics-Collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/results/06_clustering_cell_type/30PC_0.08res/gb_sb_1467_30PC_0.08res_counts_per_sample_per_cluster.csv",
                       header = TRUE)

pheno_data <- count_data %>% dplyr::select(sample_id, animal_model,total_cells_per_sample) %>% unique(.)

##with sex info from Ayushi
male_mice <- paste0("sample",c("A", "C", "D", "F", "G", "H", "I", "J", "K"))
pheno_data[["sex"]] <- "female"
pheno_data$sex[pheno_data$sample_id %in% male_mice] <- "male"

orig_meta_data <- read.csv("/Volumes/Jain-Boinformatics-Collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/data/V25-03_Samples.csv",
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

res <- logodds_optimized_norm_factors_likelihood_test_design(counts_res, pheno_data)
#res_tmm <- logodds_tmm_norm_factors_design(counts_res, pheno_data)


#res$results %>% write.csv(.,
                 # "gb_sb_1467_30PC_0.08res_normalized_log_odds_diff_comp.csv",
                  #row.names = F)


pheno_data %<>% mutate(optim_factor = res$optim_factor$optim.norm.factor)

full_count_data <- merge(count_data, pheno_data, by="sample_id")

cell_annot <- read.csv("/Volumes/Jain-Boinformatics-Collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/results/06_clustering_cell_type/30PC_0.08res/gb_sb_1467_30PC_0.08res_sctype_brain_celltype_annotation_table.csv")

cell_annot %<>% mutate(cluster_id = paste0("cluster_", cluster))
cell_annot %<>% filter(type != "Hippocampal granule cell")

res_ko <- res$results
full_results <- merge(res_ko, cell_annot)

select_results <- full_results %>%
  dplyr::select(cluster_id, type, comparison, estimates, estimates_significance)

colnames(select_results) <- c("cluster_id",
                              "cell_type",
                              "comparison_vs_KO_wo_B3",
                              "logOR",
                              "p.value")
select_results %<>% mutate(logOR = round(logOR*100)/100,
                           p.value = round(p.value*1000)/1000)

#select_results %>%
  #write.csv(., "gb_sb_1467_30PC_0.08res_w_annot_normalized_log_odds_diff_comp.csv",
            #row.names = F)

require(ggplot2)

for(c in 0:15) {
  temp_data <- full_count_data %>% filter(cluster_id == c)
  
  pdf(paste0("gb_sb_1467_proportion_of_cells_per_condition_cluster_",c,"_",30,"PC_res_","0_08.pdf"))
  print(ggplot(temp_data, 
               aes(x=animal_model.x, 
                   y=((number_of_cells_per_sample_in_cluster+0.02)/(total_cells_per_sample.x*optim_factor)),
                   )) +
          geom_boxplot(outlier.color = NA) +
          geom_jitter(position=position_jitterdodge(dodge.width = 0.001), size = 4, aes(color = Genotype, shape = sex)) +
          scale_y_log10() +
          ylab("normalized proportion of cells per condition")+
          xlab("condition") +
          ggtitle(paste0("Cluster ", c)) +
          theme_bw())
  dev.off()
}


library(paletteer)

pdf(paste0("gb_sb_1467_proportion_of_cells_per_condition_all_cluster_",30,"PC_res_","0_08.pdf"))
  print(ggplot(full_count_data, 
               aes(x=as.factor(animal_model.x), 
                   y=((number_of_cells_per_sample_in_cluster+0.02)/(total_cells_per_sample.x*optim_factor)),
                   fill=as.factor(cluster_id))) +
          geom_bar(stat = "identity") +
          ylab("normalized proportion of cells per condition")+
          xlab("cluster id") +
          ggtitle(paste0("All clusters")) +
          scale_color_manual(values = paletteer_d("ggsci::category20_d3")[1:16])+
          theme_bw())
  dev.off()
  
setwd("/Volumes/Jain-Boinformatics-Collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/results/08_between_cluster_association_genotype/30PC_0.08res/log_odds_w_normalization/gb_sb_1467_proportion_between_condition_boxplots_scatterplots")

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


#WT - B3 vs KO - B3
##WT - B3 vs WT + B3
##KO + B3 vs KO - B3
##KO + B3 vs WT - B3
cols <- c('orange','#66CCFF', '#009933', 'gold')
cond1="KO + Vitamin B3"
cond2="KO - Vitamin B3"
temp_data <- full_count_data %>% filter(animal_model.x %in% c(cond2,cond1))

old_name_celltype_abbr<- data.frame(cluster_id=unlist(strsplit(celltype_abbr, "_"))[2*(1:length(celltype_abbr))-1], renamed_cluster_id=celltype_abbr)
old_name_celltype_abbr[2,2] <- "1_HypoTh_Neur"

temp_data$cluster_id <- as.factor(temp_data$cluster_id)
temp_data_renamed <-  left_join(temp_data, old_name_celltype_abbr) 
###logOdds plots
pheno_data %<>% mutate(animal_model = as.factor(animal_model),
                       animal_model = relevel(animal_model, ref = "KO - Vitamin B3"))

res <- logodds_optimized_norm_factors_likelihood_test_design(counts_res, pheno_data)
all_results <- res$results 
comp_results <- all_results %>% dplyr::filter(comparison %in% "groupKO + Vitamin B3")

old_name_comp_results<- data.frame(cluster_id=comp_results$cluster_id, renamed_cluster_id=celltype_abbr)
comp_results_renamed <-  left_join(comp_results, old_name_comp_results) 



#plot

#pdf(paste0("gb_sb_1467_proportion_of_cells_between_condition_all_cluster_",30,"PC_res_","0_08_",cond1,"_vs_",cond2,".pdf"), width=10)
pdf(paste0("subfig_supp_B.pdf"), width=10, height=20)

   p1<- (ggplot(temp_data_renamed, 
               aes(y=factor(renamed_cluster_id, levels=celltype_abbr), 
                   x=((number_of_cells_per_sample_in_cluster+0.02)/(total_cells_per_sample.x*optim_factor)),
                   fill=as.factor(animal_model.x))) +
          geom_boxplot() +
          xlab("normalized proportion of cells per condition")+
          ylab("cell type") +
          ggtitle(paste0("")) +
          scale_fill_manual(values = cols[c(1,2)])+  
          labs(fill='Condition')+
          theme_bw())
 
   

  p2<- (ggplot(comp_results_renamed, 
               aes(y=factor(renamed_cluster_id, levels=celltype_abbr), 
                   x=estimates,
                   color=as.numeric(estimates_significance)))+
           geom_point(aes(colour = cut(estimates_significance, c(0,0.1)), size=abs(estimates), fill=ifelse(as.numeric(estimates_significance) <=0.1, "darkred","grey"))) +
           geom_vline(xintercept=0, lty="dashed", color="grey81")+
          xlab("logOR")+
          ylab("cell type") +
          ggtitle(paste0("")) +
          theme_bw()+ theme(legend.position = "none")
          )

  ### 

 
###
 
  
cond1="KO + Vitamin B3"
cond2="WT - Vitamin B3"
temp_data <- full_count_data %>% filter(animal_model.x %in% c(cond2,cond1))

old_name_celltype_abbr<- data.frame(cluster_id=unlist(strsplit(celltype_abbr, "_"))[2*(1:length(celltype_abbr))-1], renamed_cluster_id=celltype_abbr)
old_name_celltype_abbr[2,2] <- "1_HypoTh_Neur"
temp_data$cluster_id <- as.factor(temp_data$cluster_id)
temp_data_renamed <-  left_join(temp_data, old_name_celltype_abbr) 
###logOdds plots 
pheno_data %<>% mutate(animal_model = as.factor(animal_model),
                       animal_model = relevel(animal_model, ref = "WT - Vitamin B3"))

res <- logodds_optimized_norm_factors_likelihood_test_design(counts_res, pheno_data)
all_results <- res$results 
comp_results <- all_results %>% dplyr::filter(comparison %in% "groupKO + Vitamin B3")
old_name_comp_results<- data.frame(cluster_id=comp_results$cluster_id, renamed_cluster_id=celltype_abbr)
comp_results_renamed <-  left_join(comp_results, old_name_comp_results) 


#plot
#pdf(paste0("gb_sb_1467_proportion_of_cells_between_condition_all_cluster_",30,"PC_res_","0_08_",cond1,"_vs_",cond2,".pdf"), width=10)
   p3<- (ggplot(temp_data_renamed, 
               aes(y=factor(renamed_cluster_id, levels=celltype_abbr), 
                   x=((number_of_cells_per_sample_in_cluster+0.02)/(total_cells_per_sample.x*optim_factor)),
                   fill=as.factor(animal_model.x))) +
          geom_boxplot() +
          xlab("normalized proportion of cells per condition")+
          ylab("cell type") +
          ggtitle(paste0("")) +
          scale_fill_manual(values = cols[c(2,3)])+  
          labs(fill='Condition')+
          theme_bw())
 
   

  p4<- (ggplot(comp_results_renamed, 
               aes(y=factor(renamed_cluster_id, levels=celltype_abbr), 
                   x=estimates,
                   color=as.numeric(estimates_significance)))+
           geom_point(aes(colour = cut(estimates_significance, c(0,0.1)), size=abs(estimates), fill=ifelse(as.numeric(estimates_significance) <=0.1, "darkred","grey"))) +
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

###logOdds plots 
pheno_data %<>% mutate(animal_model = as.factor(animal_model),
                       animal_model = relevel(animal_model, ref = "WT - Vitamin B3"))

res <- logodds_optimized_norm_factors_likelihood_test_design(counts_res, pheno_data)
all_results <- res$results 
comp_results <- all_results %>% dplyr::filter(comparison %in% "groupWT + Vitamin B3")
old_name_comp_results<- data.frame(cluster_id=comp_results$cluster_id, renamed_cluster_id=celltype_abbr)
comp_results_renamed <-  left_join(comp_results, old_name_comp_results) 
#plot

#pdf(paste0("gb_sb_1467_proportion_of_cells_between_condition_all_cluster_",30,"PC_res_","0_08_",cond1,"_vs_",cond2,".pdf"), width=10)
   p5<- (ggplot(temp_data_renamed, 
               aes(y=factor(renamed_cluster_id, levels=celltype_abbr), 
                   x=((number_of_cells_per_sample_in_cluster+0.02)/(total_cells_per_sample.x*optim_factor)),
                   fill=as.factor(animal_model.x))) +
          geom_boxplot() +
          xlab("normalized proportion of cells per condition")+
          ylab("cell type") +
          ggtitle(paste0("")) +
          scale_fill_manual(values = cols[c(3,4)])+  
          labs(fill='Condition')+
          theme_bw())
 
   

  p6<- (ggplot(comp_results_renamed, 
               aes(y=factor(renamed_cluster_id, levels=celltype_abbr), 
                   x=estimates,
                   color=as.numeric(estimates_significance)))+
           geom_point(aes(colour = cut(estimates_significance, c(0,0.1)), size=abs(estimates), fill=ifelse(as.numeric(estimates_significance) <=0.1, "darkred","grey"))) +
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
###logOdds plots
pheno_data %<>% mutate(animal_model = as.factor(animal_model),
                       animal_model = relevel(animal_model, ref = "WT - Vitamin B3"))

res <- logodds_optimized_norm_factors_likelihood_test_design(counts_res, pheno_data)
all_results <- res$results 
comp_results <- all_results %>% dplyr::filter(comparison %in% "groupKO - Vitamin B3")
old_name_comp_results<- data.frame(cluster_id=comp_results$cluster_id, renamed_cluster_id=celltype_abbr)
comp_results_renamed <-  left_join(comp_results, old_name_comp_results) 

#pdf(paste0("gb_sb_1467_proportion_of_cells_between_condition_all_cluster_",30,"PC_res_","0_08_",cond1,"_vs_",cond2,".pdf"))
   p7<- (ggplot(temp_data_renamed, 
               aes(y=factor(renamed_cluster_id, levels=celltype_abbr), 
                   x=((number_of_cells_per_sample_in_cluster+0.02)/(total_cells_per_sample.x*optim_factor)),
                   fill=as.factor(animal_model.x))) +
          geom_boxplot() +
          xlab("normalized proportion of cells per condition")+
          ylab("cell type") +
          ggtitle(paste0("")) +
          scale_fill_manual(values = cols[c(1,3)])+  
          labs(fill='Condition')+
          theme_bw())
 
   

  p8<- (ggplot(comp_results_renamed, 
               aes(y=factor(renamed_cluster_id, levels=celltype_abbr), 
                   x=estimates,
                   color=as.numeric(estimates_significance)))+
           geom_point(aes(colour = cut(estimates_significance, c(0,0.1)), size=abs(estimates), fill=ifelse(as.numeric(estimates_significance) <=0.1, "darkred","grey"))) +
           geom_vline(xintercept=0, lty="dashed", color="grey81")+
          xlab("logOR")+
          ylab("cell type") +
          ggtitle(paste0("")) +
          theme_bw()+ theme(legend.position = "none")
          )
library(cowplot)
plot_grid(p2, p1, p4, p3, p6, p5, p8, p7, align = "h", ncol = 2, rel_heights = c(0.2, 0.2))

dev.off()