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
  

  
  return(results)
  
}                                          


# logodds_tmm_norm_factors_design <- function(counts_res, pheno_data) {
  
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
    glmerFit <- glmer(cbind((count+1), (total_cells.norm - count)) ~ (1|sampleID) + group,
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



count_data <- read.csv("/Volumes/Jain-Boinformatics-Collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/results/06_clustering_cell_type/30PC_0.08res/gb_sb_1467_30PC_0.08res_counts_per_sample_per_cluster.csv",
                       header = TRUE)

pheno_data <- count_data %>% dplyr::select(sample_id, animal_model,total_cells_per_sample) %>% unique(.)

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
res %>% write.csv(.,
                  "gb_sb_1467_30PC_0.08res_normalized_log_odds_diff_comp.csv",
                  row.names = F)
