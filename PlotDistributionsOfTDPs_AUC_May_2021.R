rm(list = ls())
setwd("~/Dropbox (Gladstone)/PFOCRInPathwayAnalyses")
require(dplyr)
require(tidyr)
require(magrittr)
require(ggplot2)
pval_res <- readRDS("./PFOCRInPathwayAnalyses/pvalue_results_human_voom.rds")
load(paste0(getwd(), "/PFOCRInPathwayAnalyses/databases.RData"))
load(paste0(getwd(), "/PFOCRInPathwayAnalyses/databases_pfocr_3intersect_v2.RData"))
load(paste0(getwd(), "/PFOCRInPathwayAnalyses/pfocr_miny_April2021.RData"))
ocr_length_3sets <- sapply(pfocr_3sets_list, length)
ocr_length_miny <- sapply(pfocr_miny_list, length)
go_length <- sapply(go_list, length)
wp_length <- sapply(wp_list, length)

GeneSetSizes <- data.frame(Length=c(ocr_length_3sets, ocr_length_miny, go_length, wp_length), Database=c(rep("ocr_3sets", length(ocr_length_3sets)), rep("ocr_miny", length(ocr_length_miny)), rep("go", length(go_length)), rep("wp", length(wp_length))))
GeneSetSizes %<>% filter(Length >= 3 & Length <= 500)
pdf("Gene_sets_size_distribution_w_2_pfocr_databases.pdf")
ggplot(GeneSetSizes, aes(x=Length, color=Database)) + stat_ecdf() +xlab("Gene set size") + ylab("Probability") +
  theme(text = element_text(size = 15))
dev.off()
ggplot(GeneSetSizes, aes(x=Length, color=Database)) + geom_density(adjust=3)

pvalue_results_human_voom=readRDS("./PFOCRInPathwayAnalyses/pvalue_results_human_voom.rds")
gene_entrez=readRDS("./PFOCRInPathwayAnalyses/gene_entrez.rds")

GetFractionAssociationsInPathwayGenes <- function(index, pvalue_results_human_voom, path_annotation, gene_entrez) {
  path_genes <- path_annotation %>%
    .$gene %>%
    unique() %>%
    as.data.frame() 
  colnames(path_genes)[1] <- "ENTREZID"
  temp <- p.adjust(pvalue_results_human_voom[,index], method = "holm")
  diff_res <- data.frame(gene_entrez, adj.pval=temp)
  diff_res %<>% merge(.,path_genes) %>%
    filter(!is.na(adj.pval))
  
  Nass <- diff_res %>%
    filter(adj.pval < 0.05) %>%
    nrow()
  
  return(Nass/nrow(diff_res))
}



GSEid <- list.files(path = "./PFOCRInPathwayAnalyses/GSE/")
minSig <- 5
minSize <- 3
theta <- 0.9
Nperm <- 1000
thetaL <- as.character(100*theta)
ChooseGSE <- 1

permute_res_dir <- paste0(getwd(),"/permuted_geneset_databases_results_updated_pfocr_gene_set_sizes_may_2021/")

PermuteResFiles  <- permute_res_dir %>%
  list.files() %>%
  as.data.frame() 

colnames(PermuteResFiles) <- "RObject"

PermuteResFiles %<>% mutate(RObject=as.character(RObject)) %>%
  filter(!grepl("pfocr", RObject)) 

GetGSEIndices <- function(x) {
  return(strsplit(x, "_")[[1]][3])
}

GSEIndices <- PermuteResFiles %>%
  .$RObject %>% 
  sapply(.,GetGSEIndices) %>%
  as.integer() %>%
  unique()

quantiles <- c(0.50, 0.75, 0.90, 0.95, 0.99)

get_area <- function(index, rNsig, quantiles) {
  temp_tdp <- rNsig[index,-1] %>% as.numeric()
  area <- 0
  for(i in 1:(length(quantiles) - 1)) {
    area <- area + 
      0.5*(quantiles[i+1] - quantiles[i])*(temp_tdp[i] + temp_tdp[i+1])
  }
  area <- area/0.49
  return(area)
}

get_area_tdp <- function(temp_tdp, quantiles) {
  area <- 0
  for(i in 1:(length(quantiles) - 1)) {
    area <- area + 
      0.5*(quantiles[i+1] - quantiles[i])*(temp_tdp[i] + temp_tdp[i+1])
  }
  area <- area/0.49
  return(area)
}

plot_GSE <- function(ChooseGSE, plot_data, observed_data, nsig_data, database, oNsig, tdp.fullL) {
  pdf(paste0("GSE_index_", ChooseGSE,"_Permuted_",database,"_plots_tdp.full_", tdp.fullL,".pdf"))
  print(ggplot(nsig_data, aes(x=Nsig)) + 
    geom_histogram() +
    geom_vline(xintercept = oNsig, lty=2, col="red") +
    ggtitle(paste0(database, ": GSE index = ", ChooseGSE, " tdp.full = ", tdp.fullL)) + 
    xlab("No of significant gene sets") +
    theme(text = element_text(size = 12)))
  
  print(ggplot(nsig_data, aes(x=Nsig)) + 
    stat_ecdf() +
    geom_vline(xintercept = oNsig, lty=2, col="red") +
    ggtitle(paste0(database, ": GSE index = ", ChooseGSE, " tdp.full = ", tdp.fullL)) +
    xlab("No of significant gene sets") +
    theme(text = element_text(size = 12)))
  
  print(ggplot(plot_data, aes(x=Value)) +
          geom_histogram() +
          facet_wrap(~Metric, nrow=3) +
          geom_vline(observed_data, mapping=aes(xintercept=Value), lty=2, col="red") +
          ggtitle(paste0(database, ": GSE index = ", ChooseGSE, " tdp.full = ", tdp.fullL)) +
          xlab("TDP.bound") +
          theme(text = element_text(size = 12)) )
  
  print(ggplot(plot_data, aes(x=Value)) +
    stat_ecdf() +
    facet_wrap(~Metric, nrow=3) +
    geom_vline(observed_data, mapping=aes(xintercept=Value), lty=2, col="red") +
    ggtitle(paste0(database, ": GSE index = ", ChooseGSE, " tdp.full = ", tdp.fullL)) + 
    xlab("TDP.bound") +
    theme(text = element_text(size = 12)) )
  
  dev.off()
}

plot_TDPQuantiles_AUC <- function(ChooseGSE, GSEid, ocr_length, go_length, wp_length, minSig, minSize, theta, quantiles) {
  print(ChooseGSE)
  q_ocr <- NA
  q_ocr_miny <- NA
  q_wp <- NA
  q_go <- NA
  
  permute_res <- readRDS(paste0(permute_res_dir, "GSE_index_", ChooseGSE,"_", Nperm, "_perm_result.rds"))
  
  rand_ocr_tdp <- sapply(1:Nperm, get_area, permute_res$rNsig_pfcor, quantiles)
  
  rand_ocr_tdp_miny <- sapply(1:Nperm, get_area, permute_res$rNsig_pfcor_miny, quantiles)
  
  rand_wp_tdp <- sapply(1:Nperm, get_area, permute_res$rNsig_wp, quantiles)
  
  rand_go_tdp <- sapply(1:Nperm, get_area, permute_res$rNsig_go, quantiles)
  
  rand_ocr_nsig_mean <- permute_res$rNsig_pfcor %>%
    as.data.frame() %>%
    .$Nsig %>%
    mean()
  
  rand_ocr_nsig_mean_miny <- permute_res$rNsig_pfcor_miny %>%
    as.data.frame() %>%
    .$Nsig %>%
    mean()
  
  rand_wp_nsig_mean <- permute_res$rNsig_wp %>%
    as.data.frame() %>%
    .$Nsig %>%
    mean()
  
  rand_go_nsig_mean <- permute_res$rNsig_go %>%
    as.data.frame() %>%
    .$Nsig %>%
    mean()
  
  rand_ocr_nsig_sd <- permute_res$rNsig_pfcor %>%
    as.data.frame() %>%
    .$Nsig %>%
    sd()
  
  rand_ocr_nsig_sd_miny <- permute_res$rNsig_pfcor_miny %>%
    as.data.frame() %>%
    .$Nsig %>%
    sd()
  
  rand_wp_nsig_sd <- permute_res$rNsig_wp %>%
    as.data.frame() %>%
    .$Nsig %>%
    sd()
  
  rand_go_nsig_sd <- permute_res$rNsig_go %>%
    as.data.frame() %>%
    .$Nsig %>%
    sd()
  
  tdp.fullL <- gsub("[.]", "_", round(1000*permute_res$TDPbound_full)/1000)
  ##plot random distribution versus observed
  ##WP-figures
  wp_plot_data <- data.frame(permute_res$rNsig_wp[,-1], area=rand_wp_tdp) %>%
    gather(.,Metric,Value)
  
  wp_observed_data <- data.frame(Value = c(get_area_tdp(permute_res$oTDP_bound_90_wp, quantiles),(permute_res$oTDP_bound_90_wp)))
  wp_observed_data %<>% mutate(Metric=c("area", paste0("quantile_", c(50, 75, 90, 95, 99))))
  
  wp_nsig_data <- data.frame(Nsig=permute_res$rNsig_wp[,1])
  
  plot_GSE(ChooseGSE, wp_plot_data, wp_observed_data, wp_nsig_data, "WP", permute_res$oNsig_wp, tdp.fullL)
  
  ##GO-figures
  go_plot_data <- data.frame(permute_res$rNsig_go[,-1], area=rand_go_tdp) %>%
    gather(.,Metric,Value)
  
  go_observed_data <- data.frame(Value = c(get_area_tdp(permute_res$oTDP_bound_90_go, quantiles),(permute_res$oTDP_bound_90_go)))
  go_observed_data %<>% mutate(Metric=c("area", paste0("quantile_", c(50, 75, 90, 95, 99))))
  
  go_nsig_data <- data.frame(Nsig=permute_res$rNsig_go[,1])
  
  plot_GSE(ChooseGSE, go_plot_data, go_observed_data, go_nsig_data, "GO", permute_res$oNsig_go, tdp.fullL)
 
  ##PFOCR-figures
  pfocr_plot_data <- data.frame(permute_res$rNsig_pfcor[,-1], area=rand_ocr_tdp) %>%
    gather(.,Metric,Value)
  
  pfocr_observed_data <- data.frame(Value = c(get_area_tdp(permute_res$oTDP_bound_90_pfocr, quantiles),(permute_res$oTDP_bound_90_pfocr)))
  pfocr_observed_data %<>% mutate(Metric=c("area", paste0("quantile_", c(50, 75, 90, 95, 99))))
  
  pfocr_nsig_data <- data.frame(Nsig=permute_res$rNsig_pfcor[,1])
  
  plot_GSE(ChooseGSE, pfocr_plot_data, pfocr_observed_data, pfocr_nsig_data, "PFOCR-3sets", permute_res$oNsig_pfocr, tdp.fullL)
 
  ##PFOCR-miny-figures
  pfocr_miny_plot_data <- data.frame(permute_res$rNsig_pfcor_miny[,-1], area=rand_ocr_tdp_miny) %>%
    gather(.,Metric,Value)
  
  pfocr_miny_observed_data <- data.frame(Value = c(get_area_tdp(permute_res$oTDP_bound_90_pfocr_miny, quantiles),(permute_res$oTDP_bound_90_pfocr_miny)))
  pfocr_miny_observed_data %<>% mutate(Metric=c("area", paste0("quantile_", c(50, 75, 90, 95, 99))))
  
  pfocr_miny_nsig_data <- data.frame(Nsig=permute_res$rNsig_pfcor_miny[,1])
  
  plot_GSE(ChooseGSE, pfocr_miny_plot_data, pfocr_miny_observed_data, pfocr_miny_nsig_data, "PFOCR-miny", permute_res$oNsig_pfocr_miny, tdp.fullL)
  
}

plot_TDPQuantiles_AUC(148, GSEid, ocr_length, go_length, wp_length, minSig, minSize, theta, quantiles)

get_random_prob <- function(s, tdp_full) {
  pbinom(min(s, floor(tdp_full*s)+1),s, prob=tdp_full, lower.tail = FALSE)
}

getTDPQuantiles_AUC <- function(ChooseGSE, GSEid, ocr_length, go_length, wp_length, minSig, minSize, theta, quantiles, pfocr_annotation_3sets, pfocr_miny_annotation, wp_annotation, go_annotation,gene_entrez,pvalue_results_human_voom) {
  print(ChooseGSE)
  q_ocr <- NA
  q_ocr_miny <- NA
  q_wp <- NA
  q_go <- NA
  
  permute_res <- readRDS(paste0(permute_res_dir, "GSE_index_", ChooseGSE,"_", Nperm, "_perm_result.rds"))
  
  rand_ocr_tdp <- sapply(1:Nperm, get_area, permute_res$rNsig_pfcor, quantiles)
 
  rand_ocr_tdp_miny <- sapply(1:Nperm, get_area, permute_res$rNsig_pfcor_miny, quantiles)
  
  rand_wp_tdp <- sapply(1:Nperm, get_area, permute_res$rNsig_wp, quantiles)
  
  rand_go_tdp <- sapply(1:Nperm, get_area, permute_res$rNsig_go, quantiles)
  
  rand_ocr_nsig_mean <- permute_res$rNsig_pfcor %>%
    as.data.frame() %>%
    .$Nsig %>%
    mean()
  
  rand_ocr_nsig_mean_miny <- permute_res$rNsig_pfcor_miny %>%
    as.data.frame() %>%
    .$Nsig %>%
    mean()
  
  rand_wp_nsig_mean <- permute_res$rNsig_wp %>%
    as.data.frame() %>%
    .$Nsig %>%
    mean()
  
  rand_go_nsig_mean <- permute_res$rNsig_go %>%
    as.data.frame() %>%
    .$Nsig %>%
    mean()
  
  rand_ocr_nsig_sd <- permute_res$rNsig_pfcor %>%
    as.data.frame() %>%
    .$Nsig %>%
    sd()
  
  rand_ocr_nsig_sd_miny <- permute_res$rNsig_pfcor_miny %>%
    as.data.frame() %>%
    .$Nsig %>%
    sd()
  
  rand_wp_nsig_sd <- permute_res$rNsig_wp %>%
    as.data.frame() %>%
    .$Nsig %>%
    sd()
  
  rand_go_nsig_sd <- permute_res$rNsig_go %>%
    as.data.frame() %>%
    .$Nsig %>%
    sd()
  
  
  # ocr_result <- read.table(paste0(getwd(), "/PFOCRInPathwayAnalyses/GSE/",GSEid[ChooseGSE],"/rSEA/PFOCR_3sets/result.txt"), header = T)
  # tempIndices <- match(ocr_result$set_id, names(ocr_length))
  # temp_ocr_length <- ocr_length[tempIndices]
  # Nocr <- sum(!is.na(ocr_result$Comp.adjP) & ocr_result$Comp.adjP < 0.05 & temp_ocr_length >= minSize, na.rm = T)
  Nocr <- permute_res$oNsig_pfocr
  if(Nocr > minSig) {
    # sig_ocr_result <- data.frame(ocr_result[!is.na(ocr_result$Comp.adjP) & ocr_result$Comp.adjP < 0.05 & temp_ocr_length >= minSize, ], database=rep("ocr", Nocr))
    # q_ocr <- quantile(sig_ocr_result$TDP.estimate, theta, na.rm = T)
    q_ocr <- get_area_tdp(permute_res$oTDP_bound_90_pfocr, quantiles)
    q_ocr <- 100*(q_ocr - mean(rand_ocr_tdp))/(mean(rand_ocr_tdp) + 0.0001)
  }
  diffocr <- (Nocr - rand_ocr_nsig_mean)
  zocr <- diffocr/(rand_ocr_nsig_sd+1)
  
  Nocr_miny <- permute_res$oNsig_pfocr_miny
  if(Nocr_miny > minSig) {
    # sig_ocr_result <- data.frame(ocr_result[!is.na(ocr_result$Comp.adjP) & ocr_result$Comp.adjP < 0.05 & temp_ocr_length >= minSize, ], database=rep("ocr", Nocr))
    # q_ocr <- quantile(sig_ocr_result$TDP.estimate, theta, na.rm = T)
    q_ocr_miny <- get_area_tdp(permute_res$oTDP_bound_90_pfocr_miny, quantiles)
    q_ocr_miny <- 100*(q_ocr_miny - mean(rand_ocr_tdp_miny))/(mean(rand_ocr_tdp_miny) + 0.0001)
  }
  diffocr_miny <- (Nocr_miny - rand_ocr_nsig_mean_miny)
  zocr_miny <- diffocr_miny/(rand_ocr_nsig_sd_miny+1)
  
  # wp_result <- read.table(paste0(getwd(), "/PFOCRInPathwayAnalyses/GSE/",GSEid[ChooseGSE],"/rSEA/WP/result.txt"), header = T)
  # tempIndices <- match(wp_result$set_id, names(wp_length))
  # temp_wp_length <- wp_length[tempIndices]
  # Nwp <- sum(!is.na(wp_result$Comp.adjP) & wp_result$Comp.adjP < 0.05 & temp_wp_length >= minSize, na.rm = T)
  Nwp <- permute_res$oNsig_wp
  if(Nwp > minSig) {
    # sig_wp_result <- data.frame(wp_result[!is.na(wp_result$Comp.adjP) & wp_result$Comp.adjP < 0.05 & temp_wp_length >= minSize, ], database=rep("wp", Nwp))
    # q_wp <- quantile(sig_wp_result$TDP.estimate, theta, na.rm = T)
    q_wp <- get_area_tdp(permute_res$oTDP_bound_90_wp, quantiles)
    q_wp <- 100*(q_wp - mean(rand_wp_tdp))/(mean(rand_wp_tdp) + 0.0001)
  }
  diffwp <- (Nwp - rand_wp_nsig_mean)
  zwp <- diffwp/(rand_wp_nsig_sd+1)
  
  # go_result <- read.table(paste0(getwd(), "/PFOCRInPathwayAnalyses/GSE/",GSEid[ChooseGSE],"/rSEA/GO/result.txt"), header = T)
  # tempIndices <- match(go_result$set_id, names(go_length))
  # temp_go_length <- go_length[tempIndices]
  # Ngo <- sum(!is.na(go_result$Comp.adjP) & go_result$Comp.adjP < 0.05 & temp_go_length >= minSize & temp_go_length < 500, na.rm = T)
  Ngo <- permute_res$oNsig_go
  if(Ngo > minSig) {
    # sig_go_result <- data.frame(go_result[!is.na(go_result$Comp.adjP) & go_result$Comp.adjP < 0.05 & temp_go_length >= minSize & temp_go_length < 500, ], database=rep("go", Ngo))
    # q_go <- quantile(sig_go_result$TDP.estimate, theta, na.rm = T)
    q_go <- get_area_tdp(permute_res$oTDP_bound_90_go, quantiles)
    q_go <- 100*(q_go - mean(rand_go_tdp))/(mean(rand_go_tdp) + 0.0001)
  }
  diffgo <- (Ngo - rand_go_nsig_mean)
  zgo <- diffgo/(rand_go_nsig_sd+1)
  
  # PlotData <- rbind(sig_ocr_result, sig_wp_result)
  # PlotData <- data.frame(rbind(PlotData, sig_go_result))
  # ggplot(PlotData, aes(x=TDP.estimate, color=database)) + stat_ecdf() + geom_vline(xintercept = q_go,lty=2)
  # ggplot(PlotData, aes(y=TDP.estimate, x=database)) + geom_boxplot()
 
  ##theoretical estimates
  
  ##ocr
  path_genes_assoc_fraction <- GetFractionAssociationsInPathwayGenes(ChooseGSE, pvalue_results_human_voom, pfocr_annotation_3sets, gene_entrez)
  ocr_gene_set_nonunique_sizes <- GeneSetSizes %>% filter(Database == "ocr_3sets") %>% .$Length
  ocr_gene_set_sizes <- GeneSetSizes %>% filter(Database == "ocr_3sets") %>% .$Length %>% unique()
  random_ocr_Nsig_sets <- sum(table(ocr_gene_set_nonunique_sizes)*sapply(ocr_gene_set_sizes, get_random_prob,path_genes_assoc_fraction))

  print(c(path_genes_assoc_fraction, permute_res$TDPbound_full))
    ##given gene set size
  ##ocr_miny
  path_genes_assoc_fraction <- GetFractionAssociationsInPathwayGenes(ChooseGSE, pvalue_results_human_voom, pfocr_miny_annotation, gene_entrez)
  ocr_miny_gene_set_nonunique_sizes <- GeneSetSizes %>% filter(Database == "ocr_miny") %>% .$Length
  ocr_miny_gene_set_sizes <- GeneSetSizes %>% filter(Database == "ocr_miny") %>% .$Length %>% unique()
  random_ocr_miny_Nsig_sets <- sum(table(ocr_miny_gene_set_nonunique_sizes)*sapply(ocr_miny_gene_set_sizes, get_random_prob,path_genes_assoc_fraction))
  print(c(path_genes_assoc_fraction, permute_res$TDPbound_full))
  
  ##wp
  path_genes_assoc_fraction <- GetFractionAssociationsInPathwayGenes(ChooseGSE, pvalue_results_human_voom, wp_annotation, gene_entrez)
  wp_gene_set_nonunique_sizes <- GeneSetSizes %>% filter(Database == "wp") %>% .$Length
  wp_gene_set_sizes <- GeneSetSizes %>% filter(Database == "wp") %>% .$Length %>% unique()
  random_wp_Nsig_sets <- sum(table(wp_gene_set_nonunique_sizes)*sapply(wp_gene_set_sizes, get_random_prob,path_genes_assoc_fraction))
  print(c(path_genes_assoc_fraction, permute_res$TDPbound_full))
  
  ##go
  path_genes_assoc_fraction <- GetFractionAssociationsInPathwayGenes(ChooseGSE, pvalue_results_human_voom, go_annotation, gene_entrez)
  go_gene_set_nonunique_sizes <- GeneSetSizes %>% filter(Database == "go") %>% .$Length
  go_gene_set_sizes <- GeneSetSizes %>% filter(Database == "go") %>% .$Length %>% unique()
  random_go_Nsig_sets <- sum(table(go_gene_set_nonunique_sizes)*sapply(go_gene_set_sizes, get_random_prob,path_genes_assoc_fraction))
  print(c(path_genes_assoc_fraction, permute_res$TDPbound_full))
  
  return(c(q_ocr, q_ocr_miny,q_wp, q_go, zocr, zocr_miny,zwp, zgo, diffocr, diffocr_miny, diffwp, diffgo, get_area_tdp(permute_res$oTDP_bound_90_pfocr, quantiles), get_area_tdp(permute_res$oTDP_bound_90_pfocr_miny, quantiles),get_area_tdp(permute_res$oTDP_bound_90_wp, quantiles), get_area_tdp(permute_res$oTDP_bound_90_go, quantiles), mean(rand_ocr_tdp), mean(rand_ocr_tdp_miny), mean(rand_wp_tdp), mean(rand_go_tdp), permute_res$TDPbound_full, rand_ocr_nsig_mean, random_ocr_Nsig_sets, rand_ocr_nsig_mean_miny, random_ocr_miny_Nsig_sets, rand_wp_nsig_mean, random_wp_Nsig_sets, rand_go_nsig_mean, random_go_Nsig_sets, ChooseGSE))
}

##Results for GSEIndices 339 and 398 missing from the results 
GSEIndices %<>% 
  setdiff(.,c(339,398,473, 475, 497)) 
q_results <- t(sapply(GSEIndices, getTDPQuantiles_AUC, GSEid, ocr_length, go_length, wp_length, minSig, minSize, theta, quantiles,pfocr_annotation_3sets, pfocr_miny_annotation, wp_annotation, go_annotation,gene_entrez,pvalue_results_human_voom))
q_results <- data.frame(q_results)
colnames(q_results) <- c("median_centered_tdp_ocr_3sets", "median_centered_tdp_ocr_miny", "median_centered_tdp_wp", "median_centered_tdp_go", "z_ocr", "z_ocr_miny","z_wp", "z_go", "diff_ocr", "diff_ocr_miny","diff_wp", "diff_go","observed_area_tdp_ocr_3sets", "observed_area_tdp_ocr_miny", "observed_area_tdp_wp", "observed_area_tdp_go", "mean_rand_area_tdp_ocr_3sets", "mean_rand_area_tdp_ocr_miny", "mean_rand_area_tdp_wp", "mean_rand_area_tdp_go", "TDP.full", "perm_Nocr","rand_Nsig_ocr","perm_Nocr_miny","rand_Nsig_ocr_miny","perm_Nwp","rand_Nsig_wp", "perm_Ngo","rand_Nsig_go", "GSE_index")

ggplot(q_results, aes(x=perm_Ngo, y=rand_Nsig_go)) + geom_point() + geom_abline(intercept = 0, slope = 1, lty=2, col="red")
ggplot(q_results, aes(x=perm_Nwp, y=rand_Nsig_wp)) + geom_point() + geom_abline(intercept = 0, slope = 1, lty=2, col="red")
ggplot(q_results, aes(x=perm_Nocr, y=rand_Nsig_ocr)) + geom_point() + geom_abline(intercept = 0, slope = 1, lty=2, col="red")
ggplot(q_results, aes(x=perm_Nocr_miny, y=rand_Nsig_ocr_miny)) + geom_point() + geom_abline(intercept = 0, slope = 1, lty=2, col="red")

tdp_thresh <- 0.05
tdp_threshL <- gsub("[.]","_",as.character(tdp_thresh))
q_results_thresh <- q_results %>%
  filter(TDP.full > tdp_thresh)

pdf(paste0("All_databases_observed_vs_random_area_nsig_thresh_",tdp_threshL,".pdf"))
print(ggplot(q_results_thresh, aes(x=mean_rand_area_tdp_wp, y=observed_area_tdp_wp, color=TDP.full)) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0, lty=2, col="red")  +
  ggtitle("WP") +
  theme(text = element_text(size = 12))) 

print(ggplot(q_results_thresh, aes(x=-mean_rand_area_tdp_wp+observed_area_tdp_wp, color=TDP.full)) + 
        geom_histogram() +
        geom_vline(xintercept = 0, lty=2, col="red")  +
        ggtitle("WP") +
        xlab("TDP: Observed - Random") +
        xlim(-0.2,0.2) +
        theme(text = element_text(size = 12))) 

print(ggplot(q_results_thresh, aes(x=-mean_rand_area_tdp_wp+observed_area_tdp_wp, color=TDP.full)) + 
        stat_ecdf() +
        geom_vline(xintercept = 0, lty=2, col="red")  +
        ggtitle("WP") +
        xlab("TDP: Observed - Random") +
        xlim(-0.2,0.2) +
        theme(text = element_text(size = 12))) 

print(ggplot(q_results_thresh, aes(x=diff_wp, color=TDP.full)) + 
        geom_histogram() +
        geom_vline(xintercept = 0, lty=2, col="red")  +
        ggtitle("WP") +
        xlab("Nsig: Observed - Random") +
        theme(text = element_text(size = 12))) 

print(ggplot(q_results_thresh, aes(x=diff_wp, color=TDP.full)) + 
        stat_ecdf() +
        geom_vline(xintercept = 0, lty=2, col="red")  +
        ggtitle("WP") +
        xlab("Nsig: Observed - Random") +
        theme(text = element_text(size = 12))) 

print(ggplot(q_results_thresh, aes(x=mean_rand_area_tdp_ocr_3sets, y=observed_area_tdp_ocr_3sets, color=TDP.full)) + 
  geom_point() +
  geom_abline(slope=1, intercept = 0, lty=2, col="red")  +
  ggtitle("PFOCR-3sets") +
  theme(text = element_text(size = 12))) 

print(ggplot(q_results_thresh, aes(x=-mean_rand_area_tdp_ocr_3sets+observed_area_tdp_ocr_3sets, color=TDP.full)) + 
        geom_histogram() +
        geom_vline(xintercept = 0, lty=2, col="red")  +
        ggtitle("PFOCR-3sets") +
        xlab("Observed - Random") +
        xlim(-0.2,0.2) +
        theme(text = element_text(size = 12))) 

print(ggplot(q_results_thresh, aes(x=-mean_rand_area_tdp_ocr_3sets+observed_area_tdp_ocr_3sets, color=TDP.full)) + 
        stat_ecdf() +
        geom_vline(xintercept = 0, lty=2, col="red")  +
        ggtitle("PFOCR-3sets") +
        xlab("Observed - Random") +
        xlim(-0.2,0.2) +
        theme(text = element_text(size = 12))) 
print(ggplot(q_results_thresh, aes(x=diff_ocr, color=TDP.full)) + 
        geom_histogram() +
        geom_vline(xintercept = 0, lty=2, col="red")  +
        ggtitle("PFOCR-3sets") +
        xlab("Nsig: Observed - Random") +
        theme(text = element_text(size = 12))) 

print(ggplot(q_results_thresh, aes(x=diff_ocr, color=TDP.full)) + 
        stat_ecdf() +
        geom_vline(xintercept = 0, lty=2, col="red")  +
        ggtitle("PFOCR-3sets") +
        xlab("Nsig: Observed - Random") +
        theme(text = element_text(size = 12))) 

print(ggplot(q_results_thresh, aes(x=mean_rand_area_tdp_ocr_miny, y=observed_area_tdp_ocr_miny, color=TDP.full)) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0, lty=2, col="red")  +
  ggtitle("PFOCR-miny") +
  theme(text = element_text(size = 12))) 

print(ggplot(q_results_thresh, aes(x=-mean_rand_area_tdp_ocr_miny+observed_area_tdp_ocr_miny, color=TDP.full)) + 
        geom_histogram() +
        geom_vline(xintercept = 0, lty=2, col="red")  +
        ggtitle("PFOCR-miny") +
        xlab("TDP: Observed - Random") +
        xlim(-0.2,0.2) +
        theme(text = element_text(size = 12))) 

print(ggplot(q_results_thresh, aes(x=-mean_rand_area_tdp_ocr_miny+observed_area_tdp_ocr_miny, color=TDP.full)) + 
        stat_ecdf() +
        geom_vline(xintercept = 0, lty=2, col="red")  +
        ggtitle("PFOCR-miny") +
        xlab("TDP: Observed - Random") +
        xlim(-0.2,0.2) +
        theme(text = element_text(size = 12))) 
print(ggplot(q_results_thresh, aes(x=diff_ocr_miny, color=TDP.full)) + 
        geom_histogram() +
        geom_vline(xintercept = 0, lty=2, col="red")  +
        ggtitle("PFOCR-miny") +
        xlab("Nsig: Observed - Random") +
        theme(text = element_text(size = 12))) 

print(ggplot(q_results_thresh, aes(x=diff_ocr_miny, color=TDP.full)) + 
        stat_ecdf() +
        geom_vline(xintercept = 0, lty=2, col="red")  +
        ggtitle("PFOCR-miny") +
        xlab("Nsig: Observed - Random") +
        theme(text = element_text(size = 12))) 

print(ggplot(q_results_thresh, aes(x=mean_rand_area_tdp_go, y=observed_area_tdp_go, color=TDP.full)) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0, lty=2, col="red")  +
  ggtitle("GO") +
  theme(text = element_text(size = 12))) 

print(ggplot(q_results_thresh, aes(x=-mean_rand_area_tdp_go+observed_area_tdp_go, color=TDP.full)) + 
        geom_histogram() +
        geom_vline(xintercept = 0, lty=2, col="red")  +
        ggtitle("GO") +
        xlab("Observed - Random") +
        xlim(-0.2,0.2) +
        theme(text = element_text(size = 12))) 

print(ggplot(q_results_thresh, aes(x=-mean_rand_area_tdp_go+observed_area_tdp_go, color=TDP.full)) + 
        stat_ecdf() +
        geom_vline(xintercept = 0, lty=2, col="red")  +
        ggtitle("GO") +
        xlab("Observed - Random") +
        xlim(-0.2,0.2) +
        theme(text = element_text(size = 12))) 
print(ggplot(q_results_thresh, aes(x=diff_go, color=TDP.full)) + 
        geom_histogram() +
        geom_vline(xintercept = 0, lty=2, col="red")  +
        ggtitle("GO") +
        xlab("Nsig: Observed - Random") +
        theme(text = element_text(size = 12))) 

print(ggplot(q_results_thresh, aes(x=diff_go, color=TDP.full)) + 
        stat_ecdf() +
        geom_vline(xintercept = 0, lty=2, col="red")  +
        ggtitle("GO") +
        xlab("Nsig: Observed - Random") +
        theme(text = element_text(size = 12))) 

dev.off()

ggplot(q_results, aes(x=diff_wp, y=diff_go)) +
  geom_point() +
  geom_smooth(method = "lm", se=FALSE) +
  geom_hline(yintercept = 0, lty=2) +
  geom_vline(xintercept = 0, lty=2) +
  xlim(-100,100) +
  ylim(-100,100) +
  geom_abline(slope = 1, intercept = 0, lty=2, col="red")

ggplot(q_results, aes(x=diff_wp, y=diff_ocr)) +
  geom_point() +
  geom_smooth(method = "lm", se=FALSE) +
  geom_hline(yintercept = 0, lty=2) +
  geom_vline(xintercept = 0, lty=2) +
  xlim(-100,100) +
  ylim(-1000,1000) +
  geom_abline(slope = 1, intercept = 0, lty=2, col="red")

ggplot(q_results, aes(x=diff_go, y=diff_ocr)) +
  geom_point() +
  geom_smooth(method = "lm", se=FALSE) +
  geom_hline(yintercept = 0, lty=2) +
  geom_vline(xintercept = 0, lty=2) +
  xlim(-100,100) +
  ylim(-10000,10000) +
  geom_abline(slope = 1, intercept = 0, lty=2, col="red")

ggplot(q_results, aes(x=z_wp, y=z_go)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, lty=2, col="red") +
  geom_hline(yintercept = 0, lty=2) +
  geom_vline(xintercept = 0, lty=2)

ggplot(q_results, aes(x=z_wp, y=z_ocr)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, lty=2, col="red") +
  geom_hline(yintercept = 0, lty=2) +
  geom_vline(xintercept = 0, lty=2)

ggplot(q_results, aes(x=z_go, y=z_ocr)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, lty=2, col="red") +
  geom_hline(yintercept = 0, lty=2) +
  geom_vline(xintercept = 0, lty=2)

z_results <- q_results_thresh %>%
  select(z_ocr, z_ocr_miny, z_wp, z_go) %>%
  gather(Database, Zstat, z_ocr:z_go)


ggplot(z_results, aes(x=Zstat, color=Database)) + 
  geom_density(adjust=5) + 
  xlim(-10,10) +
  geom_vline(xintercept = c(0,2), lty=2)
pdf(paste0("pathway_databases_w_pfocr_3sets_Zstat_random_pathways_comparison_tdp_thresh_level_", tdp_threshL,"_percent.pdf"))
ggplot(z_results, aes(x=Zstat, color=Database)) + 
  stat_ecdf() + 
  xlim(-10,10) +
  geom_vline(xintercept = c(0,2), lty=2) +
  theme(text = element_text(size = 15))
dev.off()

z_results <- q_results_thresh %>%
  select(z_ocr, z_ocr_miny, z_wp, z_go)
for(i in 1:(4 - 1)) {
  for(j in (i+1):4) {
    databases <- names(z_results)
    x <- databases[i]
    y <- databases[j]
    tempRes <- t.test(z_results[,i], z_results[,j], paired = T)
    tempRes_w <- wilcox.test(z_results[,i], z_results[,j], paired = T)
    temp_diff <- median((z_results[,i] - z_results[,j]), na.rm = T)
    print(ggplot(z_results, aes(x=z_results[,i], y=z_results[,j])) 
          + geom_point() 
          + geom_abline(slope = 1, intercept = 0, lty=2, col="red") 
          + xlab(x) + ylab(y) 
          + geom_hline(yintercept = 0, lty=2) 
          + geom_vline(xintercept = 0, lty=2)
          + theme(text = element_text(size = 12))  
          + ggtitle(paste0("Mean difference = ", round(1e4*temp_diff)/1e4, ";", "pvalue = ", tempRes_w$p.value)))
    
    print(ggplot(z_results, aes(x=z_results[,i]-z_results[,j])) 
          + geom_histogram() 
          + geom_vline(xintercept = 0, lty=2, col="red") 
          + xlab(paste0(x, "-", y)) 
          + theme(text = element_text(size = 12))  
          + ggtitle(paste0("Mean difference = ", round(1e4*temp_diff)/1e4, ";", "pvalue = ", tempRes_w$p.value)))
    
    print(ggplot(z_results, aes(x=z_results[,i]-z_results[,j])) 
          + stat_ecdf() 
          + geom_vline(xintercept = 0, lty=2, col="red") 
          + xlab(paste0(x, "-", y)) 
          + theme(text = element_text(size = 12))  
          + ggtitle(paste0("Mean difference = ", round(1e4*temp_diff)/1e4, ";", "pvalue = ", tempRes_w$p.value)))
    
  }
}

require(ggplot2)
# q_results %<>% filter(median_centered_tdp_ocr_3sets > 0 & median_centered_tdp_wp > 0 & median_centered_tdp_wp > 0)
pdf(paste0("pathway_databases_w_pfocr_3sets_miny_centered_random_pathways_comparison_tdp_thresh_level_", tdp_threshL,".pdf"))
for(i in 1:(4 - 1)) {
  for(j in (i+1):4) {
    databases <- names(q_results)
    x <- databases[i]
    y <- databases[j]
    tempRes <- t.test(q_results_thresh[,i], q_results_thresh[,j], paired = T)
    tempRes_w <- wilcox.test(q_results_thresh[,i], q_results_thresh[,j], paired = T)
    temp_diff <- median((q_results_thresh[,i] - q_results_thresh[,j]), na.rm = T)
    print(ggplot(q_results_thresh, aes(x=q_results_thresh[,i], y=q_results_thresh[,j])) 
          + geom_point() 
          + geom_abline(slope = 1, intercept = 0, lty=2, col="red") 
          + xlab(x) + ylab(y) 
          + geom_hline(yintercept = 0, lty=2) 
          + geom_vline(xintercept = 0, lty=2)
          + theme(text = element_text(size = 12))  
          + ggtitle(paste0("Mean difference = ", round(1e4*temp_diff)/1e4, ";", "pvalue = ", tempRes_w$p.value)))

    print(ggplot(q_results_thresh, aes(x=q_results_thresh[,i]-q_results_thresh[,j])) 
          + geom_histogram() 
          + geom_vline(xintercept = 0, lty=2, col="red") 
          + xlab(paste0(x,"-" ,y) )
          + xlim(-40,30)
          + theme(text = element_text(size = 12))  
          + ggtitle(paste0("Mean difference = ", round(1e4*temp_diff)/1e4, ";", "pvalue = ", tempRes_w$p.value)))
    print(ggplot(q_results_thresh, aes(x=q_results_thresh[,i]-q_results_thresh[,j])) 
          + stat_ecdf() 
          + geom_vline(xintercept = 0, lty=2, col="red") 
          + xlab(paste0(x,"-" ,y) )
          + xlim(-40,30)
          + theme(text = element_text(size = 12))  
          + ggtitle(paste0("Mean difference = ", round(1e4*temp_diff)/1e4, ";", "pvalue = ", tempRes_w$p.value)))
    
  }
}
dev.off()

pdf(paste0("pathway_databases_w_pfocr_3sets_raw_comparison_tdp_thresh_level_", tdp_threshL,".pdf"))
for(i in 1:(4 - 1)) {
  for(j in (i+1):4) {
    databases <- names(q_results_thresh)
    x <- databases[(12+i)]
    y <- databases[(12+j)]
    tempRes <- t.test(q_results_thresh[,(12+i)], q_results_thresh[,(12+j)], paired = T)
    print(ggplot(q_results_thresh, aes(x=q_results_thresh[,(12+i)], y=q_results_thresh[,(12+j)])) 
          + geom_point() 
          + geom_abline(slope = 1, intercept = 0, lty=2, col="red") 
          + xlab(x) + ylab(y) 
          + geom_hline(yintercept = 0, lty=2) 
          + geom_vline(xintercept = 0, lty=2)
          + theme(text = element_text(size = 12))  
          + ggtitle(paste0("Mean difference = ", round(1e4*tempRes$estimate)/1e4, ";", "pvalue = ", round(1e4*tempRes$p.value)/1e4)))
  }
}
dev.off()

get_pfocr_TDPQuantiles <- function(ChooseGSE, GSEid, ocr_length, go_length, wp_length, minSig, minSize, theta) {
  print(ChooseGSE)
  q_ocr_orig <- NA
  q_ocr_noalias <- NA
  q_ocr_nobe0 <- NA
  q_ocr_nobe2 <- NA
  q_ocr_nobe3 <- NA
  q_ocr_noprev <- NA
  
  permute_res <- readRDS(paste0(getwd(), ))
  ocr_result <- read.table(paste0(getwd(), "/GSE/",GSEid[ChooseGSE],"/rSEA/PFOCR/result.txt"), header = T)
  tempIndices <- match(ocr_result$set_id, names(ocr_length))
  temp_ocr_length <- ocr_length[tempIndices]
  Nocr <- sum(!is.na(ocr_result$Comp.adjP) & ocr_result$Comp.adjP < 0.05 , na.rm = T)
  if(Nocr > minSig) {
    sig_ocr_result <- data.frame(ocr_result[!is.na(ocr_result$Comp.adjP) & ocr_result$Comp.adjP < 0.05 , ], database=rep("ocr_noalias", Nocr))
    q_ocr_orig <- quantile(sig_ocr_result$TDP.estimate, theta, na.rm = T)
  }
  
  ocr_result <- read.table(paste0(getwd(), "/GSE/",GSEid[ChooseGSE],"/rSEA/PFOCR_noalias/result.txt"), header = T)
  tempIndices <- match(ocr_result$set_id, names(ocr_length))
  temp_ocr_length <- ocr_length[tempIndices]
  Nocr <- sum(!is.na(ocr_result$Comp.adjP) & ocr_result$Comp.adjP < 0.05 , na.rm = T)
  if(Nocr > minSig) {
    sig_ocr_result <- data.frame(ocr_result[!is.na(ocr_result$Comp.adjP) & ocr_result$Comp.adjP < 0.05 , ], database=rep("ocr_noalias", Nocr))
    q_ocr_noalias <- quantile(sig_ocr_result$TDP.estimate, theta, na.rm = T)
  }
  
  ocr_result <- read.table(paste0(getwd(), "/GSE/",GSEid[ChooseGSE],"/rSEA/PFOCR_nobe0/result.txt"), header = T)
  tempIndices <- match(ocr_result$set_id, names(ocr_length))
  temp_ocr_length <- ocr_length[tempIndices]
  Nocr <- sum(!is.na(ocr_result$Comp.adjP) & ocr_result$Comp.adjP < 0.05 , na.rm = T)
  if(Nocr > minSig) {
    sig_ocr_result <- data.frame(ocr_result[!is.na(ocr_result$Comp.adjP) & ocr_result$Comp.adjP < 0.05 , ], database=rep("ocr_nobe0", Nocr))
    q_ocr_nobe0 <- quantile(sig_ocr_result$TDP.estimate, theta, na.rm = T)
  }
  
  
  ocr_result <- read.table(paste0(getwd(), "/GSE/",GSEid[ChooseGSE],"/rSEA/PFOCR_nobe2/result.txt"), header = T)
  tempIndices <- match(ocr_result$set_id, names(ocr_length))
  temp_ocr_length <- ocr_length[tempIndices]
  Nocr <- sum(!is.na(ocr_result$Comp.adjP) & ocr_result$Comp.adjP < 0.05 , na.rm = T)
  if(Nocr > minSig) {
    sig_ocr_result <- data.frame(ocr_result[!is.na(ocr_result$Comp.adjP) & ocr_result$Comp.adjP < 0.05 , ], database=rep("ocr_nobe2", Nocr))
    q_ocr_nobe2 <- quantile(sig_ocr_result$TDP.estimate, theta, na.rm = T)
  }
  
  ocr_result <- read.table(paste0(getwd(), "/GSE/",GSEid[ChooseGSE],"/rSEA/PFOCR_nobe3/result.txt"), header = T)
  tempIndices <- match(ocr_result$set_id, names(ocr_length))
  temp_ocr_length <- ocr_length[tempIndices]
  Nocr <- sum(!is.na(ocr_result$Comp.adjP) & ocr_result$Comp.adjP < 0.05 , na.rm = T)
  if(Nocr > minSig) {
    sig_ocr_result <- data.frame(ocr_result[!is.na(ocr_result$Comp.adjP) & ocr_result$Comp.adjP < 0.05 , ], database=rep("ocr_nobe3", Nocr))
    q_ocr_nobe3 <- quantile(sig_ocr_result$TDP.estimate, theta, na.rm = T)
  }
  
  ocr_result <- read.table(paste0(getwd(), "/GSE/",GSEid[ChooseGSE],"/rSEA/PFOCR_noprev/result.txt"), header = T)
  tempIndices <- match(ocr_result$set_id, names(ocr_length))
  temp_ocr_length <- ocr_length[tempIndices]
  Nocr <- sum(!is.na(ocr_result$Comp.adjP) & ocr_result$Comp.adjP < 0.05 , na.rm = T)
  if(Nocr > minSig) {
    sig_ocr_result <- data.frame(ocr_result[!is.na(ocr_result$Comp.adjP) & ocr_result$Comp.adjP < 0.05 , ], database=rep("ocr_noprev", Nocr))
    q_ocr_noprev <- quantile(sig_ocr_result$TDP.estimate, theta, na.rm = T)
  }
  
  # PlotData <- rbind(sig_ocr_result, sig_wp_result)
  # PlotData <- data.frame(rbind(PlotData, sig_go_result))
  # ggplot(PlotData, aes(x=TDP.estimate, color=database)) + stat_ecdf() + geom_vline(xintercept = q_go,lty=2)
  # ggplot(PlotData, aes(y=TDP.estimate, x=database)) + geom_boxplot()
  
  return(c(q_ocr_orig, q_ocr_noalias, q_ocr_nobe0, q_ocr_nobe2, q_ocr_nobe3, q_ocr_noprev))
}

ocr_q_results <- t(sapply(1:length(GSEid), get_pfocr_TDPQuantiles, GSEid, ocr_length, go_length, wp_length, minSig, minSize, theta, quantiles))
ocr_q_results <- data.frame(ocr_q_results)
colnames(ocr_q_results) <- c("median_tdp_ocr_orig","median_tdp_ocr_noalias", "median_tdp_ocr_nobe0", "median_tdp_ocr_nobe2", "median_tdp_ocr_nobe3", "median_tdp_ocr_noprev")

pdf(paste0("ocr_databases_comparison_theta_level_", thetaL,"_percent.pdf"))
for(i in 1:(ncol(ocr_q_results) - 1)) {
  for(j in (i+1):ncol(ocr_q_results)) {
    databases <- names(ocr_q_results)
    x <- databases[i]
    y <- databases[j]
    tempRes <- t.test(ocr_q_results[,i], ocr_q_results[,j], paired = T)
    print(ggplot(ocr_q_results, aes(x=ocr_q_results[,i], y=ocr_q_results[,j])) 
          + geom_point() 
          + geom_abline(slope = 1, intercept = 0, lty=2, col="red") 
          + xlab(x) + ylab(y)
          + ggtitle(paste0("Mean difference = ", round(1e4*tempRes$estimate)/1e4, ";", "pvalue = ", round(1e4*tempRes$p.value)/1e4)))
  }
}
dev.off()
