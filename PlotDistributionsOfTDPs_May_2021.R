rm(list = ls())
setwd("~/Dropbox (Gladstone)/PFOCRInPathwayAnalyses")
require(dplyr)
require(tidyr)
require(magrittr)
require(ggplot2)

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


GSEid <- list.files(path = "./PFOCRInPathwayAnalyses/GSE/")
minSig <- 5
minSize <- 3
theta <- 0.9
Nperm <- 10
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
  as.integer()


getTDPQuantiles <- function(ChooseGSE, GSEid, ocr_length, go_length, wp_length, minSig, minSize, theta) {
  print(ChooseGSE)
  q_ocr <- NA
  q_ocr_miny <- NA
  q_wp <- NA
  q_go <- NA
  
  permute_res <- readRDS(paste0(permute_res_dir, "GSE_index_", ChooseGSE,"_", Nperm, "_perm_result.rds"))
  rand_ocr_tdp <- permute_res$rNsig_pfcor %>%
    as.data.frame() %>%
    .$quantile_90 %>%
    mean()
 
  rand_ocr_tdp_miny <- permute_res$rNsig_pfcor_miny %>%
    as.data.frame() %>%
    .$quantile_90 %>%
    mean()
  
  rand_wp_tdp <- permute_res$rNsig_wp %>%
    as.data.frame() %>%
    .$quantile_90 %>%
    mean()
  
  rand_go_tdp <- permute_res$rNsig_go %>%
    as.data.frame() %>%
    .$quantile_90 %>%
    mean()

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
    q_ocr <- permute_res$oTDP_bound_90_pfocr[3]
    q_ocr <- q_ocr - rand_ocr_tdp
  }
  diffocr <- (Nocr - rand_ocr_nsig_mean)
  zocr <- diffocr/(rand_ocr_nsig_sd+1)
  
  Nocr_miny <- permute_res$oNsig_pfocr_miny
  if(Nocr_miny > minSig) {
    # sig_ocr_result <- data.frame(ocr_result[!is.na(ocr_result$Comp.adjP) & ocr_result$Comp.adjP < 0.05 & temp_ocr_length >= minSize, ], database=rep("ocr", Nocr))
    # q_ocr <- quantile(sig_ocr_result$TDP.estimate, theta, na.rm = T)
    q_ocr_miny <- permute_res$oTDP_bound_90_pfocr_miny[3]
    q_ocr_miny <- q_ocr_miny - rand_ocr_tdp_miny
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
    q_wp <- permute_res$oTDP_bound_90_wp[3]
    q_wp <- q_wp - rand_wp_tdp
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
    q_go <- permute_res$oTDP_bound_90_go[3]
    q_go <- q_go - rand_go_tdp
  }
  diffgo <- (Ngo - rand_go_nsig_mean)
  zgo <- diffgo/(rand_go_nsig_sd+1)
  
  # PlotData <- rbind(sig_ocr_result, sig_wp_result)
  # PlotData <- data.frame(rbind(PlotData, sig_go_result))
  # ggplot(PlotData, aes(x=TDP.estimate, color=database)) + stat_ecdf() + geom_vline(xintercept = q_go,lty=2)
  # ggplot(PlotData, aes(y=TDP.estimate, x=database)) + geom_boxplot()
  
  return(c(q_ocr, q_ocr_miny,q_wp, q_go, zocr, zocr_miny,zwp, zgo, diffocr, diffocr_miny, diffwp, diffgo, permute_res$oTDP_bound_90_pfocr[3], permute_res$oTDP_bound_90_pfocr_miny[3],permute_res$oTDP_bound_90_wp[3], permute_res$oTDP_bound_90_go[3]))
}

q_results <- t(sapply(GSEIndices, getTDPQuantiles, GSEid, ocr_length, go_length, wp_length, minSig, minSize, theta))
q_results <- data.frame(q_results)
colnames(q_results) <- c("median_centered_tdp_ocr_3sets", "median_centered_tdp_ocr_miny", "median_centered_tdp_wp", "median_centered_tdp_go", "z_ocr", "z_ocr_miny","z_wp", "z_go", "diff_ocr", "diff_ocr_miny","diff_wp", "diff_go","median_raw_tdp_ocr_3sets", "median_raw_tdp_ocr_miny", "median_raw_tdp_wp", "median_raw_tdp_go")

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

z_results <- q_results %>%
  select(z_ocr, z_ocr_miny, z_wp, z_go) %>%
  gather(Database, Zstat, z_ocr:z_go)

ggplot(z_results, aes(x=Zstat, color=Database)) + 
  geom_density(adjust=5) + 
  xlim(-10,10) +
  geom_vline(xintercept = c(0,2), lty=2)
pdf(paste0("pathway_databases_w_pfocr_3sets_Zstat_random_pathways_comparison_theta_level_", thetaL,"_percent.pdf"))
ggplot(z_results, aes(x=Zstat, color=Database)) + 
  stat_ecdf() + 
  xlim(-10,10) +
  geom_vline(xintercept = c(0,2), lty=2) +
  theme(text = element_text(size = 15))
dev.off()

require(ggplot2)
# q_results %<>% filter(median_centered_tdp_ocr_3sets > 0 & median_centered_tdp_wp > 0 & median_centered_tdp_wp > 0)
pdf(paste0("pathway_databases_w_pfocr_3sets_miny_centered_random_pathways_comparison_theta_level_", thetaL,"_percent.pdf"))
for(i in 1:(4 - 1)) {
  for(j in (i+1):4) {
    databases <- names(q_results)
    x <- databases[i]
    y <- databases[j]
    tempRes <- t.test(q_results[,i], q_results[,j], paired = T)
    print(ggplot(q_results, aes(x=q_results[,i], y=q_results[,j])) 
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

pdf(paste0("pathway_databases_w_pfocr_3sets_raw_comparison_theta_level_", thetaL,"_percent.pdf"))
for(i in 1:(3 - 1)) {
  for(j in (i+1):3) {
    databases <- names(q_results)
    x <- databases[(9+i)]
    y <- databases[(9+j)]
    tempRes <- t.test(q_results[,(9+i)], q_results[,(9+j)], paired = T)
    print(ggplot(q_results, aes(x=q_results[,(9+i)], y=q_results[,(9+j)])) 
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

ocr_q_results <- t(sapply(1:length(GSEid), get_pfocr_TDPQuantiles, GSEid, ocr_length, go_length, wp_length, minSig, minSize, theta))
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
