# scriptDir=/wynton/group/gladstone/biocore/projects/pfocr_pathway_enrichment_evaluation/PFOCRInPathwayAnalyses/
# containerDir=/wynton/group/gladstone/biocore/containers
# dataDir=/wynton/group/gladstone/biocore/projects/PFOCR/
# export SINGULARITY_BINDPATH="$containerDir,$scriptDir,$dataDir"
# singularity exec $containerDir/pathway_enrichment_pathway_databases_evaluation_latest.sif R

##################################required libraries
require("rhdf5")
require("tools")
require('edgeR')
require(GEOquery)
require(rSEA)
require("clusterProfiler")
require("org.Hs.eg.db")
require(dplyr)
require(GO.db)
require(bayesbio)
require(readxl)
require(pROC)
args <- commandArgs(trailingOnly=TRUE)

GSE_index <- as.integer(args[1])
print(GSE_index)
outdir <- "/wynton/group/gladstone/biocore/projects/pfocr_pathway_enrichment_evaluation/permuted_geneset_databases_results_updated_pfocr_gene_set_sizes_may_2021/"

min_set_size <- 3
max_set_size <- 500 
dataDir <- "/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/"
database_lists1 <- load("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/databases_pfocr_3intersect_v2.RData")#pfocr_3sets
database_lists2 <- load("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/databases.RData")#has wp, pfocr, go
database_lists3 <- load("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/pfocr_miny_April2021.RData")#has wp, pfocr, go

database_lists <- append(database_lists1, database_lists2)
database_lists <- append(database_lists, database_lists3)

database_lists <- unname(unlist(sapply(database_lists, grep, pattern="_list$", value = T, perl = T)))
for (db in database_lists) {
  eval(call("<-", as.name(db),  Filter(Negate(is.null), lapply(get(db), function(x){
    if(length(x) < min_set_size | length(x) > max_set_size)
      NULL
    else
      x
  }))
  ))
}



##################################required data sets
pvalue_results_human_voom=readRDS("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/pvalue_results_human_voom.rds")
logFC_results_human_voom=readRDS("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/logFC_results_human_voom.rds")
merged_filtered_gses=readRDS("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/merged_filtered_gses.rds")

load("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/databases_pfocr_3sets.RData")#has wp, pfocr, go
gene_entrez=readRDS("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/gene_entrez.rds")

gene_entrez=readRDS("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/gene_entrez.rds")

##################################required functions
source("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/run_rSEA3.r")
source("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/run_gsea2.r")
source("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/run_ORA3.r")

rsea_results_human_voom_pfocr <- vector("list", 9)
names(rsea_results_human_voom_pfocr) <- c("set_id", "ID", "Coverage", "TDP.bound","TDP.estimate", "SC.adjP", "Comp.adjP")

rsea_results_human_voom_go <- vector("list", 9)
names(rsea_results_human_voom_go) <- c("set_id", "ID",  "Coverage", "TDP.bound","TDP.estimate", "SC.adjP", "Comp.adjP")

rsea_results_human_voom_wp <- vector("list", 9)
names(rsea_results_human_voom_wp) <- c("set_id", "ID", "Coverage", "TDP.bound","TDP.estimate", "SC.adjP", "Comp.adjP")

set.seed(1234)

Nperm <- 10
qLevels <- c(0.5, 0.75, 0.9, 0.95, 0.99)
names_qLevels <- paste("quantile", round(100*qLevels), sep="_")
PermuteDatabase <- function(p, GSE_index, path_list, path_annotation,pvalue_results_human_voom, gene_entrez) {
  #print(p)
  temp_list <- path_list
  temp_annotation <- path_annotation
  temp_annotation$gene <- path_annotation$gene[sample(1:nrow(path_annotation), nrow(path_annotation), replace = FALSE)]
  for(index in 1:length(temp_list)) {
    temp_list[[index]] <- temp_annotation %>%
      subset(set_id == names(temp_list)[index]) %>%
      pull(gene) %>%
      as.character()
  }
  
  data_m=data.frame(rownames(pvalue_results_human_voom),pvalue_results_human_voom[,GSE_index])
  colnames(data_m)=c("Gene","pvalue")
  merged=merge(data_m,gene_entrez,by="Gene",all = "L",all.x=TRUE,all.y=FALSE)
  merged <- merged[!is.na(merged$pvalue) & !is.na(merged$ENTREZID),]

  res <- SEA(merged$pvalue, merged$ENTREZID, pathlist = temp_list)
  # res=apply(as.matrix(GSE_index),1,function(x){run_rSEA3(as.matrix(pvalue_results_human_voom[,x]),rsea_results_human_voom_go,temp_list,temp_annotation)})
  Nsig <- sum(res$Comp.adjP < 0.05, na.rm = TRUE)
  TempRes <- data.frame(Comp.adjP=res$Comp.adjP, TDP.bound=res$TDP.bound)
  TempRes %<>% filter(Comp.adjP < 0.05)
  TDP_bound_90 <- TempRes %>%
  	.$TDP.bound %>%
	quantile(.,c(0.5, 0.75, 0.9, 0.95, 0.99), na.rm=TRUE)
  Res <- append(Nsig, TDP_bound_90)
  names(Res) <- append("Nsig", names_qLevels)
  print(Res)
  return(Res)
}

data_m=data.frame(rownames(pvalue_results_human_voom),pvalue_results_human_voom[,GSE_index])
colnames(data_m)=c("Gene","pvalue")
merged=merge(data_m,gene_entrez,by="Gene",all = "L",all.x=TRUE,all.y=FALSE)
merged <- merged[!is.na(merged$pvalue) & !is.na(merged$ENTREZID),]


TDPestimate_full <- setTDP(merged$pvalue, merged$ENTREZID, alpha = 0.05)$TDP.estimate
TDPbound_full <- setTDP(merged$pvalue, merged$ENTREZID, alpha = 0.05)$TDP.bound



# rsea_results_human_voom_wp=apply(as.matrix(GSE_index),1,function(x){run_rSEA3(as.matrix(pvalue_results_human_voom[,x]),rsea_results_human_voom_wp,wp_list,wp_annotation)})
rsea_results_human_voom_wp=SEA(merged$pvalue, merged$ENTREZID, pathlist = wp_list)
TempRes <- data.frame(Comp.adjP=rsea_results_human_voom_wp$Comp.adjP, TDP.bound=rsea_results_human_voom_wp$TDP.bound)
TempRes %<>% filter(Comp.adjP < 0.05)
TDP_bound_90 <- TempRes %>%
  .$TDP.bound %>%
  quantile(.,c(0.5, 0.75, 0.9, 0.95, 0.99), na.rm=TRUE)
oTDP_bound_90_wp <- TDP_bound_90
oNsig_wp <- sum(rsea_results_human_voom_wp$Comp.adjP < 0.05, na.rm = TRUE)
oTDP_bound_90_wp <- TDP_bound_90
print(oNsig_wp)
rNsig_wp <- t(sapply(1:Nperm, PermuteDatabase, GSE_index, wp_list, wp_annotation,pvalue_results_human_voom, gene_entrez))


#rsea_results_human_voom_go=apply(as.matrix(GSE_index),1,function(x){run_rSEA3(as.matrix(pvalue_results_human_voom[,x]),rsea_results_human_voom_go,go_list,go_annotation)})
rsea_results_human_voom_go=SEA(merged$pvalue, merged$ENTREZID, pathlist = go_list)
TempRes <- data.frame(Comp.adjP=rsea_results_human_voom_go$Comp.adjP, TDP.bound=rsea_results_human_voom_go$TDP.bound)
TempRes %<>% filter(Comp.adjP < 0.05)
TDP_bound_90 <- TempRes %>%
  .$TDP.bound %>%
  quantile(.,c(0.5, 0.75, 0.9, 0.95, 0.99), na.rm=TRUE)
oTDP_bound_90_go <- TDP_bound_90

oNsig_go <- sum(rsea_results_human_voom_go$Comp.adjP < 0.05, na.rm = TRUE)
print(oNsig_go)
rNsig_go <- t(sapply(1:Nperm, PermuteDatabase, GSE_index, go_list, go_annotation,pvalue_results_human_voom, gene_entrez))

# rsea_results_human_voom_pfocr=apply(as.matrix(GSE_index),1,function(x){run_rSEA3(as.matrix(pvalue_results_human_voom[,x]),rsea_results_human_voom_pfocr,pfocr_3sets_list,pfocr_annotation_3sets)})
rsea_results_human_voom_pfocr=SEA(merged$pvalue, merged$ENTREZID, pathlist = pfocr_3sets_list)
TempRes <- data.frame(Comp.adjP=rsea_results_human_voom_pfocr$Comp.adjP, TDP.bound=rsea_results_human_voom_pfocr$TDP.bound)
TempRes %<>% filter(Comp.adjP < 0.05)
TDP_bound_90 <- TempRes %>%
  .$TDP.bound %>%
  quantile(.,c(0.5, 0.75, 0.9, 0.95, 0.99), na.rm=TRUE)
oTDP_bound_90_pfocr <- TDP_bound_90

oNsig_pfocr <- sum(rsea_results_human_voom_pfocr$Comp.adjP < 0.05, na.rm = TRUE)
print(oNsig_pfocr)
rNsig_pfcor <- t(sapply(1:Nperm, PermuteDatabase, GSE_index, pfocr_3sets_list, pfocr_annotation_3sets,pvalue_results_human_voom, gene_entrez))


rsea_results_human_voom_pfocr_miny=SEA(merged$pvalue, merged$ENTREZID, pathlist = pfocr_miny_list)
TempRes <- data.frame(Comp.adjP=rsea_results_human_voom_pfocr_miny$Comp.adjP, TDP.bound=rsea_results_human_voom_pfocr_miny$TDP.bound)
TempRes %<>% filter(Comp.adjP < 0.05)
TDP_bound_90 <- TempRes %>%
  .$TDP.bound %>%
  quantile(.,c(0.5, 0.75, 0.9, 0.95, 0.99), na.rm=TRUE)
oTDP_bound_90_pfocr_miny <- TDP_bound_90

oNsig_pfocr_miny <- sum(rsea_results_human_voom_pfocr_miny$Comp.adjP < 0.05, na.rm = TRUE)
print(oNsig_pfocr_miny)
rNsig_pfcor_miny <- t(sapply(1:Nperm, PermuteDatabase, GSE_index, pfocr_miny_list, pfocr_annotation_miny,pvalue_results_human_voom, gene_entrez))

GSE_index_2_100_perm_res <- list(TDPbound_full=TDPbound_full, oNsig_wp=oNsig_wp, rNsig_wp=rNsig_wp, oNsig_go=oNsig_go, rNsig_go=rNsig_go, oNsig_pfocr=oNsig_pfocr, rNsig_pfcor=rNsig_pfcor, oNsig_pfocr_miny=oNsig_pfocr_miny, rNsig_pfcor_miny=rNsig_pfcor_miny, oTDP_bound_90_wp=oTDP_bound_90_wp, oTDP_bound_90_go=oTDP_bound_90_go, oTDP_bound_90_pfocr=oTDP_bound_90_pfocr, , oTDP_bound_90_pfocr_miny=oTDP_bound_90_pfocr_miny)
saveRDS(GSE_index_2_100_perm_res, file=paste0(outdir, "GSE_index_",GSE_index,"_",Nperm,"_perm_result.rds"))

