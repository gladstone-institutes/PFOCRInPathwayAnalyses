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
require(dplyr)
require(magrittr)

args <- commandArgs(trailingOnly=TRUE)

GSE_index <- as.integer(args[1])
print(GSE_index)
outdir <- "/wynton/group/gladstone/biocore/projects/pfocr_pathway_enrichment_evaluation/permuted_geneset_databases_gsea_results_july_2021/"

min_set_size <- 3
max_set_size <- 500 
dataDir <- "/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/"
database_lists1 <- load("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/databases_pfocr_3intersect_v2.RData")#pfocr_3sets
database_lists2 <- load("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/databases.RData")#has wp, pfocr, go
database_lists3 <- load("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/pfocr_miny_April2021.RData")#has wp, pfocr, go

database_lists <- append(database_lists1, database_lists2)
database_lists <- append(database_lists, database_lists3)

database_annotation <- unname(unlist(sapply(database_lists, grep, pattern="annotation", value = T, perl = T)))
database_lists <- unname(unlist(sapply(database_lists, grep, pattern="list$", value = T, perl = T)))
for (db in database_lists) {
  eval(call("<-", as.name(db),  Filter(Negate(is.null), lapply(get(db), function(x){
    if(length(x) < min_set_size | length(x) > max_set_size)
      NULL
    else
      x
  }))
  ))
}

##need to filter annotation data as well
wp_annotation %<>% filter(set_id %in% names(wp_list))
go_annotation %<>% filter(set_id %in% names(go_list))
pfocr_annotation %<>% filter(set_id %in% names(pfocr_list))
pfocr_annotation_3sets %<>% filter(set_id %in% names(pfocr_3sets_list))
pfocr_miny_annotation %<>% filter(set_id %in% names(pfocr_miny_list))

##################################required data sets
pvalue_results_human_voom=readRDS("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/pvalue_results_human_voom.rds")
logFC_results_human_voom=readRDS("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/logFC_results_human_voom.rds")
merged_filtered_gses=readRDS("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/merged_filtered_gses.rds")

load("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/databases_pfocr_3sets.RData")#has wp, pfocr, go
gene_entrez=readRDS("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/gene_entrez.rds")

gene_entrez=readRDS("/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses/gene_entrez.rds")


set.seed(1234)

Nperm <- 100
qLevels <- c(0.5, 0.75, 0.9, 0.95, 0.99)
names_qLevels <- paste("quantile", round(100*qLevels), sep="_")


choose_random_genes_in_path <- function(path, path_annotation, gene_prob) {
  
  random_genes <- sample(names(gene_prob), size=length(path), replace = FALSE, prob = gene_prob)
  
  return(random_genes)
}

PermuteDatabase <- function(p, GSE_index, path_list, path_annotation,pvalue_results_human_voom, gene_entrez,path_gene_prob, logFC_results_human_voom) {
  #print(p)
  temp_annotation <- path_annotation
  
  temp_list <- lapply(path_list, choose_random_genes_in_path, path_annotation, path_gene_prob)
  temp_annotation$gene <- unlist(temp_list)
  
  data_m=data.frame(rownames(pvalue_results_human_voom),pvalue_results_human_voom[,GSE_index])
  colnames(data_m)=c("Gene","pvalue")
  merged=merge(data_m,gene_entrez,by="Gene",all = "L",all.x=TRUE,all.y=FALSE)
  merged <- merged[!is.na(merged$pvalue) & !is.na(merged$ENTREZID),]
  
  logFC_data=data.frame(rownames(logFC_results_human_voom),logFC_results_human_voom[,GSE_index])
  ### geneList prep
  colnames(logFC_data) <-"FC"
  logFC_data <- add_rownames(as.data.frame(logFC_data), "Gene")
  gene_list <- merged %>%
    left_join(logFC_data, by="Gene") %>%
    mutate(Score = sign(as.numeric(FC)) * - log10(as.numeric(as.character(pvalue)))) %>%
    dplyr::select(c("Score","ENTREZID")) %>%
    drop_na(ENTREZID) %>%
    dplyr::arrange(desc(Score))
  gene_list <- unlist(split(gene_list[, 1], gene_list[, 2]))
  
  gene_list = sort(gene_list[unique(names(gene_list))], decreasing = TRUE)
  
  res <- GSEA(
    gene_list,
    pAdjustMethod="BH",
    TERM2GENE = temp_annotation[,c("set_id","gene")],
    TERM2NAME = temp_annotation[,c("set_id","set_id")]    ,
    minGSSize = 1,
    maxGSSize = 100000,
    pvalueCutoff = 1,
    verbose=FALSE)
  # res=apply(as.matrix(GSE_index),1,function(x){run_rSEA3(as.matrix(pvalue_results_human_voom[,x]),rsea_results_human_voom_go,temp_list,temp_annotation)})
  Nsig <- sum(res$qvalues < 0.05, na.rm = TRUE)
  TempRes <- data.frame(qvalues=res$qvalues, aNES=abs(res$NES))
  TempRes %<>% filter(qvalues < 0.05)
  TDP_bound_90 <- TempRes %>%
    .$aNES %>%
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

logFC_data=data.frame(rownames(logFC_results_human_voom),logFC_results_human_voom[,GSE_index])
### geneList prep
colnames(logFC_data) <-"FC"
logFC_data <- add_rownames(as.data.frame(logFC_data), "Gene")
gene_list <- merged %>%
  left_join(logFC_data, by="Gene") %>%
  mutate(Score = sign(as.numeric(FC)) * - log10(as.numeric(as.character(pvalue)))) %>%
  dplyr::select(c("Score","ENTREZID")) %>%
  drop_na(ENTREZID) %>%
  dplyr::arrange(desc(Score))
gene_list <- unlist(split(gene_list[, 1], gene_list[, 2]))

gene_list = sort(gene_list[unique(names(gene_list))], decreasing = TRUE)

##WP
rsea_results_human_voom_wp <- GSEA(
  gene_list,
  pAdjustMethod="BH",
  TERM2GENE = wp_annotation[,c("set_id","gene")],
  TERM2NAME = wp_annotation[,c("set_id","set_id")]    ,
  minGSSize = 1,
  maxGSSize = 100000,
  pvalueCutoff = 1,
  verbose=FALSE)

TempRes <- data.frame(qvalues=rsea_results_human_voom_wp$qvalues, aNES=rsea_results_human_voom_wp$NES)
TempRes %<>% filter(qvalues < 0.05)
TDP_bound_90 <- TempRes %>%
  .$aNES %>%
  quantile(.,c(0.5, 0.75, 0.9, 0.95, 0.99), na.rm=TRUE)
oTDP_bound_90_wp <- TDP_bound_90
oNsig_wp <- sum(rsea_results_human_voom_wp$qvalues < 0.05, na.rm = TRUE)
oTDP_bound_90_wp <- TDP_bound_90
print(oNsig_wp)

unique_genes <- wp_annotation %>%
  .$gene %>%
  unique()

wp_gene_prob <- wp_annotation %>%
  .$gene %>%
  table() %>%
  divide_by(length(unique_genes))

rNsig_wp <- t(sapply(1:Nperm, PermuteDatabase, GSE_index, wp_list, wp_annotation,pvalue_results_human_voom, gene_entrez, wp_gene_prob, logFC_results_human_voom))

##GO
rsea_results_human_voom_go <- GSEA(
  gene_list,
  pAdjustMethod="BH",
  TERM2GENE = go_annotation[,c("set_id","gene")],
  TERM2NAME = go_annotation[,c("set_id","set_id")]    ,
  minGSSize = 1,
  maxGSSize = 100000,
  pvalueCutoff = 1,
  verbose=FALSE)

TempRes <- data.frame(qvalues=rsea_results_human_voom_go$qvalues, aNES=abs(rsea_results_human_voom_go$NES))
TempRes %<>% filter(qvalues < 0.05)
TDP_bound_90 <- TempRes %>%
  .$aNES %>%
  quantile(.,c(0.5, 0.75, 0.9, 0.95, 0.99), na.rm=TRUE)
oTDP_bound_90_go <- TDP_bound_90

oNsig_go <- sum(rsea_results_human_voom_go$qvalues < 0.05, na.rm = TRUE)
print(oNsig_go)

unique_genes <- go_annotation %>%
  .$gene %>%
  unique()

go_gene_prob <- go_annotation %>%
  .$gene %>%
  table() %>%
  divide_by(length(unique_genes))

rNsig_go <- t(sapply(1:Nperm, PermuteDatabase, GSE_index, go_list, go_annotation,pvalue_results_human_voom, gene_entrez, go_gene_prob, logFC_results_human_voom))

##PFCOR-3sets
rsea_results_human_voom_pfocr <- GSEA(
  gene_list,
  pAdjustMethod="BH",
  TERM2GENE = pfocr_annotation_3sets[,c("set_id","gene")],
  TERM2NAME = pfocr_annotation_3sets[,c("set_id","set_id")]    ,
  minGSSize = 1,
  maxGSSize = 100000,
  pvalueCutoff = 1,
  verbose=FALSE)

TempRes <- data.frame(qvalues=rsea_results_human_voom_pfocr$qvalues, aNES=abs(rsea_results_human_voom_pfocr$NES))
TempRes %<>% filter(qvalues < 0.05)
TDP_bound_90 <- TempRes %>%
  .$aNES %>%
  quantile(.,c(0.5, 0.75, 0.9, 0.95, 0.99), na.rm=TRUE)
oTDP_bound_90_pfocr <- TDP_bound_90

unique_genes <- pfocr_annotation_3sets %>%
  .$gene %>%
  unique()

pfocr_3sets_gene_prob <- pfocr_annotation_3sets %>%
  .$gene %>%
  table() %>%
  divide_by(length(unique_genes))

oNsig_pfocr <- sum(rsea_results_human_voom_pfocr$qvalues < 0.05, na.rm = TRUE)
print(oNsig_pfocr)
rNsig_pfcor <- t(sapply(1:Nperm, PermuteDatabase, GSE_index, pfocr_3sets_list, pfocr_annotation_3sets,pvalue_results_human_voom, gene_entrez, pfocr_3sets_gene_prob, logFC_results_human_voom))

##PFCOR-miny
rsea_results_human_voom_pfocr_miny <- GSEA(
  gene_list,
  pAdjustMethod="BH",
  TERM2GENE = pfocr_miny_annotation[,c("set_id","gene")],
  TERM2NAME = pfocr_miny_annotation[,c("set_id","set_id")]    ,
  minGSSize = 1,
  maxGSSize = 100000,
  pvalueCutoff = 1,
  verbose=FALSE)

TempRes <- data.frame(qvalues=rsea_results_human_voom_pfocr_miny$qvalues, aNES=abs(rsea_results_human_voom_pfocr_miny$NES))
TempRes %<>% filter(qvalues < 0.05)
TDP_bound_90 <- TempRes %>%
  .$aNES %>%
  quantile(.,c(0.5, 0.75, 0.9, 0.95, 0.99), na.rm=TRUE)
oTDP_bound_90_pfocr_miny <- TDP_bound_90

oNsig_pfocr_miny <- sum(rsea_results_human_voom_pfocr_miny$qvalues < 0.05, na.rm = TRUE)
print(oNsig_pfocr_miny)

unique_genes <- pfocr_miny_annotation %>%
  .$gene %>%
  unique()

pfocr_miny_gene_prob <- pfocr_miny_annotation %>%
  .$gene %>%
  table() %>%
  divide_by(length(unique_genes))

rNsig_pfcor_miny <- t(sapply(1:Nperm, PermuteDatabase, GSE_index, pfocr_miny_list, pfocr_miny_annotation,pvalue_results_human_voom, gene_entrez, pfocr_miny_gene_prob, logFC_results_human_voom))

GSE_index_2_100_perm_res <- list(oNsig_wp=oNsig_wp, rNsig_wp=rNsig_wp, oNsig_go=oNsig_go, rNsig_go=rNsig_go, oNsig_pfocr=oNsig_pfocr, rNsig_pfcor=rNsig_pfcor, oNsig_pfocr_miny=oNsig_pfocr_miny, rNsig_pfcor_miny=rNsig_pfcor_miny, oTDP_bound_90_wp=oTDP_bound_90_wp, oTDP_bound_90_go=oTDP_bound_90_go, oTDP_bound_90_pfocr=oTDP_bound_90_pfocr, oTDP_bound_90_pfocr_miny=oTDP_bound_90_pfocr_miny)
saveRDS(GSE_index_2_100_perm_res, file=paste0(outdir, "GSE_index_",GSE_index,"_",Nperm,"_gsea_update_perm_result.rds"))

