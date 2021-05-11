##########################################################################################knockout benchmark
require(pROC)

setwd("/Users/mingyoungshin/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/PFOCRInPathwayAnalyses/")
pvalue_results_human_voom<-readRDS("pvalue_results_human_voom.rds")

library(dplyr)
library(magrittr)

rsea_results_human_voom_go=list()

## process result.txt files for counts
for(i in 1:ncol(pvalue_results_human_voom)){
  rsea_results_human_voom_go[[i]]<-read.table(paste0("GSE/",colnames(pvalue_results_human_voom)[i],"/rSEA/GO/result.txt"),header=TRUE)
}

knockout_target_gse<-read.table("knockout_target_gse",header=FALSE)

database_lists2<-load("/Users/mingyoungshin/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/PFOCRInPathwayAnalyses/pfocr_miny_April2021.RData")
database_lists1<-load("/Users/mingyoungshin/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/PFOCRInPathwayAnalyses/databases.RData")#has wp, pfocr, go
database_lists3<-load("/Users/mingyoungshin/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/PFOCRInPathwayAnalyses/databases_pfocr_3sets.RData")
min_set_size <- 3
max_set_size <- 500 

database_lists2<- unname(unlist(sapply(database_lists2, grep, pattern="_list$", value = T, perl = T)))
for (db in database_lists2) {
  eval(call("<-", as.name(db),  Filter(Negate(is.null), lapply(get(db), function(x){
    if(length(x) < min_set_size | length(x) > max_set_size)
      NULL
    else
      x
  }))
  ))
}

database_lists1 <- unname(unlist(sapply(database_lists1, grep, pattern="_list$", value = T, perl = T)))
for (db in database_lists1) {
  eval(call("<-", as.name(db),  Filter(Negate(is.null), lapply(get(db), function(x){
    if(length(x) < min_set_size | length(x) > max_set_size)
      NULL
    else
      x
  }))
  ))
}

database_lists3<- unname(unlist(sapply(database_lists3, grep, pattern="_list$", value = T, perl = T)))
for (db in database_lists3) {
  eval(call("<-", as.name(db),  Filter(Negate(is.null), lapply(get(db), function(x){
    if(length(x) < min_set_size | length(x) > max_set_size)
      NULL
    else
      x
  }))
  ))
}

########GO
merged_go<-read.table("/Users/mingyoungshin/git/PFOCRInPathwayAnalyses/human_knockout_go.txt",header=TRUE)
merged_go<-read.table("/Users/mingyoungshin/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/PFOCRInPathwayAnalyses/human_knockout_go_April2021version.txt",header=TRUE)
targets=which(is.na(merged_go$GO)==FALSE)#get GSEs matched with knock genes
target_gses=merged_go[targets,2]
length(target_gses)
target_gses=target_gses[-which(target_gses%in%c("GSE75291","GSE78506","GSE141126","GSE76421"))]
length(target_gses)
target_index=match(target_gses,colnames(pvalue_results_human_voom))#retrieve their index
target_index_go=target_index

rsea_go_sig_gene_count=NULL
rsea_go_aucs=NULL
rsea_go_aucs_ci=NULL
rsea_go_50lower=NULL
rsea_go_TRUEs=NULL
for(i in unique(target_index)){
  print(i)
  rsea_go_sig_gene_count=c(rsea_go_sig_gene_count,sum(rsea_results_human_voom_go[[i]]$Com.adjP <0.05,na.rm=TRUE))
  
  if( sum(is.na(rsea_results_human_voom_go[[i]])) ==sum(lengths(rsea_results_human_voom_go[[i]])) ){next}#pass null cases
  gse_match=which(as.character(merged_go$GSE)==colnames(pvalue_results_human_voom)[i])
  for(j in gse_match){
    print(merged_go$GO[j])
    if(is.na(merged_go$GO[j])){next}
    
    true_sig=merged_go$GO[j]
    true_sig=as.character(strsplit(as.character(true_sig),";")[[1]])
    true_sig=true_sig[true_sig%in%(names(go_list))]
    print(length(true_sig))
    response=apply(as.matrix(1:length(rsea_results_human_voom_go[[i]]$set_id)),1,function(x) rsea_results_human_voom_go[[i]]$set_id[x]%in%true_sig)
    print(sum(response))
    print(table(response))
    rsea_go_TRUEs=rbind(rsea_go_TRUEs,c(i,merged_go$Symbol[j],as.numeric(table(response)['TRUE'])))
    
    if(length(table(response))==1) next
    
   
    
    roc_r=roc(response,as.numeric(rsea_results_human_voom_go[[i]]$Comp.adjP))
    auc_r=auc(roc_r)
    
    if(auc_r<=0.5){rsea_go_50lower=rbind(rsea_go_50lower,c(i,j))}
    
    auc_r_ci=ci(auc_r)
    if(!is.na(auc_r_ci[1])){
      if(auc_r_ci[1] < roc_r$auc[1] & auc_r_ci[2] > roc_r$auc[1]){
        rsea_go_aucs_ci=c(rsea_go_aucs_ci,auc_r[1])
      }
    }

    
    rsea_go_aucs=c(rsea_go_aucs,auc_r[1])
  }
  
  
}


boxplot(rsea_go_aucs)




########wp

rsea_results_human_voom_wp=list()

## process result.txt files for counts
for(i in 1:ncol(pvalue_results_human_voom)){
  rsea_results_human_voom_wp[[i]]<-read.table(paste0("GSE/",colnames(pvalue_results_human_voom)[i],"/rSEA/WP/result.txt"),header=TRUE)
}


merged_wp<-read.table("/Users/mingyoungshin/git/PFOCRInPathwayAnalyses/human_knockout_wp.txt",header=TRUE)
merged_wp<-read.table("/Users/mingyoungshin/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/PFOCRInPathwayAnalyses/human_knockout_wp_April2021version.txt",header=TRUE)

targets=which(is.na(merged_wp$WP)==FALSE)#get GSEs matched with knock genes
target_gses=merged_wp[targets,2]
length(target_gses)
target_gses=target_gses[-which(target_gses%in%c("GSE75291","GSE78506","GSE141126","GSE76421"))]
length(target_gses)

target_index=match(target_gses,colnames(pvalue_results_human_voom))#retrieve their index
target_index_wp=target_index

rsea_wp_sig_gene_count=NULL
rsea_wp_aucs=NULL
rsea_wp_aucs_ci=NULL
rsea_wp_50lower=NULL
rsea_wp_TRUEs=NULL
for(i in unique(target_index)){
  print(i)
  rsea_wp_sig_gene_count=c(rsea_wp_sig_gene_count,sum(rsea_results_human_voom_wp[[i]]$Com.adjP <0.05,na.rm=TRUE))
  
  if( sum(is.na(rsea_results_human_voom_wp[[i]])) ==sum(lengths(rsea_results_human_voom_wp[[i]])) ){next}#pass null cases
  gse_match=which(as.character(merged_wp$GSE)==colnames(pvalue_results_human_voom)[i])
  for(j in gse_match){
    print(merged_wp$WP[j])
    if(is.na(merged_wp$WP[j])){next}
    
    true_sig=merged_wp$WP[j]
    true_sig=as.character(strsplit(as.character(true_sig),";")[[1]])
    true_sig=true_sig[true_sig%in%(names(wp_list))]
    print(length(true_sig))
    response=apply(as.matrix(1:length(rsea_results_human_voom_wp[[i]]$set_id)),1,function(x) rsea_results_human_voom_wp[[i]]$set_id[x]%in%true_sig)
    print(sum(response))
    print(table(response))
    
    rsea_wp_TRUEs=rbind(rsea_wp_TRUEs,c(i,merged_wp$Symbol[j],as.numeric(table(response)['TRUE'])))
    
    if(length(table(response))==1) next
    
    roc_r=roc(response,as.numeric(rsea_results_human_voom_wp[[i]]$Comp.adjP))
    auc_r=auc(roc_r)
    if(auc_r<=0.5){rsea_wp_50lower=rbind(rsea_wp_50lower,c(i,j))}
    
    auc_r_ci=ci(auc_r)
    if(!is.na(auc_r_ci[1])){    if(auc_r_ci[1] < roc_r$auc[1] & auc_r_ci[2] > roc_r$auc[1]){
      rsea_wp_aucs_ci=c(rsea_wp_aucs_ci,auc_r[1])
    }}

    
    rsea_wp_aucs=c(rsea_wp_aucs,auc_r[1])
  }
  
  
}


boxplot(rsea_wp_aucs)



########pfocr

rsea_results_human_voom_pfocr=list()

## process result.txt files for counts
for(i in 1:ncol(pvalue_results_human_voom)){
  rsea_results_human_voom_pfocr[[i]]<-read.table(paste0("GSE/",colnames(pvalue_results_human_voom)[i],"/rSEA/PFOCR/result.txt"),header=TRUE)
}


merged_pfocr<-read.table("/Users/mingyoungshin/git/PFOCRInPathwayAnalyses/human_knockout_pfocr.txt",header=TRUE)
merged_pfocr<-read.table("/Users/mingyoungshin/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/PFOCRInPathwayAnalyses/human_knockout_pfocr_April2021version.txt",header=TRUE)

targets=which(is.na(merged_pfocr$PFOCR)==FALSE)#get GSEs matched with knock genes
target_gses=merged_pfocr[targets,2]
length(target_gses)
target_gses=target_gses[-which(target_gses%in%c("GSE75291","GSE78506","GSE141126","GSE76421"))]
length(target_gses)

target_index=match(target_gses,colnames(pvalue_results_human_voom))#retrieve their index
target_index_pfocr=target_index

rsea_pfocr_sig_gene_count=NULL
rsea_pfocr_aucs=NULL
rsea_pfocr_aucs_ci=NULL
rsea_pfocr_50lower=NULL
rsea_pfocr_TRUEs=NULL
for(i in unique(target_index)){
  print(i)
  rsea_pfocr_sig_gene_count=c(rsea_pfocr_sig_gene_count,sum(rsea_results_human_voom_pfocr[[i]]$Com.adjP <0.05,na.rm=TRUE))
  
  if( sum(is.na(rsea_results_human_voom_pfocr[[i]])) ==sum(lengths(rsea_results_human_voom_pfocr[[i]])) ){next}#pass null cases
  gse_match=which(as.character(merged_pfocr$GSE)==colnames(pvalue_results_human_voom)[i])
  for(j in gse_match){
    print(merged_pfocr$PFOCR[j])
    if(is.na(merged_pfocr$PFOCR[j])){next}
    
    true_sig=merged_pfocr$PFOCR[j]
    true_sig=as.character(strsplit(as.character(true_sig),";")[[1]])
    true_sig=true_sig[true_sig%in%(names(pfocr_list))]
    print(length(true_sig))
    response=apply(as.matrix(1:length(rsea_results_human_voom_pfocr[[i]]$set_id)),1,function(x) rsea_results_human_voom_pfocr[[i]]$set_id[x]%in%true_sig)
    print(sum(response))
    print(table(response))
    rsea_pfocr_TRUEs=rbind(rsea_pfocr_TRUEs,c(i,merged_pfocr$Symbol[j],as.numeric(table(response)['TRUE'])))
    
    if(length(table(response))==1) next
    
    roc_r=roc(response,as.numeric(rsea_results_human_voom_pfocr[[i]]$Comp.adjP))
    auc_r=auc(roc_r)
    if(auc_r<=0.5){rsea_pfocr_50lower=rbind(rsea_pfocr_50lower,c(i,j))}
    
    auc_r_ci=ci(auc_r)
    if(!is.na(auc_r_ci[1])){    if(auc_r_ci[1] < roc_r$auc[1] & auc_r_ci[2] > roc_r$auc[1]){
      rsea_pfocr_aucs_ci=c(rsea_pfocr_aucs_ci,auc_r[1])
    }}

    
    rsea_pfocr_aucs=c(rsea_pfocr_aucs,auc_r[1])
  }
  
  
}


boxplot(rsea_pfocr_aucs)



library(readxl)
library(org.Hs.eg.db)

data<-read_excel("/Users/mingyoungshin/git/PFOCRInPathwayAnalyses/human_knockout_genelist_alias_updated.xlsx", col_names = F)
data_gene=as.character(unlist(data[,2]))

gene_entrez <- bitr(data_gene, fromType = "SYMBOL",
                    toType = c("ENTREZID"),
                    OrgDb = org.Hs.eg.db, drop=FALSE)


colnames(gene_entrez)[1]="Gene"
gene_entrez=unique(gene_entrez)


tb=table(gene_entrez$Gene)
tb[which(tb>1)]


result_pfocr=matrix(NA,nrow=nrow(gene_entrez),ncol=3)
for(i in 1:nrow(result_pfocr)){
  result_pfocr[i,1]=gene_entrez$Gene[i]
  result_pfocr[i,2]=gene_entrez$ENTREZID[i]
  
  for(j in 1:length(pfocr_3sets_list)){
    temp=NULL
    if( result_pfocr[i,2]%in%pfocr_3sets_list[[j]]){
      temp=names(pfocr_3sets_list)[j]
    }
    if(length(temp)>0){result_pfocr[i,3]=paste(result_pfocr[i,3],temp,sep=";")}
    
  }
}

result_pfocr=gsub("NA;","",result_pfocr)

nrow(result_pfocr)-sum(is.na(result_pfocr[,3]))
64
colnames(result_pfocr)=c("Symbol","Entrez","PFOCR")
colnames(data)=c("GSE","Symbol")
merged_pfocr=merge(data,result_pfocr,by="Symbol")

nrow(merged_pfocr)-sum(is.na(merged_pfocr[,3]))
85
write.table(merged_pfocr,"human_knockout_pfocr_3sets.txt",col.names = colnames(merged_pfocr),quote=FALSE)




library(readxl)
library(org.Hs.eg.db)

data<-read_excel("/Users/mingyoungshin/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/PFOCRInPathwayAnalyses/human_knockout_genelist_alias_updated_April2021.xlsx", col_names = F)
data_gene=as.character(unlist(data[,2]))

gene_entrez <- bitr(data_gene, fromType = "SYMBOL",
                    toType = c("ENTREZID"),
                    OrgDb = org.Hs.eg.db, drop=FALSE)


colnames(gene_entrez)[1]="Gene"
gene_entrez=unique(gene_entrez)


tb=table(gene_entrez$Gene)
tb[which(tb>1)]


result_pfocr=matrix(NA,nrow=nrow(gene_entrez),ncol=3)
for(i in 1:nrow(result_pfocr)){
  result_pfocr[i,1]=gene_entrez$Gene[i]
  result_pfocr[i,2]=gene_entrez$ENTREZID[i]
  
  for(j in 1:length(pfocr_3sets_list)){
    temp=NULL
    if( result_pfocr[i,2]%in%pfocr_3sets_list[[j]]){
      temp=names(pfocr_3sets_list)[j]
    }
    if(length(temp)>0){result_pfocr[i,3]=paste(result_pfocr[i,3],temp,sep=";")}
    
  }
}

result_pfocr=gsub("NA;","",result_pfocr)

nrow(result_pfocr)-sum(is.na(result_pfocr[,3]))
40
colnames(result_pfocr)=c("Symbol","Entrez","PFOCR")
colnames(data)=c("GSE","Symbol")
merged_pfocr=merge(data,result_pfocr,by="Symbol")

nrow(merged_pfocr)-sum(is.na(merged_pfocr[,3]))
53
write.table(merged_pfocr,"human_knockout_pfocr_3sets_April2021version.txt",col.names = colnames(merged_pfocr),quote=FALSE)



rsea_results_human_voom_pfocr_3sets=list()

## process result.txt files for counts
for(i in 1:ncol(pvalue_results_human_voom)){
  rsea_results_human_voom_pfocr_3sets[[i]]<-read.table(paste0("GSE/",colnames(pvalue_results_human_voom)[i],"/rSEA/PFOCR_3sets/result.txt"),header=TRUE)
}



merged_pfocr_3sets<-read.table("/Users/mingyoungshin/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/PFOCRInPathwayAnalyses/human_knockout_pfocr_3sets.txt",header=TRUE)
merged_pfocr_3sets<-read.table("/Users/mingyoungshin/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/PFOCRInPathwayAnalyses/human_knockout_pfocr_3sets_April2021version.txt",header=TRUE)

targets=which(is.na(merged_pfocr_3sets$PFOCR)==FALSE)#get GSEs matched with knock genes
target_gses=merged_pfocr_3sets[targets,2]
length(target_gses)
target_gses=target_gses[-which(target_gses%in%c("GSE75291","GSE78506","GSE141126","GSE76421"))]
length(target_gses)

target_index=match(target_gses,colnames(pvalue_results_human_voom))#retrieve their index
target_index_pfocr_3sets=target_index

rsea_pfocr_3sets_sig_gene_count=NULL
rsea_pfocr_3sets_aucs=NULL
rsea_pfocr_3sets_aucs_ci=NULL
rsea_pfocr_3set_50lower=NULL
rsea_pfocr_3rsets_TRUEs=NULL
for(i in unique(target_index)){
  print(i)
  rsea_pfocr_3sets_sig_gene_count=c(rsea_pfocr_3sets_sig_gene_count,sum(rsea_results_human_voom_pfocr_3sets[[i]]$Com.adjP <0.05,na.rm=TRUE))
  
  if( sum(is.na(rsea_results_human_voom_pfocr_3sets[[i]])) ==sum(lengths(rsea_results_human_voom_pfocr_3sets[[i]])) ){next}#pass null cases
  gse_match=which(as.character(merged_pfocr_3sets$GSE)==colnames(pvalue_results_human_voom)[i])
  for(j in gse_match){
    print(merged_pfocr_3sets$PFOCR[j])
    if(is.na(merged_pfocr_3sets$PFOCR[j])){next}
    
    true_sig=merged_pfocr_3sets$PFOCR[j]
    true_sig=as.character(strsplit(as.character(true_sig),";")[[1]])
    true_sig=true_sig[true_sig%in%(names(pfocr_3sets_list))]
    print(length(true_sig))
    response=apply(as.matrix(1:length(rsea_results_human_voom_pfocr_3sets[[i]]$set_id)),1,function(x) rsea_results_human_voom_pfocr_3sets[[i]]$set_id[x]%in%true_sig)
    print(sum(response))
    print(table(response))
    rsea_pfocr_3rsets_TRUEs=rbind(rsea_pfocr_3rsets_TRUEs,c(i,merged_pfocr_3sets$Symbol[j],as.numeric(table(response)['TRUE'])))
    
    if(length(table(response))==1) next
    
    roc_r=roc(response,as.numeric(rsea_results_human_voom_pfocr_3sets[[i]]$Comp.adjP))
    auc_r=auc(roc_r)
    if(auc_r<=0.5){rsea_pfocr_3set_50lower=rbind(rsea_pfocr_3set_50lower,c(i,j))}
    
    auc_r_ci=ci(auc_r)
    if(!is.na(auc_r_ci[1])){
      if(auc_r_ci[1] < roc_r$auc[1] & auc_r_ci[2] > roc_r$auc[1]){
        rsea_pfocr_3sets_aucs_ci=c(rsea_pfocr_3sets_aucs_ci,auc_r[1])
      }
    }

    
    rsea_pfocr_3sets_aucs=c(rsea_pfocr_3sets_aucs,auc_r[1])
  }
  
  
}




library(ggplot2)

data=rbind( cbind(rsea_pfocr_3sets_aucs,rep("PFOCR_3sets",length(rsea_pfocr_3sets_aucs))),
  cbind(rsea_pfocr_aucs,rep("PFOCR",length(rsea_pfocr_aucs))),
            cbind(rsea_wp_aucs,rep("WP",length(rsea_wp_aucs))),
            cbind(rsea_go_aucs,rep("GO",length(rsea_go_aucs))) )
colnames(data)=c("AUC","Data")

data=data.frame(data)
data$AUC=as.numeric(as.character(data$AUC))

# Change violin plot colors by groups
p<-ggplot(data, aes(x=Data, y=AUC, fill=Data)) +
  geom_boxplot() + scale_fill_manual(values=c("antiquewhite3", "darkseagreen", "darksalmon","lightblue")) + theme_classic()
p

data_ci=rbind( cbind(rsea_pfocr_3sets_aucs_ci,rep("PFOCR_3sets",length(rsea_pfocr_3sets_aucs_ci))),
            cbind(rsea_pfocr_aucs_ci,rep("PFOCR",length(rsea_pfocr_aucs_ci))),
            cbind(rsea_wp_aucs_ci,rep("WP",length(rsea_wp_aucs_ci))),
            cbind(rsea_go_aucs_ci,rep("GO",length(rsea_go_aucs_ci))) )
colnames(data_ci)=c("AUC","Data")
data_ci=data.frame(data_ci)
data_ci$AUC=as.numeric(as.character(data_ci$AUC))

p<-ggplot(data_ci, aes(x=Data, y=AUC, fill=Data)) +
  geom_boxplot() + scale_fill_manual(values=c("antiquewhite3", "darkseagreen", "darksalmon","lightblue")) + theme_classic()
p

saveRDS(data_ci,"data_ci.rds")

#################assess sets with auc <0.5
Reduce(intersect, list(rsea_go_50lower[,1],rsea_wp_50lower[,1],rsea_pfocr_3set_50lower[,1],rsea_pfocr_50lower[,1]))#0
library(ggvenn)
ggvenn(
  list(rsea_go_50lower=rsea_go_50lower[,1],rsea_wp_50lower=rsea_wp_50lower[,1],rsea_pfocr_3set_50lower=rsea_pfocr_3set_50lower[,1],rsea_pfocr_50lower=rsea_pfocr_50lower[,1]), 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

table(c(rsea_go_50lower[,1],rsea_wp_50lower[,1],rsea_pfocr_3set_50lower[,1],rsea_pfocr_50lower[,1]))

#11  30  61  77  81  93 110 120 121 131 156 195 211 234 244 317 332 
#1   2   3   2   1   1   2   1   3   3   2   1   1   1   2   2   1 





#####################################################################################PCA plot
#################data 1
i=61
merged_filtered_gses[[i]]$GSE#GSE75215

geo_data<-GEOquery::getGEO(merged_filtered_gses[[i]]$GSE,getGPL = FALSE)
show(geo_data)

pheno=pData(geo_data[[1]])
dim(pheno)
colnames(pheno)
gsm_title=cbind(rownames(pheno),as.character(pheno$title))
head(gsm_title)

# Selected samples to be extracted
samp = gsm_title[,1]

# Identify columns to be extracted
sample_locations = which(samples %in% samp)
sample_locations

# extract gene expression from compressed data
expression = rhdf5::h5read(destination_file, "data/expression", index=list(1:length(genes), sample_locations))

rownames(expression) = genes
colnames(expression) = samples[sample_locations]

###DEG

filter <- apply(expression, 1, function(x) length(x[x>5])>=2) #minimum count (4) is satified in at least two conditions
FilterGeneCountData <- expression[filter,]

###group GSMs
Group=gsm_title[,2]
group_index=grep("group",names(merged_filtered_gses[[i]]))
selected_all=NULL
for(gi in group_index){
  selected=which( (Group%in%unlist(merged_filtered_gses[[i]][gi])) ==TRUE)
  Group[selected]=paste0("group",gi)
  selected_all=c(selected_all,selected)
}
gsm_title1=cbind(gsm_title,Group)
gsm_title1=gsm_title1[selected_all,]
colnames(gsm_title1)=c("GSM","Pheno","Group")

gsm_title1=as.data.frame(gsm_title1)
gsm_title1$Group <- as.factor(gsm_title1$Group)

TempIndices <- match(gsm_title1[,1], colnames(FilterGeneCountData))

colN=colnames(FilterGeneCountData)

FilterGeneCountData=FilterGeneCountData[,TempIndices]

y <- DGEList(counts=FilterGeneCountData, group=NULL)
y <- calcNormFactors(y)


cpm <- cpm(y, log = TRUE, prior.count = 0.01)
rv <- apply(cpm,1,var) 

#Select genes with highest variance.
keep <- order(rv, decreasing = TRUE)[1:500]
selected <- cpm[keep, ] %>% t()

#Transpose is needed to ensure that each row is a vector.
pca <- prcomp(selected, scale=T, center = T)

stddev <- pca$sdev
pc1_var <- round(100*stddev[1]^2/sum(stddev^2))
pc2_var <- round(100*stddev[2]^2/sum(stddev^2))
pc3_var <- round(100*stddev[3]^2/sum(stddev^2))
PlotData <- data.frame(cbind(PC1 = pca$x[,1], PC2 = pca$x[,2],PC3 = pca$x[,3]))
PlotData <- gsm_title1[, c("Group")] %>%
  cbind(PlotData, .)
colnames(PlotData)[4]="Group"
gp1=ggplot(PlotData, aes(x=PC1, y=PC2, color=Group, shape=Group)) + 
  geom_point(size=4.5) +  
  xlab(paste("PC1:", pc1_var, "% variance")) + 
  ylab(paste("PC2:", pc2_var, "% variance"))

gp1
jpeg(paste0("PCA_",merged_filtered_gses[[i]]$GSE,"_index",i,".jpeg"))
print(gp1)
dev.off()

####################### 
GSE68925 ROCK1
GSE75215 PRMT5
GSE59598 ZNF423
GSE99596 SIRT7
GSE96637 ROR1
GSE76421 JAG1
GSE108500 FOXP1
GSE77140 TRIM24
GSE118532 RNF168
#######################check true gene set count

tb=table(c(rsea_go_50lower[,1],rsea_wp_50lower[,1],rsea_pfocr_3set_50lower[,1],rsea_pfocr_50lower[,1]))
tb=cbind(names(tb),tb)
colnames(tb)=c("GSE_index","Count_database_low_performance")
colnames(rsea_pfocr_3rsets_TRUEs)=c("GSE_index","knockgene","Count_True_geneset")
colnames(rsea_pfocr_TRUEs)=c("GSE_index","knockgene","Count_True_geneset")
colnames(rsea_wp_TRUEs)=c("GSE_index","knockgene","Count_True_geneset")
colnames(rsea_go_TRUEs)=c("GSE_index","knockgene","Count_True_geneset")

rsea_pfocr_3rsets_TRUEs_merged=merge(rsea_pfocr_3rsets_TRUEs,tb,by="GSE_index")
rsea_pfocr_3rsets_TRUEs_merged=rsea_pfocr_3rsets_TRUEs_merged[order(rsea_pfocr_3rsets_TRUEs_merged[,4]),]

rsea_pfocr_TRUEs_merged=merge(rsea_pfocr_TRUEs,tb,by="GSE_index")
rsea_pfocr_TRUEs_merged=rsea_pfocr_TRUEs_merged[order(rsea_pfocr_TRUEs_merged[,4]),]

rsea_wp_TRUEs_merged=merge(rsea_wp_TRUEs,tb,by="GSE_index")
rsea_wp_TRUEs_merged=rsea_wp_TRUEs_merged[order(rsea_wp_TRUEs_merged[,4]),]

rsea_go_TRUEs_merged=merge(rsea_go_TRUEs,tb,by="GSE_index")
rsea_go_TRUEs_merged=rsea_go_TRUEs_merged[order(rsea_go_TRUEs_merged[,4]),]

save(rsea_pfocr_3rsets_TRUEs_merged,rsea_pfocr_TRUEs_merged,rsea_wp_TRUEs_merged,rsea_go_TRUEs_merged,file="rsea_GSE_knockoutGene_TRUEgeneset_lowPerformanceDBcount.RData")


permute_databases.sh
ob - mean(null)/sd(null)

