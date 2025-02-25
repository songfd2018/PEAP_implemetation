# Loading libraries and data
rm(list=ls())
library(ggplot2)
library(survival)
library(xtable)
library(fgsea)
library(ggpubr)
library(gbm)
library(caret)
library(limma)
library(edgeR)
library(glmnet)
library(modelr)
library(purrr)
library(caret)
library(dplyr)
library(corrplot)
library(nnet)
library(DESeq2)
library(readr)
library(readxl)
library(genefu)

## Input parameters##
subtype_limit = 1000
lambda_fix = 0.02
set.seed(1234)
setwd("D:/PEAP/PEAP_02/New_data/")
# TCGA 
load("TCGA_PEAP_preprocess.Rdata")
pheno_sample_tcga = pheno_sample

# GSE 
load("GSE5327_PEAP_preprocess.Rdata")
pheno_sample_GSE5327 = pheno_sample
load("GSE2034_PEAP_preprocess.Rdata")
pheno_sample_GSE2034 = pheno_sample
load("SCANB_PEAP_preprocess.Rdata")
GeneExp_BRCA_SCANB=log2(GeneExp_BRCA_SCANB+1)


# fold change
load("PEAP_Top1500_genes_eachpheno.Rdata")

rm(pheno_sample)
# phenotype
name_pheno <- c("Basal", "Her2","LumA", "LumB", "Normal")
clinical_brca_raw=read_tsv("TCGA_survival.tsv")

cindex_list2034 = list()
cindex_list5327 = list()
cindex_listSCANB = list()
cox_list = list()

# Extract Survival data

# TCGA
surv_dat = clinical_brca_raw[,c("case_submitter_id",
                                "vital_status", 
                                "days_to_last_follow_up",
                                "days_to_death")]
surv_dat$surv_time = ifelse(surv_dat$vital_status=="Alive",
                            surv_dat$days_to_last_follow_up,
                            surv_dat$days_to_death)
surv_dat$vital_num = as.numeric(surv_dat$vital_status=="Dead")
surv_dat=surv_dat[,c("case_submitter_id","vital_num","surv_time")]
surv_dat =surv_dat[surv_dat$surv_time!= "'--",]
surv_dat$surv_time= as.numeric(surv_dat$surv_time)/30
surv_dat=unique(surv_dat[surv_dat$surv_time>0,])
surv_dat$vital_num= ifelse(surv_dat$surv_time > 60,
                           0,
                           surv_dat$vital_num)
surv_dat$surv_time= ifelse(surv_dat$surv_time > 60,
                           60,
                           surv_dat$surv_time)



print("limma start")
# DEG limma
## Feature Selection for DEG limma

# Normalization
DGE_obj <- DGEList(counts = 2^GeneExp_BRCA - 1, group = pheno_sample_tcga)
keep.exprs <- filterByExpr(DGE_obj)
DGE_obj <- DGE_obj[keep.exprs,, keep.lib.sizes=FALSE]
DGE_obj <- calcNormFactors(DGE_obj)
logCPM <- cpm(DGE_obj, log=TRUE, prior.count=3)

DEG_df_list = list()
limmaG_list = list()
selecG_limma_list = list()
top_gene_table_list = list()

# Main: one vs rest
for(i in 1:5){ 
  pheno = name_pheno[i]
  # limma
  pheno_sample_other <- pheno_sample_tcga
  for(j in 1:length(pheno_sample_tcga)){
    if(pheno!=pheno_sample_other[j]){
      pheno_sample_other[j] <- "other"
    }
  }
  design <- model.matrix(~pheno_sample_other) 
  colnames(design) <- gsub("pheno_sample_other", "", colnames(design))
  fit <- lmFit(logCPM, design) 
  fit <- eBayes(fit, trend=TRUE) 
  top_gene_table <- topTable(fit, coef=2,
                             number = subtype_limit,
                             adjust.method = "BH",p.value = 0.005)
  top_gene_table_list = append(top_gene_table_list, list(top_gene_table))
  top_genes <- rownames(top_gene_table)
  N <- ncol(logCPM)
  P <- length(top_genes)
  DEG_df_eachPheno <- matrix(NA, nrow = N, ncol = P)
  gene_list <-rownames(logCPM)
  
  for(p in 1:P){
    g <- which(gene_list == top_genes[p])
    DEG_df_eachPheno[,p] <- logCPM[g,]
  }
  DEG_df_eachPheno <- as.data.frame(DEG_df_eachPheno)
  rownames(DEG_df_eachPheno) <- colnames(logCPM)
  colnames(DEG_df_eachPheno) <- top_genes
  DEG_df_eachPheno$phenotype = pheno_sample_tcga
  DEG_df_list = append(DEG_df_list, list(DEG_df_eachPheno))
  limmaGene = colnames(DEG_df_eachPheno)[-length(colnames(DEG_df_eachPheno))]
  limmaG_list = append(limmaG_list, list(limmaGene))
}
limmaG = unique(unlist(limmaG_list))
DEG_limma = as.data.frame(t(logCPM)[,limmaG])
DEG_limma$case_submitter_id = gsub("\\.", "-",substr(row.names(DEG_limma),1,12))
DEG_limma$subtype = pheno_sample_tcga


## Cox
surv_limma = merge(surv_dat,DEG_limma, by="case_submitter_id")
x=as.matrix(surv_limma[,!colnames(surv_limma)%in%c("case_submitter_id","vital_num","surv_time","subtype")])
y=Surv(surv_limma$surv_time,surv_limma$vital_num)
limma_fit <- glmnet(x, y, family = "cox", alpha = 1, lambda = lambda_fix)
print("limma end")


# DEG edgeR
print("edgeR start")
## Feature Selection for DEG edgeR

# Normalization
DGE_obj <- DGEList(counts = 2^GeneExp_BRCA - 1, group = pheno_sample_tcga)
keep.exprs <- filterByExpr(DGE_obj)
DGE_obj <- DGE_obj[keep.exprs,, keep.lib.sizes=FALSE]
DGE_obj <- calcNormFactors(DGE_obj)
DGE_obj$group <- pheno_sample_tcga

# calculate preparations
design <- model.matrix(~pheno_sample_tcga)
colnames(design)=levels(as.factor(pheno_sample_tcga))
DGE_obj <- estimateDisp(DGE_obj, design)
logCPM <- cpm(DGE_obj, log=TRUE, prior.count=3)
N <- ncol(logCPM)
fit <- glmQLFit(DGE_obj, design)

DEG_df_list = list()
edgeRG_list = list()
selecG_edgeR_list = list()
top_gene_table_list = list()
# extract top genes and do lasso
for (i in 1:5) {
  contrast_vec = c(0,0,0,0,0)
  contrast_vec[i]=1
  result <- glmQLFTest(fit, contrast = contrast_vec)
  # output target gene
  top_gene_table <- topTags(result, n = subtype_limit, p.value = 0.005,adjust.method = "BH")
  top_gene_table_list = append(top_gene_table_list, list(top_gene_table))
  top_genes=rownames(top_gene_table)
  
  P <- length(top_genes)
  
  DEG_df_eachPheno <- matrix(NA, nrow = N, ncol = P)
  gene_list <-rownames(logCPM)
  
  for(p in 1:P){
    g <- which(gene_list == top_genes[p])
    DEG_df_eachPheno[,p] <- logCPM[g,]
  }
  DEG_df_eachPheno <- as.data.frame(DEG_df_eachPheno)
  rownames(DEG_df_eachPheno) <- colnames(logCPM)
  colnames(DEG_df_eachPheno) <- top_genes
  DEG_df_eachPheno$phenotype = pheno_sample_tcga
  DEG_df_list = append(DEG_df_list, list(DEG_df_eachPheno))
  edgeRGene = colnames(DEG_df_eachPheno)[-length(colnames(DEG_df_eachPheno))]
  edgeRG_list = append(edgeRG_list, list(edgeRGene))
}
edgeRG = unique(unlist(edgeRG_list))
DEG_edgeR = as.data.frame(t(logCPM)[,edgeRG])
DEG_edgeR$case_submitter_id = gsub("\\.", "-",substr(row.names(DEG_edgeR),1,12))
DEG_edgeR$subtype = pheno_sample_tcga


## Cox
surv_edgeR = merge(surv_dat,DEG_edgeR, by="case_submitter_id")
x=as.matrix(surv_edgeR[,!colnames(surv_edgeR)%in%c("case_submitter_id","vital_num","surv_time","subtype")])
y=Surv(surv_edgeR$surv_time,surv_edgeR$vital_num)
edgeR_fit <- glmnet(x, y, family = "cox", alpha = 1, lambda = lambda_fix)
print("edgeR end")




# Wilcoxon
print("Wilcoxon start")
## Feature Selection for DEG edgeR

# Normalization
DGE_obj <- DGEList(counts = 2^GeneExp_BRCA - 1, group = pheno_sample_tcga)
keep.exprs <- filterByExpr(DGE_obj)
DGE_obj <- DGE_obj[keep.exprs,, keep.lib.sizes=FALSE]
DGE_obj <- calcNormFactors(DGE_obj)
DGE_obj$group <- pheno_sample_tcga

# calculate preparations
design <- model.matrix(~pheno_sample_tcga)
colnames(design)=levels(as.factor(pheno_sample_tcga))
DGE_obj <- estimateDisp(DGE_obj, design)
logCPM <- cpm(DGE_obj, log=TRUE, prior.count=3)

pheno_name=unique(pheno_sample_tcga)
all_top_genes = data.frame()
wilcoxon_pvalue_table=data.frame()
for(i in 1:5){
  for(j in 1:nrow(logCPM)){
    gene_for_test = logCPM[j,]
    x=unlist(gene_for_test[pheno_sample_tcga==pheno_name[i]])
    y=unlist(gene_for_test[pheno_sample_tcga!=pheno_name[i]])
    wil=wilcox.test(x,y)
    wil_result = data.frame(genename=rownames(logCPM)[j], pval=wil$p.value)
    wilcoxon_pvalue_table<-rbind(wilcoxon_pvalue_table, wil_result)
  }
}

wilcoxon_pvalue_table.adj=wilcoxon_pvalue_table
for(i in 1:5){
  wilcoxon_pvalue_table.adj[(nrow(logCPM)*(i-1)+1):(nrow(logCPM)*i),2]=p.adjust(wilcoxon_pvalue_table[(nrow(logCPM)*(i-1)+1):(nrow(logCPM)*i),2], method = "BH")
}



DEG_df_list = list()
wilcoxonG_list = list()
selecG_edgeR_list = list()
top_gene_table_list = list()
# extract top genes and do lasso
for (i in 1:5) {
  contrast_vec = c(0,0,0,0,0)
  contrast_vec[i]=1
  # output target gene
  
  top_genes=wilcoxon_pvalue_table[(nrow(logCPM)*(i-1)+1):(nrow(logCPM)*i),][which(wilcoxon_pvalue_table[(nrow(logCPM)*(i-1)+1):(nrow(logCPM)*i),2]<=5e-10),1]
  
  P <- length(top_genes)
  
  DEG_df_eachPheno <- matrix(NA, nrow = N, ncol = P)
  gene_list <-rownames(logCPM)
  
  for(p in 1:P){
    g <- which(gene_list == top_genes[p])
    DEG_df_eachPheno[,p] <- logCPM[g,]
  }
  DEG_df_eachPheno <- as.data.frame(DEG_df_eachPheno)
  rownames(DEG_df_eachPheno) <- colnames(logCPM)
  colnames(DEG_df_eachPheno) <- top_genes
  DEG_df_eachPheno$phenotype = pheno_sample_tcga
  DEG_df_list = append(DEG_df_list, list(DEG_df_eachPheno))
  wilcoxonGene = colnames(DEG_df_eachPheno)[-length(colnames(DEG_df_eachPheno))]
  wilcoxonG_list = append(wilcoxonG_list, list(wilcoxonGene))
}
wilcoxonG = unique(unlist(wilcoxonG_list))
DEG_wilcoxon = as.data.frame(t(logCPM)[,wilcoxonG])
DEG_wilcoxon$case_submitter_id = gsub("\\.", "-",substr(row.names(DEG_wilcoxon),1,12))
DEG_wilcoxon$subtype = pheno_sample_tcga


## Cox
surv_wilcoxon = merge(surv_dat,DEG_wilcoxon, by="case_submitter_id")
x=as.matrix(surv_wilcoxon[,!colnames(surv_wilcoxon)%in%c("case_submitter_id","vital_num","surv_time","subtype")])
y=Surv(surv_wilcoxon$surv_time,surv_edgeR$vital_num)
wilcoxon_fit <- glmnet(x, y, family = "cox", alpha = 1, lambda = lambda_fix)
print("Wilcoxon end")






# DEG DESeq2
print("DESeq2 start")
## Feature Selection for DEG DESeq2
DEG_df_list = list()
DESeq2G_list = list()
selecG_DESeq2_list = list()
top_gene_table_list = list()

# extract top genes and do lasso
colData <- data.frame(condition = pheno_sample_tcga)
DGE_obj <- DGEList(counts = round(2^GeneExp_BRCA - 1), group = pheno_sample_tcga)
logCPM <- cpm(DGE_obj, log=TRUE, prior.count=3)
dim(logCPM)
for (i in 1:5) {
  pheno_sample_tcga_copy = pheno_sample_tcga
  pheno_sample_tcga_copy = ifelse(pheno_sample_tcga_copy==name_pheno[i],
                                  name_pheno[i],"other")
  colData <- data.frame(condition = pheno_sample_tcga)
  dds <- DESeqDataSetFromMatrix(countData = round(2^GeneExp_BRCA - 1,0),
                                colData = colData,
                                design = model.matrix(~pheno_sample_tcga_copy))
  dds <- DESeq(dds)
  res <- results(dds)
  N <- ncol(logCPM)
  res_adj <- res[which(res$padj < 0.005), ]
  # output target gene
  top_gene_table <- res_adj[order(res_adj$padj), ]
  top_gene_table_list = append(top_gene_table_list, list(top_gene_table))
  top_genes <- rownames(res_adj)[1:subtype_limit]
  top_genes = top_genes[!is.na(top_genes)]
  P <- length(top_genes)
  
  DEG_df_eachPheno <- matrix(NA, nrow = N, ncol = P)
  gene_list <-rownames(logCPM)
  for(p in 1:P){
    g <- which(gene_list == top_genes[p])
    DEG_df_eachPheno[,p] <- logCPM[g,]
  }
  DEG_df_eachPheno <- as.data.frame(DEG_df_eachPheno)
  rownames(DEG_df_eachPheno) <- colnames(logCPM)
  colnames(DEG_df_eachPheno) <- top_genes
  DEG_df_eachPheno$phenotype = pheno_sample_tcga
  DEG_df_list = append(DEG_df_list, list(DEG_df_eachPheno))
  DESeq2Gene = colnames(DEG_df_eachPheno)[-length(colnames(DEG_df_eachPheno))]
  DESeq2G_list = append(DESeq2G_list, list(DESeq2Gene))
}
DESeq2G = unique(unlist(DESeq2G_list))
DEG_DESeq2= as.data.frame(t(logCPM)[,DESeq2G])
DEG_DESeq2$case_submitter_id = gsub("\\.", "-",substr(row.names(DEG_DESeq2),1,12))
DEG_DESeq2$subtype = pheno_sample_tcga

## Cox

surv_DESeq2 = merge(surv_dat,DEG_DESeq2, by="case_submitter_id")
x=as.matrix(surv_DESeq2[,!colnames(surv_DESeq2)%in%c("case_submitter_id","vital_num","surv_time","subtype")])
y=Surv(surv_DESeq2$surv_time,surv_DESeq2$vital_num)
DESeq2_fit <- glmnet(x, y, family = "cox", alpha = 1, lambda = lambda_fix)
print("DESeq2 end")


# DEP
print("peap start")
## DEP fold change
Selected_genepairs = paste0(pmin(top1500_Genes[,1],top1500_Genes[,2]),
                            "_",
                            pmax(top1500_Genes[,1],top1500_Genes[,2]))
genePair = unique(names(table(Selected_genepairs)))
genePairMatrix =  matrix(unlist(strsplit(genePair,"_")), ncol = 2, byrow = TRUE)  
gene1 = GeneExp_BRCA[as.vector(genePairMatrix[,1]),]
row.names(gene1) = NULL
gene2 = GeneExp_BRCA[as.vector(genePairMatrix[,2]),]
row.names(gene2) = NULL
pheno_row_foldchange = gene1 - gene2
row.names(pheno_row_foldchange) = genePair
DEP_df = as.data.frame(t(pheno_row_foldchange))
peapP_list = list()
for (i in c(1:5)){
  start = 1+(1500*(i-1))
  end = subtype_limit+(1500*(i-1))
  peapP = paste(pmin(top1500_Genes[start:end, 1],
                     top1500_Genes[start:end, 2]),
                "_",
                pmax(top1500_Genes[start:end, 1],
                     top1500_Genes[start:end, 2]),
                sep="")
  peapP_list = append(peapP_list, list(peapP))
}
peapP = unique(unlist(peapP_list))
DEP_peap = DEP_df[,peapP]
DEP_peap$case_submitter_id = gsub("\\.", "-",substr(row.names(DEP_peap),1,12))
DEP_peap$subtype = pheno_sample_tcga

## Cox

surv_peap = merge(surv_dat,DEP_peap, by="case_submitter_id")
x=as.matrix(surv_peap[,!colnames(surv_peap)%in%c("case_submitter_id","vital_num","surv_time","subtype")])
y=Surv(surv_peap$surv_time,surv_peap$vital_num)
peap_fit <- glmnet(x, y, family = "cox", alpha = 1, lambda = lambda_fix)
print("peap end")

# PAM 50
print("PAM50 start")
data(pam50)
pam50_genes = row.names(pam50$centroids)

# extract pam50 genes and do lasso
colData <- data.frame(condition = pheno_sample_tcga)
DGE_obj <- DGEList(counts = round(2^GeneExp_BRCA - 1), group = pheno_sample_tcga)
logCPM <- cpm(DGE_obj, log=TRUE, prior.count=3)
row.names(logCPM)
pam50_g = pam50_genes[pam50_genes %in% row.names(logCPM)]
DEG_pam50= as.data.frame(t(logCPM)[,pam50_g])
DEG_pam50$case_submitter_id = gsub("\\.", "-",substr(row.names(DEG_pam50),1,12))
DEG_pam50$subtype = pheno_sample_tcga

## Cox

surv_pam50 = merge(surv_dat,DEG_pam50, by="case_submitter_id")
x=as.matrix(surv_pam50[,!colnames(surv_pam50)%in%c("case_submitter_id","vital_num","surv_time","subtype")])
y=Surv(surv_pam50$surv_time,surv_pam50$vital_num)
pam50_fit <- glmnet(x, y, family = "cox", alpha = 1, lambda = 0.02)
print("PAM50 end")

# 2034
print("2034")
for (year_limit in 1:5) {
  print(c("Year Limit: ", year_limit))
  ## GSE2034
  
  GSE2034_survival <- read_excel("GSE2034_survival.xlsx")
  surv_dat2034 = GSE2034_survival[,2:7]
  colnames(surv_dat2034) = c("GEO_id","node","surv_time",
                             "vital_num","ER_status","Brain_relapses")
  surv_dat2034$vital_num= ifelse(surv_dat2034$surv_time>year_limit*12,
                                 0,
                                 surv_dat2034$vital_num)
  surv_dat2034$surv_time= ifelse(surv_dat2034$surv_time>year_limit*12,
                                 year_limit*12,
                                 surv_dat2034$surv_time)
  ## limma
  # Normalization
  DGE_obj_2034 <- DGEList(counts = 2^GeneExp_BRCA_GSE2034 - 1,
                          group = pheno_sample_GSE2034)
  keep.exprs.2034 <- filterByExpr(DGE_obj_2034)
  DGE_obj_2034 <- calcNormFactors(DGE_obj_2034)
  logCPM_2034 <- as.data.frame(t(cpm(DGE_obj_2034, log=TRUE, prior.count=3)))
  logCPM_2034$subtype=pheno_sample_GSE2034
  surv_limma2034 = merge(surv_dat2034, 
                         logCPM_2034[,limmaG],
                         by.x="GEO_id",by.y="row.names")
  
  
  pred = predict(limma_fit,as.matrix(surv_limma2034[,limmaG]),type="response")
  tmp=Cindex(pred,Surv(surv_limma2034$surv_time,surv_limma2034$vital_num))
  cindex_list2034 = c(cindex_list2034, tmp)
  
  ## edgeR
  # Normalization
  DGE_obj_2034 <- DGEList(counts = 2^GeneExp_BRCA_GSE2034 - 1,
                          group = pheno_sample_GSE2034)
  keep.exprs.2034 <- filterByExpr(DGE_obj_2034)
  DGE_obj_2034 <- calcNormFactors(DGE_obj_2034)
  logCPM_2034 <- as.data.frame(t(cpm(DGE_obj_2034, log=TRUE, prior.count=3)))
  logCPM_2034$subtype=pheno_sample_GSE2034
  surv_edgeR2034 = merge(surv_dat2034, 
                         logCPM_2034[,edgeRG],
                         by.x="GEO_id",by.y="row.names")
  pred = predict(edgeR_fit,as.matrix(surv_edgeR2034[,edgeRG]),type="response")
  tmp=Cindex(pred,Surv(surv_edgeR2034$surv_time,surv_edgeR2034$vital_num))
  cindex_list2034 = c(cindex_list2034, tmp)
  
  
  ## wilcoxon
  # Normalization
  DGE_obj_2034 <- DGEList(counts = 2^GeneExp_BRCA_GSE2034 - 1,
                          group = pheno_sample_GSE2034)
  keep.exprs.2034 <- filterByExpr(DGE_obj_2034)
  DGE_obj_2034 <- calcNormFactors(DGE_obj_2034)
  logCPM_2034 <- as.data.frame(t(cpm(DGE_obj_2034, log=TRUE, prior.count=3)))
  logCPM_2034$subtype=pheno_sample_GSE2034
  surv_wilcoxon2034 = merge(surv_dat2034, 
                         logCPM_2034[,wilcoxonG],
                         by.x="GEO_id",by.y="row.names")
  pred = predict(wilcoxon_fit,as.matrix(surv_wilcoxon2034[,wilcoxonG]),type="response")
  tmp=Cindex(pred,Surv(surv_wilcoxon2034$surv_time,surv_wilcoxon2034$vital_num))
  cindex_list2034 = c(cindex_list2034, tmp)
  
  
  ## DESeq2
  # Normalization
  DGE_obj_2034 <- DGEList(counts = 2^GeneExp_BRCA_GSE2034 - 1,
                          group = pheno_sample_GSE2034)
  keep.exprs.2034 <- filterByExpr(DGE_obj_2034)
  DGE_obj_2034 <- calcNormFactors(DGE_obj_2034)
  logCPM_2034 <- as.data.frame(t(cpm(DGE_obj_2034, log=TRUE, prior.count=3)))
  logCPM_2034$subtype=pheno_sample_GSE2034
  surv_DESeq22034 = merge(surv_dat2034,
                          logCPM_2034[,DESeq2G],
                          by.x="GEO_id",by.y="row.names")
  pred = predict(DESeq2_fit,as.matrix(surv_DESeq22034[,DESeq2G]),type="response")
  tmp=Cindex(pred,Surv(surv_DESeq22034$surv_time,surv_DESeq22034$vital_num))
  cindex_list2034 = c(cindex_list2034, tmp)
  
  # peap
  genePairMatrix_2034 =  matrix(unlist(strsplit(peapP,"_")), ncol = 2, byrow = TRUE)
  gene1 = GeneExp_BRCA_GSE2034[as.vector(genePairMatrix_2034[,1]),]
  row.names(gene1) = NULL
  gene2 = GeneExp_BRCA_GSE2034[as.vector(genePairMatrix_2034[,2]),]
  row.names(gene2) = NULL
  pheno_row_foldchange_GSE2034 = gene1 - gene2
  row.names(pheno_row_foldchange_GSE2034) = peapP
  DEP_df = as.data.frame(t(pheno_row_foldchange_GSE2034))
  DEP_df$subtype = pheno_sample_GSE2034
  surv_peap2034 = merge(surv_dat2034, 
                        DEP_df[,c(peapP,"subtype")],
                        by.x="GEO_id",by.y="row.names")
  
  
  pred = predict(peap_fit,as.matrix(surv_peap2034[,peapP]),type="response")
  tmp=Cindex(pred,Surv(surv_peap2034$surv_time,surv_peap2034$vital_num))
  cindex_list2034 = c(cindex_list2034, tmp)
  
  # pam50
  DGE_obj_2034 <- DGEList(counts = 2^GeneExp_BRCA_GSE2034 - 1,
                          group = pheno_sample_GSE2034)
  keep.exprs.2034 <- filterByExpr(DGE_obj_2034)
  DGE_obj_2034 <- calcNormFactors(DGE_obj_2034)
  logCPM_2034 <- as.data.frame(t(cpm(DGE_obj_2034, log=TRUE, prior.count=3)))
  logCPM_2034$subtype=pheno_sample_GSE2034
  surv_pam502034 = merge(surv_dat2034,
                         logCPM_2034[,pam50_g],
                         by.x="GEO_id",by.y="row.names")
  pred = predict(pam50_fit,as.matrix(surv_pam502034[,pam50_g]),type="response")
  tmp=Cindex(pred,Surv(surv_pam502034$surv_time,surv_pam502034$vital_num))
  cindex_list2034 = c(cindex_list2034, tmp)
}

# 5327

print("5327")
for (year_limit in 1:5) {
  print(c("Year Limit: ", year_limit))
  ## GSE2034
  surv_dat5327 <- read_excel("GSE5327_survival.xlsx")
  surv_dat5327$vital_num= ifelse(surv_dat5327$surv_time>year_limit*12,
                                 0,
                                 surv_dat5327$vital_num)
  surv_dat5327$surv_time= ifelse(surv_dat5327$surv_time>year_limit*12,
                                 year_limit*12,
                                 surv_dat5327$surv_time)
  ## limma
  # Normalization
  DGE_obj_5327 <- DGEList(counts = 2^GeneExp_BRCA_GSE5327 - 1,
                          group = pheno_sample_GSE5327)
  keep.exprs.2034 <- filterByExpr(DGE_obj_5327)
  DGE_obj_5327 <- calcNormFactors(DGE_obj_5327)
  logCPM_5327 <- as.data.frame(t(cpm(DGE_obj_5327, log=TRUE, prior.count=3)))
  logCPM_5327$subtype=pheno_sample_GSE5327
  surv_dat5327$GEO_id
  surv_limma5327 = merge(surv_dat5327, 
                         logCPM_5327[,limmaG],
                         by.x="GEO_id",by.y="row.names")
  
  
  pred = predict(limma_fit,as.matrix(surv_limma5327[,limmaG]),type="response")
  tmp=Cindex(pred,Surv(surv_limma5327$surv_time,surv_limma5327$vital_num))
  cindex_list5327 = c(cindex_list5327, tmp)
  
  ## edgeR
  # Normalization
  DGE_obj_5327 <- DGEList(counts = 2^GeneExp_BRCA_GSE5327 - 1,
                          group = pheno_sample_GSE5327)
  keep.exprs.5327 <- filterByExpr(DGE_obj_5327)
  DGE_obj_5327 <- calcNormFactors(DGE_obj_5327)
  logCPM_5327 <- as.data.frame(t(cpm(DGE_obj_5327, log=TRUE, prior.count=3)))
  logCPM_5327$subtype=pheno_sample_GSE5327
  surv_edgeR5327 = merge(surv_dat5327, 
                         logCPM_5327[,edgeRG],
                         by.x="GEO_id",by.y="row.names")
  pred = predict(edgeR_fit,as.matrix(surv_edgeR5327[,edgeRG]),type="response")
  tmp=Cindex(pred,Surv(surv_edgeR5327$surv_time,surv_edgeR5327$vital_num))
  cindex_list5327 = c(cindex_list5327, tmp)
  
  
  ## wilcoxon
  # Normalization
  DGE_obj_5327 <- DGEList(counts = 2^GeneExp_BRCA_GSE5327 - 1,
                          group = pheno_sample_GSE5327)
  keep.exprs.5327 <- filterByExpr(DGE_obj_5327)
  DGE_obj_5327 <- calcNormFactors(DGE_obj_5327)
  logCPM_5327 <- as.data.frame(t(cpm(DGE_obj_5327, log=TRUE, prior.count=3)))
  logCPM_5327$subtype=pheno_sample_GSE5327
  surv_wilcoxon5327 = merge(surv_dat5327, 
                            logCPM_5327[,wilcoxonG],
                            by.x="GEO_id",by.y="row.names")
  pred = predict(wilcoxon_fit,as.matrix(surv_wilcoxon5327[,wilcoxonG]),type="response")
  tmp=Cindex(pred,Surv(surv_wilcoxon5327$surv_time,surv_wilcoxon5327$vital_num))
  cindex_list5327 = c(cindex_list5327, tmp)
  
  ## DESeq2
  # Normalization
  DGE_obj_5327 <- DGEList(counts = 2^GeneExp_BRCA_GSE5327 - 1,
                          group = pheno_sample_GSE5327)
  keep.exprs.5327 <- filterByExpr(DGE_obj_5327)
  DGE_obj_5327 <- calcNormFactors(DGE_obj_5327)
  logCPM_5327 <- as.data.frame(t(cpm(DGE_obj_5327, log=TRUE, prior.count=3)))
  logCPM_5327$subtype=pheno_sample_GSE5327
  surv_DESeq25327 = merge(surv_dat5327,
                          logCPM_5327[,DESeq2G],
                          by.x="GEO_id",by.y="row.names")
  pred = predict(DESeq2_fit,as.matrix(surv_DESeq25327[,DESeq2G]),type="response")
  tmp=Cindex(pred,Surv(surv_DESeq25327$surv_time,surv_DESeq25327$vital_num))
  cindex_list5327 = c(cindex_list5327, tmp)
  
  # peap
  genePairMatrix_5327 =  matrix(unlist(strsplit(peapP,"_")), ncol = 2, byrow = TRUE)
  gene1 = GeneExp_BRCA_GSE5327[as.vector(genePairMatrix_5327[,1]),]
  row.names(gene1) = NULL
  gene2 = GeneExp_BRCA_GSE5327[as.vector(genePairMatrix_5327[,2]),]
  row.names(gene2) = NULL
  pheno_row_foldchange_GSE5327 = gene1 - gene2
  row.names(pheno_row_foldchange_GSE5327) = peapP
  DEP_df = as.data.frame(t(pheno_row_foldchange_GSE5327))
  DEP_df$subtype = pheno_sample_GSE5327
  surv_peap5327 = merge(surv_dat5327, 
                        DEP_df[,c(peapP,"subtype")],
                        by.x="GEO_id",by.y="row.names")
  
  
  pred = predict(peap_fit,as.matrix(surv_peap5327[,peapP]),type="response")
  tmp=Cindex(pred,Surv(surv_peap5327$surv_time,surv_peap5327$vital_num))
  cindex_list5327 = c(cindex_list5327, tmp)
  
  # pam50
  DGE_obj_5327 <- DGEList(counts = 2^GeneExp_BRCA_GSE5327 - 1,
                          group = pheno_sample_GSE5327)
  keep.exprs.5327 <- filterByExpr(DGE_obj_5327)
  DGE_obj_5327 <- calcNormFactors(DGE_obj_5327)
  logCPM_5327 <- as.data.frame(t(cpm(DGE_obj_5327, log=TRUE, prior.count=3)))
  logCPM_5327$subtype=pheno_sample_GSE5327
  surv_pam505327 = merge(surv_dat5327,
                         logCPM_5327[,pam50_g],
                         by.x="GEO_id",by.y="row.names")
  pred = predict(pam50_fit,as.matrix(surv_pam505327[,pam50_g]),type="response")
  tmp=Cindex(pred,Surv(surv_pam505327$surv_time,surv_pam505327$vital_num))
  cindex_list5327 = c(cindex_list5327, tmp) 
  
  # wilcoxon rank sum test
  
}


print("Combined")
surv_datCOMBINED <- rbind(surv_dat2034[,c("GEO_id","vital_num", "surv_time")], surv_dat5327)
GeneExp_BRCA_COMBINED <- cbind(GeneExp_BRCA_GSE2034,GeneExp_BRCA_GSE5327)
pheno_sample_COMBINED <- c(pheno_sample_GSE2034, pheno_sample_GSE5327)
cindex_listCOMBINED = list()
for (year_limit in 1:5) {
  print(c("Year Limit: ", year_limit))
  ## GSE2034
  GSE2034_survival <- read_excel("GSE2034_survival.xlsx")
  surv_dat2034 = GSE2034_survival[,2:7]
  colnames(surv_dat2034) = c("GEO_id","node","surv_time",
                             "vital_num","ER_status","Brain_relapses")
  surv_dat5327 <- read_excel("GSE5327_survival.xlsx")
  surv_datCOMBINED <- rbind(surv_dat2034[,c("GEO_id","vital_num", "surv_time")], surv_dat5327)
  surv_datCOMBINED$vital_num= ifelse(surv_datCOMBINED$surv_time>year_limit*12,
                                 0,
                                 surv_datCOMBINED$vital_num)
  surv_datCOMBINED$surv_time= ifelse(surv_datCOMBINED$surv_time>year_limit*12,
                                 year_limit*12,
                                 surv_datCOMBINED$surv_time)
  ## limma
  # Normalization
  DGE_obj_COMBINED <- DGEList(counts = 2^GeneExp_BRCA_COMBINED - 1,
                          group = pheno_sample_COMBINED)
  keep.exprs.COMBINED <- filterByExpr(DGE_obj_COMBINED)
  DGE_obj_COMBINED <- calcNormFactors(DGE_obj_COMBINED)
  logCPM_COMBINED <- as.data.frame(t(cpm(DGE_obj_COMBINED, log=TRUE, prior.count=3)))
  logCPM_COMBINED$subtype=pheno_sample_COMBINED
  surv_datCOMBINED$GEO_id
  surv_limmaCOMBINED = merge(surv_datCOMBINED, 
                         logCPM_COMBINED[,limmaG],
                         by.x="GEO_id",by.y="row.names")
  
  
  pred = predict(limma_fit,as.matrix(surv_limmaCOMBINED[,limmaG]),type="response")
  tmp=Cindex(pred,Surv(surv_limmaCOMBINED$surv_time,surv_limmaCOMBINED$vital_num))
  cindex_listCOMBINED = c(cindex_listCOMBINED, tmp)
  
  ## edgeR
  # Normalization
  DGE_obj_COMBINED <- DGEList(counts = 2^GeneExp_BRCA_COMBINED - 1,
                          group = pheno_sample_COMBINED)
  keep.exprs.COMBINED <- filterByExpr(DGE_obj_COMBINED)
  DGE_obj_COMBINED <- calcNormFactors(DGE_obj_COMBINED)
  logCPM_COMBINED <- as.data.frame(t(cpm(DGE_obj_COMBINED, log=TRUE, prior.count=3)))
  logCPM_COMBINED$subtype=pheno_sample_COMBINED
  surv_edgeRCOMBINED = merge(surv_datCOMBINED, 
                         logCPM_COMBINED[,edgeRG],
                         by.x="GEO_id",by.y="row.names")
  pred = predict(edgeR_fit,as.matrix(surv_edgeRCOMBINED[,edgeRG]),type="response")
  tmp=Cindex(pred,Surv(surv_edgeRCOMBINED$surv_time,surv_edgeRCOMBINED$vital_num))
  cindex_listCOMBINED = c(cindex_listCOMBINED, tmp)
  
  
  ## wilcoxon
  # Normalization
  DGE_obj_COMBINED <- DGEList(counts = 2^GeneExp_BRCA_COMBINED - 1,
                          group = pheno_sample_COMBINED)
  keep.exprs.COMBINED <- filterByExpr(DGE_obj_COMBINED)
  DGE_obj_COMBINED <- calcNormFactors(DGE_obj_COMBINED)
  logCPM_COMBINED <- as.data.frame(t(cpm(DGE_obj_COMBINED, log=TRUE, prior.count=3)))
  logCPM_COMBINED$subtype=pheno_sample_COMBINED
  surv_wilcoxonCOMBINED = merge(surv_datCOMBINED, 
                            logCPM_COMBINED[,wilcoxonG],
                            by.x="GEO_id",by.y="row.names")
  pred = predict(wilcoxon_fit,as.matrix(surv_wilcoxonCOMBINED[,wilcoxonG]),type="response")
  tmp=Cindex(pred,Surv(surv_wilcoxonCOMBINED$surv_time,surv_wilcoxonCOMBINED$vital_num))
  cindex_listCOMBINED = c(cindex_listCOMBINED, tmp)
  
  ## DESeq2
  # Normalization
  DGE_obj_COMBINED <- DGEList(counts = 2^GeneExp_BRCA_COMBINED - 1,
                          group = pheno_sample_COMBINED)
  keep.exprs.COMBINED <- filterByExpr(DGE_obj_COMBINED)
  DGE_obj_COMBINED <- calcNormFactors(DGE_obj_COMBINED)
  logCPM_COMBINED <- as.data.frame(t(cpm(DGE_obj_COMBINED, log=TRUE, prior.count=3)))
  logCPM_COMBINED$subtype=pheno_sample_COMBINED
  surv_DESeq2COMBINED = merge(surv_datCOMBINED,
                          logCPM_COMBINED[,DESeq2G],
                          by.x="GEO_id",by.y="row.names")
  pred = predict(DESeq2_fit,as.matrix(surv_DESeq2COMBINED[,DESeq2G]),type="response")
  tmp=Cindex(pred,Surv(surv_DESeq2COMBINED$surv_time,surv_DESeq2COMBINED$vital_num))
  cindex_listCOMBINED = c(cindex_listCOMBINED, tmp)
  
  # peap
  genePairMatrix_COMBINED =  matrix(unlist(strsplit(peapP,"_")), ncol = 2, byrow = TRUE)
  gene1 = GeneExp_BRCA_COMBINED[as.vector(genePairMatrix_COMBINED[,1]),]
  row.names(gene1) = NULL
  gene2 = GeneExp_BRCA_COMBINED[as.vector(genePairMatrix_COMBINED[,2]),]
  row.names(gene2) = NULL
  pheno_row_foldchange_COMBINED = gene1 - gene2
  row.names(pheno_row_foldchange_COMBINED) = peapP
  DEP_df = as.data.frame(t(pheno_row_foldchange_COMBINED))
  DEP_df$subtype = pheno_sample_COMBINED
  surv_peapCOMBINED = merge(surv_datCOMBINED, 
                        DEP_df[,c(peapP,"subtype")],
                        by.x="GEO_id",by.y="row.names")
  
  
  pred = predict(peap_fit,as.matrix(surv_peapCOMBINED[,peapP]),type="response")
  tmp=Cindex(pred,Surv(surv_peapCOMBINED$surv_time,surv_peapCOMBINED$vital_num))
  cindex_listCOMBINED = c(cindex_listCOMBINED, tmp)
  
  # pam50
  DGE_obj_COMBINED <- DGEList(counts = 2^GeneExp_BRCA_COMBINED - 1,
                          group = pheno_sample_COMBINED)
  keep.exprs.COMBINED <- filterByExpr(DGE_obj_COMBINED)
  DGE_obj_COMBINED <- calcNormFactors(DGE_obj_COMBINED)
  logCPM_COMBINED <- as.data.frame(t(cpm(DGE_obj_COMBINED, log=TRUE, prior.count=3)))
  logCPM_COMBINED$subtype=pheno_sample_COMBINED
  surv_pam50COMBINED = merge(surv_datCOMBINED,
                         logCPM_COMBINED[,pam50_g],
                         by.x="GEO_id",by.y="row.names")
  pred = predict(pam50_fit,as.matrix(surv_pam50COMBINED[,pam50_g]),type="response")
  tmp=Cindex(pred,Surv(surv_pam50COMBINED$surv_time,surv_pam50COMBINED$vital_num))
  cindex_listCOMBINED = c(cindex_listCOMBINED, tmp)
}


# print("SCANB")
# load("surv_SCANB.Rdata")
# surv_datSCANB$surv_time = surv_datSCANB$surv_time/30
# cindex_listSCANB=list()
# for (year_limit in 1:10) {
#   print(c("Year Limit: ", year_limit))
#   ## GSE2034
#   load("surv_SCANB.Rdata")
# 
#   surv_datSCANB$vital_num= ifelse(surv_datSCANB$surv_time>year_limit*12,
#                                  0,
#                                  1)
#   surv_datSCANB$surv_time= ifelse(surv_datSCANB$surv_time>year_limit*12,
#                                      year_limit*12,
#                                      surv_datSCANB$surv_time)
# 
#   ## limma
#   # Normalization
#   DGE_obj_SCANB <- DGEList(counts = 2^GeneExp_BRCA_SCANB - 1,
#                           group = pheno_sample_SCANB)
#   keep.exprs.SCANB <- filterByExpr(DGE_obj_SCANB)
#   DGE_obj_SCANB <- calcNormFactors(DGE_obj_SCANB)
#   logCPM_SCANB <- as.data.frame(t(cpm(DGE_obj_SCANB, log=TRUE, prior.count=3)))
#   logCPM_SCANB$subtype=pheno_sample_SCANB
#   surv_datSCANB$Patient_id
#   surv_limmaSCANB = merge(surv_datSCANB, 
#                          logCPM_SCANB[,limmaG],
#                          by.x="Patient_id",by.y="row.names")
#   
#   
#   pred = predict(limma_fit,as.matrix(surv_limmaSCANB[,limmaG]),type="response")
#   tmp=Cindex(pred,Surv(surv_limmaSCANB$surv_time,surv_limmaSCANB$vital_num))
#   cindex_listSCANB = c(cindex_listSCANB, tmp)
#   
#   ## edgeR
#   # Normalization
#   DGE_obj_SCANB <- DGEList(counts = 2^GeneExp_BRCA_SCANB - 1,
#                           group = pheno_sample_SCANB)
#   keep.exprs.SCANB <- filterByExpr(DGE_obj_SCANB)
#   DGE_obj_SCANB <- calcNormFactors(DGE_obj_SCANB)
#   logCPM_SCANB <- as.data.frame(t(cpm(DGE_obj_SCANB, log=TRUE, prior.count=3)))
#   logCPM_SCANB$subtype=pheno_sample_SCANB
#   surv_edgeRSCANB = merge(surv_datSCANB, 
#                          logCPM_SCANB[,edgeRG],
#                          by.x="Patient_id",by.y="row.names")
#   pred = predict(edgeR_fit,as.matrix(surv_edgeRSCANB[,edgeRG]),type="response")
#   tmp=Cindex(pred,Surv(surv_edgeRSCANB$surv_time,surv_edgeRSCANB$vital_num))
#   cindex_listSCANB = c(cindex_listSCANB, tmp)
#   
#   ## DESeq2
#   # Normalization
#   DGE_obj_SCANB <- DGEList(counts = 2^GeneExp_BRCA_SCANB - 1,
#                           group = pheno_sample_SCANB)
#   keep.exprs.SCANB <- filterByExpr(DGE_obj_SCANB)
#   DGE_obj_SCANB <- calcNormFactors(DGE_obj_SCANB)
#   logCPM_SCANB <- as.data.frame(t(cpm(DGE_obj_SCANB, log=TRUE, prior.count=3)))
#   logCPM_SCANB$subtype=pheno_sample_SCANB
#   surv_DESeq2SCANB = merge(surv_datSCANB,
#                           logCPM_SCANB[,DESeq2G],
#                           by.x="Patient_id",by.y="row.names")
#   pred = predict(DESeq2_fit,as.matrix(surv_DESeq2SCANB[,DESeq2G]),type="response")
#   tmp=Cindex(pred,Surv(surv_DESeq2SCANB$surv_time,surv_DESeq2SCANB$vital_num))
#   cindex_listSCANB = c(cindex_listSCANB, tmp)
#   
#   # peap
#   genePairMatrix_SCANB =  matrix(unlist(strsplit(peapP,"_")), ncol = 2, byrow = TRUE)
#   gene1 = GeneExp_BRCA_SCANB[as.vector(genePairMatrix_SCANB[,1]),]
#   row.names(gene1) = NULL
#   gene2 = GeneExp_BRCA_SCANB[as.vector(genePairMatrix_SCANB[,2]),]
#   row.names(gene2) = NULL
#   pheno_row_foldchange_SCANB = gene1 - gene2
#   row.names(pheno_row_foldchange_SCANB) = peapP
#   DEP_df = as.data.frame(t(pheno_row_foldchange_SCANB))
#   DEP_df$subtype = pheno_sample_SCANB
#   surv_peapSCANB = merge(surv_datSCANB, 
#                         DEP_df[,c(peapP,"subtype")],
#                         by.x="Patient_id",by.y="row.names")
#   
#   
#   pred = predict(peap_fit,as.matrix(surv_peapSCANB[,peapP]),type="response")
#   tmp=Cindex(pred,Surv(surv_peapSCANB$surv_time,surv_peapSCANB$vital_num))
#   cindex_listSCANB = c(cindex_listSCANB, tmp)
#   
#   # pam50
#   DGE_obj_SCANB <- DGEList(counts = 2^GeneExp_BRCA_SCANB - 1,
#                           group = pheno_sample_SCANB) 
#   keep.exprs.SCANB <- filterByExpr(DGE_obj_SCANB)
#   DGE_obj_SCANB <- calcNormFactors(DGE_obj_SCANB)
#   logCPM_SCANB <- as.data.frame(t(cpm(DGE_obj_SCANB, log=TRUE, prior.count=3)))
#   logCPM_SCANB$subtype=pheno_sample_SCANB
#   surv_pam50SCANB = merge(surv_datSCANB,
#                          logCPM_SCANB[,pam50_g],
#                          by.x="Patient_id",by.y="row.names")
#   pred = predict(pam50_fit,as.matrix(surv_pam50SCANB[,pam50_g]),type="response")
#   tmp=Cindex(pred,Surv(surv_pam50SCANB$surv_time,surv_pam50SCANB$vital_num))
#   cindex_listSCANB = c(cindex_listSCANB, tmp)
# }


## Outputing
save.image(paste0("result20241103.RData"))
method_vec = c("limma", "edgeR","Wilcoxon", "DESeq2", "PEAP", "PAM50")
df = data.frame()
for (year in 1:5) {
  for (m in 1:6) {
    if (m==1){
      differential_dim = length(unique(unlist(limmaG_list)))
      coe = limma_fit$beta
      lasso_dim = length(coe[coe!=0])
    } else if (m==2){
      differential_dim = length(unique(unlist(edgeRG_list)))
      coe = edgeR_fit$beta
      lasso_dim = length(coe[coe!=0])
    } else if (m==3){
      differential_dim = length(unique(unlist(wilcoxonG_list)))
      coe = wilcoxon_fit$beta
      lasso_dim = length(coe[coe!=0])
    } else if (m==4){
      differential_dim = length(unique(unlist(DESeq2G_list)))
      coe = DESeq2_fit$beta
      lasso_dim = length(coe[coe!=0])
    } else if (m==5){
      differential_dim = length(unique(unlist(peapP_list)))
      coe = peap_fit$beta
      lasso_dim = length(coe[coe!=0])
    } else {
      differential_dim = length(pam50_g)
      coe = pam50_fit$beta
      lasso_dim = length(coe[coe!=0])
    }
    
    df = rbind(df, data.frame(
      subtype_limit = lambda_fix,
      Algo = method_vec[m],
      year_lim = year,
      cindex2034_cox = cindex_list2034[[m+6*(year-1)]],
      cindex5327_cox = cindex_list5327[[m+6*(year-1)]],
      cindexCOMBINED_cox = cindex_listCOMBINED[[m+6*(year-1)]],
      differential_dim = differential_dim,
      lasso_dim = lasso_dim
    ))
  }
}

write.csv(df, file = "result_20241103.csv")



