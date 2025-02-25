rm(list=ls())
library(ggplot2)
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
library(pROC)
library(sva)

set.seed(1234)

setwd("D:/PEAP/Infection/")
load("Infection_PEAP03_output_012.Rdata")
name_pheno <- c("0", "1", "2")


## ---------------------------------------------------------------------------------------------
lassoPreprocessBinary = function(df, target){
  # Input a data frame with a phenotype column (including all 5 phenotypes)
  # Transfer it into "binary" phenotype according to "target"
  # Output x, y suitable for glmnet
  df_copy = df
  df_copy$phenoNumeric = as.numeric(df_copy$phenotype==target)
  x = as.matrix(df_copy[,-c(which(colnames(df_copy)=="phenotype" | colnames(df_copy)=="phenoNumeric"))])
  y = unlist(df_copy[,"phenoNumeric"])
  return(list(x,y))
}

matching = function(x, matchMode){
  # If matchMode is 0: return numerical value representing phenotype
  # If matchMode is 1: return phenotype according to numerical value
  phenoV = c("0", "2", "11", "12")
  ans = ifelse(matchMode==0,
               sum(which(phenoV==x)), phenoV[x])
  return(ans)
} 


lassoPreprocessMulti = function(df){
  # Input a data frame with a phenotype column (including all 5 phenotypes)
  # Output x, y suitable for glmnet
  df_copy = df
  df_copy$phenoNumeric = lapply(df_copy$phenotype,FUN = matching, matchMode=0)
  x = as.matrix(df_copy[,-c(which(colnames(df_copy)=="phenotype" | colnames(df_copy)=="phenoNumeric"))])
  y = unlist(df_copy[,"phenoNumeric"])
  return(list(x,y))
}


## ---------------------------------------------------------------------------------------------
DEP_DEG_result <- as.data.frame(matrix(rep(0,20), nrow=5))
colnames(DEP_DEG_result) <- c("DEP","Limma","edgeR","DESeq2")
rownames(DEP_DEG_result) <- c("ACC","AUC0","AUC1","AUC2","feature_number")

all_results <- list()

lambda_seq <- c(10^(-6),10^(-4),10^(-2),10^(-1))
lambda_lasso <- lambda_seq[2]
lambda_lasso = 0.1
lasso_alpha = 0.5


## ---------------------------------------------------------------------------------------------
Selected_genepairs = paste(top100_Genes[,1],"_", top100_Genes[,2],sep="")

genePair = unique(names(table(Selected_genepairs)))
genePairMatrix =  matrix(unlist(strsplit(genePair,"_")), ncol = 2, byrow = TRUE)

gene1 = AllGeneExp[as.vector(genePairMatrix[,1]),]
row.names(gene1) = NULL
gene2 = AllGeneExp[as.vector(genePairMatrix[,2]),]
row.names(gene2) = NULL
pheno_row_foldchange_1 = gene1 - gene2
row.names(pheno_row_foldchange_1) = genePair
DEP_df = as.data.frame(t(pheno_row_foldchange_1))
DEP_df$phenotype = pheno_sample


## ---------------------------------------------------------------------------------------------
set.seed(1234)  # Set a seed for reproducibility
train_indices <- createDataPartition(DEP_df$phenotype, p = 1, list = FALSE)  # 70% for training, adjust as needed
train_data <- DEP_df[train_indices, ]
test_data <- DEP_df[-train_indices, ]
# DEP_df_test <- test_data


## ---------------------------------------------------------------------------------------------
DEP_lambda <- lambda_lasso

peapP_list = list()
selecP_list = list()
for (i in c(1:3)){
  peapP = paste(top100_Genes[(1000*i-999):(1000*i), 1],"_",
                             top100_Genes[(1000*i-999):(1000*i), 2],sep="")
  peapP_list = append(peapP_list, list(peapP))
  
  # lasso
  DEPxy = lassoPreprocessBinary(train_data[,c(peapP,"phenotype")],
                                name_pheno[i])
  # cvLassoDEP = cv.glmnet(DEPxy[[1]], DEPxy[[2]],
  #                        alpha = 1, nfolds = 10)
  # plot(cvLassoDEP)
  # DEP_lambda = cvLassoDEP$lambda.min
  finalLassoDEP = glmnet(DEPxy[[1]], DEPxy[[2]],
                       alpha = lasso_alpha, lambda=DEP_lambda) 
  result = as.matrix(coef(finalLassoDEP))
  selecP = names(result[result[,1]!=0,])[-1]
  selecP_list = append(selecP_list, list(selecP))
}
selecP = unique(unlist(selecP_list))
DEP_lasso = train_data[,selecP]
DEP_lasso$phenotype = pheno_sample[train_indices]



## ---------------------------------------------------------------------------------------------
test_set_Gene <- read.csv("data21802.csv", header = T)
test_set_Meta <- read.csv(file = "label_GSE21802.csv", header = T)

test_set_Gene <- test_set_Gene[,-1]

rownames(test_set_Gene) <- test_set_Gene[,1]
test_set_Gene <- test_set_Gene[,-1]

test_set_Gene <- log2(test_set_Gene - min(test_set_Gene) + 1)

match(test_set_Meta[,1], colnames(test_set_Gene))

genePairMatrix =  matrix(unlist(strsplit(colnames(DEP_lasso)[-ncol(DEP_lasso)],"_")), ncol = 2, byrow = TRUE)
gene1 = test_set_Gene[as.vector(genePairMatrix[,1]),]
row.names(gene1) = NULL
gene2 = test_set_Gene[as.vector(genePairMatrix[,2]),]
row.names(gene2) = NULL
pheno_row_foldchange_1 = gene1 - gene2
row.names(pheno_row_foldchange_1) = colnames(DEP_lasso)[-ncol(DEP_lasso)]

DEP_df_test = as.data.frame(t(pheno_row_foldchange_1))
DEP_df_test$phenotype = test_set_Meta$label

DEP_df_test_1 <- DEP_df_test

test_set_Gene <- read.csv("data57065.csv", header = T)
test_set_Meta <- read.csv(file = "label_GSE57065.csv", header = T)

test_set_Gene <- test_set_Gene[,-1]


rownames(test_set_Gene) <- test_set_Gene[,1]
test_set_Gene <- test_set_Gene[,-1]

match(test_set_Meta[,1], colnames(test_set_Gene))

genePairMatrix =  matrix(unlist(strsplit(colnames(DEP_lasso)[-ncol(DEP_lasso)],"_")), ncol = 2, byrow = TRUE)
gene1 = test_set_Gene[as.vector(genePairMatrix[,1]),]
row.names(gene1) = NULL
gene2 = test_set_Gene[as.vector(genePairMatrix[,2]),]
row.names(gene2) = NULL
pheno_row_foldchange_1 = gene1 - gene2
row.names(pheno_row_foldchange_1) = colnames(DEP_lasso)[-ncol(DEP_lasso)]

DEP_df_test = as.data.frame(t(pheno_row_foldchange_1))
DEP_df_test$phenotype = test_set_Meta$label

DEP_df_test_2 <- DEP_df_test
DEP_df_test_combined <- rbind(DEP_df_test_1, DEP_df_test_2)


## ---------------------------------------------------------------------------------------------
set.seed(123)
gbm_mod = gbm(phenotype ~.,
              data = DEP_lasso,
              distribution = "multinomial",
              cv.folds = 10,
              shrinkage = .01,
              n.trees = 1000,
              verbose = T)
gbm.perf(gbm_mod)


pred = predict.gbm(object = gbm_mod,
                   newdata = DEP_df_test_combined,
                   n.trees = 1000,
                   type = "response")



labels = colnames(pred)[apply(pred, 1, which.max)]

result = data.frame(DEP_df_test_combined$phenotype, b = labels)

cm_DEP = confusionMatrix(data = as.factor(as.numeric(labels)), reference= as.factor(DEP_df_test_combined$phenotype))
cm_DEP

result <- as.data.frame(result)


DEP_DEG_result[1,1] <- cm_DEP$overall[1]
DEP_DEG_result[2,1] <- roc(result[,1]==0, pred[,,1][,1])$auc[1]
DEP_DEG_result[3,1] <- roc(result[,1]==1, pred[,,1][,2])$auc[1]
DEP_DEG_result[4,1] <- roc(result[,1]==2, pred[,,1][,3])$auc[1]
DEP_DEG_result[5,1] <- ncol(DEP_lasso)-1

save(DEP_lasso, DEP_df_test_combined, gbm_mod, pred, cm_DEP, result, file = paste0("DEP_lambda_",lambda_lasso,".Rdata"))


## ---------------------------------------------------------------------------------------------
DEG_lambda <- lambda_lasso
# Normalization
set.seed(123)
# custom_lambda <- seq((1/12), 1, by = (1/12))

DGE_obj <- DGEList(counts = 2^AllGeneExp-1, group = pheno_sample)
keep.exprs <- filterByExpr(DGE_obj)
DGE_obj <- DGE_obj[keep.exprs,, keep.lib.sizes=FALSE]
DGE_obj <- calcNormFactors(DGE_obj)
logCPM <- cpm(DGE_obj, log=TRUE, prior.count=3)

DEG_df_list = list()
limmaG_list = list()
selecG_list = list()
top_gene_table_list = list()

# Main: one vs rest
for(i in 1:3){ 
  pheno = name_pheno[i]
  print(pheno)
  # limma
  pheno_sample_other <- pheno_sample
  for(j in 1:length(pheno_sample)){
    if(pheno!=pheno_sample_other[j]){
      pheno_sample_other[j] <- "other"
    }
  }
  design <- model.matrix(~pheno_sample_other) 
  colnames(design) <- gsub("pheno_sample_other", "", colnames(design))
  fit <- lmFit(logCPM, design) 
  fit <- eBayes(fit, trend=TRUE) 
  top_gene_table <- topTable(fit, coef=2,
                             number = 6283,adjust.method = "BH",p.value = 0.005)
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
  DEG_df_eachPheno$phenotype = pheno_sample
  DEG_df_list = append(DEG_df_list, list(DEG_df_eachPheno))
  limmaGene = colnames(DEG_df_eachPheno)
  limmaG_list = append(limmaG_list, list(limmaGene))
  
  # lasso
  DEGxy = lassoPreprocessBinary(DEG_df_eachPheno, name_pheno[i])
  # cvLassoDEG = cv.glmnet(DEGxy[[1]], DEGxy[[2]],
  #                        alpha = 1, nfolds = 10, lambda = custom_lambda)
  # plot(cvLassoDEG)
  # DEG_lambda = cvLassoDEG$lambda.min
  # print(DEG_lambda)
  finalLassoDEG = glmnet(DEGxy[[1]], DEGxy[[2]],
                       alpha = lasso_alpha, lambda=DEG_lambda) 
  result = as.matrix(coef(finalLassoDEG))
  selecG = names(result[result[,1]!=0,])[-1]
  selecG_list = append(selecG_list, list(selecG))
}
limmaG = unique(unlist(limmaG_list))
selecG = unique(unlist(selecG_list))
DEG_lasso = as.data.frame(t(logCPM)[,selecG])
DEG_lasso$phenotype = pheno_sample


## ---------------------------------------------------------------------------------------------
test_set_Gene <- read.csv("data21802.csv", header = T)
test_set_Meta <- read.csv(file = "label_GSE21802.csv", header = T)

test_set_Gene <- test_set_Gene[,-1]


rownames(test_set_Gene) <- test_set_Gene[,1]
test_set_Gene <- test_set_Gene[,-1]

test_set_Gene <- log2(test_set_Gene - min(test_set_Gene) + 1)


match(test_set_Meta[,1], colnames(test_set_Gene))

DEG_df_test <- as.data.frame(matrix(nrow = ncol(test_set_Gene), ncol = ncol(DEG_lasso)))
colnames(DEG_df_test) <- colnames(DEG_lasso)
rownames(DEG_df_test) <- colnames(test_set_Gene)

for(i in 1:ncol(DEG_df_test)){
  if(colnames(DEG_df_test)[i] %in% rownames(test_set_Gene)){
    DEG_df_test[,i] <- t(test_set_Gene[which(rownames(test_set_Gene)==colnames(DEG_lasso)[i]),])
  }
}
colnames(DEG_df_test) <- as.character(colnames(DEG_lasso)) 

DEG_df_test$phenotype <-  test_set_Meta$label
DEG_df_test_1 <- DEG_df_test


test_set_Gene <- read.csv("data57065.csv", header = T)
test_set_Meta <- read.csv(file = "label_GSE57065.csv", header = T)

test_set_Gene <- test_set_Gene[,-1]

rownames(test_set_Gene) <- test_set_Gene[,1]
test_set_Gene <- test_set_Gene[,-1]

DEG_df_test <- as.data.frame(matrix(nrow = ncol(test_set_Gene), ncol = ncol(DEG_lasso)))
colnames(DEG_df_test) <- colnames(DEG_lasso)
rownames(DEG_df_test) <- colnames(test_set_Gene)

for(i in 1:ncol(DEG_df_test)){
  if(colnames(DEG_df_test)[i] %in% rownames(test_set_Gene)){
    DEG_df_test[,i] <- t(test_set_Gene[which(rownames(test_set_Gene)==colnames(DEG_lasso)[i]),])
  }
}
colnames(DEG_df_test) <- as.character(colnames(DEG_lasso)) 

DEG_df_test$phenotype <-  test_set_Meta$label
DEG_df_test_2 <- DEG_df_test

DEG_df_test_combined <- rbind(DEG_df_test_1, DEG_df_test_2)


## ---------------------------------------------------------------------------------------------
gbm_mod = gbm(phenotype ~.,
              data = DEG_lasso,
              distribution = "multinomial",
              cv.folds = 10,
              shrinkage = .01,
              n.trees = 1000,
              verbose = T)
gbm.perf(gbm_mod)


pred = predict.gbm(object = gbm_mod,
                   newdata = as.data.frame(DEG_df_test_combined),
                   n.trees = 1000,
                   type = "response")
  
labels = colnames(pred)[apply(pred, 1, which.max)]
result = data.frame(DEG_df_test_combined$phenotype, b = labels)
  
cm_DEG = confusionMatrix(data = as.factor(as.numeric(labels)), reference= as.factor(DEG_df_test_combined$phenotype))
cm_DEG


DEP_DEG_result[1,2] <- cm_DEG$overall[1]
DEP_DEG_result[2,2] <- roc(result[,1]==0, pred[,,1][,1])$auc[1]
DEP_DEG_result[3,2] <- roc(result[,1]==1, pred[,,1][,2])$auc[1]
DEP_DEG_result[4,2] <- roc(result[,1]==2, pred[,,1][,3])$auc[1]
DEP_DEG_result[5,2] <- ncol(DEG_lasso)-1


save(DEG_lasso, DEG_df_test_combined, gbm_mod, pred, cm_DEG, result, file = paste0("DEP_limma_",lambda_lasso,".Rdata"))


## ---------------------------------------------------------------------------------------------
DEG_lambda <- lambda_lasso

set.seed(123)
custom_lambda <- seq((1/12), 1, by = (1/12))

# Normalization
pheno_sample <- as.character(pheno_sample)
subtype_limit = 500
DGE_obj <- DGEList(counts = 2^AllGeneExp - 1, group = pheno_sample)
keep.exprs <- filterByExpr(DGE_obj)
DGE_obj <- DGE_obj[keep.exprs,, keep.lib.sizes=FALSE]
DGE_obj <- calcNormFactors(DGE_obj)
DGE_obj$group <- pheno_sample

# calculate preparations
design <- model.matrix(~pheno_sample)
colnames(design)=levels(as.factor(pheno_sample))
DGE_obj <- estimateDisp(DGE_obj, design)
logCPM <- cpm(DGE_obj, log=TRUE, prior.count=3)
N <- ncol(logCPM)
fit <- glmQLFit(DGE_obj, design)

DEG_df_list = list()
edgeRG_list = list()
selecG_edgeR_list = list()
top_gene_table_list = list()
# extract top genes and do lasso
for (i in 1:3) {
  contrast_vec = c(0,0,0)
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
  DEG_df_eachPheno$phenotype = pheno_sample
  DEG_df_list = append(DEG_df_list, list(DEG_df_eachPheno))
  edgeRGene = colnames(DEG_df_eachPheno)[-length(colnames(DEG_df_eachPheno))]
  edgeRG_list = append(edgeRG_list, list(edgeRGene))
  
  # lasso
  DEGxy = lassoPreprocessBinary(DEG_df_eachPheno, name_pheno[i])
  # cvLassoDEG = cv.glmnet(DEGxy[[1]], DEGxy[[2]],
  #                        alpha = 1, nfolds = 10, lambda = custom_lambda)
  # plot(cvLassoDEG)
  # DEG_lambda = cvLassoDEG$lambda.min
  # print(DEG_lambda)
  finalLassoDEG = glmnet(DEGxy[[1]], DEGxy[[2]],
                         alpha = lasso_alpha, lambda=DEG_lambda) 
  result = as.matrix(coef(finalLassoDEG))
  selecG = names(result[result[,1]!=0,])[-1]
  selecG_edgeR_list = append(selecG_edgeR_list, list(selecG))
}
edgeRG = unique(unlist(edgeRG_list))
selecG = unique(unlist(selecG_edgeR_list))
DEG_lasso = as.data.frame(t(logCPM)[,selecG])
DEG_lasso$phenotype = pheno_sample


## ---------------------------------------------------------------------------------------------
test_set_Gene <- read.csv("data21802.csv", header = T)
test_set_Meta <- read.csv(file = "label_GSE21802.csv", header = T)

test_set_Gene <- test_set_Gene[,-1]


rownames(test_set_Gene) <- test_set_Gene[,1]
test_set_Gene <- test_set_Gene[,-1]

test_set_Gene <- log2(test_set_Gene - min(test_set_Gene) + 1)


match(test_set_Meta[,1], colnames(test_set_Gene))

DEG_df_test <- as.data.frame(matrix(nrow = ncol(test_set_Gene), ncol = ncol(DEG_lasso)))
colnames(DEG_df_test) <- colnames(DEG_lasso)
rownames(DEG_df_test) <- colnames(test_set_Gene)

for(i in 1:ncol(DEG_df_test)){
  if(colnames(DEG_df_test)[i] %in% rownames(test_set_Gene)){
    DEG_df_test[,i] <- t(test_set_Gene[which(rownames(test_set_Gene)==colnames(DEG_lasso)[i]),])
  }
}
colnames(DEG_df_test) <- as.character(colnames(DEG_lasso)) 

DEG_df_test$phenotype <-  test_set_Meta$label
DEG_df_test_1 <- DEG_df_test


test_set_Gene <- read.csv("data57065.csv", header = T)
test_set_Meta <- read.csv(file = "label_GSE57065.csv", header = T)

test_set_Gene <- test_set_Gene[,-1]

# test_set_Gene <- test_set_Gene[-which(test_set_Gene[,1]=="1-Mar"),]
# test_set_Gene <- test_set_Gene[-which(test_set_Gene[,1]=="2-Mar"),]


rownames(test_set_Gene) <- test_set_Gene[,1]
test_set_Gene <- test_set_Gene[,-1]

DEG_df_test <- as.data.frame(matrix(nrow = ncol(test_set_Gene), ncol = ncol(DEG_lasso)))
colnames(DEG_df_test) <- colnames(DEG_lasso)
rownames(DEG_df_test) <- colnames(test_set_Gene)

for(i in 1:ncol(DEG_df_test)){
  if(colnames(DEG_df_test)[i] %in% rownames(test_set_Gene)){
    DEG_df_test[,i] <- t(test_set_Gene[which(rownames(test_set_Gene)==colnames(DEG_lasso)[i]),])
  }
}
colnames(DEG_df_test) <- as.character(colnames(DEG_lasso)) 

DEG_df_test$phenotype <-  test_set_Meta$label
DEG_df_test_2 <- DEG_df_test

DEG_df_test_combined <- rbind(DEG_df_test_1, DEG_df_test_2)
colnames(DEG_df_test_combined) <- colnames(DEG_lasso)


## ---------------------------------------------------------------------------------------------
gbm_mod = gbm(phenotype ~.,
              data = DEG_lasso,
              distribution = "multinomial",
              cv.folds = 10,
              shrinkage = .01,
              n.trees = 1000,
              verbose = T)
gbm.perf(gbm_mod)


pred = predict.gbm(object = gbm_mod,
                   newdata = as.data.frame(DEG_df_test_combined),
                   n.trees = 1000,
                   type = "response")
  
labels = colnames(pred)[apply(pred, 1, which.max)]
result = data.frame(DEG_df_test_combined$phenotype, b = labels)
  
cm_DEG = confusionMatrix(data = as.factor(as.numeric(labels)), reference= as.factor(DEG_df_test_combined$phenotype))
cm_DEG


DEP_DEG_result[1,3] <- cm_DEG$overall[1]
DEP_DEG_result[2,3] <- roc(result[,1]==0, pred[,,1][,1])$auc[1]
DEP_DEG_result[3,3] <- roc(result[,1]==1, pred[,,1][,2])$auc[1]
DEP_DEG_result[4,3] <- roc(result[,1]==2, pred[,,1][,3])$auc[1]
DEP_DEG_result[5,3] <- ncol(DEG_lasso)-1


save(DEG_lasso, DEG_df_test_combined, gbm_mod, pred, cm_DEG, result, file = paste0("DEP_edgeR_",lambda_lasso,".Rdata"))


## ---------------------------------------------------------------------------------------------
DEG_lambda <- lambda_lasso

set.seed(123)
# custom_lambda <- seq((1/11), 1, by = (1/11))

DEG_df_list = list()
DESeq2G_list = list()
selecG_DESeq2_list = list()
top_gene_table_list = list()
# extract top genes and do lasso
colData <- data.frame(condition = pheno_sample)
DGE_obj <- DGEList(counts = round(2^AllGeneExp - 1), group = pheno_sample)
logCPM <- cpm(DGE_obj, log=TRUE, prior.count=3)
dim(logCPM)
for (i in 1:3) {
  pheno_sample_tcga_copy = pheno_sample
  pheno_sample_tcga_copy = ifelse(pheno_sample_tcga_copy==name_pheno[i],
                                  name_pheno[i],"other")
  colData <- data.frame(condition = pheno_sample)
  dds <- DESeqDataSetFromMatrix(countData = round(2^AllGeneExp - 1,0),
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
  DEG_df_eachPheno$phenotype = pheno_sample
  DEG_df_list = append(DEG_df_list, list(DEG_df_eachPheno))
  DESeq2Gene = colnames(DEG_df_eachPheno)[-length(colnames(DEG_df_eachPheno))]
  DESeq2G_list = append(DESeq2G_list, list(DESeq2Gene))
  
  # lasso
  DEGxy = lassoPreprocessBinary(DEG_df_eachPheno, name_pheno[i])
  # cvLassoDEG = cv.glmnet(DEGxy[[1]], DEGxy[[2]],
  #                        alpha = 1, nfolds = 10, lambda = custom_lambda)
  # plot(cvLassoDEG)
  # DEG_lambda = cvLassoDEG$lambda.min
  # print(DEG_lambda)
  finalLassoDEG = glmnet(DEGxy[[1]], DEGxy[[2]],
                         alpha = lasso_alpha, lambda=DEG_lambda) 
  result = as.matrix(coef(finalLassoDEG))
  selecG = names(result[result[,1]!=0,])[-1]
  selecG_DESeq2_list = append(selecG_DESeq2_list, list(selecG))
  print(length(selecG))
}
DESeq2G = unique(unlist(DESeq2G_list))
selecG = unique(unlist(selecG_DESeq2_list))
DEG_lasso = as.data.frame(t(logCPM)[,selecG])
DEG_lasso$phenotype = pheno_sample


## ---------------------------------------------------------------------------------------------
test_set_Gene <- read.csv("data21802.csv", header = T)
test_set_Meta <- read.csv(file = "label_GSE21802.csv", header = T)

test_set_Gene <- test_set_Gene[,-1]


rownames(test_set_Gene) <- test_set_Gene[,1]
test_set_Gene <- test_set_Gene[,-1]

test_set_Gene <- log2(test_set_Gene - min(test_set_Gene) + 1)


match(test_set_Meta[,1], colnames(test_set_Gene))

DEG_df_test <- as.data.frame(matrix(nrow = ncol(test_set_Gene), ncol = ncol(DEG_lasso)))
colnames(DEG_df_test) <- colnames(DEG_lasso)
rownames(DEG_df_test) <- colnames(test_set_Gene)

for(i in 1:ncol(DEG_df_test)){
  if(colnames(DEG_df_test)[i] %in% rownames(test_set_Gene)){
    DEG_df_test[,i] <- t(test_set_Gene[which(rownames(test_set_Gene)==colnames(DEG_lasso)[i]),])
  }
}
colnames(DEG_df_test) <- as.character(colnames(DEG_lasso)) 

DEG_df_test$phenotype <-  test_set_Meta$label
DEG_df_test_1 <- DEG_df_test


test_set_Gene <- read.csv("data57065.csv", header = T)
test_set_Meta <- read.csv(file = "label_GSE57065.csv", header = T)

test_set_Gene <- test_set_Gene[,-1]

# test_set_Gene <- test_set_Gene[-which(test_set_Gene[,1]=="1-Mar"),]
# test_set_Gene <- test_set_Gene[-which(test_set_Gene[,1]=="2-Mar"),]


rownames(test_set_Gene) <- test_set_Gene[,1]
test_set_Gene <- test_set_Gene[,-1]

DEG_df_test <- as.data.frame(matrix(nrow = ncol(test_set_Gene), ncol = ncol(DEG_lasso)))
colnames(DEG_df_test) <- colnames(DEG_lasso)
rownames(DEG_df_test) <- colnames(test_set_Gene)

for(i in 1:ncol(DEG_df_test)){
  if(colnames(DEG_df_test)[i] %in% rownames(test_set_Gene)){
    DEG_df_test[,i] <- t(test_set_Gene[which(rownames(test_set_Gene)==colnames(DEG_lasso)[i]),])
  }
}
colnames(DEG_df_test) <- as.character(colnames(DEG_lasso)) 

DEG_df_test$phenotype <-  test_set_Meta$label
DEG_df_test_2 <- DEG_df_test

DEG_df_test_combined <- rbind(DEG_df_test_1, DEG_df_test_2)


## ---------------------------------------------------------------------------------------------
gbm_mod = gbm(phenotype ~.,
              data = DEG_lasso,
              distribution = "multinomial",
              cv.folds = 10,
              shrinkage = .01,
              n.trees = 1000,
              verbose = T)
gbm.perf(gbm_mod)


pred = predict.gbm(object = gbm_mod,
                   newdata = as.data.frame(DEG_df_test_combined),
                   n.trees = 1000,
                   type = "response")
  
labels = colnames(pred)[apply(pred, 1, which.max)]
result = data.frame(DEG_df_test_combined$phenotype, b = labels)
  
cm_DEG = confusionMatrix(data=as.factor(as.numeric(labels)), reference=as.factor(DEG_df_test_combined$phenotype))
cm_DEG

DEP_DEG_result[1,4] <- cm_DEG$overall[1]
DEP_DEG_result[2,4] <- roc(result[,1]==0, pred[,,1][,1])$auc[1]
DEP_DEG_result[3,4] <- roc(result[,1]==1, pred[,,1][,2])$auc[1]
DEP_DEG_result[4,4] <- roc(result[,1]==2, pred[,,1][,3])$auc[1]
DEP_DEG_result[5,4] <- ncol(DEG_lasso)-1

save(DEG_lasso, DEG_df_test_combined, gbm_mod, pred, cm_DEG, result, file = paste0("DEP_DESeq2_",lambda_lasso,".Rdata"))


## ---------------------------------------------------------------------------------------------
save.image(paste0("DEP_DEG_result_lambda_", lambda_lasso, "_alpha_", lasso_alpha, ".Rdata"))
write.csv(DEP_DEG_result, file = paste0("result_lambda_", lambda_lasso,  "_alpha_", lasso_alpha, ".csv"))

