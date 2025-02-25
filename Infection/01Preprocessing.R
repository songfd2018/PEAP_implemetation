rm(list=ls())
library(ggplot2)
library(clusterProfiler) # Download KEGG pathway
library(biomaRt) # gene annotation conversion
library(dplyr)
library(edgeR)
library(qusage) # for read.gmt
library(plyr)
ver <- "Infection"

setwd("D:/PEAP/Infection")

Infection_Gene_exp_data <- read.csv("Gene_df_infection_big.csv", header = T)
Infection_meta_data <- read.csv("meta_data_infection_big.csv", header = T, row.names = 1)

Gene_matrix <- Infection_Gene_exp_data
colnames(Gene_matrix)[1] <- "Genes"
Gene_matrix <- aggregate(. ~ Genes, data = Gene_matrix, mean)

Infection_Gene_exp_data <- Gene_matrix


# Infection_Gene_exp_data <- Infection_Gene_exp_data[-1,]
rownames(Infection_Gene_exp_data) <- Infection_Gene_exp_data[,1]
Infection_Gene_exp_data <- Infection_Gene_exp_data[,-1]
Infection_meta_data$label[which(Infection_meta_data$label == 12)] = 1 

Infection_meta_data <- Infection_meta_data[match(colnames(Infection_Gene_exp_data), rownames(Infection_meta_data)),]


pheno_sample <- Infection_meta_data$label[which(colnames(Infection_Gene_exp_data)==rownames(Infection_meta_data))]

setwd("D:/PEAP/Infection/")

KEGG_list <- read.gmt("c2.cp.kegg.v2022.1.Hs.symbols.gmt")
gene_symbol <- rownames(Infection_Gene_exp_data)
num_KEGG_pathway <- length(KEGG_list)

nrow(Infection_Gene_exp_data)
length(unique(unlist(KEGG_list)))


KEGG_genes <- list()
KEGG_num <- NULL
for(i in 1:num_KEGG_pathway){
  KEGG_genes[[i]] <- which(gene_symbol %in% KEGG_list[[i]])
  KEGG_num <- c(KEGG_num, length(KEGG_genes[[i]]))
}

length(unique(unlist(KEGG_genes)))
hist(KEGG_num)

save.image("Preprocessing_wksp_Infection_20240404.Rdata")

setwd("D:/PEAP/Infection/data_for_02_Infection_final_ver/")

for(i in 1:num_KEGG_pathway){
  cur_pathway <- names(KEGG_list)[i]
  common_genes <- KEGG_genes[[i]]
  
  out_df <- data.frame(GeneIndex = common_genes, 
                       GeneSymbol = gene_symbol[common_genes])
  
  print(paste0("In the ",i,"-th pathway---",cur_pathway,"---there are ", length(common_genes), " genes in gene expression profile."))
  
  if(length(common_genes) > 10){
    GeneExp_pathway <- Infection_Gene_exp_data[common_genes, ]
    
    # exclude non-expressed genes
    expressed_genes <- rowSums(GeneExp_pathway) > 0
    GeneExp <- GeneExp_pathway[expressed_genes, ]
    
    save(GeneExp, pheno_sample, file = paste0("Gene_expression_KEGG",i,"_v",ver,".RData"))
  }
}

