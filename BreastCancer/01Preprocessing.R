rm(list=ls())
library(ggplot2)
#library(clusterProfiler) # Download KEGG pathway
library(biomaRt) # gene annotation conversion
library(dplyr)
library(edgeR)
library(qusage) # for read.gmt
ver <- 8
setwd("D:/PEAP/PEAP_02/")
# load("Preprocessing_wksp.Rdata")
load("shared_genes.Rdata") # Common gene names of training set and testing set. A char vector 

## Load RPKM counts and identify the subtypes of samples

# We download the gene expression profiles and metadata of GDC TCGA breast cancer database from UCSC Xena  <https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Breast%20Cancer%20(BRCA)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443>. 
# There are 60,483 features and 1,217 samples in the gene expression matrix, while 1,247 samples are included in the metadata. We need to align samples in the gene expression matrix with those in the metadata and obtain the corresponding phenotype.



GeneExp_BRCA <- read.table(file = 'TCGA-BRCA.htseq_fpkm.tsv.gz', sep = '\t', header = TRUE)
rownames(GeneExp_BRCA) <- GeneExp_BRCA$Ensembl_ID
GeneExp_BRCA <- GeneExp_BRCA[,-1]
sampleID_dat <- colnames(GeneExp_BRCA)

metadata_BRCA <- read.table(file = '../data/TCGA.BRCA.sampleMap_BRCA_clinicalMatrix', sep = '\t',
                           header = TRUE)


## ----mapping----------------------------------------------------------------------------------
gene_annotation <- read.table("gencode.v22.annotation.gene.probeMap", header = TRUE, sep = "\t")

gene_exp <- rowSums(GeneExp_BRCA)
order_genes <- order(gene_exp, decreasing = TRUE)
GeneExp_BRCA <- GeneExp_BRCA[order_genes,]

gene_matching <- match(rownames(GeneExp_BRCA), gene_annotation$id)
gene_annotation <- gene_annotation[gene_matching, ]

duplicated_genes <- duplicated(gene_annotation$gene)
GeneExp_BRCA <- GeneExp_BRCA[!duplicated_genes,]
gene_annotation <- gene_annotation[!duplicated_genes,]

rownames(GeneExp_BRCA) <- gene_annotation$gene
share_id <- which(rownames(GeneExp_BRCA) %in% GSE5327_genes)
GeneExp_BRCA <- GeneExp_BRCA[share_id,]


sampleID_dat <- gsub("\\.", "-", sampleID_dat)
patientID_dat <- substr(sampleID_dat, 1,15)
# split_sampleID <- matrix(unlist(strsplit(sampleID_dat,split='-')), ncol = 4, byrow = TRUE)

matched_index_meta <- match(patientID_dat, metadata_BRCA$sampleID)

## check whether there are some repeated sample IDs, the samples with the sample ID in the metadata are distinguished by -01A or -01B in the gene expression matrix
length(matched_index_meta)
length(unique(matched_index_meta))

phenotype <- metadata_BRCA$PAM50Call_RNAseq
pheno_sample <- phenotype[matched_index_meta]

unknown <- which(pheno_sample == "")
pheno_sample <- pheno_sample[-unknown]
GeneExp_BRCA <- GeneExp_BRCA[,-unknown]


# load("Preprocessing_wksp.Rdata")


## ---- HVG-------------------------------------------------------------------------------------
Total_expression_samples <- colSums(GeneExp_BRCA)
summary(Total_expression_samples)

N <- ncol(GeneExp_BRCA)
G <- nrow(GeneExp_BRCA)
mean_expression <- rowSums(GeneExp_BRCA)/N
sd_expression <- apply(GeneExp_BRCA, 1, sd)

mead_sd_curve <- data.frame(GeneName = rownames(GeneExp_BRCA),
                            MeanLevel = mean_expression,
                            SDLevel = sd_expression)

ggplot(mead_sd_curve, aes(x=MeanLevel, y=SDLevel)) +
  geom_point()


selected_genes <- sd_expression > 2
mead_sd_curve$Selected <- factor(selected_genes,levels = c("TRUE","FALSE"))

ggplot(mead_sd_curve, aes(x=MeanLevel, y=SDLevel)) +
  geom_point(aes(color = Selected))

GeneExp <-  GeneExp_BRCA[selected_genes, ]
save(GeneExp, pheno_sample, file = paste0("Gene_expression_HVGs_v",ver,".RData"))


## ----KEGG-------------------------------------------------------------------------------------
KEGG_list <- read.gmt("c2.cp.kegg.v2022.1.Hs.symbols.gmt")
gene_symbol <- rownames(GeneExp_BRCA)
num_KEGG_pathway <- length(KEGG_list)

nrow(GeneExp_BRCA)
length(unique(unlist(KEGG_list)))

save.image("BRCA_overlap_preprocessing")


## ----common_genes-----------------------------------------------------------------------------
KEGG_genes <- list()
KEGG_num <- NULL
for(i in 1:num_KEGG_pathway){
  KEGG_genes[[i]] <- which(gene_symbol %in% KEGG_list[[i]])
  KEGG_num <- c(KEGG_num, length(KEGG_genes[[i]]))
}

length(unique(unlist(KEGG_genes)))


## ----pathway----------------------------------------------------------------------------------

for(i in 1:num_KEGG_pathway){
  cur_pathway <- names(KEGG_list)[i]
  common_genes <- KEGG_genes[[i]]

  out_df <- data.frame(GeneIndex = common_genes, 
                       GeneSymbol = gene_symbol[common_genes])
  
  print(paste0("In the ",i,"-th pathway---",cur_pathway,"---there are ", length(common_genes), " genes in gene expression profile."))
  
  if(length(common_genes) > 10){
    GeneExp_pathway <- GeneExp_BRCA[common_genes, ]
    
    # exclude non-expressed genes
    expressed_genes <- rowSums(GeneExp_pathway) > 0
    GeneExp <- GeneExp_pathway[expressed_genes, ]
    
    save(GeneExp, pheno_sample, file = paste0("gene_expression/Gene_expression_KEGG",i,"_v",ver,".RData"))
  }
}



