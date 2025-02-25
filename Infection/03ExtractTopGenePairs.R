##############################################################
# Parse the enrichment scores and p-value for all gene pairs #
##############################################################
rm(list=ls())
ver <- "Infection"
library(ggplot2)
library(xtable)
library(fgsea)
library(ggpubr)
library(qusage) # for read.gmt

setwd("D:/PEAP/Infection/")


# Load pathway
KEGG_list <- read.gmt("c2.cp.kegg.v2022.1.Hs.symbols.gmt")
KEGG_names <- names(KEGG_list)
num_KEGG_pathway <- length(KEGG_list)
ESres_collection <- NULL
AllGeneExp <- data.frame()

for(p in 1:num_KEGG_pathway){
  res_file <- paste0("output_infection_0406/PEAP_KEGG",p,"_v",ver,".RData")
  if(file.exists(res_file)){
    load(res_file)
    PEAP_res$pathway <- KEGG_names[p]
    ESres_collection <- rbind(ESres_collection, PEAP_res)
    GeneExp$Genename <- rownames(GeneExp)
    
    AllGeneExp <- rbind(AllGeneExp, GeneExp)
  }
}

AllGeneExp <- AllGeneExp[!duplicated(AllGeneExp$Genename), ]
AllGeneExp$Genename <- NULL

# adjust for hypothesis testing by FDR across all gene pairs and all phenotype
adjusted_pval <- matrix(p.adjust(unlist(ESres_collection[, 2 + N_pheno + 1:N_pheno]), method = "BH"),ncol = N_pheno)

phenotype_names <- names(table(pheno_sample)) 
colnames(adjusted_pval) <- phenotype_names



# load count data matrix
load("Preprocessing_wksp_Infection_20240404.Rdata")


gene_set <- rownames(AllGeneExp)
top100_Genes <- data.frame()



for(ph in 1:N_pheno){
  
  ph_name <- phenotype_names[ph]
  adjp_infor <- cbind(ESres_collection[,c(1:2,2+ph,5+ph,3 + N_pheno * 2, 12)])
  adjp_Ph <- adjp_infor[order(adjusted_pval[,ph]),]
  adjp_Ph[,4] <- log10(adjp_Ph[,4])
  colnames(adjp_Ph) <- c("Gene 1", "Gene 2", "ES", "log10(pval)")
  rownames(adjp_Ph) <- NULL
  xtable(adjp_Ph[1:20,])

  # pdf(file = paste0("Top5_Enrich_",ph_name,".pdf"), width = 20, height = 20)
  # Draw Enrichment plot
  # ESplot <- NULL
  # for(g in 1:5){
  # 
  #   g1 <- which(gene_set == adjp_Ph[g,1])
  #   g2 <- which(gene_set == adjp_Ph[g,2])
  # 
  #   fc_genepair <- unlist(AllGeneExp[g1,] - AllGeneExp[g2,])
  # 
  #   names(fc_genepair) <- sample_names
  #   ESplot[[g]] <- plotEnrichment(cell_name_list[[ph_name]],fc_genepair) # + labs(title="Basal subtype")
  # }
  # 
  # ggarrange(ESplot[[1]], ESplot[[2]], ESplot[[3]], ESplot[[4]], ESplot[[5]],
  #           labels = paste0(adjp_Ph[1:5, 1],"-",adjp_Ph[1:5, 2]),
  #           ncol = 1, nrow = 5)
  # dev.off()

  # Draw Heatmap for top 100 gene pairs
  # Delete same gene pairs
  one_pheno_top100 <- data.frame()
  used_genes <- c()
  t <- 1
  for(i in 1:nrow(adjp_Ph)){
    if(adjp_Ph[i,1] %in% used_genes == FALSE & adjp_Ph[i,2] %in% used_genes == FALSE){
      used_genes <- append(used_genes, adjp_Ph[i,1])
      used_genes <- append(used_genes, adjp_Ph[i,2])
      one_pheno_top100 <- rbind(one_pheno_top100, adjp_Ph[i,])
      t <- t+1
    }
    if(t==1001){
      break;
    }
  }
  
  #Store top 100 gene pairs of each phenotype
  top100_Genes <- rbind(top100_Genes, one_pheno_top100)
}
save.image("Nonoverlap_top100GenesOfEachGenepairs_Infection_20240414.Rdata")

save(pheno_sample, AllGeneExp, top100_Genes, file = "Infection_PEAP03_output_012.Rdata")


load("Preprocessing_wksp_Infection_20240404.Rdata")
load("Nonoverlap_top100GenesOfEachGenepairs_Infection_20240406.Rdata")
n_pheno = N_pheno
##################################################################################
# Fisher's Exact test to find out the association between phenotype and pathways #
##################################################################################


# Combine fold change into the data frame
ESres_noduplicated <- ESres_collection[!duplicated(ESres_collection[,1:2]),] 

# Adjust for multiple hypothesis testing
adjusted_pval_BH <- matrix(p.adjust(unlist(ESres_noduplicated[, 2 + n_pheno + 1:n_pheno]), 
                                    method = "BH"),ncol = n_pheno)
adjusted_pval_BH[is.na(adjusted_pval_BH)] <- 1
fc_noduplicated <- as.matrix(ESres_noduplicated[, 9:11])

# Add DE indicators to the data frame of ES results
thres <- -70
thres_fc <- 1
DE_indicators <- adjusted_pval_BH < 10^thres & abs(fc_noduplicated) > thres_fc
colnames(DE_indicators) <- paste0("DE.", phenotype_names)

ESres_noduplicated <- data.frame(ESres_noduplicated, DE_indicators)

# match the names of gene pairs back to the original data frame with replicated gene pairs
matching_genepairs <- match(paste0(ESres_collection[,1], "_", ESres_collection[,2]), 
                            paste0(ESres_noduplicated[,1], "_", ESres_noduplicated[,2]))

ESres_collection <- data.frame(ESres_collection, ESres_noduplicated[matching_genepairs, 13:15])

#####################
# Fisher exact test #
#####################
pathway_names <- unique(ESres_collection$pathway)
n_pathway <- length(pathway_names)
Fisher_odds <- matrix(NA, n_pathway, n_pheno)
Fisher_pval <- matrix(NA, n_pathway, n_pheno)

rownames(Fisher_odds) <- pathway_names
colnames(Fisher_odds) <- phenotype_names

rownames(Fisher_pval) <- pathway_names
colnames(Fisher_pval) <- phenotype_names

for(i in 1:n_pathway){
  for(j in 1:n_pheno){
    phenotype_cur <- pathway_names[i]
    fisher_res <- fisher.test(ESres_collection$pathway == phenotype_cur,
                              ESres_collection[, j + 12], alternative = "greater")
    
    Fisher_odds[i,j] <- fisher_res$estimate
    Fisher_pval[i,j] <- fisher_res$p.value
  }
}

Fisher_adjpval <- matrix(p.adjust(Fisher_pval, method = "BH"), nrow = n_pathway) 
colnames(Fisher_adjpval) <- colnames(Fisher_pval)
rownames(Fisher_adjpval) <- rownames(Fisher_pval)

Fisher_log_adjpval <- ifelse(log10(Fisher_adjpval) < -50, 50, -log10(Fisher_adjpval))

#######################################
# Heatmap of all significant pathways #
#######################################
library(pheatmap)
thres_fisher <- 10
apply(Fisher_log_adjpval > thres_fisher, 2, sum)

apply(Fisher_log_adjpval > 5, 2, sum)

FDR_10_pathway <- NULL
for(j in 1:n_pheno){
  order_path <- order(Fisher_log_adjpval[,j], decreasing = TRUE) 
  FDR_10_pathway <- c(FDR_10_pathway, order_path[Fisher_log_adjpval[order_path,j] > thres_fisher])
}

FDR_10_pathway <- FDR_10_pathway[!duplicated(FDR_10_pathway)]

# Col_annotation
anno <- data.frame(Phenotype = phenotype_names)
rownames(anno) <- phenotype_names

pheno_color <- c("#34A853", "#FBBC05", "#EA4335")
names(pheno_color) <- phenotype_names


annoCol <- list(Phenotype = pheno_color)

logpvalue_FDR10 <- Fisher_log_adjpval[FDR_10_pathway,]

library(RColorBrewer)
PEAP_FDR <- 70
FC <- 1
Fisher_FDR <- 10
save(ESres_collection, file="Infection_ESres.Rdata")

pdf(file = paste0("Heatmap_association_PEAPFDR",PEAP_FDR,
                  "_FC",FC,"_FisherFDR",Fisher_FDR, ".pdf"), 
    width = 10, height = 4)


pheatmap(logpvalue_FDR10,
         breaks = seq(0,20,0.2),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_col = anno,
         annotation_colors = annoCol,
         color = colorRampPalette(brewer.pal(n = 7, name = "OrRd"))(100),
         fontsize_row = 12, 
         cellwidth = 30, cellheight = 20)
dev.off()

save.image("0406pval_FisherExactTest.RData")
















