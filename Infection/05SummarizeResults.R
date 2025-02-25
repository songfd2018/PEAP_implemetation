rm(list = ls())
library(getopt)
library(fgsea)
setwd("D:/PEAP/Infection/")
load("marker_genepairs.Rdata")
load("Infection_PEAP03_output_012.Rdata")
load("Preprocessing_wksp_Infection_20240404.Rdata")

marker_gene_pairs = unlist(selecP_list)
gene_pair_ind = paste0(top100_Genes[,1], "_", top100_Genes[,2])

marker_gene_pairs_ind = match(marker_gene_pairs, gene_pair_ind)

setwd("D:/PEAP/Infection/data_for_02_Infection_final_ver/")

index <- 1

N_pheno = 3
pheno_names=c("healthy", "bacterial", "viral")
gene1 <- rep(NA, length(marker_gene_pairs))
gene2 <- rep(NA, length(marker_gene_pairs))
ES <- matrix(NA, length(marker_gene_pairs), N_pheno)
pvalue <- matrix(NA, length(marker_gene_pairs), N_pheno)
colnames(ES) <- pheno_names
colnames(pvalue) <- pheno_names


start_loop <- Sys.time()

for(j in marker_gene_pairs_ind){
  filename = which(names(KEGG_list)==top100_Genes[j,6])
  load(paste0("Gene_expression_KEGG",filename,"_vInfection.RData"))
  gene_set <- rownames(GeneExp)
  sample_names <- colnames(GeneExp)
  
  # define a list of cell labels
  pheno_names <- names(table(pheno_sample))
  
  cell_name_list <- list()
  N_pheno <- length(pheno_names)
  for(i in 1:N_pheno){
    cell_name_list[[i]] <- sample_names[which(pheno_sample == pheno_names[i])]
  }
  names(cell_name_list) <- pheno_names
  
  set.seed(42)



  g1 = which(rownames(GeneExp)==top100_Genes[j,1])
  g2 = which(rownames(GeneExp)==top100_Genes[j,2])
  
  if(g1 %% 100 == 0){print(g1)}
  gene1[index] <- gene_set[g1]
  gene2[index] <- gene_set[g2]
  
  fc_genepair <- unlist(GeneExp[g1,] - GeneExp[g2,])
  names(fc_genepair) <- sample_names
  
  
  # head(fgseaRes[order(pval), ])
  if(sum(is.na(fc_genepair)) == 0){
    if(diff(range(fc_genepair)) != 0){
      fgseaRes <- fgseaMultilevel(pathways = cell_name_list, 
                        stats = fc_genepair,
                        eps = 0, scoreType = "std")
      
      ES[index, ] <- fgseaRes$ES
      pvalue[index, ] <- fgseaRes$pval
      
    }else{
      
      cat("The fold changes of gene pair (",g1,",",g2,") are all zeros.\n", sep = "")
    }
    
  }else{
    
    cat("There exist missing values in the fold change of gene pair (",g1,",",g2,").\n",sep = "")
    
  }
  # cat("Finish running on gene pair (",g1,",",g2,")...\n",sep = "")
  index <- index + 1
  
  if(index %% 100 == 0){
    mid_loop <- Sys.time()
    cat("It takes ",difftime(mid_loop,start_loop,units = "mins"),
        " mins to finish the calculation of ",index," gene pairs.\n")
  }
  cat(index)
}


end_loop <- Sys.time()
cat("It takes ",difftime(end_loop,start_loop,units = "mins"),
    " mins to calculate ES and nominal p-value for all gene pairs.\n")


PEAP_res <- data.frame(gene1 = gene1, gene2 = gene2, 
                       ES = ES, pvalue = pvalue)

save.image("PEAP05_output_std.Rdata")



load("PEAP05_output_std.Rdata")

PEAP_res_std = PEAP_res

load("PEAP05_output_pos.Rdata")

PEAP_res_pos = PEAP_res

load("PEAP05_output_neg.Rdata")

PEAP_res_neg = PEAP_res

colnames(PEAP_res_std) = paste0(colnames(PEAP_res_std), ".std")
colnames(PEAP_res_pos) = paste0(colnames(PEAP_res_pos), ".pos")
colnames(PEAP_res_neg) = paste0(colnames(PEAP_res_neg), ".neg")


PEAP_all = cbind(PEAP_res_std, PEAP_res_pos[,3:8], PEAP_res_neg[,3:8])

for(i in 1:nrow(PEAP_all)){
  for(j in 1:3){
    if(is.na(PEAP_all[i,j+5])){
      if(PEAP_all[i,j+2]>0){
        PEAP_all[i,j+5] = PEAP_all[i,j+11]
      }
      else(PEAP_all[i,j+5] = PEAP_all[i,j+17])
    }
  }

}

PEAP_all = PEAP_all[,1:8]
colnames(PEAP_all) = c("gene1", "gene2", "ES.healthy","ES.bacterial","ES.viral", "pvalue.healthy", "pvalue.bacterial", "pvalue.viral")
write.csv(PEAP_all, file = "PEAP_res_infection_20240422.csv",row.names = FALSE)









