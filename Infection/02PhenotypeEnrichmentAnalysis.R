######################################################
# Calculate the Enrichment score and nominal p-value #
######################################################
library(getopt)
library(fgsea)

# run by Rscript 02PEAP_v7.R -d Gene_expression_hsa00190_v7
spec = matrix(c('dataset', 'd', 1, "character",
                'seed', 's', 2, "integer"), byrow=TRUE, ncol=4)
opt = getopt(spec)

filename <- opt$dataset
if (is.null(opt$seed)) { 
  seed <- 42
}else{
  seed <- opt$seed
}

split_name <- unlist(strsplit(filename,split='_'))
set_name <- split_name[3]
ver_infor <- split_name[4]

# setwd("D:/CUHKSZ/Collaboration/GSEA/code/v6")
# source("wks.r")
# setwd("/home/songfangda/PEAP/01BRCA/v7/")
setwd("D:/PEAP/Infection/data_for_02_3/")


load(paste0("gene_expression/",filename,".RData"))
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

set.seed(seed)
n_gene <- nrow(GeneExp)
N_pairs <- n_gene * (n_gene - 1)/2
gene1 <- rep(NA, N_pairs)
gene2 <- rep(NA, N_pairs)
ES <- matrix(NA, N_pairs, N_pheno)
pvalue <- matrix(NA, N_pairs, N_pheno)
colnames(ES) <- pheno_names
colnames(pvalue) <- pheno_names

start_loop <- Sys.time()
index <- 1
for(g1 in 1:(n_gene-1)){
  for(g2 in (g1+1):n_gene){
    
    gene1[index] <- gene_set[g1]
    gene2[index] <- gene_set[g2]
    
    fc_genepair <- unlist(GeneExp[g1,] - GeneExp[g2,])
    names(fc_genepair) <- sample_names

    
    # head(fgseaRes[order(pval), ])
    if(sum(is.na(fc_genepair)) == 0){
      if(diff(range(fc_genepair)) != 0){
        fgseaRes <- fgsea(pathways = cell_name_list, 
                          stats = fc_genepair,
                          eps = 0)
        
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
  }
}

end_loop <- Sys.time()
cat("It takes ",difftime(end_loop,start_loop,units = "mins"),
    " mins to calculate ES and nominal p-value for all gene pairs.\n")




PEAP_res <- data.frame(gene1 = gene1, gene2 = gene2, 
                       ES = ES, pvalue = pvalue)



n_pairs <- nrow(PEAP_res)
phenotype_names <- c("0", "1", "2")
pathway_names <- unique(PEAP_res$pathway)
n_pheno <- 3

list_phenotype <- list()
for(j in 1:n_pheno){
  list_phenotype[[j]] <- which(pheno_sample == phenotype_names[j])
}

fc_gp <- matrix(NA, n_pairs, n_pheno)
for(i in 1:n_pairs){
  gene_pairs <- unlist(PEAP_res[i,1:2])
  cur_logratio <- unlist(GeneExp[gene_pairs[1],] - GeneExp[gene_pairs[2],])
  for(j in 1:n_pheno){
    fc_gp[i,j] <- mean(cur_logratio[list_phenotype[[j]]]) - mean(cur_logratio[-list_phenotype[[j]]])
  }
}

colnames(fc_gp) <- c("fc.healthy", "fc.bacterial", "fc.viral")

PEAP_res <- cbind(PEAP_res, fc_gp)






# plotEnrichment(cell_name_list[["Basal"]],fc_genepair) + labs(title="Basal subtype")
if(!dir.exists("output")){
  dir.create("output")
}
save.image(file = paste0("output/PEAP_",set_name,"_",ver_infor,".RData"))
