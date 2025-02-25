### ROC plot

rm(list = ls())
library(pROC)
library(RColorBrewer)
library(ggplot2)
setwd("D:/PEAP/Infection/")

# The input here comes from PEAP04
load("DEP_lambda_0.1.Rdata")
pred_peap <- pred

load("DEP_limma_0.1.Rdata")
pred_limma <- pred

load("DEP_edgeR_0.1.Rdata")
pred_edgeR <- pred

load("DEP_DESeq2_0.1.Rdata")
pred_DESeq <- pred

# Sample data (replace this with your actual data)
set.seed(123)

method_names <- c("Healthy", "Bacterial", "Viral")
actual_class <- DEP_df_test_combined$phenotype

auc_all <- list()

pdf(paste0("ROCp.pdf"),width = 15, height = 5)
par(mfrow=c(1,3))
for(i in 0:2){
  predicted_scores_peap <- pred_peap[,,1][,i+1]
  predicted_scores_limma <- pred_limma[,,1][,i+1]
  predicted_scores_edgeR <- pred_edgeR[,,1][,i+1]
  predicted_scores_DESeq <- pred_DESeq[,,1][,i+1]
  
  roc_peap <- roc(as.numeric(actual_class == as.character(i)), predicted_scores_peap, smooth=F)
  roc_limma <- roc(as.numeric(actual_class == as.character(i)), predicted_scores_limma, smooth=T)
  roc_edgeR <- roc(as.numeric(actual_class == as.character(i)), predicted_scores_edgeR, smooth=T)
  roc_DESeq <- roc(as.numeric(actual_class == as.character(i)), predicted_scores_DESeq, smooth=T)
  
  auc_all[[i+1]] <- c(roc_limma$auc, roc_edgeR$auc, roc_DESeq$auc, roc_peap$auc)
  
  
  # pdf(paste0("ROC_",method_names[i+1],".pdf"),width = 5, height = 5)
  # ROC plot 1
  roc(as.numeric(actual_class == as.character(i)), predicted_scores_peap, plot=TRUE, legacy.axes=FALSE, percent=TRUE, col=rgb(248,230,32, maxColorValue = 255), lwd=4, print.auc=FALSE, cex.lab = 2, cex.axis = 2, main = method_names[i+1], cex.main=2)
  
  plot.roc(as.numeric(actual_class == as.character(i)), predicted_scores_limma, percent=TRUE, col=rgb(145,213,66, maxColorValue = 255), lwd=4, print.auc=FALSE, add=TRUE, smooth = TRUE, lty = 6)
  
  plot.roc(as.numeric(actual_class == as.character(i)), predicted_scores_edgeR, percent=TRUE, col=rgb(53,183,119, maxColorValue = 255), lwd=4, print.auc=FALSE, add=TRUE, smooth = TRUE, lty = 5)
  
  plot.roc(as.numeric(actual_class == as.character(i)), predicted_scores_DESeq, percent=TRUE, col=rgb(38,104,141, maxColorValue = 255), lwd=4, print.auc=FALSE, add=TRUE, smooth = TRUE, lty = 4)
  
  # Add legend
  # legend("bottom",
  #        legend=c("PEAP", "Limma","edgeR", "DESeq2"),
  #        lty = c(1,6,5,4),
  #        col=c(rgb(248,230,32, maxColorValue = 255), rgb(145,213,66, maxColorValue = 255), rgb(53,183,119, maxColorValue = 255), rgb(38,104,141, maxColorValue = 255)),
  #        lwd=6, xpd = TRUE, bty = "n", seg.len=5)
  # dev.off()
}
dev.off()



### Enrechment Score plot
plotEnrichment_plus <- function (pathway, stats, gseaParam = 1, ticksSize = 0.2) 
{
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
                          returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  diff <- (max(tops) - min(bottoms))/8
  x = y = NULL
  g <- ggplot(toPlot, aes(x = x, y = y)) + 
    geom_point(color = "black", size = 0.1) + 
    geom_hline(yintercept = max(tops), colour = "red", linetype = "dashed") + 
    geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed") + geom_hline(yintercept = 0, colour = "black") + geom_line(color = "black") + theme_bw() + 
    geom_segment(data = data.frame(x = pathway), mapping = aes(x = x, y = 0.8, xend = x, yend = 0.95), size = ticksSize) + 
    theme(panel.border = element_blank(), panel.grid.minor = element_blank()) + 
    ylim(c(-1,1)) +
    xlim(c(0,2500)) +
    labs(x = "Rank", y = "")
  g
}

rm(list = ls())
setwd("D:/PEAP/Infection/")
library(fgsea)
library(ggplot2)


load("Sepsis_PEAP03_output_012.Rdata")
load("ES.Rdata")

##########################################
# Heatmap to visualize enrichment scores #
##########################################
library(fgsea)
load("ES.Rdata")

ESres_top5 <- NULL
for(i in 1:3){
  order_pval <- order(ESres_collection[,i + 5])
  ESres_cur <- ESres_collection[order_pval[1:5], ]
  for(j in 1:5){
    # Ensure all enrichment scores are positive
    if(ESres_cur[j, i + 2] < 0){
      ESres_cur[j, 1:2] <- ESres_cur[j, 2:1]
      ESres_cur[j, 3:5] <- -ESres_cur[j, 3:5]
    }
  }
  ESres_top5 <- rbind(ESres_top5, ESres_cur)
}

# Select the top 5 gene pairs with the smallest p-value for each phenotype
library(pheatmap)

ES_matrix <- ESres_top5[,3:5]
ES_matrix <- rbind(ES_matrix, ESres_collection[order_pval[100001:100005], 3:5])

colnames(ES_matrix) <- c("Healthy", "Baterial", "Viral")
rownames(ES_matrix) <- paste0("Gene Pairs ",1:20)

pheatmap(ES_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE, 
         color = colorRampPalette(brewer.pal(9, "OrRd"))(100))

heatmap_sample <- ES_matrix
annotation_label_col <- as.factor(c(rep("phenotype1", 1), rep("phenotype2", 1), rep("phenotype3", 1)))
main_groups_col <- data.frame(row.names = colnames(heatmap_sample), phenotype=annotation_label_col)

annotation_label_row <- as.factor(c(rep("pathway1", 5), rep("pathway2", 5), rep("pathway3", 10)))
main_groups_row <- data.frame(row.names = rownames(heatmap_sample), pathway = annotation_label_row)

ann_colors=list(phenotype=c(phenotype1='#55AF7B',phenotype2='#FAC230',phenotype3='#EB4537'), pathway=c(pathway1=brewer.pal(3,"Set3")[1], pathway2=brewer.pal(3,"Set3")[2], pathway3=brewer.pal(3,"Set3")[3]))

pdf("heatmap_pathwayAnalysis.pdf", height = 30, width = 10)
pheatmap(heatmap_sample, cluster_rows = F, cluster_cols = F, annotation_colors = ann_colors,cellwidt = 45,cellheight = 15, color = colorRampPalette(brewer.pal(9, "OrRd"))(100), annotation_row = main_groups_row, annotation_col = main_groups_col)
dev.off()


### Enrichment Score Plot
library(gridExtra)
library(fgsea)
library(ggplot2)

load("Infection_PEAP03_output_012.Rdata")
load(paste0("DEP_lambda_",0.1,".Rdata"))

gene_pair_name <- colnames(DEP_df_test_combined)[-ncol(DEP_df_test_combined)]

sample_name <- colnames(AllGeneExp)
sample_pheno <- pheno_sample

cell_name_list <- list("0"=sample_name[which(sample_pheno==0)], "1"=sample_name[which(sample_pheno==1)],"2"=sample_name[which(sample_pheno==2)])

esdata_all <- as.data.frame(t(DEP_lasso))

gene_pair_ref <- paste0(top100_Genes[,1], "_", top100_Genes[,2])


# All Enrichment Score plot

for(i in 1:87){
  gene_pair_class <- -1
  esdata <- as.numeric(esdata_all[i,])
  esdata <- as.vector(esdata)
  names(esdata) <- sample_name
  current_gene_pair_name <- rownames(esdata_all[i,])
  
  if(current_gene_pair_name %in% gene_pair_ref[1:1000]){
    gene_pair_class = 0
  }else if(current_gene_pair_name %in% gene_pair_ref[1001:2000]){
    gene_pair_class = 1
  }else if(current_gene_pair_name %in% gene_pair_ref[2001:3000]){
    gene_pair_class = 2
  }
  
  p_index=0
  ESplot_0 <- plotEnrichment_plus(cell_name_list[[p_index+1]], esdata) #+
  # labs(title=paste0(current_gene_pair_name, "s",gene_pair_class, "a", p_index)) + theme(plot.title = element_text(size = 10))
  
  p_index=1
  ESplot_1 <- plotEnrichment_plus(cell_name_list[[p_index+1]], esdata)# +
  # labs(title=paste0(current_gene_pair_name, "s",gene_pair_class, "a", p_index)) + theme(plot.title = element_text(size = 10))
  
  p_index=2
  ESplot_2 <- plotEnrichment_plus(cell_name_list[[p_index+1]], esdata)# +
  # labs(title=paste0(current_gene_pair_name, "s",gene_pair_class, "a", p_index)) + theme(plot.title = element_text(size = 10))
  
  
  
  ESplot_all <- grid.arrange(ESplot_0, ESplot_1, ESplot_2, nrow=1, ncol=3)
  ggsave(ESplot_all, file = paste0("ESplot_genepair",i,".jpg"))
}

# Sample Enrichment Score plot

i=1
gene_pair_class <- -1
esdata <- as.numeric(esdata_all[i,])
esdata <- as.vector(esdata)
names(esdata) <- sample_name
current_gene_pair_name <- rownames(esdata_all[i,])

if(current_gene_pair_name %in% gene_pair_ref[1:1000]){
  gene_pair_class = 0
}else if(current_gene_pair_name %in% gene_pair_ref[1001:2000]){
  gene_pair_class = 1
}else if(current_gene_pair_name %in% gene_pair_ref[2001:3000]){
  gene_pair_class = 2
}

p_index=0
ESplot_0 <- plotEnrichment_plus(cell_name_list[[p_index+1]], esdata)
p_index=1
ESplot_1 <- plotEnrichment_plus(cell_name_list[[p_index+1]], esdata)
p_index=2
ESplot_2 <- plotEnrichment_plus(cell_name_list[[p_index+1]], esdata)


i=40
gene_pair_class <- -1
esdata <- as.numeric(esdata_all[i,])
esdata <- as.vector(esdata)
names(esdata) <- sample_name
current_gene_pair_name <- rownames(esdata_all[i,])

if(current_gene_pair_name %in% gene_pair_ref[1:1000]){
  gene_pair_class = 0
}else if(current_gene_pair_name %in% gene_pair_ref[1001:2000]){
  gene_pair_class = 1
}else if(current_gene_pair_name %in% gene_pair_ref[2001:3000]){
  gene_pair_class = 2
}

p_index=0
ESplot_3 <- plotEnrichment_plus(cell_name_list[[p_index+1]], esdata)
p_index=1
ESplot_4 <- plotEnrichment_plus(cell_name_list[[p_index+1]], esdata)
p_index=2
ESplot_5 <- plotEnrichment_plus(cell_name_list[[p_index+1]], esdata)



i=72
gene_pair_class <- -1
esdata <- as.numeric(esdata_all[i,])
esdata <- as.vector(esdata)
names(esdata) <- sample_name
current_gene_pair_name <- rownames(esdata_all[i,])

if(current_gene_pair_name %in% gene_pair_ref[1:1000]){
  gene_pair_class = 0
}else if(current_gene_pair_name %in% gene_pair_ref[1001:2000]){
  gene_pair_class = 1
}else if(current_gene_pair_name %in% gene_pair_ref[2001:3000]){
  gene_pair_class = 2
}

p_index=0
ESplot_6 <- plotEnrichment_plus(cell_name_list[[p_index+1]], esdata)
p_index=1
ESplot_7 <- plotEnrichment_plus(cell_name_list[[p_index+1]], esdata)
p_index=2
ESplot_8 <- plotEnrichment_plus(cell_name_list[[p_index+1]], esdata)



ESplot_all <- grid.arrange(ESplot_3, ESplot_4, ESplot_5, ESplot_6, ESplot_7, ESplot_8, nrow=2, ncol=3)
ggsave(ESplot_all, file = paste0("ESplot_3by3.jpg"))


# Accuracy, Sensitivity, Specificity, AUC

library(RColorBrewer)

# Inputs from PEAP04
load("DEP_lambda_0.1.Rdata")
cm_peap <- cm_DEP

load("DEP_limma_0.1.Rdata")
cm_limma <- cm_DEG

load("DEP_edgeR_0.1.Rdata")
cm_edgeR <- cm_DEG

load("DEP_DESeq2_0.1.Rdata")
cm_DESeq <- cm_DEG




PhenotypeName <- c("Healthy", "Bacterial", "Viral")

for(i in 0:2){
  df = as.data.frame(matrix(data = NA, nrow = 12, ncol=3))
  colnames(df) <- c("Method", "DataName", "Data")
  phenotype_index = i+1
  df$Method <- c(rep("Limma", 3), rep("edgeR", 3), rep("DESeq2", 3), rep("PEAP", 3))
  df$DataName <- rep(c("Sensitivity", "Sepcificity", "AUC"), 4)
  data_onephenotype <- c(cm_limma$byClass[phenotype_index, 1:2], auc_all[[phenotype_index]][1], cm_edgeR$byClass[phenotype_index, 1:2], auc_all[[phenotype_index]][2],cm_DESeq$byClass[phenotype_index, 1:2], auc_all[[phenotype_index]][3],cm_peap$byClass[phenotype_index, 1:2], auc_all[[phenotype_index]][4])
  df$Data <- data_onephenotype
  df$Method <- factor(df$Method,levels = c("PEAP", "Limma", "edgeR", "DESeq2"))
  
  
  df_plot <- ggplot(df, aes(x=DataName,y=Data,fill=Method)) + 
    geom_col(position = 'dodge') +
    ggtitle(paste0(PhenotypeName[phenotype_index])) +
    scale_fill_manual(values=c(rgb(248,230,32, maxColorValue = 255), rgb(145,213,66, maxColorValue = 255), rgb(53,183,119, maxColorValue = 255), rgb(38,104,141, maxColorValue = 255))) +
    xlab("") +
    ylab("Value") +
    theme(plot.title = element_text(size = 30), axis.text = element_text(size = 20), axis.title.y = element_text(size = 20), legend.text = element_text(size = 18), legend.title = element_text(size = 20), legend.position = "bottom")
  ggsave(df_plot, filename = paste0("result_", PhenotypeName[phenotype_index],".jpg"))
}


df <- data.frame(Method = c("Limma", "edgeR","DESeq2","PEAP"), acc = c(cm_limma$overall[1],cm_edgeR$overall[1],cm_DESeq$overall[1],cm_peap$overall[1])) 
df$Method <- factor(df$Method,levels = c("PEAP", "Limma", "edgeR", "DESeq2"))

df_plot <- ggplot(df, aes(x=Method,y=acc, fill=Method)) + 
  geom_col() +
  ggtitle("Accuracy") +
  scale_fill_manual(values=c(rgb(248,230,32, maxColorValue = 255), rgb(145,213,66, maxColorValue = 255), rgb(53,183,119, maxColorValue = 255), rgb(38,104,141, maxColorValue = 255))) +
  xlab("") +
  ylab("Value") +
  theme(plot.title = element_text(size = 30), axis.text = element_text(size = 20), axis.title.y = element_text(size = 20), legend.text = element_text(size = 18), legend.title = element_text(size = 20), legend.position = "bottom")
ggsave(df_plot, filename = "Accuracy_result.jpg")





rm(list = ls())

library(ggplot2)
library(pheatmap)
library(RColorBrewer)

set.seed(1)
setwd("D:/PEAP/Infection/")

heatmap_sample <- matrix(data=rnorm(15*5, 1, 2), nrow=5, ncol = 15)

heatmap_sample[1,1:5] = heatmap_sample[1,1:5] + 15

heatmap_sample[2:3,6:10] = heatmap_sample[2:3,6:10] + 15

heatmap_sample[4:5,11:15] = heatmap_sample[4:5,11:15] + 15



colnames(heatmap_sample) <- paste0("sample", seq(1,15,by=1))
rownames(heatmap_sample) <- paste0("gene", seq(1,5,by=1))
pheatmap(heatmap_sample, cluster_rows = F, cluster_cols = F)
annotation_label_col <- as.factor(c(rep("phenotype1", 5), rep("phenotype2", 5), rep("phenotype3", 5)))
main_groups_col <- data.frame(row.names = colnames(heatmap_sample), phenotype=annotation_label_col)

annotation_label_row <- as.factor(c(rep("pathway1", 1), rep("pathway2", 2), rep("pathway3", 2)))
main_groups_row <- data.frame(row.names = rownames(heatmap_sample), pathway = annotation_label_row)

ann_colors=list(phenotype=c(phenotype1='#55AF7B',phenotype2='#FAC230',phenotype3='#EB4537'), pathway=c(pathway1=brewer.pal(3,"Set3")[1], pathway2=brewer.pal(3,"Set3")[2], pathway3=brewer.pal(3,"Set3")[3]))
pheatmap(heatmap_sample, cluster_rows = F, cluster_cols = F, annotation_row = main_groups_row, annotation_col = main_groups_col, annotation_colors = ann_colors, cellwidt = 20,cellheight = 20, color = rev(colorRampPalette(brewer.pal(9, "PiYG"))(100)))


heatmap_sample <- matrix(data=rnorm(15*10, 1, 2), nrow=10, ncol = 15)

heatmap_sample[1:3,1:5] = heatmap_sample[1:3,1:5] + 15

heatmap_sample[4:6,6:10] = heatmap_sample[4:6,6:10] + 15

heatmap_sample[7:10,11:15] = heatmap_sample[7:10,11:15] + 15


colnames(heatmap_sample) <- paste0("sample", seq(1,15,by=1))
rownames(heatmap_sample) <- paste0("genepairs", seq(1,10,by=1))
pheatmap(heatmap_sample, cluster_rows = F, cluster_cols = F)
annotation_label_col <- as.factor(c(rep("phenotype1", 5), rep("phenotype2", 5), rep("phenotype3", 5)))
main_groups_col <- data.frame(row.names = colnames(heatmap_sample), phenotype=annotation_label_col)

annotation_label_row <- as.factor(c(rep("pathway1", 3), rep("pathway2", 3), rep("pathway3", 4)))
main_groups_row <- data.frame(row.names = rownames(heatmap_sample), pathway = annotation_label_row)

ann_colors=list(phenotype=c(phenotype1='#55AF7B',phenotype2='#FAC230',phenotype3='#EB4537'), pathway=c(pathway1=brewer.pal(3,"Set3")[1], pathway2=brewer.pal(3,"Set3")[2], pathway3=brewer.pal(3,"Set3")[3]))
pheatmap(heatmap_sample, cluster_rows = F, cluster_cols = F, annotation_row = main_groups_row, annotation_col = main_groups_col, annotation_colors = ann_colors, cellwidt = 20,cellheight = 20, color = rev(colorRampPalette(brewer.pal(9, "PiYG"))(100)))
















### for pathway
set.seed(123)
heatmap_sample <- matrix(data=rnorm(3*3, 0, 1.5), nrow=3, ncol = 3)

heatmap_sample[1,1] = heatmap_sample[1,1] + 30
heatmap_sample[2,1] = heatmap_sample[2,1] + 25

heatmap_sample[3,2] = heatmap_sample[4,2] + 29


colnames(heatmap_sample) <- paste0("phenotype", seq(1,3,by=1))
rownames(heatmap_sample) <- paste0("pathway", seq(1,3,by=1))
annotation_label_col <- as.factor(c(rep("phenotype1", 1), rep("phenotype2", 1), rep("phenotype3", 1)))
main_groups_col <- data.frame(row.names = colnames(heatmap_sample), phenotype=annotation_label_col)

annotation_label_row <- as.factor(c(rep("pathway1", 1), rep("pathway2", 1), rep("pathway3", 1)))
main_groups_row <- data.frame(row.names = rownames(heatmap_sample), pathway = annotation_label_row)

ann_colors=list(phenotype=c(phenotype1='#55AF7B',phenotype2='#FAC230',phenotype3='#EB4537'), pathway=c(pathway1=brewer.pal(3,"Set3")[1], pathway2=brewer.pal(3,"Set3")[2], pathway3=brewer.pal(3,"Set3")[3]))

pheatmap(heatmap_sample, cluster_rows = F, cluster_cols = F, annotation_colors = ann_colors,cellwidt = 80,cellheight = 60, color = colorRampPalette(brewer.pal(9, "OrRd"))(100), annotation_row = main_groups_row, annotation_col = main_groups_col)





rm(list = ls())
library(fgsea)
load("Nonoverlap_top100GenesOfEachGenepairs_Infection_20240406.Rdata")

ESres_top5 <- NULL
for(i in 1:3){
  order_pval <- order(ESres_collection[,i + 5])
  ESres_cur <- ESres_collection[order_pval[1:5], ]
  for(j in 1:5){
    # Ensure all enrichment scores are positive
    if(ESres_cur[j, i + 2] < 0){
      ESres_cur[j, 1:2] <- ESres_cur[j, 2:1]
      ESres_cur[j, 3:5] <- -ESres_cur[j, 3:5]
    }
  }
  ESres_top5 <- rbind(ESres_top5, ESres_cur)
}




################################################ 
# Draw scatter plot with marignal distribution #
################################################
library(ggplot2)
library(ggExtra)

used_gene_list <- c(ESres_top5$gene1, ESres_top5$gene2)



plotEnrichment_plus <- function (pathway, stats, gseaParam = 1, ticksSize = 0.2) 
{
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
                          returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  diff <- (max(tops) - min(bottoms))/8
  x = y = NULL
  g <- ggplot(toPlot, aes(x = x, y = y)) + 
    geom_point(color = "black", size = 0.1) + 
    geom_hline(yintercept = max(tops), colour = "red", linetype = "dashed") + 
    geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed") + geom_hline(yintercept = 0, colour = "black") + geom_line(color = "black") + theme_bw() + 
    geom_segment(data = data.frame(x = pathway), mapping = aes(x = x, y = 0.8, xend = x, yend = 0.95), size = ticksSize) + 
    theme(panel.border = element_blank(), panel.grid.minor = element_blank()) + 
    ylim(c(-1,1)) +
    xlim(c(0,2500)) +
    labs(x = "Rank", y = "")
  g
}


################################################ 
# Draw scatter plot with marignal distribution #
################################################
library(ggplot2)
library(ggExtra)
example_df <- iris[iris$Species != "versicolor", ]

p <- ggplot(example_df) +
  geom_point(aes(x = Sepal.Length, y = Sepal.Width, color = Species), 
             alpha = 0.6, shape = 16, size = 3) +  
  scale_color_brewer(palette = "Dark2") + 
  theme_bw() +
  theme(legend.position = "bottom") + 
  labs(x = "Gene 1", y = "Gene 2") 

ggMarginal(p, type = "density", groupColour = TRUE, groupFill = TRUE)

example_df <- iris
example_df$new_label <- factor(sample(1:2, size = nrow(iris), replace = TRUE))


######################################
# Scatter plot with marginal density #
######################################
p <- ggplot(example_df) +
  geom_point(aes(x = Sepal.Length, y = Sepal.Width, color = new_label), 
             alpha = 0.6, shape = 16, size = 3) +  
  scale_color_brewer(palette = "Dark2") + 
  theme_bw() +
  theme(legend.position = "bottom") + 
  labs(x = "Gene 1", y = "Gene 2") 

ggMarginal(p, type = "density", groupColour = TRUE, groupFill = TRUE)

ggsave(p, file = "density_plot_lowES_eg.jpg")



library(fgsea)

sample_name_dep <- colnames(AllGeneExp)
sample_pheno <- pheno_sample

cell_name_list <- list("0"=sample_name_dep[which(sample_pheno==0)], "1"=sample_name_dep[which(sample_pheno==1)],"2"=sample_name_dep[which(sample_pheno==2)])

sample_genepair_ind = 4
index2=which(rownames(AllGeneExp)==ESres_top5$gene1[sample_genepair_ind])
index1=which(rownames(AllGeneExp)==ESres_top5$gene2[sample_genepair_ind])
esdata <- AllGeneExp[index1,]-AllGeneExp[index2,]

esdata <- as.numeric(esdata)
esdata <- as.vector(esdata)
names(esdata) <- sample_name_dep

p_index=2

plotEnrichment_plus(cell_name_list[[p_index+1]], esdata)+labs(title=paste0(rownames(AllGeneExp)[index1], "_", rownames(AllGeneExp)[index2], "_pheotype_", p_index))

ESplot <- plotEnrichment_plus(cell_name_list[[p_index+1]], esdata)+labs(title=paste0(rownames(AllGeneExp)[index1], "_", rownames(AllGeneExp)[index2], "_pheotype_", p_index))
ggsave(ESplot, file="ES_example_esplot.jpg")



############################
# supplementary table      #
############################


rm(list = ls())
setwd("D:/PEAP/Infection/")
load("Infection_ESres.Rdata")
load("marker_genepairs.Rdata")

marker_gp_ind <- paste0(ESres_collection$gene1,"_", ESres_collection$gene2)


ESres_marker <- data.frame()

for(i in 1:3){
  match_ind <- match(selecP_list[[i]], marker_gp_ind)
  ESres_marker <- rbind(ESres_marker, ESres_collection[match_ind,])
}
colnames(ESres_marker)[c(3:8, 13:15)] <- c("ES.healthy", "ES.bacterial", "ES.viral", "pvalue.healthy", "pvalue.bacterial", "pvalue.viral", "DE.healthy", "DE.bacterial", "DE.viral")
ESres_marker = ESres_marker[,c(1,2,3,6,9,13,4,7,10,14,5,8,11,15,12)]

selected_phenotype_ind = data.frame(Selected_phenotype = c(rep("healthy", length(selecP_list[[1]])),rep("bacterial", length(selecP_list[[2]])) , rep("viral", length(selecP_list[[3]]))))
ESres_marker = cbind(selected_phenotype_ind, ESres_marker)

write.csv(ESres_marker, "supplementary_table1.csv", row.names = F)








