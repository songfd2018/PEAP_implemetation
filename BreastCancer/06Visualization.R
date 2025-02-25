rm(list = ls())
library(ggplot2)
library(tidyr)
library(dplyr)
library(readxl)
library(gridExtra)
library(fgsea)
library(tidyverse)
library(broom)
library(glmnet)
library(survival)
library(survminer)
library(easyTCGA)
setwd("D:/PEAP/PEAP_02")
load("result.RData")
C_index <- read_excel("D:/PEAP/PEAP_02/New_data/C_index_20241030.xlsx")
load("ESres_BRCA.Rdata")

PEAP_pvalue <- read.csv("ESres_marker_BRCA_v4.csv")
PEAP_pvalue <- PEAP_pvalue[,c(2,3,5,9,13,17,21)]
PEAP_pvalue[,-c(1:2)]=-log10(PEAP_pvalue[,-c(1:2)])

### hazard ratio
surv_obj <- Surv(surv_peap$surv_time, surv_peap$vital_num)
cox = as.vector(peap_fit$beta)
names(cox) = as.vector(row.names(peap_fit$beta))
cox = cox[cox!=0]
covariates <- sapply(names(cox), function(x) {
  if (grepl("-", x, fixed = TRUE)) {
    return(paste0("`", x, "`"))
  } else {
    return(x)
  }
})

formula_string <- paste("surv_obj", "~", paste(covariates, collapse = " + "))
hr = coxph(as.formula(formula_string), data = surv_peap)%>%
  tidy() %>%
  mutate(upper = estimate + 1.96 * std.error,
         lower = estimate - 1.96 * std.error)

hr$term[38] = "C1QB_HLA-DRA"
hr_genepairs_label = match(hr$term, Selected_genepairs)
t = 1
hr_phenotype = c(0)
for(i in 1:nrow(hr)){
  if(hr_genepairs_label[i]>=1500*t){
    hr_phenotype = c(hr_phenotype, i-1)
    t=t+1}
}
hr_phenotype = c(hr_phenotype, nrow(hr))

PEAP_pvalue_matrix <- PEAP_pvalue[match(hr$term, paste0(PEAP_pvalue$gene2,"_",PEAP_pvalue$gene1)),]
PEAP_pvalue_matrix[is.na(PEAP_pvalue_matrix$gene1),] = PEAP_pvalue[match(hr$term[is.na(PEAP_pvalue_matrix$gene1)], paste0(PEAP_pvalue$gene1,"_",PEAP_pvalue$gene2)),]


for(i in 1:5){
  hr[(hr_phenotype[i]+1):hr_phenotype[i+1],] = hr[(hr_phenotype[i]+1):hr_phenotype[i+1],][order(PEAP_pvalue_matrix[(hr_phenotype[i]+1):hr_phenotype[i+1],i+2], decreasing = T),]
}

PEAP_pvalue_matrix <- PEAP_pvalue[match(hr$term, paste0(PEAP_pvalue$gene2,"_",PEAP_pvalue$gene1)),]
PEAP_pvalue_matrix[is.na(PEAP_pvalue_matrix$gene1),] = PEAP_pvalue[match(hr$term[is.na(PEAP_pvalue_matrix$gene1)], paste0(PEAP_pvalue$gene1,"_",PEAP_pvalue$gene2)),]


hr %>%
  mutate(across(all_of(c("estimate", "lower", "upper")), exp)) %>%
  filter(estimate > 1) %>%
  ggplot(aes(x = estimate, y = reorder(term, c(25:1)))) +
  geom_linerange(aes(xmin = lower, xmax = upper), size = 1, alpha = 0.5) +
  geom_point(size = 2,color="red",alpha=0.5) +
  theme_minimal(base_size = 10) +
  scale_color_manual(values = c("green4", "red3"), guide = "none") +
  xlim(c(0, 6.5)) +
  labs(title = "Hazard ratio for Gene Pairs", y = NULL,
       x = "Hazard ratio estimate (95% C.I.)") +
  theme(axis.text.y = element_text(hjust = 0, size = 10),
        axis.text.x = element_text(size = 10),
        plot.title = element_text(hjust = 0.4),
        axis.title.x = element_text(hjust = 0.4))



### Supplementary
hr %>%
  mutate(across(all_of(c("estimate", "lower", "upper")), exp)) %>%
  filter(estimate <= 1) %>%
  ggplot(aes(x = estimate, y = reorder(term, c(29:1)))) +
  geom_linerange(aes(xmin = lower, xmax = upper), size = 1, alpha = 0.5) +
  geom_point(size = 2,color="red",alpha=0.5) +
  theme_minimal(base_size = 10) +
  scale_color_manual(values = c("green4", "red3"), guide = "none") +
  xlim(c(0, 6.5)) +
  labs(title = "Hazard ratio for Gene Pairs", y = NULL,
       x = "Hazard ratio estimate (95% C.I.)") +
  theme(axis.text.y = element_text(hjust = 0, size = 10),
        axis.text.x = element_text(size = 10),
        plot.title = element_text(hjust = 0.4),
        axis.title.x = element_text(hjust = 0.4))

ggsave("Hazard Ratio Supp.pdf")





hr_1 = hr %>%
  mutate(across(all_of(c("estimate", "lower", "upper")), exp)) %>%
  filter(estimate > 1)

PEAP_pvalue_matrix <- PEAP_pvalue[match(hr_1$term, paste0(PEAP_pvalue$gene2,"_",PEAP_pvalue$gene1)),]
PEAP_pvalue_matrix[is.na(PEAP_pvalue_matrix$gene1),] = PEAP_pvalue[match(hr_1$term[is.na(PEAP_pvalue_matrix$gene1)], paste0(PEAP_pvalue$gene1,"_",PEAP_pvalue$gene2)),]

rownames(PEAP_pvalue_matrix) = hr_1$term
PEAP_pvalue_matrix <- PEAP_pvalue_matrix[3:7]

colnames(PEAP_pvalue_matrix) <- name_pheno
rownames(PEAP_pvalue_matrix) <- paste0("Gene Pairs ",1:25)

library(pheatmap)
library(RColorBrewer)
pheatmap(PEAP_pvalue_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE, 
         color = colorRampPalette(brewer.pal(9, "OrRd"))(100))

heatmap_sample <- PEAP_pvalue_matrix
annotation_label_col <- as.factor(name_pheno)
main_groups_col <- data.frame(row.names = colnames(heatmap_sample), phenotype=annotation_label_col)

annotation_label_row <- as.factor(c(rep(name_pheno[1], hr_phenotype[2]-hr_phenotype[1]), rep(name_pheno[2], hr_phenotype[3]-hr_phenotype[2]), rep(name_pheno[3], hr_phenotype[4]-hr_phenotype[3]), rep(name_pheno[4], hr_phenotype[5]-hr_phenotype[4]), rep(name_pheno[5], hr_phenotype[6]-hr_phenotype[5])))
annotation_label_row <- annotation_label_row[match(hr_1$term, hr$term)]
main_groups_row <- data.frame(row.names = rownames(heatmap_sample), pathway = annotation_label_row)

ann_colors=list(phenotype=c(Basal=brewer.pal(5,"Set1")[1],Her2=brewer.pal(5,"Set1")[2], LumA=brewer.pal(5,"Set1")[3], LumB=brewer.pal(5,"Set1")[4], Normal=brewer.pal(5,"Set1")[5]), pathway=c(Basal=brewer.pal(5,"Set1")[1],Her2=brewer.pal(5,"Set1")[2], LumA=brewer.pal(5,"Set1")[3], LumB=brewer.pal(5,"Set1")[4], Normal=brewer.pal(5,"Set1")[5]))

bk <- c(seq(0,15,by=0.1))
p <- pheatmap(heatmap_sample, cluster_rows = T, cluster_cols = F, annotation_colors = ann_colors, color = colorRampPalette(brewer.pal(9, "OrRd"))(150),annotation_col = main_groups_col, breaks = bk)
# ,cellwidt = 80,cellheight = 60

row_cluster_order = p$tree_row$order
hr_1 = hr_1[row_cluster_order, ]
hr_1 %>%
  ggplot(aes(x = estimate, y = reorder(term, c(25:1)))) +
  geom_linerange(aes(xmin = lower, xmax = upper), size = 1, alpha = 0.5) +
  geom_point(size = 2,color="red",alpha=0.5) +
  theme_minimal(base_size = 10) +
  scale_color_manual(values = c("green4", "red3"), guide = "none") +
  xlim(c(0, 6.5)) +
  labs(title = "Hazard ratio for Gene Pairs", y = NULL,
       x = "Hazard ratio estimate (95% C.I.)") +
  theme(axis.text.y = element_text(hjust = 0, size = 10),
        axis.text.x = element_text(size = 10),
        plot.title = element_text(hjust = 0.4),
        axis.title.x = element_text(hjust = 0.4))



### p value heatmap supp
hr_2 = hr %>%
  mutate(across(all_of(c("estimate", "lower", "upper")), exp)) %>%
  filter(estimate <= 1)




PEAP_pvalue_matrix <- PEAP_pvalue[match(hr_2$term, paste0(PEAP_pvalue$gene2,"_",PEAP_pvalue$gene1)),]
PEAP_pvalue_matrix[is.na(PEAP_pvalue_matrix$gene1),] = PEAP_pvalue[match(hr_2$term[is.na(PEAP_pvalue_matrix$gene1)], paste0(PEAP_pvalue$gene1,"_",PEAP_pvalue$gene2)),]

rownames(PEAP_pvalue_matrix) = hr_2$term
PEAP_pvalue_matrix <- PEAP_pvalue_matrix[3:7]

colnames(PEAP_pvalue_matrix) <- name_pheno
rownames(PEAP_pvalue_matrix) <- paste0("Gene Pairs ",1:29)

library(pheatmap)
library(RColorBrewer)
pheatmap(PEAP_pvalue_matrix,
         cluster_rows = T,
         cluster_cols = FALSE, 
         color = colorRampPalette(brewer.pal(9, "OrRd"))(100))

heatmap_sample <- PEAP_pvalue_matrix
rownames(heatmap_sample) = hr_2$term
annotation_label_col <- as.factor(name_pheno)
main_groups_col <- data.frame(row.names = colnames(heatmap_sample), phenotype=annotation_label_col)

annotation_label_row <- as.factor(c(rep(name_pheno[1], hr_phenotype[2]-hr_phenotype[1]), rep(name_pheno[2], hr_phenotype[3]-hr_phenotype[2]), rep(name_pheno[3], hr_phenotype[4]-hr_phenotype[3]), rep(name_pheno[4], hr_phenotype[5]-hr_phenotype[4]), rep(name_pheno[5], hr_phenotype[6]-hr_phenotype[5])))
annotation_label_row <- annotation_label_row[match(hr_2$term, hr$term)]
main_groups_row <- data.frame(row.names = hr_2$term, enriched_phenotype = annotation_label_row)

ann_colors=list(phenotype=c(Basal=brewer.pal(5,"Set1")[1],Her2=brewer.pal(5,"Set1")[2], LumA=brewer.pal(5,"Set1")[3], LumB=brewer.pal(5,"Set1")[4], Normal=brewer.pal(5,"Set1")[5]), enriched_phenotype=c(Basal=brewer.pal(5,"Set1")[1],Her2=brewer.pal(5,"Set1")[2], LumA=brewer.pal(5,"Set1")[3], LumB=brewer.pal(5,"Set1")[4], Normal=brewer.pal(5,"Set1")[5]))

bk <- c(seq(0,15,by=0.1))
pheatmap(heatmap_sample, cluster_rows = F,filename = "heatmap_kstest_pval_Supp.pdf", cluster_cols = F, annotation_colors = ann_colors, color = colorRampPalette(brewer.pal(9, "OrRd"))(150), annotation_row = main_groups_row, annotation_col = main_groups_col, breaks = bk)











### ES plot
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
    geom_hline(yintercept = max(tops),
               colour = "red",
               linetype = "dashed") + 
    geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed") + geom_hline(yintercept = 0, colour = "black") + geom_line(color = "black") + theme_bw() + 
    geom_segment(data = data.frame(x = pathway), mapping = aes(x = x, y = 0.8, xend = x, yend = 0.95), size = ticksSize) + 
    theme(panel.border = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 15)) + 
    ylim(c(-1,1)) +
    xlim(c(0,1000)) +
    labs(x = "Rank", y = "Enrichment Score")
  g
}

gene_pair_name <- colnames(DEP_peap)[-c(ncol(DEP_peap)-1,ncol(DEP_peap))]

sample_name <- colnames(GeneExp_BRCA)
sample_pheno <- pheno_sample_tcga


set.seed(1234)

cell_name_list <- list("0"=sample_name[which(sample_pheno=="Basal")],
                       "1"=sample_name[which(sample_pheno=="Her2")],
                       "2"=sample_name[which(sample_pheno=="LumA")],
                       "3"=sample_name[which(sample_pheno=="LumB")],
                       "4"=sample_name[which(sample_pheno=="Normal")])
esdata_all <- as.data.frame(t(DEP_peap))

gene_pair_ref <- paste0(pmin(top1500_Genes[,1],top1500_Genes[,2]),
                        "_",
                        pmax(top1500_Genes[,1],top1500_Genes[,2]))




i=which(gene_pair_name=="ALOX15B_EPHX2")
print(i)
gene_pair_class <- -1
esdata <- as.numeric(esdata_all[i,])
esdata <- as.vector(esdata)
names(esdata) <- sample_name
current_gene_pair_name <- rownames(esdata_all[i,])

if(current_gene_pair_name %in% gene_pair_ref[1:1500]){
  gene_pair_class = 0
}else if(current_gene_pair_name %in% gene_pair_ref[1501:3000]){
  gene_pair_class = 1
}else if(current_gene_pair_name %in% gene_pair_ref[3001:4500]){
  gene_pair_class = 2
}else if(current_gene_pair_name %in% gene_pair_ref[4501:6000]){
  gene_pair_class = 3
}else if(current_gene_pair_name %in% gene_pair_ref[6001:7500]){
  gene_pair_class = 4
}

p_index=0
ESplot_0 <- plotEnrichment_plus(cell_name_list[[p_index+1]], esdata)+labs(title=paste0(current_gene_pair_name, "s",gene_pair_class, "a", p_index)) + theme(plot.title = element_text(size = 10))

p_index=1
ESplot_1 <- plotEnrichment_plus(cell_name_list[[p_index+1]], esdata)+labs(title=paste0(current_gene_pair_name, "s",gene_pair_class, "a", p_index)) + theme(plot.title = element_text(size = 10))

p_index=2
ESplot_2 <- plotEnrichment_plus(cell_name_list[[p_index+1]], esdata)+labs(title=paste0(current_gene_pair_name, "s",gene_pair_class, "a", p_index)) + theme(plot.title = element_text(size = 10))

p_index=3
ESplot_3 <- plotEnrichment_plus(cell_name_list[[p_index+1]], esdata)+labs(title=paste0(current_gene_pair_name, "s",gene_pair_class, "a", p_index)) + theme(plot.title = element_text(size = 10))

p_index=4
ESplot_4 <- plotEnrichment_plus(cell_name_list[[p_index+1]], esdata)+labs(title=paste0(current_gene_pair_name, "s",gene_pair_class, "a", p_index)) + theme(plot.title = element_text(size = 10))


ESplot_all <- grid.arrange(ESplot_0,ESplot_3,nrow=1, ncol=2)


library(GseaVis)
library(clusterProfiler)

m_t2g = tibble(gs_name = pheno_sample_tcga, entrez_gene = names(esdata))
genelist <- sort(esdata,decreasing = T)

x = GSEA(geneList = genelist, TERM2GENE = m_t2g)
x <- GSEA(genelist, 
          TERM2GENE = m_t2g,
          minGSSize = 1,
          maxGSSize = 500,
          pvalueCutoff = 1,
          pAdjustMethod = "BH",
          seed = 1)

library(RColorBrewer)
gseaNb_plus <- function (object = NULL, subPlot = 3, lineSize = 0.8, geneSetID = NULL, 
          rmSegment = FALSE, termWidth = 40, segCol = "red", addGene = NULL, 
          geneCol = NULL, arrowAngle = 20, arrowLength = 0.2, arrowEnd = "first", 
          arrowType = "closed", curveCol = c("#76BA99", "#EB4747", 
                                             "#996699"), htCol = c("#08519C", "#A50F15"), rankCol = c("#08519C", 
                                                                                                      "white", "#A50F15"), rankSeq = 5000, htHeight = 0.3, 
          force = 20, max.overlaps = 50, geneSize = 4, newGsea = FALSE, 
          addPoint = TRUE, newCurveCol = c("#336699", "white", "#993399"), 
          newHtCol = c("#336699", "white", "#993399"), rmHt = FALSE, 
          addPval = FALSE, pvalX = 0.9, pvalY = 0.9, pvalSize = 4, 
          pCol = "grey30", pHjust = 1, rmPrefix = TRUE, nesDigit = 2, 
          pDigit = 2, markTopgene = FALSE, topGeneN = 5, kegg = FALSE, 
          legend.position = "right", add.geneExpHt = FALSE, exp = NULL, 
          scale.exp = TRUE, sample.order = NULL, exp.col = c("blue", 
                                                             "white", "red"), ht.legend = TRUE, ght.relHight = 0.4, 
          ght.geneText.size = 6, ght.facet = FALSE, ght.facet.scale = "free", 
          termID.order = NULL, rank.gene = NULL, rank.gene.nudgey = 2) 
{
  gsdata <- purrr::map_df(geneSetID, function(setid) {
    gsInfo(object, geneSetID = setid) %>% dplyr::mutate(id = setid)
  })
  if (kegg == FALSE) {
    gsdata1 <- purrr::map_df(unique(gsdata$Description), 
                             function(setid) {
                               tmp <- gsdata %>% dplyr::filter(Description == 
                                                                 setid) %>% dplyr::mutate(gene_name = names(object@geneList)) %>% 
                                 dplyr::filter(position == 1)
                             })
  }
  else {
    gene2Symbol <- object@gene2Symbol %>% data.frame()
    gsdata1 <- purrr::map_df(unique(gsdata$Description), 
                             function(setid) {
                               tmp <- gsdata %>% dplyr::filter(Description == 
                                                                 setid) %>% dplyr::mutate(gene_name = gene2Symbol$.) %>% 
                                 dplyr::filter(position == 1)
                             })
  }
  data_ga <- data.frame(object) %>% dplyr::filter(ID %in% geneSetID)
  data_ga <- data_ga[unique(gsdata$id), ]
  niceTit <- purrr::map_chr(unique(gsdata$Description), function(x) {
    tit <- unlist(strsplit(x, split = "_"))
    if (length(tit) == 1) {
      niceTit <- paste(stringr::str_to_title(tit[1:length(tit)]), 
                       collapse = " ") %>% stringr::str_wrap(., width = termWidth)
    }
    else {
      if (rmPrefix == TRUE) {
        niceTit <- paste(stringr::str_to_title(tit[2:length(tit)]), 
                         collapse = " ") %>% stringr::str_wrap(., width = termWidth)
      }
      else {
        niceTit <- paste(stringr::str_to_title(tit[1:length(tit)]), 
                         collapse = " ") %>% stringr::str_wrap(., width = termWidth)
      }
    }
  })
  if (length(geneSetID) != 1) {
    ledend.t <- niceTit
    niceTit <- ""
  }
  if (length(geneSetID) == 1) {
    line <- ggplot2::geom_line(ggplot2::aes_(color = ~runningScore), 
                               size = lineSize)
    line.col <- ggplot2::scale_color_gradient(low = curveCol[1], 
                                              high = curveCol[2])
    legend.position = "none"
  }
  else {
    mulcol <- curveCol
    names(mulcol) <- unique(gsdata$Description)
    line <- ggplot2::geom_line(ggplot2::aes_(color = ~Description), 
                               size = lineSize)
    line.col <- ggplot2::scale_color_manual(values = mulcol, 
                                            labels = ledend.t, name = "Phenotype")
    legend.position = legend.position
  }
  pcurve <- ggplot2::ggplot(gsdata, ggplot2::aes_(x = ~x, y = ~runningScore)) + 
    line + line.col + ggplot2::geom_hline(yintercept = 0, 
                                          size = lineSize, color = "black", lty = "dashed") + ggplot2::theme_bw(base_size = 14) + 
    ggplot2::scale_x_continuous(expand = c(0, 0)) + ggplot2::theme(legend.position = legend.position, 
                                                                   legend.box.background = ggplot2::element_blank(), plot.title = ggplot2::element_text(hjust = 0.5), 
                                                                   panel.grid = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), 
                                                                   axis.text.x = ggplot2::element_blank(), axis.line.x = ggplot2::element_blank(), 
                                                                   axis.title.x = ggplot2::element_blank(), legend.background = ggplot2::element_rect(fill = "transparent"), 
                                                                   plot.margin = ggplot2::margin(t = 0.2, r = 0.2, b = 0, 
                                                                                                 l = 0.2, unit = "cm")) + ggplot2::ylab("Running Sum") + 
    ggplot2::ggtitle(niceTit)
  midpoint <- sum(range(gsdata$runningScore))/2
  pnew <- ggplot2::ggplot(gsdata, ggplot2::aes_(x = ~x, y = ~runningScore, 
                                                color = ~runningScore)) + ggplot2::geom_hline(yintercept = 0, 
                                                                                              size = lineSize, color = "black", lty = "dashed") + ggplot2::geom_line(size = lineSize) + 
    ggplot2::geom_segment(data = gsdata1, ggplot2::aes_(xend = ~x, 
                                                        yend = 0)) + ggplot2::theme_bw(base_size = 14) + 
    ggplot2::scale_color_gradient2(low = newCurveCol[1], 
                                   mid = newCurveCol[2], high = newCurveCol[3], midpoint = midpoint) + 
    ggplot2::scale_x_continuous(expand = c(0, 0)) + ggplot2::theme(legend.position = "none", 
                                                                   plot.title = ggplot2::element_text(hjust = 0.5), axis.ticks.x = ggplot2::element_blank(), 
                                                                   axis.text.x = ggplot2::element_blank(), axis.line.x = ggplot2::element_blank(), 
                                                                   axis.title.x = ggplot2::element_blank(), legend.background = ggplot2::element_rect(fill = "transparent"), 
                                                                   plot.margin = ggplot2::margin(t = 0.2, r = 0.2, b = 0, 
                                                                                                 l = 0.2, unit = "cm")) + ggplot2::ylab("Running Enrichment Score") + 
    ggplot2::ggtitle(niceTit)
  if (addPoint == TRUE) {
    panother <- pnew + ggplot2::geom_point()
  }
  else {
    panother <- pnew
  }
  if (newGsea == FALSE) {
    pcurveRes <- pcurve
  }
  else {
    pcurveRes <- panother
  }
  if (is.null(addGene)) {
    plabel <- pcurveRes
  }
  else {
    if (markTopgene == TRUE) {
      geneLabel <- gsdata1 %>% dplyr::arrange(x) %>% dplyr::slice_head(n = topGeneN)
    }
    else {
      geneLabel <- gsdata1 %>% dplyr::filter(gene_name %in% 
                                               addGene)
    }
    if (nrow(geneLabel) == 0) {
      message("Your gene is not in this pathway! Please choose again!")
    }
    else {
      if (rmSegment == TRUE) {
        if (is.null(geneCol)) {
          plabel <- pcurveRes + ggrepel::geom_text_repel(data = geneLabel, 
                                                         ggplot2::aes_(label = ~gene_name), force = force, 
                                                         max.overlaps = max.overlaps, size = geneSize, 
                                                         fontface = "italic", arrow = ggplot2::arrow(angle = arrowAngle, 
                                                                                                     length = ggplot2::unit(arrowLength, "cm"), 
                                                                                                     ends = arrowEnd, type = arrowType))
        }
        else {
          plabel <- pcurveRes + ggrepel::geom_text_repel(data = geneLabel, 
                                                         ggplot2::aes_(label = ~gene_name), force = force, 
                                                         max.overlaps = max.overlaps, size = geneSize, 
                                                         fontface = "italic", color = geneCol, arrow = ggplot2::arrow(angle = arrowAngle, 
                                                                                                                      length = ggplot2::unit(arrowLength, "cm"), 
                                                                                                                      ends = arrowEnd, type = arrowType))
        }
      }
      else {
        if (is.null(geneCol)) {
          plabel <- pcurveRes + ggplot2::geom_segment(data = geneLabel, 
                                                      ggplot2::aes_(xend = ~x, yend = 0), color = segCol) + 
            ggrepel::geom_text_repel(data = geneLabel, 
                                     ggplot2::aes_(label = ~gene_name), force = force, 
                                     max.overlaps = max.overlaps, size = geneSize, 
                                     fontface = "italic", arrow = ggplot2::arrow(angle = arrowAngle, 
                                                                                 length = ggplot2::unit(arrowLength, "cm"), 
                                                                                 ends = arrowEnd, type = arrowType))
        }
        else {
          plabel <- pcurveRes + ggplot2::geom_segment(data = geneLabel, 
                                                      ggplot2::aes_(xend = ~x, yend = 0), color = segCol) + 
            ggrepel::geom_text_repel(data = geneLabel, 
                                     ggplot2::aes_(label = ~gene_name), force = force, 
                                     max.overlaps = max.overlaps, size = geneSize, 
                                     fontface = "italic", color = geneCol, arrow = ggplot2::arrow(angle = arrowAngle, 
                                                                                                  length = ggplot2::unit(arrowLength, "cm"), 
                                                                                                  ends = arrowEnd, type = arrowType))
        }
      }
    }
  }
  if (addPval == TRUE) {
    pLabel <- paste0("NES: ", round(data_ga$NES, digits = nesDigit), 
                     "\n", "Pvalue: ", ifelse(data_ga$pvalue < 0.001, 
                                              "< 0.001", round(data_ga$pvalue, digits = pDigit)), 
                     "\n", "Ajusted Pvalue: ", ifelse(data_ga$p.adjust < 
                                                        0.001, "< 0.001", round(data_ga$p.adjust, digits = pDigit)), 
                     "\n", sep = " ")
    px <- pvalX * nrow(gsdata[which(gsdata$id == geneSetID[1]), 
    ])
    py <- pvalY * sum(abs(range(gsdata$runningScore))) + 
      min(gsdata$runningScore)
    if (length(geneSetID) == 1) {
      pLabelOut <- plabel + ggplot2::annotate(geom = "text", 
                                              x = px, y = py, label = pLabel, size = pvalSize, 
                                              color = pCol, fontface = "italic", hjust = pHjust)
    }
    else {
      mytable <- tibble::tibble(x = px, y = py, table = list(tibble::tibble(NES = round(data_ga$NES, 
                                                                                        digits = nesDigit), Pvalue = ifelse(data_ga$pvalue < 
                                                                                                                              0.001, "< 0.001", round(data_ga$pvalue, digits = pDigit)), 
                                                                            `Ajusted Pvalue` = ifelse(data_ga$p.adjust < 
                                                                                                        0.001, "< 0.001", round(data_ga$p.adjust, digits = pDigit)))))
      pLabelOut <- plabel + ggpp::geom_table(data = mytable, 
                                             ggplot2::aes(px, py, label = table))
    }
  }
  else {
    pLabelOut <- plabel
  }
  if (length(geneSetID) == 1) {
    line.col <- ggplot2::scale_color_manual(values = "black")
  }
  if (add.geneExpHt == TRUE) {
    pseg.b = 0
  }
  else {
    pseg.b = 0.2
  }
  pseg <- ggplot2::ggplot(gsdata, ggplot2::aes_(x = ~x, y = ~runningScore, 
                                                color = ~Description)) + ggplot2::geom_segment(data = gsdata1, 
                                                                                               ggplot2::aes_(x = ~x, xend = ~x, y = 0, yend = 1), show.legend = F) + 
    line.col + ggplot2::scale_x_continuous(expand = c(0, 
                                                      0)) + ggplot2::scale_y_continuous(expand = c(0, 0)) + 
    ggplot2::theme_bw(base_size = 14) + ggplot2::theme(axis.ticks = ggplot2::element_blank(), 
                                                       axis.text = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(), 
                                                       panel.grid = ggplot2::element_blank(), axis.line.x = ggplot2::element_blank(), 
                                                       strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_blank(), 
                                                       panel.spacing = ggplot2::unit(0.1, "cm"), plot.margin = ggplot2::margin(t = 0, 
                                                                                                                               r = 0.2, b = pseg.b, l = 0.2, unit = "cm")) + ggplot2::xlab("Rank in Ordered Dataset") + 
    ggplot2::facet_wrap(~Description, ncol = 1)
  if (subPlot > 2) {
    pseg <- pseg + ggplot2::theme(axis.title.x = ggplot2::element_blank())
  }
  else {
    pseg <- pseg
  }
  d <- purrr::map_df(unique(gsdata$Description), function(setid) {
    tmp <- gsdata %>% dplyr::filter(Description == setid)
    v <- seq(1, sum(tmp$position), length.out = 9)
    inv <- findInterval(rev(cumsum(tmp$position)), v)
    if (min(inv) == 0) {
      inv <- inv + 1
    }
    color <- (grDevices::colorRampPalette(c(htCol[1], "white", 
                                            htCol[2])))(10)
    ymin <- 0
    yy <- htHeight
    xmin <- which(!duplicated(inv))
    xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
    d <- data.frame(ymin = ymin, ymax = yy, xmin = xmin, 
                    xmax = xmax, col = color[unique(inv)], Description = setid)
  })
  pseg_ht <- pseg + ggplot2::geom_rect(ggplot2::aes_(xmin = ~xmin, 
                                                     xmax = ~xmax, ymin = ~ymin, ymax = ~ymax, fill = ~I(col)), 
                                       data = d, alpha = 0.8, inherit.aes = FALSE)
  pseg_ht1 <- pseg_ht + ggplot2::xlab("") + ggplot2::theme(axis.title.x = ggplot2::element_blank(), 
                                                           plot.margin = ggplot2::margin(t = -0.1, r = 0.2, b = 0, 
                                                                                         l = 0.2, unit = "cm"))
  if (add.geneExpHt == TRUE) {
    prank.b = 0
  }
  else {
    prank.b = 0.2
  }
  prank <- ggplot2::ggplot(gsdata[which(gsdata$Description == 
                                          unique(gsdata$Description)[1]), ], ggplot2::aes_(x = ~x, 
                                                                                           y = ~geneList)) + ggplot2::geom_col(ggplot2::aes_(fill = ~geneList), 
                                                                                                                               width = 1, color = NA, show.legend = F) + ggplot2::scale_fill_gradient2(low = rankCol[1], 
                                                                                                                                                                                                       mid = rankCol[2], high = rankCol[3], midpoint = 0) + 
    ggplot2::geom_hline(yintercept = 0, size = 0.8, color = "black", 
                        lty = "dashed") + ggplot2::scale_x_continuous(breaks = seq(0, 
                                                                                   nrow(gsdata), rankSeq)) + ggplot2::theme_bw(base_size = 14) + 
    ggplot2::theme(panel.grid = ggplot2::element_blank(), 
                   plot.margin = ggplot2::margin(t = -0.1, r = 0.2, 
                                                 b = prank.b, l = 0.2, unit = "cm")) + ggplot2::coord_cartesian(expand = 0) + 
    ggplot2::ylab("Ranked List") + ggplot2::xlab("REL Rank of Cells")
  if (add.geneExpHt == TRUE) {
    prank <- prank + ggplot2::theme(axis.title.x = ggplot2::element_blank())
  }
  else {
    prank <- prank
  }
  if (kegg == FALSE) {
    rank.g <- data.frame(logfc = object@geneList, gene_name = names(object@geneList)) %>% 
      dplyr::mutate(x = 1:length(object@geneList))
  }
  else {
    rank.g <- data.frame(logfc = object@geneList, gene_name = object@gene2Symbol) %>% 
      dplyr::mutate(x = 1:length(object@geneList))
  }
  if (!is.null(rank.gene)) {
    target.rank.g <- rank.g %>% dplyr::filter(gene_name %in% 
                                                rank.gene) %>% dplyr::mutate(vjust = ifelse(logfc > 
                                                                                              0, "bottom", "top"), nudge_y = ifelse(logfc > 0, 
                                                                                                                                    -rank.gene.nudgey, rank.gene.nudgey))
    prank <- prank + ggrepel::geom_text_repel(data = target.rank.g, 
                                              ggplot2::aes(x = as.numeric(x), y = 0, label = gene_name, 
                                                           vjust = vjust, nudge_y = nudge_y), max.overlaps = 200, 
                                              direction = "x", angle = 90, fontface = "italic", 
                                              size = geneSize)
  }
  d <- purrr::map_df(unique(d$Description), function(x) {
    tmp <- d %>% dplyr::filter(Description == x)
    htcolor <- (grDevices::colorRampPalette(newHtCol))(nrow(tmp))
    tmp <- tmp %>% dplyr::mutate(htcol = htcolor)
  })
  ht <- ggplot2::ggplot(gsdata, ggplot2::aes_(x = ~x, y = ~runningScore)) + 
    ggplot2::geom_rect(ggplot2::aes_(xmin = ~xmin, xmax = ~xmax, 
                                     ymin = ~ymin, ymax = ~ymax, fill = ~I(htcol)), data = d, 
                       alpha = 0.8, inherit.aes = FALSE) + ggplot2::scale_x_continuous(expand = c(0, 
                                                                                                  0)) + ggplot2::scale_y_continuous(expand = c(0, 0)) + 
    ggplot2::theme_bw(base_size = 14) + ggplot2::theme(panel.grid = ggplot2::element_blank(), 
                                                       axis.ticks = ggplot2::element_blank(), axis.text = ggplot2::element_blank(), 
                                                       axis.title = ggplot2::element_blank(), plot.margin = ggplot2::margin(t = 0, 
                                                                                                                            r = 0.2, b = 0.2, l = 0.2, unit = "cm"))
  if (add.geneExpHt == TRUE) {
    target.g <- purrr::map_df(data_ga$ID, function(x) {
      tmp <- data_ga %>% dplyr::filter(ID == x)
      coregene <- unique(unlist(strsplit(tmp$core_enrichment, 
                                         split = "\\/")))
      output <- data.frame(gene_name = coregene, ID = x, 
                           Description = tmp$Description) %>% dplyr::distinct(., 
                                                                              gene_name, .keep_all = TRUE)
    })
    gpos <- if (kegg == TRUE) {
      match(target.g$gene_name, object@gene2Symbol)
    }
    else {
      match(target.g$gene_name, names(object@geneList))
    }
    ginfo <- target.g %>% dplyr::mutate(gpos = gpos) %>% 
      dplyr::arrange(gpos)
    if (scale.exp == TRUE) {
      gexp <- t(scale(t(exp[, 2:ncol(exp)]), scale = TRUE, 
                      center = TRUE)) %>% data.frame()
      gexp$gene_name <- exp[, 1]
    }
    else {
      gexp <- exp
      colnames(gexp)[1] <- "gene_name"
    }
    exp.long <- gexp %>% dplyr::filter(gene_name %in% unique(ginfo$gene_name)) %>% 
      dplyr::left_join(., ginfo[, 1:2], by = "gene_name") %>% 
      reshape2::melt(., id.vars = c("gene_name", "ID"))
    exp.long$gene_name <- factor(exp.long$gene_name, levels = unique(ginfo$gene_name))
    if (!is.null(sample.order)) {
      exp.long$variable <- factor(exp.long$variable, levels = sample.order)
    }
    if (!is.null(termID.order)) {
      exp.long$ID <- factor(exp.long$ID, levels = termID.order)
    }
    ght <- ggplot2::ggplot(exp.long) + ggplot2::geom_tile(ggplot2::aes(x = gene_name, 
                                                                       y = variable, fill = value), color = NA, show.legend = ht.legend) + 
      ggplot2::theme_bw(base_size = 14) + ggplot2::coord_cartesian(expand = 0) + 
      ggplot2::scale_fill_gradient2(low = exp.col[1], mid = exp.col[2], 
                                    high = exp.col[3], midpoint = 0, name = "Z-Score") + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, 
                                                         vjust = 0.5, hjust = 1, size = ght.geneText.size), 
                     axis.text = ggplot2::element_text(color = "black"), 
                     axis.ticks.x = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(), 
                     plot.margin = ggplot2::margin(t = -0.1, r = 0.2, 
                                                   b = 0.2, l = 0.2, unit = "cm")) + ggplot2::scale_y_discrete(position = "right") + 
      ggplot2::xlab("") + ggplot2::ylab("")
    if (ght.facet == TRUE) {
      fght <- ght + ggplot2::facet_wrap(~ID, ncol = 1, 
                                        scales = ght.facet.scale, strip.position = "left") + 
        ggplot2::theme(strip.background = ggplot2::element_rect(color = NA, 
                                                                fill = "grey90"), strip.placement = "outside")
    }
    else {
      fght <- ght
    }
  }
  if (newGsea == FALSE) {
    if (subPlot == 1) {
      pres <- pLabelOut
    }
    else if (subPlot == 2) {
      if (rmHt == FALSE) {
        pres <- aplot::plot_list(gglist = list(pLabelOut, 
                                               pseg_ht), ncol = 1, heights = c(0.8, 0.2))
      }
      else {
        pres <- aplot::plot_list(gglist = list(pLabelOut, 
                                               pseg), ncol = 1, heights = c(0.8, 0.2))
      }
    }
    else if (subPlot == 3) {
      if (rmHt == FALSE) {
        pres <- aplot::plot_list(gglist = list(pLabelOut, 
                                               pseg_ht1, prank), ncol = 1, heights = c(0.5, 
                                                                                       0.2, 0.3))
      }
      else {
        pres <- aplot::plot_list(gglist = list(pLabelOut, 
                                               pseg, prank), ncol = 1, heights = c(0.5, 0.2, 
                                                                                   0.3))
      }
    }
    else {
      message("Please give 1/2/3 parameters!")
    }
  }
  else {
    if (rmHt == FALSE) {
      pres <- aplot::plot_list(gglist = list(pLabelOut, 
                                             ht), ncol = 1, heights = c(0.9, 0.1))
    }
    else {
      pres <- pLabelOut
    }
  }
  if (add.geneExpHt == TRUE) {
    pfinal <- aplot::plot_list(gglist = list(pres + ggplot2::xlab(""), 
                                             fght), ncol = 1, heights = c(1 - ght.relHight, ght.relHight))
  }
  else {
    pfinal <- pres
  }
  return(pfinal)
}



gsea_plot <- gseaNb_plus(object = x,
       geneSetID = name_pheno,
       subPlot = 2,
       termWidth = 35,
       legend.position = "right",
       addPval = F,
       pvalX = 0.05,pvalY = 0.05,
       curveCol = c(brewer.pal(5,"Set1")), htCol = c("#08519C", "#A50F15"), rankCol = c("#08519C","white", "#A50F15"),
       )
gsea_plot

ggsave("gsea_plot_3.pdf", height = 8, width = 8)





### C-index
desired_order <- c("PEAP", "LIMMA", "edgeR", "DESeq2","PAM50")
# Reshape the data from wide to long format
long_df <- C_index %>%
  pivot_longer(cols = -Survival_time, names_to = "Method", values_to = "C_index")
long_df$Method <- factor(long_df$Method, levels = desired_order)
# Create the bar plot
ggplot(long_df, aes(x = factor(Survival_time,
                               labels = c("1-year survival",
                                          "2-year survival",
                                          "3-year survival",
                                          "4-year survival",
                                          "5-year survival")), y = C_index, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x="",y = "C-index", title = "C-index") +
  theme_minimal() +
  geom_text(aes(label = sprintf("%.2f", C_index)), # Format the label to 2 decimal places
            position = position_dodge(width = 0.9), vjust = -0.25, # Adjust vertical position
            size =5, # Adjust text size as needed
            color = "black") +
  scale_y_continuous(breaks = seq(0.1, 0.8, by = 0.1), limits = c(0, 0.8))+
  scale_fill_manual(values=c(rgb(248,230,32,maxColorValue = 255),
                             rgb(145,213,66,maxColorValue = 255),
                             rgb(53,183,119,maxColorValue = 255),
                             rgb(38,104,141,maxColorValue = 255),
                             rgb(68,4,90, maxColorValue =255)))+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))

### KM_plot
plotdata=surv_peap2034
risk_scores <- predict(peap_fit, as.matrix(plotdata[,peapP]),type="link")
risk_threshold <- median(risk_scores)
plotdata$risk_group <- ifelse(risk_scores > risk_threshold, "High risk", "Low risk")
surv_obj <- Surv(plotdata$surv_time, plotdata$vital_num)
km_fit <- survfit(surv_obj ~ risk_group, data = plotdata)
surv_plot= ggsurvplot(
  km_fit,
  data = plotdata,
  size = 1,                 # change line size
  palette =
    c("#E7B800", "#2E9FDF"),# custom color palettes
  conf.int = TRUE,          # Add confidence interval            # Add p-value
  legend.labs =
    c("High Risk", "Low Risk"),    # Change legend labels
  ggtheme = theme(axis.text.y = element_text(size = 15),
                  axis.text.x = element_text(size = 15),
                  legend.direction = "vertical"),      # Change ggplot2 theme
  ylim = c(0.5, 1),
  title="All sample"
)
log_rank_test = survdiff(surv_obj ~ risk_group, data = plotdata)

surv_plot$plot <- surv_plot$plot + 
  annotate("text", x = 15, y = 0.65, label = paste("P-value =", format.pval(log_rank_test$pvalue)), 
           size = 5, color = "black")
print(surv_plot)
plotdata=surv_peap2034[surv_peap2034$subtype=="Basal",]
risk_scores <- predict(peap_fit, as.matrix(plotdata[,peapP]),type="link")
risk_threshold <- median(risk_scores)
plotdata$risk_group <- ifelse(risk_scores > risk_threshold, "High risk", "Low risk")
surv_obj <- Surv(plotdata$surv_time, plotdata$vital_num)
km_fit <- survfit(surv_obj ~ risk_group, data = plotdata)
surv_plot= ggsurvplot(
  km_fit,
  data = plotdata,
  size = 1,                 # change line size
  palette =
    c("#E7B800", "#2E9FDF"),# custom color palettes
  conf.int = TRUE,          # Add confidence interval            # Add p-value
  legend.labs =
    c("High Risk", "Low Risk"),    # Change legend labels
  ggtheme = theme(axis.text.y = element_text(size = 15),
                  axis.text.x = element_text(size = 15),
                  legend.direction = "vertical"),      # Change ggplot2 theme
  ylim = c(0.5, 1),
  title="All sample"
)
log_rank_test = survdiff(surv_obj ~ risk_group, data = plotdata)

surv_plot$plot <- surv_plot$plot + 
  annotate("text", x = 15, y = 0.65, label = paste("P-value =", format.pval(log_rank_test$pvalue)), 
           size = 5, color = "black")
print(surv_plot)

## Supple table 1
cox = as.vector(peap_fit$beta)
names(cox) = as.vector(row.names(peap_fit$beta))
cox = cox[cox!=0]
peap_selected = data.frame(
  pair_name=names(cox)
)
ESres_noduplicated$pair_name = paste0(pmin(ESres_noduplicated$gene1,
                                           ESres_noduplicated$gene2),
                                      "_",
                                      pmax(ESres_noduplicated$gene1,
                                           ESres_noduplicated$gene2))
ESres_noduplicated_v2 = merge(ESres_noduplicated,peap_selected,by="pair_name")
covariates <- sapply(ESres_noduplicated_v2[,1], function(x) {
  if (grepl("-", x, fixed = TRUE)) {
    return(paste0("`", x, "`"))
  } else {
    return(x)
  }
})
ESres_noduplicated_v2[,1]=covariates
hr = hr[,c("term","estimate","upper","lower")]
names(hr)=c("term",
            "Hazard_Ratio_Estimate",
            "Hazard_Ratio_upper",
            "Hazard_Ratio_lower")
ESres_noduplicated_v2=merge(ESres_noduplicated_v2,hr,
                            by.x ="pair_name",
                            by.y = "term")
ESres_noduplicated_v2$pair_name = NULL
ESres_noduplicated_v2$term = NULL
write.csv(ESres_noduplicated_v2,"ESres_marker_BRCA_v2.csv",row.names = F)
peap_final_pairs=as.data.frame(ESres_noduplicated_v2[,1:2])
save(peap_final_pairs,file="peap_final_pairs.RDa")




