# Building gene pair signiture for the classification of infectious diseases

## Preprocessing
    
    Conduct log2 transformation on the gene expression levels and extract with-pathway gene pairs in terms of 186 KEGG pathways.
    
    - Input: Gene expression matrices and metadata downloaded from the GEO platform
    - Output: Relative expression level matrix for within-pathway genes and samples from three different phenotypes, including healthy, viral-infected and bacterial-infected samples

## Phenotype Enrichment Analysis for Gene Pairs

    Apply weighted KS testing to relative expression levels calculates the phenotype-specific enrichment score for gene pairs, which will be used to rank the gene pairs. Parallel computing is applied to compute the phenotype-enrichment scores and p-values for different pathways in parallel.

    - Output: Matrices of phenotype-enrichment scores and p-values for all within-pathway gene pairs and all phenotypes

## Extracting Top Gene Pairs

    Integrate the phenotype-enrichment scores and p-values of 186 KEGG pathways and extract the top non-overlapping gene pairs with the smallest p-values for the downstream analysis. Association analysis between phenotype and pathways is also conducted.
    
    - Output: Top non-overlapping gene pairs matrix and results of association analysis

## Building Gene Pair Signature for Phenotype Classification and Benchmarking

    Benchmark PEAP with three popular DE detection methods. In particular, we conduct a two-step approach to build the gene (pair) signitures for all the four competing methods. Evaluate the classification performance in terms of clustering accuracy, AUC, sensitivity and specificity measures.
   
    - Output: Classification performance for four different methods.

## Summarizing Phenotype-enrichment Scores and P-values for Selected Gene Pairs
    
    Summarize the phenotype-enrichment scores and p-values of the built gene pair signature as Table S1.

    - Output: Table of PE scores and p-values

## Visualization

    Visualize association analysis, running sum of test statistics and classification performance as Figure 3. 




