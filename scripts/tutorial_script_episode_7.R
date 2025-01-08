library(SummarizedExperiment)
library(DESeq2)
library(gplots)
library(microbenchmark)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(simplifyEnrichment)
library(tidyverse)

# read the example dataset which is a `RangedSummarizedExperiment` object
se <- readRDS("data/GSE96870_se.rds")

# only restrict to mRNA (protein-coding genes)
se <- se[rowData(se)$gbkey == "mRNA"]

#converting assay(se) from a data frame into a matrix
assay(se) <- as.matrix(assay(se)) 

# construct a 'DESeqDataSet' object where we also specify the experimental design
dds <- DESeqDataSet(se, design = ~ sex + time)

# perform DESeq2 analysis
dds <- DESeq(dds)

# obtain DESeq2 results, here we only want Male vs Female in the "sex" variable.
#log fold change in gene expression for males relative to females at Day0.
#Since there are no interaction terms between sex and time, the model assumes that the 
#difference between males and females (as captured by sex_Male_vs_Female) is constant across all time points.
resSex <- results(dds, contrast = c("sex", "Male", "Female"))

# extract DE genes with padj < 0.05
sexDE <- as.data.frame(subset(resSex, padj < 0.05))

# the list of DE genes
sexDEgenes <- rownames(sexDE)

# checking the number of DE genes
length(sexDEgenes)

#Getting a list of genes from sex chromosomes
geneGR <- rowRanges(se)
totalGenes <- rownames(se)
XYGeneSet <- totalGenes[as.vector(seqnames(geneGR)) %in% c("X", "Y")]

# plot a venn diagram of the gene list from the DE analysis and the XY gene list
plot(venn(list("sexDEgenes" = sexDEgenes, "XY gene set" = XYGeneSet)))
title(paste0("|universal| = ",length(totalGenes)))
