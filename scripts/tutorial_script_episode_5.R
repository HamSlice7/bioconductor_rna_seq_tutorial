suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(DESeq2)
  library(ggplot2)
  library(ExploreModelMatrix)
  library(cowplot)
  library(ComplexHeatmap)
  library(apeglm)
})

#Importing the SummarizedExperiment object
se <- readRDS("data/GSE96870_se.rds")

#Filtering out genes with a count across the samples of less than 5
se <- se[rowSums(assay(se)) > 5, ]

#Converting the counts from se to a matrix
counts_matrix <- as.matrix(assay(se))
assay(se) <- counts_matrix

#Creating a DESeqDataSet object from the SummarizedExperiment object with a design matrix of sez and time. design -> determining how the counts for each gene depend on the sex and time of the samples
dds <- DESeqDataSet(se, design = ~ sex + time)
