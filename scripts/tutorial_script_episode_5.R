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

#Creating a DESeqDataSet object from the SummarizedExperiment object with a design matrix of sez and time. 
#design -> determining how the counts for each gene depend on the sex and time of the samples
dds <- DESeqDataSet(se, design = ~ sex + time)

#Estimating the size factor for each sample to normalize the data
dds <- estimateSizeFactors(dds)

#Estimating dispersion for each gene
dds <- estimateDispersions(dds)

plotDispEsts(dds)

#Fitting a GLM model for each gene to find log2 fold change and then testing for significane using Wald test
dds <- nbinomWaldTest(dds)

##DAY 8 vs Day 0
resTime <- results(dds, contrast = c("time", "Day8", "Day0"))
summary(resTime)

#Looking at the genes with the lowest p-values
head(resTime[order(resTime$pvalue),])


##Male vs Female
resSex <- results(dds, contrast = c("sex", "Male", "Female"))
summary(resSex)
#Looking at the genes with the lowest p-values
head(resSex[order(resSex$pvalue),])


#plotting log2 fold change and gene mean count
plotMA(resTime)

#Shrinking the log2 fold change of genes with low mean count
resTimelfc <- lfcShrink(dds, coef = "time_Day8_vs_Day0", res = resTime)

plotMA(resTimelfc)


# Transform counts
vsd <- vst(dds, blind = TRUE)

# Get top DE genes
genes <- resTime[order(resTime$pvalue), ] |>
  head(10) |>
  rownames()
heatmapData <- assay(vsd)[genes, ]

# Scale counts for visualization
heatmapData <- t(scale(t(heatmapData)))

# Add annotation
heatmapColAnnot <- data.frame(colData(vsd)[, c("time", "sex")])
heatmapColAnnot <- HeatmapAnnotation(df = heatmapColAnnot)


# Plot as heatmap
ComplexHeatmap::Heatmap(heatmapData,
                        top_annotation = heatmapColAnnot,
                        cluster_rows = TRUE, cluster_columns = FALSE)

#output results
head(as.data.frame(resTime))
head(as.data.frame(rowRanges(se)))

temp <- cbind(as.data.frame(rowRanges(se)), as.data.frame(resTime))

write.csv(temp, file = "output/Day8vsDay0.csv")