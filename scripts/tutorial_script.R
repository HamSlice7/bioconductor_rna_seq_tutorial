install.packages(c("BiocManager", "remotes"))
BiocManager::install(c("tidyverse", "SummarizedExperiment",
                       "ExploreModelMatrix", "AnnotationDbi", "org.Hs.eg.db", 
                       "org.Mm.eg.db", "csoneson/ConfoundingExplorer",
                       "DESeq2", "vsn", "ComplexHeatmap", "hgu95av2.db",
                       "RColorBrewer", "hexbin", "cowplot", "iSEE",
                       "clusterProfiler", "enrichplot", "kableExtra",
                       "msigdbr", "gplots", "ggplot2", "simplifyEnrichment",
                       "apeglm", "microbenchmark", "Biostrings",
                       "SingleCellExperiment"))

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(hgu95av2.db)
  library(SummarizedExperiment)
})

#Reading in the expression (count) matrix
counts <- read.csv("data/GSE96870_counts_cerebellum.csv", 
                   row.names = 1)
dim(counts)

#Reading in the sample annotation data
coldata <- read.csv("data/GSE96870_coldata_cerebellum.csv", row.names = 1)

dim(coldata)

#Reading in the gene annotation data
rowranges <- read.delim("data/GSE96870_rowranges.tsv",
                        sep = "\t",
                        colClasses = c(ENTREZID = "character"),
                        header = TRUE,
                        quote = "",
                        row.names = 5)
dim(rowranges)

#Checking the types of gene and the number of each
table(rowranges$gbkey)

#Making sure the sample ID's are in the same order for the counts data and the sample annotation data
all.equal(colnames(counts), rownames(coldata))

#Making sure the genes are in the same order for the counts data and the gene annotation data
all.equal(rownames(counts), rownames(rowranges))

# If the first is not TRUE, you can match up the samples/columns in
# counts with the samples/rows in coldata like this (which is fine
# to run even if the first was TRUE):
tempindx <- match(colnames(counts), rownames(coldata))
coldata <- coldata[tempindx,]

#Check again:
all.equal(colnames(counts), rownames(coldata))

#One final check:
stopifnot(colnames(count) == rownames(coldata),
          rownames(counts) == rownames(rowranges))

#Creating the SummarizedExperiment object
se <- SummarizedExperiment(assays = list(counts),
                           rowRanges = as(rowranges, "GRanges"),
                           colData = coldata)
