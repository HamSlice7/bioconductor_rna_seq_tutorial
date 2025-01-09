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


##Fisher's exact test

#Constructing 2x2 contingency table
n <- nrow(se) #total number of genes
n_01 <- length(XYGeneSet) #size of the XY gene set
n_10 <- length(sexDEgenes) #total number of DE genes
n_11 <- length(intersect(sexDEgenes, XYGeneSet)) #number of DE genes in the XY gene set
n_12 <- n_10 - n_11 #total number of genes - number of DE genes in the XY gene set = number of DE genes not in the XY gene set
n_21 <- n_01 - n_11 #size of the XY gene set - number of DE genes in the XY gene set = number of genes in XY gene set not DE
n_20 <- n - n_10 # total number of genes - total number of DE genes = total number of non DE genes
n_02 <- n - n_01 # total number of gens - size of the XY gene set = total number of genes not in the XY gene set
n_22 <- n_02 - n_12 # total number of genes not in the XY gene set - number of DE genes not in the XY gene set= not total numbr of genes not in XY gene set that are not DE

matrix(c(n_11, n_12, n_10, n_21, n_22, n_20, n_01, n_02, n), nrow=3, byrow = TRUE)

#performing Fisher's exact test
t <- fisher.test(matrix(c(n_11, n_12, n_21, n_22), nrow = 2, byrow = TRUE), alternative = "greater")
t$p.value

#Using the hypergeometric distribution to calculate p-value of the enrichment 
#(probability of picking DE genes equal to or larger than n11 (13) under the 
#assumption of independence) from the hypergeometric distribution from n_01 draws
#(number of genes in XY gene set)
1 - phyper(n_11 -1, n_10, n_20, n_01)

#Switching n_01 (number of genes in XY data set) with n_10 (number of DE genes in total). Draw 54 (n_10), 
#probability of picking gene in XY gene set equal to or larger than observation n11 (13)
1 - phyper(n_11 -1, n_01, n_02, n_10)


#comparing speed of fisher.test() and phyper()
library(microbenchmark)
microbenchmark(
  fisher = fisher.test(matrix(c(n_11, n_12, n_21, n_22), nrow = 2, byrow = TRUE),
                       alternative = "greater"),
  hyper = 1 - phyper(n_11 - 1, n_10, n_20, n_01)
)

#List of of three gene sets as vectors
lt <- list(gene_set_1 = c("gene_1", "gene_2", "gene_3"),
           gene_set_2 = c("gene_1", "gene_3", "gene_4", "gene_5", "gene_6"),
           gene_set_3 = c("gene_4", "gene_7")
)

#converting lt to a data frame
df = data.frame(gene_set = rep(names(lt), times = sapply(lt, length)), gene = unname(unlist(lt)))

#converting df back to the list
split(df$gene, df$gene_set)


##ORA with clusterProfiler
resTime <- results(dds, contrast = c("time", "Day8", "Day0"))
summary(resTime)

#Filtering genes in resTime and then converting to a data frame
timeDE <- as.data.frame(subset(resTime, padj < 0.05 & abs(log2FoldChange) > log2(1.5)))

#Getting a list of filtered DE genes
timeDEgenes <- rownames(timeDE)

#Performing ORA on GO gene set
resTimeGO <- enrichGO(gene = timeDEgenes, keyType = "SYMBOL", ont = "BP", OrgDb = org.Mm.eg.db)

resTimeGoTable <- as.data.frame(resTimeGO)
head(resTimeGoTable)


#KEGG pathway enrichment
EntrezIDs <- mapIds(org.Mm.eg.db, keys = timeDEgenes, keytype = "SYMBOL", column = "ENTREZID", multiVals = "first")
EntrezIDs <- EntrezIDs[!is.na(EntrezIDs)]
head(EntrezIDs)

resTimeKEGG <- enrichKEGG(gene = EntrezIDs, organism = "mmu", pvalueCutoff = 1, qvalueCutoff = 1)
resTimeKEGGTable <- as.data.frame(resTimeKEGG)
head(resTimeKEGGTable)

#MSigDB enrichment
gene_sets = msigdbr(category = "H", species = "mouse")
head(gene_sets)

resTimeHallmark <- enricher(gene = timeDEgenes, TERM2GENE = gene_sets[,c("gs_name", "gene_symbol")], pvalueCutoff = 1, qvalueCutoff = 1)
resTimeHallmarkTable <- as.data.frame(resTimeHallmark)
head(resTimeHallmarkTable)
