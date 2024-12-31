library(SummarizedExperiment)
library(DESeq2)
library(vsn)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(hexbin)
library(iSEE)

#Loading in the SummarizedExperiment object created last episode
se <- readRDS("data/GSE96870_se.rds")

#Looking at the number of genes in se
nrow(se)
#Looking at the number of samples in se
ncol(se)

# Removing genes/rows that don't have > 5 total counts across samples
se <- se[rowSums(assay(se)) > 5,]

#Seeing how many genes would be present if the threshold is changed
length(which(rowSums(assay(se)) > 20))

#Looking at the number of each type of gene that survived filtering
table(rowData(se)$gbkey)

#Creating a new column of the sum of all counts across the genes for each sample
se$LibSize <- colSums(assay(se))

# Plot the libSize by using R's native pipe |>
# to extract the colData, turn it into a regular
# data frame then send to ggplot:
colData(se) |>
  as.data.frame() |>
  ggplot(aes(x = Label, y = LibSize/1e6, fill = Group)) +
    geom_col() +
    theme_bw() +
    labs(x = "Sample", y = "Total counts in millions") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
