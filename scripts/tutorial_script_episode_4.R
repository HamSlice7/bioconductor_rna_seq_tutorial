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

#Converting the counts from se to a matrix
counts_matrix <- as.matrix(assay(se))
assay(se) <- counts_matrix

#Convert the SummarizedExperiment object to a DEseqDataSet
dds <- DESeqDataSet(se, design = ~ sex + time)

#Getting the size factors for each sample
dds <- estimateSizeFactors(dds)

# Plot the size factors against library size
# and look for any patterns by group:
ggplot(data.frame(libSize = colSums(assay(dds)),
                  sizeFactor = sizeFactors(dds),
                  Group = dds$Group),
                  aes(x = libSize, y = sizeFactor, col = Group)) +
                    geom_point(size = 5) +
                    theme_bw() +
                    labs(x = "Library Size", y = "Size Factor")

#Plotting row standard deviation vs row means (row = counts of a gene across the samples)
meanSdPlot(assay(dds), ranks = FALSE)

#Applying variance stabilizing transformation to the count data
vsd <- vst(dds, blind = TRUE)

#Plotting row standard deviation vs row means for the transformed count data
meanSdPlot(assay(vsd), blind = TRUE)

#Creating a distance matrix from the counts matrix. Calculates Euclidean distance between the rows.
dst <- dist(t(assay(vsd)))

colours <- colorRampPalette(brewer.pal(9, "Blues"))(225)

Heatmap(
  as.matrix(dst),
  col = colours,
  name = "Euclidean\ndistance",
  cluster_rows = hclust(dst),
  cluster_columns = hclust(dst),
  bottom_annotation = columnAnnotation(
    sex = vsd$sex,
    time = vsd$time,
    col = list(sex = c(Female = "red", Male = "blue"),
               time = c(Day0 = "yellow", Day4 = "forestgreen", Day8 = "purple"))
  ))

#Creating a PCA plot
pcaData <- plotPCA(vsd, intgroup = c("sex", "time"),returnData = TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = sex, shape = time), size = 5) +
  theme_minimal() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  scale_color_manual(values = c(Male = "blue", Female = "red"))
