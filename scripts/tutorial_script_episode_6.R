suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(ExploreModelMatrix)
  library(dplyr)
  library(DESeq2)
})

meta <- read.csv("data/GSE96870_coldata_all.csv", row.names = 1)

#All mice on day0 are not infected
table(meta$time, meta$infection)

#All mice are 8 weeks
table(meta$age)

#visualizing the number of observations for each combination of the three predictor variables
vd <- VisualizeDesign(sampleData = meta, designFormula = ~  tissue + time + sex)

vd$cooccurrenceplots


