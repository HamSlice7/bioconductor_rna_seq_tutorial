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

vd$designmatrix

#Subsetting data to only include non-infected (Day 0) and spinal cord tissue 
meta_noninf_spc <- meta %>%
  filter(time == "Day0" & tissue == "Spinalcord")

## Use ExploreModelMatrix to create a design matrix and visualizations, given 
## the desired design formula
vd <- VisualizeDesign(sampleData = meta_noninf_spc, designFormula = ~ sex)

vd$designmatrix

vd$plotlist

#We can also generate a model matrix like this
model.matrix(~sex, data = meta_noninf_spc)

#Subsetting data to include all day0 rows
meta_noninf <- meta %>%
  filter(time == "Day0")


vd <- VisualizeDesign(sampleData = meta_noninf, designFormula = ~ sex + tissue)

vd$designmatrix

vd$plotlist

## Define a design including an interaction term
## Note that ~ sex * tissue is equivalent to 
## ~ sex + tissue + sex:tissue
vd <- VisualizeDesign(sampleData = meta_noninf, designFormula = ~ sex * tissue)

vd$designmatrix

vd$plotlist


#Combining sex and tissue variables
meta_noninf <- meta %>% 
  filter(time == "Day0")
meta_noninf$sex_tissue <- paste0(meta_noninf$sex, "_", meta_noninf$tissue)
meta_noninf

vd <- VisualizeDesign(sampleData = meta_noninf, designFormula = ~ 0 + sex_tissue)
vd$designmatrix
vd$plotlist
