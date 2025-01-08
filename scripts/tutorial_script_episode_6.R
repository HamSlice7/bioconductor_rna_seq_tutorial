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

#testing for differences between tissues after accounting for baseline differences between mice for the female non-infected mice
meta_fem_day0 <- meta %>% 
  filter(sex == "Female" & time == "Day0")

#ensure that mouse is treated as a categorical variable 
meta_fem_day0$mouse <- factor(meta_fem_day0$mouse)

vd <- VisualizeDesign(sampleData = meta_fem_day0, designFormula = ~ mouse + tissue)

vd$designmatrix

vd$plotlist

#Subsetting for female mice and time (day0 and day4 mice)
meta_fem_day_04 <- meta %>% 
  filter(sex == "Female" & time %in% c("Day0", "Day4")) %>% 
  droplevels()

# ensure that mouse is treated as a categorical variable
meta_fem_day_04$mouse <- factor(meta_fem_day_04$mouse)

#Creating a design matrix with mouse number base line and the features
design <- model.matrix(~ mouse, data = meta_fem_day_04)

#Logical vectors are converted to numeric (TRUE = 1, FALSE = 0) when added to a matrix.the = operator is used within the cbind() function to name the new columns being added to the matrix. Adding spinalcord day 0 and day 4 features to the design matrix
design <- cbind(design, Spc.Day0 = meta_fem_day_04$tissue == "Spinalcord" & meta_fem_day_04$time == "Day0", Spc.Day4 = meta_fem_day_04$tissue == "Spinalcord" & meta_fem_day_04$time == "Day4")

#Validating the experimental design for comparing the two tissues in each infection group (Day0 and Day4), and contrast the tissue differences between the infection groups
vd <- VisualizeDesign(sampleData = meta_fem_day_04 %>% 
                        select(time, tissue, mouse), designMatrix = design)
vd$designmatrix
vd$plotlist

#Repeating DE analysis from last tutorial
se <- readRDS("data/GSE96870_se.rds")
se <- se[rowSums(assay(se)) > 5, ]
assay(se) <- as.matrix(assay(se))
dds <- DESeqDataSet(se, design = ~ sex + time)
dds <- DESeq(dds)

#Listing the attributes of the dds object
attributes(dds)

#Accessing the model matrix of the DE analysis
attr(dds, "modelMatrix")

#accessing column names (independent variables) of the design matrix
resultsNames(dds)

#Visualizing the experimental design
vd <- VisualizeDesign(sampleData = colData(dds)[,c("sex", "time")], designMatrix = attr(dds, "modelMatrix"), flipCoordFitted = TRUE)
vd$plotlist

#Extracting DE analysis for day0 vs day8
resultsNames(dds) #coefficients from the design matrix
resTime <- results(dds, contrast = c("time", "Day8", "Day0"))
resTimeNum <- results(dds, contrast = c(0,0,0,1)) #alternate way of specifying day0 vs day 8 (log chnage represented by fourth coefficient)

#looking at results and seeing that resTime and resTimeNum are the same
summary(resTime)
summary(resTimeNum)

##LogFc
plot(resTime$log2FoldChange, resTimeNum$log2FoldChange)
abline

##-log10(p-value)
plot(-log10(resTime$pvalue), -log10(resTimeNum$pvalue))
abline(0,1)

#Redoing the analysis but allow an interaction between sex and time
se <- readRDS("data/GSE96870_se.rds")
se <- se[rowSums(assay(se)) > 5, ]
assay(se) <- as.matrix(assay(se))
dds <- DESeqDataSet(se, design = ~ sex * time)
dds <- DESeq(dds)

attr(dds, "modelMatrix")

#looking at the experimental design
vd <- VisualizeDesign(sampleData = colData(dds)[,c("sex", "time")], designMatrix = attr(dds, "modelMatrix"), flipCoordFitted = TRUE)
vd$plotlist

## Day8 vs Day0, female
resTimeFemale <- results(dds, contrast = c("time", "Day8", "Day0"))
summary(resTimeFemale)

## Interaction effect (difference in Day8-Day0 effect between Male and Female)
resTimeInt <- results(dds, name = "sexMale.timeDay8")
summary(resTimeInt)


###Let’s try to fit this model with the second approach mentioned above, namely to create a single factor.
se <- readRDS("data/GSE96870_se.rds")
se <- se[rowSums(assay(se)) > 5, ]
se$sex_time <- paste0(se$sex, "_", se$time) #creating a sex+time column for the sample data
assay(se) <- as.matrix(assay(se))
dds <- DESeqDataSet(se, design = ~ 0 + sex_time)
dds <- DESeq(dds)

#looking at the model matrix to look st the independent variables and the samples
attr(dds, "modelMatrix")

#looking at the experimental design
vd <- VisualizeDesign(sampleData = colData(dds)[,c("sex", "time")], designMatrix = attr(dds, "modelMatrix"), flipCoordFitted = TRUE)
vd$plotlist

##Day8 vs Day0, female
resTimeFemaleSingle <- results(dds, contrast = c("sex_time", "Female_Day8", "Female_Day0"))

##Interaction effect(difference in Day8-Day0 effect between Male and Female)
resultsNames(dds)
resTimeIntSingle <- results(dds, contrast = c(1,0,-1,-1,0,1))

#Same results - female day 0 vs day 8
summary(resTimeFemale)
summary(resTimeFemaleSingle)

#Same results - test for the difference in the Day8-Day0 effect between Male and Female
summary(resTimeInt)
summary(resTimeIntSingle)
