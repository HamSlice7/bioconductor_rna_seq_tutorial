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
