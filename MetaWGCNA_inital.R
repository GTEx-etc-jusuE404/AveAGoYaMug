#### this is a script to attempt a meta-analysis of 
#### 2 seperate datasets using WGCNA
#### Datasets are TCGA & Metabric (no lncs)


setwd("~/Bioinformatics Work/Meth & RNA/Meta-analysis WGCNA")

library(WGCNA)
library(impute)
library(dynamicTreeCut)
library(qvalue)
library(flashClust)
library(Hmisc)


#### load in the exp files
### ROWS = SAMPLES, COLUMNS = GENES

exp1 <- read.table("genomicMatrix", sep = "\t", header = TRUE, row.names = 1) #TCGA
  exp1 <- as.data.frame(t(exp1))
  
exp2 <- read.table("Exp_final_METABRIC.txt", sep = "\t", header = TRUE, row.names = 1) #METABRIC

#### preprocess the data using goodgenes (WCGNA)

gsg = goodSamplesGenes(exp1,verbose = 5);
gsg$allOK


## if return is true, all good to go
### Otherwise

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(exp1)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(exp1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  exp1 = exp1[gsg$goodSamples, gsg$goodGenes]
}


## Create a dendro of the samples to look for outliers

sampleTree = flashClust(dist(exp1), method = "average");


# Plot the sample tree

sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, 
     main = "Sample clustering to detect outliers", 
     sub="", xlab="", 
     cex.lab = 1.5,
     cex.axis = 1.5, 
     cex.main = 2, labels = FALSE)

## Might remove the 4 outliers on the right
## Plot a line to show the cut

abline(h = 270, col = "red");

# Determine clusters under the line (use h from above)

clust = cutreeStatic(sampleTree, cutHeight = 270, minSize = 10)
table(clust)

# clust 1 contains the samples we want to keep.

keepSamples1 = (clust==1)
datExpr1 = exp1[keepSamples1, ]
nGenes1 = ncol(datExpr1)
nSamples1 = nrow(datExpr1)

############################# Again for the metabric set


gsg = goodSamplesGenes(exp2,verbose = 5);
gsg$allOK


## if return is true, all good to go
### Otherwise

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(exp2)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(exp2)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  exp2 = exp2[gsg$goodSamples, gsg$goodGenes]
}


## Create a dendro of the samples to look for outliers

sampleTree = flashClust(dist(exp2), method = "average");


# Plot the sample tree

sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, 
     main = "Sample clustering to detect outliers", 
     sub="", xlab="", 
     cex.lab = 1.5,
     cex.axis = 1.5, 
     cex.main = 2, labels = FALSE)

## Might remove the 4 outliers on the right
## Plot a line to show the cut

abline(h = 110, col = "red");

# Determine clusters under the line (use h from above)

clust = cutreeStatic(sampleTree, cutHeight = 110, minSize = 10)
table(clust)

# clust 1 contains the samples we want to keep.

keepSamples = (clust==1)
datExpr2 = exp2[keepSamples, ]
nGenes2 = ncol(datExpr2)
nSamples2 = nrow(datExpr2)



##### NOW TRANSPOSE  
### ### ROWS = genes, COLUMNS = samples

datExpr1 <- as.data.frame(t(datExpr1))
datExpr2 <- as.data.frame(t(datExpr2))

#### only common genes can be included

commongenes <- intersect(rownames(datExpr1), rownames(datExpr2))
datExpr1 <- datExpr1[commongenes,]
datExpr2 <- datExpr2[commongenes,]

###
save(datExpr1, datExpr2, file = "MetaAnalysis_trimmed_input.RData")

#### correlating general network properties
### assessing the comparability of the two data sets by correlating measures of 
## average gene expression and overall connectivity between the data sets

## determine softpower for the sets

