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
library(blockmodeling)

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

############################################################################################
######  ***** READ IN THE DATA FILE ****** #########################

load(file = "MetaAnalysis_trimmed_input.RData")

#### correlating general network properties
### assessing the comparability of the two data sets by correlating measures of 
## average gene expression and overall connectivity between the data sets

## determine softpower for the sets
################

####### Set up a multiexpr file so that you can do soft-threshold in one graph



# We work with two sets:
nSets = 2;
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
  setLabels = c("TCGA", "METABRIC")

# Form multi-set expression data: columns starting from 9 contain actual expression data.
  multiExpr = vector(mode = "list", length = nSets)
  multiExpr[[1]] = list(data = as.data.frame(t(datExpr1)));
    names(multiExpr[[1]]$data) = rownames(datExpr1);
    rownames(multiExpr[[1]]$data) = names(datExpr1);
  multiExpr[[2]] = list(data = as.data.frame(t(datExpr2)));
    names(multiExpr[[2]]$data) = rownames(datExpr2);
    rownames(multiExpr[[2]]$data) = names(datExpr2);
# Check that the data has the correct format for many functions operating on multiple sets:
  exprSize = checkSets(multiExpr)

############
# Choose a set of soft-thresholding powers
  powers = c(seq(4,10,by=1), seq(12,20, by=2));
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                     verbose = 2, networkType = "signed")[[2]]);
collectGarbage();
# Plot the results:
colors = c("black", "red")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
  for (col in 1:length(plotCols))
    
  {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
  }
}
# Plot the quantities in the chosen columns vs. the soft thresholding power
sizeGrWindow(8, 6)
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
  if (set==1)
  {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]);
    addGrid();
  }
  if (col==1)
  {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set]);
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set]);
  if (col==1)
  {
    legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
  } else
    legend("topright", legend = setLabels, col = colors, pch = 20) ;
}



#####  Well, it's not great, but let's try with 9 , 10, 12

commongenes <- intersect(rownames(datExpr1), rownames(datExpr2))

###  
softpower = 9
rankExpr1 <- rank(rowMeans(datExpr1))
rankExpr2 <- rank(rowMeans(datExpr2))
random5000 <- sample(commongenes, 5000)
rankConn1 <- rank(softConnectivity
                  (t(datExpr1[random5000,]), type = "signed", power = softpower))

rankConn2 <- rank(softConnectivity
                  (t(datExpr2[random5000,]), type = "signed", power = softpower))


pdf("generalNetworkProperties1.pdf", height=10, width=9)
par(mfrow=c(2,2))
verboseScatterplot(rankExpr1,rankExpr2, xlab="Ranked Expression (A1)", 
                   ylab="Ranked Expression (A2)")
verboseScatterplot(rankConn1,rankConn2, xlab="Ranked Connectivity (A1)", 
                   ylab="Ranked Connectivity (A2)")

dev.off()

##############  NOW TO RUN THE WGCNA
softPower = 9

adjacency1 = adjacency(t(datExpr1),power=softPower,type="signed");
diag(adjacency1)=0
dissTOM1   = 1-TOMsimilarity(adjacency1, TOMType="signed")
geneTree1  = flashClust(as.dist(dissTOM1), method="average")

adjacency2 = adjacency(t(datExpr2),power=softPower,type="signed");
diag(adjacency2)=0
dissTOM2   = 1-TOMsimilarity(adjacency2, TOMType="signed")
geneTree2  = flashClust(as.dist(dissTOM2), method="average")

save.image("MetaAn.RData")  



pdf("dendrogram.pdf",height=6,width=16)
par(mfrow=c(1,2))
plot(geneTree1,xlab="",sub="",
     main="Gene clustering on TOM-based dissimilarity (TCGA)",labels=FALSE,hang=0.04);
plot(geneTree2,xlab="",sub="",
     main="Gene clustering on TOM-based dissimilarity (METABRIC)", labels=FALSE,hang=0.04); 
dev.off()

##### Now time to do the module creation
### will determine modules based on TCGA, since that was the one used in the modules
### for the initial work

mColorh=NULL
for (ds in 0:3){
  tree = cutreeHybrid(dendro = geneTree1, pamStage=FALSE,
                      minClusterSize = (30-3*ds), cutHeight = 0.99, 
                      deepSplit = ds, distM = dissTOM1)
  mColorh=cbind(mColorh,labels2colors(tree$labels));
}
pdf("Module_choices.pdf", height=10,width=25); 
plotDendroAndColors(geneTree1, mColorh, paste("dpSplt =",0:3), main = "",dendroLabels=FALSE);
dev.off()
modules1 =  mColorh[,1] # (Chosen based on plot below)


### Deepsplit = 0 has larger modules?
###### Next looking at prinicple components
### 1st PC = module Eigengene
### therefore, if ME for module X does a thing, so will all members in that module

PCs1A    = moduleEigengenes(t(datExpr1),  colors=modules1) 
ME_1A    = PCs1A$eigengenes
distPC1A = 1-abs(cor(ME_1A,use="p"))
distPC1A = ifelse(is.na(distPC1A), 0, distPC1A)
pcTree1A = hclust(as.dist(distPC1A),method="a") 
MDS_1A   = cmdscale(as.dist(distPC1A),2)
colorsA1 = names(table(modules1))

save.image("PC_colours.RData")

pdf("ModuleEigengeneVisualizations.pdf",height=6,width=6)
par(mfrow=c(1,1), mar=c(0, 3, 1, 1) + 0.1, cex=1)

plot(pcTree1A, xlab="",ylab="",main="",sub="")

plot(MDS_1A, col= colorsA1,  main="MDS plot", cex=2, pch=19)

ordergenes = geneTree1$order
plot.mat(scale(log(datExpr1[ordergenes,])) , 
         rlabels= modules1[ordergenes], 
         clabels= colnames(datExpr1), 
         rcols=modulesA1[ordergenes])


for (which.module in names(table(modules1))){
  ME = ME_1A[, paste("ME",which.module, sep="")] 
  barplot(ME, col=which.module, main="", cex.main=2, 
          ylab="eigengene expression",xlab="array sample") 
}; 

dev.off();


########  Qual and quant measures of network preservation at module level
#### IE how well are modules of TCGA conserved in METABRIC  
### Impose the modules from TCGA over the genetree of METABRIC and examine

pdf("Final_modules.pdf",height=8,width=12)
plotDendroAndColors(geneTree1, modules1, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, 
                    guideHang=0.05, main="Gene dendrogram and module colors (TCGA)") 
plotDendroAndColors(geneTree2, modules1, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, 
                    guideHang=0.05, main="Gene dendrogram and module colors (METABRIC)") 
dev.off()

### There is definitely still some module preservation - the colours run together


### get a z-score summary
# (This step will take ~10-30 minutes)
multiExpr  = list(A1=list(data=t(datExprA1g)),A2=list(data=t(datExprA2g)))
multiColor = list(A1 = modulesA1)
mp=modulePreservation(multiExpr,multiColor,referenceNetworks=1,verbose=3,networkType="signed",
                      nPermutations=30,maxGoldModuleSize=100,maxModuleSize=400)
stats = mp$preservation$Z$ref.A1$inColumnsAlsoPresentIn.A2
stats[order(-stats[,2]),c(1:2)]


