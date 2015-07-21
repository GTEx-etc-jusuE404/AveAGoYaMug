
###########  doing the metaanalysis but using a different deepsplit

setwd("~/Bioinformatics Work/Meth & RNA/Meta-analysis WGCNA/DeepSplit2")

library(WGCNA)
library(impute)
library(dynamicTreeCut)
library(qvalue)
library(flashClust)
library(Hmisc)
library(blockmodeling)

######  ***** READ IN THE DATA FILE ****** #########################

load(file = "MetaAnalysis_trimmed_input.RData")


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

#### This is where we nominate a different module set
modules1 =  mColorh[,3] # (Chosen based on plot below)

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

#save.image("PC_colours.RData")

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
multiExpr  = list(A1=list(data=t(datExpr1)),A2=list(data=t(datExpr2)))
multiColor = list(A1 = modules1)
mp=modulePreservation(multiExpr,multiColor,referenceNetworks=1,verbose=3,networkType="signed",
                      nPermutations=30,maxGoldModuleSize=100,maxModuleSize=400)
stats = mp$preservation$Z$ref.A1$inColumnsAlsoPresentIn.A2
stats[order(-stats[,2]),c(1:2)]


### check outputs - Higher z.score = more preservation between data sets
### score between 5 & 10 is moderate preservation, >10 is high preservation
### "gold" module contains random genes and grey the uncharacterised so expect
## these to be lower
### ??about turq though?

### look for interest genes

Modules <- data.frame(rownames(datExpr1), mColorh)
lookfor <- c("DNMT1", "DNMT3A", "DNMT3B", "UHRF1", "EZH2", "DNMT3L", "MTHFR", "SOX10")

Interest_mods <- Modules[Modules$rownames.datExpr1 %in% lookfor,] 

#####

### kME - module membership values

geneModuleMembership1 = signedKME(t(datExpr1), ME_1A)
colnames(geneModuleMembership1)=paste("PC",colorsA1,".cor",sep=""); 

MMPvalue1=corPvalueStudent(as.matrix(geneModuleMembership1),dim(datExpr1)[[2]]); 
colnames(MMPvalue1)=paste("PC",colorsA1,".pval",sep="");

Gene       = rownames(datExpr1)
kMEtable1  = cbind(Gene,Gene,modules1)
for (i in 1:length(colorsA1))
  kMEtable1 = cbind(kMEtable1, geneModuleMembership1[,i], MMPvalue1[,i])
colnames(kMEtable1)=c("PSID","Gene","Module",
                      sort(c(colnames(geneModuleMembership1), colnames(MMPvalue1))))

write.csv(kMEtable1,"kMEtable1.csv",row.names=FALSE)

###  Now repeat for METABRIC, using the module assignments from TCGA to determine kME values.

# First calculate MEs for A2, since we haven't done that yet
PCs2A = moduleEigengenes(t(datExpr2),  colors=modules1) 
ME_2A = PCs2A$eigengenes

geneModuleMembership2 = signedKME(t(datExpr2), ME_2A)
colnames(geneModuleMembership1)=paste("PC",colorsA1,".cor",sep=""); 

MMPvalue2=corPvalueStudent(as.matrix(geneModuleMembership2),
                           dim(datExpr2)[[2]]); 
colnames(MMPvalue2)=paste("PC",colorsA1,".pval",sep="");

kMEtable2  = cbind(Gene,Gene,modules1)
for (i in 1:length(colorsA1))
  kMEtable2 = cbind(kMEtable2, geneModuleMembership2[,i], MMPvalue2[,i])
colnames(kMEtable2)=colnames(kMEtable1)

write.csv(kMEtable2,"kMEtable2.csv",row.names=FALSE)

pdf("all_kMEtable2_vs_kMEtable1.pdf",height=8,width=8)
for (c in 1:length(colorsA1)){
  verboseScatterplot(geneModuleMembership2[,c],geneModuleMembership1[,c],main=colorsA1[c],
                     xlab="kME in A2",ylab="kME in A1")
}; dev.off()

pdf("inModule_kMEtable2_vs_kMEtable1.pdf",height=8,width=8)
for (c in 1:length(colorsA1)){
  inMod = modules1== colorsA1[c]
  verboseScatterplot(geneModuleMembership2[inMod,c],geneModuleMembership1[inMod,c],main=colorsA1[c],
                     xlab="kME in A2",ylab="kME in A1")
}; dev.off()
save.image("tutorial.RData") #(optional line of code)

#### beautiful!!! 

####### ok now you can write out a file with the relevant information
### about the modules using deeplsplit3


####  NOW determine which genes are hubs in both networks
topGenesKME = NULL
for (c in 1:length(colorsA1)){
  kMErank1    = rank(-geneModuleMembership1[,c])
  kMErank2    = rank(-geneModuleMembership2[,c])
  maxKMErank  = rank(apply(cbind(kMErank1,kMErank2+.00001),1,max))
  topGenesKME = cbind(topGenesKME,Gene[maxKMErank<=100])
}; colnames(topGenesKME) = colorsA1
topGenesKME

write.table(topGenesKME, "Top100_all_modules_both_networks.txt", sep = "\t")

### This is the top 100 most interconnected genes by module (in both networks)

### Save out the relevants into a loadable RData file to prevent doing this shit again,

save(geneTree1, geneTree2, modules1, ME_1A, ME_2A, file = "DeepSplit2_processed.RData" )

write.csv(geneModuleMembership1, "Gene_MM_TCGA.csv")
write.csv(geneModuleMembership2, "Gene_MM_METABRIC.csv")
write.csv(MMPvalue1, "Gene_MM_Pvalue_TCGA.csv")
write.csv(MMPvalue2, "Gene_MM_Pvalue_METABRIC.csv")
