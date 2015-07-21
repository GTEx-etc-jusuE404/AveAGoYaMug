
#### for trait matrixes, outputting data also

####################################  
#####################################TCGA Start
## need to set the rownames for ME_1A
rownames(ME_1A) <- colnames(datExpr1)

datTraits1 <- read.table("Traits_nonmeth_TCGA.txt", sep = "\t",
                         header = TRUE, row.names = 1)

list <- intersect(rownames(ME_1trim), rownames(datTraits1))

datTraits1 <- datTraits1[list,]
ME_1trim <- ME_1A[list,]
####### OK, now do as for the 


nGenes = ncol(ME_1trim);
nSamples = nrow(ME_1trim);
# Recalculate MEs with color labels
#MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
#MEs = orderMEs(MEs0)

moduleTraitCor = cor(ME_1trim, datTraits1, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");

##########Will only display correlations
textMatrix = paste(signif(moduleTraitCor, 2),sep = "");
dim(textMatrix) = dim(moduleTraitCor)

par(mar = c(6, 8.5, 3, 3));

# Display the values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(datTraits1),
               yLabels = names(ME_1trim),
               ySymbols = names(ME_1trim),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               #textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.85,
               zlim = c(-1,1),
               main = paste("High module expression vs trait genes"))


#################################
## label and write out the module trait cor and the pvalue

cormat <- as.matrix(moduleTraitCor)
colnames(cormat) <- colnames(datTraits1)
rownames(cormat) <- colnames(ME_1trim)

write.table(cormat, "Correlation_matrix_TCGA_traits.txt", sep = "\t")

corp <- as.matrix(moduleTraitPvalue)
colnames(corp) <- colnames(datTraits1)
rownames(corp) <- colnames(ME_1trim)

write.table(corp, "Cor_Pvalue_matrix_TCGA_traits.txt", sep = "\t")

######  Let's have a look with the meth traits for TCGA

datTraits3 <- read.table("methTraits2_Stir200i.txt", sep = "\t",
                         header = TRUE, row.names = 1)

list <- intersect(rownames(datTraits3), rownames(ME_1trim))

datTraits3 <- datTraits3[list,]
ME_1trim <- ME_1trim[list,]


####### OK, now do as for the 


nGenes = ncol(ME_1trim);
nSamples = nrow(ME_1trim);
# Recalculate MEs with color labels
#MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
#MEs = orderMEs(MEs0)

moduleTraitCor = cor(ME_1trim, datTraits3, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");

##########Will only display correlations
textMatrix = paste(signif(moduleTraitCor, 2),sep = "");
dim(textMatrix) = dim(moduleTraitCor)

par(mar = c(6, 8.5, 3, 3));

# Display the values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = NULL,
               yLabels = names(ME_1trim),
               ySymbols = names(ME_1trim),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               #textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.85,
               zlim = c(-0.6,0.6),
               main = paste("High module expression vs trait genes"))




#### the correlation matrix for the modules and the methylation values


cormat <- as.matrix(moduleTraitCor)
colnames(cormat) <- colnames(datTraits3)
rownames(cormat) <- colnames(ME_1trim)

write.table(cormat, "Correlation_matrix_TCGA_meth.txt", sep = "\t")

corp <- as.matrix(moduleTraitPvalue)
colnames(corp) <- colnames(datTraits3)
rownames(corp) <- colnames(ME_1trim)

write.table(corp, "Cor_Pvalue_matrix_TCGA_meth.txt", sep = "\t")


########################################################################
############  METABRIC



datTraits2 <- read.table("Traits_trimmed_METABRIC.txt", sep = "\t",
                         header = TRUE, row.names = 1)

list <- intersect(rownames(ME_2A), rownames(datTraits2))
datTraits2 <- datTraits2[list,]


################# for the metabric

nGenes = ncol(ME_2A);
nSamples = nrow(ME_2A);
# Recalculate MEs with color labels
#MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
#MEs = orderMEs(MEs0)


moduleTraitCor = cor(ME_2A, datTraits2, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");

##########Will only display correlations
textMatrix = paste(signif(moduleTraitCor, 2),sep = "");
dim(textMatrix) = dim(moduleTraitCor)

par(mar = c(6, 8.5, 3, 3));

# Display the values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(datTraits2),
               yLabels = names(ME_2A),
               ySymbols = names(ME_2A),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               #textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.85,
               zlim = c(-1,1),
               main = paste("High module expression vs trait genes"))

###############  Correlation tables

cormat <- as.matrix(moduleTraitCor)
colnames(cormat) <- colnames(datTraits2)
rownames(cormat) <- colnames(ME_2A)

write.table(cormat, "Correlation_matrix_METABRIC_traits.txt", sep = "\t")

corp <- as.matrix(moduleTraitPvalue)
colnames(corp) <- colnames(datTraits2)
rownames(corp) <- colnames(ME_2A)

write.table(corp, "Cor_Pvalue_matrix_METABRIC_traits.txt", sep = "\t")

