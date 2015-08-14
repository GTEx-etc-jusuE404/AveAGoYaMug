meth <- read.table("Meth_TSS200_all.txt", sep = "\t", header = TRUE, row.names = 1)


list <- intersect(rownames(meth), colnames(datExpr1))

meth <- meth[list,]
datExpr1 <- datExpr1[,list]

Me_3 <- ME_1A[list,]

rownames(Me_3) <- colnames(datExpr1)



library(WGCNA)
####### OK, now do as for the 


nGenes = ncol(Me_3);
nSamples = nrow(Me_3);
# Recalculate MEs with color labels
#MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
#MEs = orderMEs(MEs0)

moduleTraitCor = cor(Me_3, meth, use = "p");
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
               yLabels = names(Me_3),
               ySymbols = names(Me_3),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               #textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.85,
               zlim = c(-1,1),
               main = paste("High module expression vs trait genes"))


write.table(moduleTraitCor, "All_meth_trait_corr.txt", sep = "\t")
write.table(moduleTraitPvalue, "All_meth)trait_pValue.txt", sep = "\t")

cormat2 <- read.table("All_meth_trait_corr.txt", sep = "\t", header = TRUE, row.names = 1)
#### 


pdf("All_promoter_meth_traits.pdf")
par(oma= c(3,1,1,1))

### Info about heatmap.2
### must be done on a matrix not a dataframe
### dendro can be set to none or some
### key gives you a colour key
### srtCol/ srtRow gives an angle for the text
##sepCol, etc gives a cell outline
### RowSideColors - can set colours for marking the cols or rows



### colorpanel in gplots
###redgreen(n)
###redblue(n)
###bluered(n)
###colorpanel(n, low, mid, high)
## or can set the colour panel before
##eg palette <- colorRampPalette(c('#f0f3ff','#0033BB'))(256)


heatmap.2(as.matrix(cormat2),dendrogram="none", 
          notecol="black",col=bluered(200),scale="none",
          labRow = NULL, 
          key=TRUE, keysize=1.5,key.title = NA, key.xlab = "Correlation",
          density.info="none", trace="none",
          cexRow=1.2,cexCol=1.2, srtCol = 75,
          main = "Traits expression \n by module in TCGA dataset")

library(pheatmap)
library(RColorBrewer)

pdf("Meth_traits_test.pdf",height=8,width=10)
fontsize = 10


pheatmap(t(cormat2), color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(100), 
         kmeans_k = NA, breaks = NA, border_color = "grey60",
         cellwidth = NA, cellheight = NA, scale = "none", 
         cluster_rows = TRUE, cluster_cols = FALSE, fontsize_row = 10,
         annotation_legend = TRUE, drop_levels = TRUE, 
         show_rownames = F,show_colnames = T, 
         treeheight_col = 0,
         main = "Methylation of gene promoters by module")
dev.off()




### trimmed to only the genes for which we have expression data?
### how many are NOT found in the expression set???



list <- intersect(rownames(datExpr1), colnames(cormat2))

cormat3 <- cormat2[,list]



pdf("Meth_traits_test_trimmed.pdf",height=8,width=10)
fontsize = 10


pheatmap(t(cormat3), color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(100), 
         kmeans_k = NA, breaks = NA, border_color = "grey60",
         cellwidth = NA, cellheight = NA, scale = "none", 
         cluster_rows = TRUE, cluster_cols = FALSE, fontsize_row = 10,
         annotation_legend = TRUE, drop_levels = TRUE, 
         show_rownames = F,show_colnames = T, 
         treeheight_col = 0,
         main = "Methylation of gene promoters by module")
dev.off()



############# looking at TRIM29

list <- intersect(rownames(clin1), colnames(datExpr1))

clin1 <- clin1[list,]
exp <- datExpr1[,list]
exp <- t(exp)
exp <- as.data.frame(exp)

meth1 <- meth[list,]


pam50 = as.factor(clin1$PAM50Call_RNAseq)




###################
### meth information by all / basal

boxplot(meth1$TRIM29 ~ pam50)

summary(meth1$TRIM29)
sd(meth1$TRIM29)
  
summary(datadata)
sd(datadata)


###################
### exp information by all / basal

boxplot(exp$TRIM29 ~ pam50)


####### correlations between them

allcor <- cor(exp$TRIM29, meth1$TRIM29)
allcor
corPvalueStudent(allcor, nSamples = 582)

bascor <- cor(datadata, expexp)
bascor
corPvalueStudent(bascor, nSamples = 76)
