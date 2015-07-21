
##### This is a script to look at the clinical manifestations of the modules


#### read in the clinical file for metabric
  clin <- read.table("Clinical_final_METABRIC.txt", sep = "\t",
                   header = TRUE, row.names = 1)
#### Intersect this with the trimmed sample list for METABRIC
  list <- intersect(colnames(datExpr2), rownames(clin))

  clin <- clin[list,]
#### Now that you have trimmed the list, set the factors you would like
Pam50 <- as.factor(clin$Pam50_subtype)
ER <- as.factor(clin$ESR1_status)
grade <- as.factor(clin$grade)
site <- as.factor(clin$metastasis_site)
mets <- as.factor(clin$MET_Cat)
type <- as.factor(clin$histological_type)

### The data frame ME_2A is the MEs for the METABRIC set
####MEs <- as.data.frame(PCs2A$eigengenes)

###### For the survival stuff
## make a new clin file for only the surv event = 1
sur <- as.factor(clin$Survival.Event)
clinsurv <- subset(clin, sur == "1")

### trim the MEs accordingly
rownames(ME_2A) <- colnames(datExpr2)
list <- intersect(rownames(clinsurv), rownames(ME_2A))
MEsurv <- ME_2A[list,]


pdf("Module_pam50.pdf",height=8,width=8)
for (c in 1:length(colorsA1)){
  inMod = modules1 == colorsA1[c]
  boxplot(ME_2A[,c] ~ Pam50,
          main=colorsA1[c],
          xlab="Pam 50",ylab="ME in METABRIC")
}; dev.off()


pdf("Module_grade.pdf",height=8,width=8)
for (c in 1:length(colorsA1)){
  inMod = modules1 == colorsA1[c]
  boxplot(ME_2A[,c] ~ grade,
          main=colorsA1[c],
          xlab="Grade",ylab="ME in METABRIC")
}; dev.off()


pdf("Module_MetSite.pdf",height=8,width=8)
for (c in 1:length(colorsA1)){
  inMod = modules1 == colorsA1[c]
  boxplot(ME_2A[,c] ~ site,
          main=colorsA1[c],
          xlab="Met Site",ylab="ME in METABRIC")
}; dev.off()

###### For MET_CAT
pdf("Module_MetCat.pdf",height=8,width=8)
for (c in 1:length(colorsA1)){
  inMod = modules1 == colorsA1[c]
  boxplot(ME_2A[,c] ~ mets,
          main=colorsA1[c],
          xlab="Met Cat",ylab="ME in METABRIC")
}; dev.off()




#### Note that this code now contains rounding of the cor and p.value
### need to change this if that bothers you.
corr_out <- matrix(nrow = 35, ncol = 2)

for (c in 1:ncol(MEsurv)){
  cor <- cor.test(clinsurv$Survival.Time, MEsurv[,c])
  corr_out[c,1] <- round(cor$estimate, digit = 2)
  corr_out[c,2] <- round(cor$p.value, digit = 2)
}

rownames(corr_out) <- colnames(MEsurv)
write.table(corr_out, "Survival_exp_correlations.txt", sep ="\t")


##### Survival v ME plots
pdf("Module_Surv.pdf",height =8,width=8)
for (c in 1:length(colorsA1)){
  inMod = modules1 == colorsA1[c]
  plot(clinsurv$Survival.Time, MEsurv[,c],
       main=paste(colorsA1[c],"\nCor = ",corr_out[c,1], "P.value = ",corr_out[c,2]),
       xlab="Survival Time (Days)",ylab="ME in METABRIC")
}; dev.off()


