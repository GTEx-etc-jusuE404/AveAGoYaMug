



setwd("~/Bioinformatics Work/Meth & RNA/Meta-analysis WGCNA")

library(WGCNA)
library(plyr)
library(flashClust)
library(ggplot2)
library(gplots)
library(lattice)
library(extrafont)
loadfonts()
##############


load(file = "MetaAnalysis_trimmed_input.RData")
load(file = "Modules_DS0.RData")
colorsA1 = names(table(modules1))

### add in the rownames for the ME_ data frames

rownames(ME_1A) <- colnames(datExpr1)
rownames(ME_2A) <- colnames(datExpr2)

#### Mean centre the expression of the top100 from blue/green/yellow 
# for basals, in the METABRIC set
#then bin them??? ANOVAs / eyeball??



clin2 <- read.table("Clinical_final_METABRIC.txt", sep = "\t",
                    header = TRUE, row.names = 1)
#### Intersect this with the trimmed sample list for METABRIC
list <- intersect(colnames(datExpr2), rownames(clin2))

clin2 <- clin2[list,]
datExpr2 <- datExpr2[,list]
### okay, let's just do the basals for now

pam50 <- as.factor(clin2$Pam50_subtype)
clin_bas <- subset(clin2, pam50 == "Basal")

list <- intersect(colnames(datExpr2), rownames(clin_bas))

datExpr3 <- datExpr2[,list]


############################## now the mod list

exp_mod <- read.table("Top100_all_modules_both_networks.txt", sep = "\t",
                      header = TRUE)

blue <- as.data.frame(exp_mod$blue) 
green <- as.data.frame(exp_mod$green)
yellow <- as.data.frame(exp_mod$yellow)

####

list <- intersect(rownames(datExpr3), blue[,1])
subblue <- datExpr3[list,]


list <- intersect(rownames(datExpr3), yellow[,1])
subyellow <- datExpr3[list,]


list <- intersect(rownames(datExpr3), green[,1])
subgreen <- datExpr3[list,]


subyellow <- datExpr3[yellow,]



# do  mean centred and then bin

# create function to center : 'colMeans()'
center_colmeans <- function(x) {
  xcenter = colMeans(x)
  x - rep(xcenter, rep.int(nrow(x), ncol(x)))
}


centre_t <- center_colmeans(t(subgreen))
meta_ave <- as.data.frame(rowMeans(centre_t))

rownames(meta_ave) <- colnames(subgreen) 
colnames(meta_ave)[1] <- "Green_signature"


centre_t <- center_colmeans(t(subblue))
meta_ave$Blue_signature <- rowMeans(centre_t)


### yellow
centre_t <- center_colmeans(t(subyellow))

meta_ave$Yellow_signature <- rowMeans(centre_t)

meta_ave$Survival_time <- clin_bas$Survival.Time
meta_ave$Survival_event <- clin_bas$Survival.Event


write.table(meta_ave, "Mean_centred_modules_with clin.txt", sep = "\t")

####### to determine the cut off points for the signature

summary(meta_ave$Green_signature)
summary(meta_ave$Blue_signature)
summary(meta_ave$Yellow_signature)

bloc <- read.table("Mean_centred_modules_with clin.txt", sep = "\t", header = TRUE, row.names = 1)

list <- intersect(rownames(bloc), rownames(clin_bas))

bloc <- bloc[list,]
clin_bas <- clin_bas[list,]

bloc$Treatment <- clin_bas$Treatment


write.table(bloc, "Mean_centre_clinical_treatment.txt", sep = "\t")
