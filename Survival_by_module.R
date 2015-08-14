####

######### This is a script to see about surivival just using the 
### modules in both TCGA and metabric

### load in the files as needed

load(file = "MetaAnalysis_trimmed_input.RData")
load(file = "Modules_DS0.RData")
colorsA1 = names(table(modules1))

### add in the rownames for the ME_ data frames

rownames(ME_1A) <- colnames(datExpr1)
rownames(ME_2A) <- colnames(datExpr2)
#### do binning of surivival by modules
### just for now, just do within the basals and use a 1/2 (??1/3)


### you're going to need the clinical info for both sets

######## start with TCGA

clin1 <- read.table("Clinical_final_TCGA.txt", sep = "\t",
                    header = TRUE, row.names = 1)
#### Intersect this with the trimmed sample list for METABRIC
list <- intersect(colnames(datExpr1), rownames(clin1))

clin1 <- clin1[list,]
datExpr1 <- datExpr1[list,]
#### Now that you have trimmed the list, set the factors you would like
Pam50_1 <- as.factor(clin1$PAM50Call_RNAseq)


clin_bas_tcga <- clin1[Pam50_1 == "Basal",]
datExpr1_bas <- datExpr1[,Pam50_1 == "Basal"]

list <- intersect(rownames(clin_bas_tcga), rownames(ME_1A))
mod_bas_tcga <- ME_1A[list,]


###now bin

##############  Trying to make a for loop over a data frame

iterations = 76  ## number of samples (rows) in the data frame
variables = 13  ## number of variables/genes (cols) in the data frame

output <- matrix(ncol=variables, nrow=iterations)


for(i in 1:ncol(mod_bas_tcga))
{
  col = mod_bas_tcga[,i]
  brk = unique(quantile(col, probs=c(0, 0.25, 0.75, 1)))
  
  if ( length(brk) < 3)
  { 
    for(j in 1:nrow(mod_bas_tcga))
      output[,i] = 'NA'
  }
  else
  {
    output[,i] =  as.factor(cut(col,
                                breaks=brk,
                                labels = FALSE, 
                                include.lowest=TRUE))
  }
}


output <- data.frame(output)
rownames(output) <- rownames(mod_bas_tcga)
colnames(output) <- colnames(mod_bas_tcga)

bads <- names(which(sapply(output, function(x) any(x == "NA"))))

for(i in 1:length(bads))
{
  output[,bads[i]] <- NULL
}


#write.table(output, "Bins_.txt", sep = "\t")


##################  Then do survival testing (Anova / cox) and order by which is sig. 
output <- read.table("Bins_exp_basal.txt", sep = "\t", row.names =1, header = TRUE)
#check if clinical file is there.  Toggle values
### To do the tops and tails, will need to exclude all the samples with a "level2" 
### ie, the middle group - they are unhelpful. 
### or do the anova one again??

library(survival)

SurvTest <- output

attach(SurvTest)

aov.out <- matrix(ncol = 5, nrow = 13)

## How to use only the top and bottom factors, and exlclude the others
#call <- factor(SurvTest$SOX10, levels = c(1, 3))

for(i in 1:ncol(SurvTest)) {
  
  #level = factor(SurvTest[,i])
  level = factor(SurvTest[,i], levels = c(1,3))
  survival_time = clin_bas_tcga$OS_Time_nature2012
  survival_event = clin_bas_tcga$OS_event_nature2012
  survival = Surv(survival_time,survival_event == 1) ~ level
  cox = coxph(survival, data = SurvTest)
  coeffs <- as.matrix(coef(summary(cox)))  ## coeffs has more info
  #aov = anova(cox)  ### anovas, for eg more than 2 groups
  aov.out[i,] <- coeffs
  #aov.out[i,] = (aov$`Pr(>|Chi|)`)[2]  ### output anovas from above
  
}

#### PUT IN GENE NAMES AS ROWNAMES HERE
rownames(aov.out) <- colnames(SurvTest)
colnames(aov.out) <- c("log Hazard Ratio", "Hazard Ratio", "SE Hazard Ratio", "Z score", "p.value" )

write.table(aov.out, "Survival_cox_all_tops_and_tails.txt", sep = "\t")










########################################################

#### read in the clinical file for metabric
clin2 <- read.table("Clinical_final_METABRIC.txt", sep = "\t",
                   header = TRUE, row.names = 1)
#### Intersect this with the trimmed sample list for METABRIC
list <- intersect(colnames(datExpr2), rownames(clin2))

clin2 <- clin2[list,]
#### Now that you have trimmed the list, set the factors you would like
Pam50_2 <- as.factor(clin2$Pam50_subtype)

clin_bas_meta <- clin2[Pam50_2 == "Basal",]
datExpr2_bas <- datExpr2[,Pam50_2 == "Basal"]

rownames(ME_2A) <- colnames(datExpr2)

list <- intersect(rownames(clin_bas_meta), rownames(ME_2A))
mod_bas_meta <- ME_2A[list,]


iterations = 177  ## number of samples (rows) in the data frame
variables = 13  ## number of variables/genes (cols) in the data frame

output <- matrix(ncol=variables, nrow=iterations)


for(i in 1:ncol(mod_bas_meta))
{
  col = mod_bas_meta[,i]
  brk = unique(quantile(col, probs=c(0, 0.25, 0.75, 1)))
  
  if ( length(brk) < 3)
  { 
    for(j in 1:nrow(mod_bas_meta))
      output[,i] = 'NA'
  }
  else
  {
    output[,i] =  as.factor(cut(col,
                                breaks=brk,
                                labels = FALSE, 
                                include.lowest=TRUE))
  }
}


output <- data.frame(output)
rownames(output) <- rownames(mod_bas_meta)
colnames(output) <- colnames(mod_bas_meta)

bads <- names(which(sapply(output, function(x) any(x == "NA"))))

for(i in 1:length(bads))
{
  output[,bads[i]] <- NULL
}


#write.table(output, "Bins_.txt", sep = "\t")


##################  Then do survival testing (Anova / cox) and order by which is sig. 
output <- read.table("Bins_exp_basal.txt", sep = "\t", row.names =1, header = TRUE)
#check if clinical file is there.  Toggle values
### To do the tops and tails, will need to exclude all the samples with a "level2" 
### ie, the middle group - they are unhelpful. 
### or do the anova one again??

library(survival)

SurvTest <- output

attach(SurvTest)

aov.out <- matrix(ncol = 5, nrow = 13)

## How to use only the top and bottom factors, and exlclude the others
#call <- factor(SurvTest$SOX10, levels = c(1, 3))

for(i in 1:ncol(SurvTest)) {
  
  #level = factor(SurvTest[,i])
  level = factor(SurvTest[,i], levels = c(1,3))
  survival_time = clin_bas_meta$Survival.Time
  survival_event = clin_bas_meta$Survival.Event
  survival = Surv(survival_time,survival_event == 1) ~ level
  cox = coxph(survival, data = SurvTest)
  coeffs <- as.matrix(coef(summary(cox)))  ## coeffs has more info
  #aov = anova(cox)  ### anovas, for eg more than 2 groups
  aov.out[i,] <- coeffs
  #aov.out[i,] = (aov$`Pr(>|Chi|)`)[2]  ### output anovas from above
  
}

#### PUT IN GENE NAMES AS ROWNAMES HERE
rownames(aov.out) <- colnames(SurvTest)
colnames(aov.out) <- c("log Hazard Ratio", "Hazard Ratio", "SE Hazard Ratio", "Z score", "p.value" )

##write.table(aov.out, "Survival_cox_all_tops_and_tails.txt", sep = "\t")


#### none of these are significant - probably the dilution effect of too
## many useless genes at the end of the module


### Survivial curves
library(survival)


mods <- c("MEblue", "MEgreen", "MEyellow")

mods_s <- output[,mods]
attach(mods_s)


survival_time = clin_bas_tcga$OS_Time_nature2012
survival_event = clin_bas_tcga$OS_event_nature2012
level = factor(mods_s$MEblue, levels = c(1,3))


fit.diff = survdiff(Surv(survival_time,survival_event == 1) ~ level) 

chisq2 = signif(1-pchisq(fit.diff$chisq,length(levels(level))-1),3) 
fit1 = survfit(Surv(survival_time,survival_event == 1)~level,
               conf.type="log-log") 
kmcolours <- c("black", "red", "green", "blue")
plot(fit1, conf.int=F,col=kmcolours,xlab="Time to death (days)", 
     ylab="Survival",main=c("All subtype survival by\n Module Expression"), 
     lwd=4,mark.time=TRUE)
legend("bottomleft",legend=levels(factor(mods_s)),
       fill = kmcolours, cex = 1)

###do another anova for signifigance
anova(coxph(Surv(OS_Time_nature2012,OS_event_nature2012 == 1)~factor(Brown), data = datClin))


### add to the end and do a prism one instead


mods <- c("MEblue", "MEgreen", "MEyellow")

mods_m <- output[,mods]

survival_time = clin_bas_meta$Survival.Time
survival_event = clin_bas_meta$Survival.Event

mods_m$survival_time <- survival_time
mods_m$survival_event <- survival_event


write.table(mods_s, "Modules_survival_TCGA.txt", sep = "\t")
write.table(mods_m, "Modules_survival_METABRIC.txt", sep = "\t")
