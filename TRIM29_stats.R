


##### TRIM 29 stats info


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

summary(exp$TRIM29)
sd(exp$TRIM29)

summary(expexp)
sd(expexp)


####### correlations between them

allcor <- cor(exp$TRIM29, meth1$TRIM29)
allcor
corPvalueStudent(allcor, nSamples = 582)

bascor <- cor(datadata, expexp)
bascor
corPvalueStudent(bascor, nSamples = 76)

