####  

#####  THIS IS A SCRIPT FOR VISANT VISULAIZATION OF MODULES

### note that you will need to datExpr files
## the modules 
## the module colour etc &
## the softpower from previous code.

## first gotta do the functions urg'
# ===================================================
#The function TOMdist1 computes a dissimilarity 
# based on the topological overlap matrix (Ravasz et al)
# Input: an Adjacency matrix with entries in [0,1]
if(exists("TOMdist1")) rm(TOMdist1);
TOMdist1=function(adjmat1, maxADJ=FALSE) {
  diag(adjmat1)=0;
  adjmat1[is.na(adjmat1)]=0;
  maxh1=max(as.dist(adjmat1) ); minh1=min(as.dist(adjmat1) ); 
  if (maxh1>1 | minh1 < 0 ) {print(paste("ERROR: the adjacency matrix contains entries that are larger than 1 or smaller than 0!!!, max=",maxh1,", min=",minh1)) } else { 
    if (  max(c(as.dist(abs(adjmat1-t(adjmat1)))))>10^(-12)   ) {print("ERROR: non-symmetric adjacency matrix!!!") } else { 
      adjmat1= (adjmat1+ t(adjmat1) )/2
      kk=apply(adjmat1,2,sum)
      maxADJconst=1
      if (maxADJ==TRUE) maxADJconst=max(c(as.dist(adjmat1 ))) 
      Dhelp1=matrix(kk,ncol=length(kk),nrow=length(kk))
      denomTOM= pmin(as.dist(Dhelp1),as.dist(t(Dhelp1)))   +as.dist(maxADJconst-adjmat1); 
      gc();gc();
      numTOM=as.dist(adjmat1 %*% adjmat1 +adjmat1);
      #TOMmatrix=numTOM/denomTOM
      # this turns the TOM matrix into a dissimilarity 
      out1=1-as.matrix(numTOM/denomTOM) 
      diag(out1)=1 
      # setting the diagonal to 1 is unconventional (it should be 0)
      # but it leads to nicer looking TOM plots... 
      out1
    }}}




visantPrepOverall <- function(colorFinal2, moduleColor, datExprrest, genes, numInt, power, signed=FALSE)
{
  ## This file returns the overall strongest connections within a module for one group of subjects
  
  ## USER inputs
  # colorFinal2 = vector with module association for each probe
  # moduleColor = color of the module... should be a member of colorFinal2
  # datExprrest = expression data for the genes corresponding to cIndex
  # genes = list of genes that correspond to the probes
  # numInt = number of interactions to output to the visant file
  
  cIndex = (colorFinal2==moduleColor)
  datExpr=datExprrest[,cIndex]
  if (signed){AdjMat =((1+cor(datExpr, use="p"))/2)^power
  } else {AdjMat =abs(cor(datExpr, use="p"))^power}
  diag(AdjMat)=0
  Degree <- apply(AdjMat,1,sum)
  Degree = Degree/max(Degree)
  meanExpr=apply(datExpr,2,mean)
  ProbeSet=colnames(datExpr)
  GeneSet=genes[cIndex];
  Module=rep(moduleColor,length(meanExpr)) 
  
  ## This file summarizes intramodular connectivity and expression for each gene in each group:
  fn=paste(moduleColor,"_connectivityOverall.csv",sep="")
  datConn = cbind(ProbeSet,GeneSet,meanExpr,Degree,Module)
  datConn = datConn[order(Degree, decreasing=TRUE),]
  write.table(datConn,file=fn,sep=",",
              row.names=F, col.names= c("probes","genes","meanExpr","kin","module"))
  write(paste(fn, "written."),"")
  
  ## TOM matrix
  distTOM <- TOMdist1(AdjMat)
  simTOM = 1-distTOM
  diag(simTOM)=0
  simTOMcutoff = simTOM
  
  ## Correlation matrix
  Pearson = cor(datExpr ,use="p")
  diag(Pearson)=0
  
  ## Dynamically determine the appropriate cutoff
  cutoff = 0.24
  len    = 10000
  dir    = "increase"
  loops  = 0
  split  = 0.01
  numInt = numInt+100
  while(len>100){
    loops = loops+1
    if (dir == "increase") { cutoff = cutoff+split; }
    if (dir == "decrease") { cutoff = cutoff-split; }
    indices = (simTOMcutoff>cutoff)
    len = sum(sum(indices))
    if (len < numInt) {dir = "decrease";}
    if (len >=numInt) {dir = "increase";}
    len = abs(len-numInt)
    if (loops>500){ len=0;}  
    if ((loops%%100)==0){ split = split/2; }
  }
  write(c(loops,cutoff,len),"")
  
  ## Output using cutoffs:
  indices = (simTOMcutoff[1,]>cutoff)
  datout=cbind(
    rep(GeneSet[1], length(ProbeSet[indices])),
    GeneSet[indices],rep(0, length(ProbeSet[indices])),
    rep("M1002", length(ProbeSet[indices])),
    simTOM[1,][indices], Pearson[1,][indices])
  colnames(datout) = c("gene1","gene2","zero","color","TO","Correlation")
  for(i in seq(2,length(ProbeSet),by=1)){
    indices = (simTOMcutoff[i,]>cutoff)
    datout=rbind(datout,cbind(
      rep(GeneSet[i], length(ProbeSet[indices])),
      GeneSet[indices],rep(0, length(ProbeSet[indices])),
      rep("M1002", length(ProbeSet[indices])),
      simTOM[i,][indices], Pearson[i,][indices]))
  }
  datout = datout[order(datout[,5],decreasing=TRUE),]
  fn = paste(moduleColor,"_visantOverall.csv",sep="")
  write.table(datout,file=fn,sep=",",row.names=F, col.names=c("gene1","gene2",
                                                              "zero","color","TO","Correlation"))
  write(paste(fn, "written."),"")
}


visantPrep <- function(colorFinal2, moduleColor, ciIndex, msIndex, datExprrest, genes, cutoff, power=power, signed=FALSE)
{
  ## This file compares two groups of subjects and returns the connections that are differential between the two groups
  
  ## USER inputs
  # colorFinal2 = vector with module association for each probe
  # moduleColor = color of the module... should be a member of colorFinal2
  # ciIndex = index of control subjects
  # msIndex = index of AD subjects
  # datExprrest = expression data for the genes corresponding to cIndex
  # genes = list of genes that correspond to the probes
  # cutoff = NOT A CUTOFF!!! NUMBER OF INTEGERS TO INCLUDE IN THE ANALYSIS!!!
  # power = power of the network
  # signed = is this a signed network TRUE/FALSE
  
  cIndex = (colorFinal2==moduleColor)
  datExprCI=datExprrest[ciIndex,cIndex]
  datExprMS=datExprrest[msIndex,cIndex]
  if (signed){AdjMatCI =((1+cor(datExprCI, use="p"))/2)^power
  } else {AdjMatCI =abs(cor(datExprCI, use="p"))^power}
  if (signed){AdjMatMS =((1+cor(datExprMS, use="p"))/2)^power
  } else {AdjMatMS =abs(cor(datExprMS, use="p"))^power}
  diag(AdjMatCI)=0
  diag(AdjMatMS)=0
  DegreeCI <- apply(AdjMatCI,1,sum)
  DegreeMS <- apply(AdjMatMS,1,sum)
  DegreeCI = DegreeCI/max(DegreeCI)
  DegreeMS = DegreeMS/max(DegreeMS)
  meanExprCI=apply(datExprCI,2,mean)
  meanExprMS=apply(datExprMS,2,mean)
  ProbeSet1=colnames(datExprMS)
  GeneSet1=genes[cIndex];
  Module=rep(moduleColor,length(meanExprMS)) 
  
  if (cutoff>=0){
    ## This file summarizes intramodular connectivity and expression for each gene in each group:
    
    write.table( cbind(ProbeSet1,GeneSet1,meanExprCI,meanExprMS,DegreeCI,DegreeMS,Module), 
                 file= paste(moduleColor,"_connectivity.csv",sep=""),sep=",",row.names=F, 
                 col.names= c("Probes","genes","meanExprCI","meanExprMS","kin_CI","kin_MS","Module"))
    
    ## Output correlation between kIn_CI and kIn_MS
    write("kIn_CI vs. kIN_MS correlation","")
    corr = cor(DegreeCI, DegreeMS, method="s")
    write(corr,"")
  }
  

}


################################################################



for(co in colorsA1[colorsA1!="grey"])
  visantPrepOverall(modules1, co, t(datExpr1), 
                    rownames(datExpr1), 500, softPower, TRUE)



for(co in colorsA1[colorsA1!="grey"])
  visantPrepOverall(modules1, co, t(datExpr2), 
                    rownames(datExpr2), 500, softPower, TRUE)


#### HUB genes specific to one network?

datExprA12g = t(cbind(datExpr1,datExpr2))
i1 = 1:dim(datExpr1)[[2]];
i2 = (1:dim(datExpr2)[[2]])+length(i1)
for (co in colorsA1[colorsA1!="grey"])
  visantPrep(modules1, co, i1, i2, datExprA12g, 
             rownames(datExpr1), 500, softPower, TRUE)


#####
### genes specific to TCGA can be found:
  ### in connectivity files, those with high kin_CI and low kin_MS are only hubs in TCGA
#### in viasnt files, genes with many connections ot TO-Ratio are hubs only in TCGA

### These last visant files need to be modified- the "0" and the colour columns
### are missing, and the columns A,C, E, F, H, and I are extraneous