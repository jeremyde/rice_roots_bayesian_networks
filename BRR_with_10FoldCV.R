####Script was used to compute genomic predictive ability (GPA) of agronomic (18) and RSA (10) traits### 
####Script were modified from R pakcage "BGLR" (Pérez,  P., and G. de los Campos, 2014)###############  
####and their 10 fold cross validations.#################################################################
setwd("C:/Users/santosh.sharma/Documents/Paper 1/GENOMICPREDICT/BayesianABC/BB/BRRAG")
#Loading and preparing the input data
library(BGLR)
library(dplyr)
GT1 <- read.csv("~.csv",header=T)
GT<-(GT1[,2:ncol(GT1)])
PT1 <-read.csv("~.csv",header=T)

N=PT1[,6]
Z=scale(GT,scale=TRUE,center=TRUE) 
Z[1:5,1:5]
############Ten Fold Cross Validation########################
## Fitting models 
nIter<-45000 #For real data sets more samples are needed
burnIn<-5000
thin<-10
folds<-10
set.seed(123) #Set seed for the random number generator
sets<-rep(1:10,60)[-1]
sets<-sets[order(runif(nrow(Z)))]
################Creating FOlds####################3
COR.CV<-rep(10,times=(folds+1))
names(COR.CV)<-c(paste('fold=',1:folds,sep=''),'Pooled')
ETA<-list(MRK=list(X=Z,model="BRR")) 
w<-rep(1/nrow(Z),folds) ## weights for pooled correlations and MSE
yHatCV<-numeric()
for(fold in 1:folds)
{
  yNa<-N
  whichNa<-which(sets==fold)
  yNa[whichNa]<-NA
  prefix<-paste('PM_BL','_fold_',fold,'_',sep='')
  fmBRR<-BGLR(y=yNa,ETA=ETA, nIter=nIter, burnIn=burnIn,thin=thin)
  yHatCV[whichNa]<-fmBRR$yHat[fmBRR$whichNa]
  w[fold]<-w[fold]*length(fmBRR$whichNa)
  COR.CV[fold]<-cor(fmBRR$yHat[whichNa],N[whichNa],use="complete")
}
COR.CV[11]<-mean(COR.CV[1:10])
COR.CV
write.csv(COR.CV,file= "~.csv")
