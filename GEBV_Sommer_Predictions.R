########Scripts were modified from R package "Sommer" (Covarrubias-Pazaran, 2016)###########
setwd("C:/Users/~")
library(snow)
library(doSNOW)
library(parallel)
detectCores(logical=FALSE)
C1<-makeCluster(4,type="SOCK")
registerDoSNOW(C1)

library(sommer)
library(dplyr)
############################################
## ====================================== ##
############################################
# GEBV in the Rice population
############################################
## ====================================== ##
############################################
GT <- read.csv("GENOMAT_IMPUTED.csv",header=T)
GT[1:5,1:5]
attach(GT)
head(GT)
PT <-read.csv("PT_Rice1.csv",header=T)
attach(PT)
head(PT)
############################################################
## Numeric marker matrices and A.mats 
############################################################
rownames(PT)<- PT$id
PT[which(PT == "--", arr.ind = TRUE)] <- NA # To make missing value NA
GENOMAT<-(GT[,2:ncol(GT)])
GENOMAT[1:5,1:5]
GT1<-GENOMAT
colnames(PT) <- paste0("X",1:ncol(PT))
PT <- as.data.frame(PT);PT$id <- as.factor(rownames(PT))
# select environment 1
rownames(GT1) <- rownames(PT)
K <- A.mat(GT1) # additive relationship matrix
colnames(K) <- rownames(K) <- rownames(PT)
# GBLUP pedigree-based approach
y.trn <- PT
head(y.trn)
## GBLUP
ans <- mmer(X17~1,
            random=~vs(id,Gu=K),
            rcov=~units,
            data=y.trn) # kinship based
ans$U$`u:id`$X17 <- as.data.frame(ans$U$`u:id`$X17) # Genomic Estimated Breeding Value(GEBV as data frame)
rownames(ans$U$`u:id`$X17) <- gsub("id","",rownames(ans$U$`u:id`$X17))
result_list<-ans$U$`u:id`$X17
write.csv(result_list,file="GP_CORD.csv")