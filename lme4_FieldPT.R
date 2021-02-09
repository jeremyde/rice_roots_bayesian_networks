#The script was used for computing variances, and BLUPs and their analysis of phenotypic data colelcted in flooded field###
####Script were modified from lme4 R pakcage (Bates et al., 2011)###################
setwd("C:/Users/~")
data1=read.csv("~",header=T)
attach(data1)
head(data1)
library(Matrix)
library(lme4)
options("scipen"=100, "digits"=4)
options(lmerControl=list(check.nobs.vs.rankZ = "warning",
                          check.nobs.vs.nlev = "warning",
                          check.nobs.vs.nRE = "warning",
                          check.nlev.gtreq.5 = "warning",
                          check.nlev.gtr.1 = "warning"))
#Random model to calculate BLUPS and variances
 PHTvarcomp <- lmer(data=data1, PHT ~ (1|RIL) + (1|Year) + (1|Rep%in%Year) + (1|RIL:Year), REML=T)
View(PHTvarcomp)
summary (PHTvarcomp)
sink("summaryPHT.csv")
summary (PHTvarcomp)
Sink()
PHTMODEL=lmer(data=data1, PHT ~ (1|RIL) + (1|Year) + (1|Rep%in%Year) + (1|RIL:Year))
PHTBLUP=ranef(PHTMODEL)
str(PHTBLUP)
PHTBLUPRIL = PHTBLUP$RIL
write.csv(PHTBLUP, file = "PHT BLUPS.csv")
RILBLUP = PHTBLUPRIL[,1]
hist(RILBLUP, col="brown")
lmean=tapply(PHT, RIL, na.rm =T, mean)
plot(RILBLUP, lmean, col ="blue")