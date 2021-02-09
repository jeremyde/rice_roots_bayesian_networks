#The script was used to compute and analyze variances and BLUPS for root system architecture (RSA) traits in Augmetned Design with replicated checks##############
####Script were modified from lme4 R pakcage (Bates et al., 2011)###################
setwd("C:/Users/~")
inputData <- read.csv("~.csv",header=T)
attach(inputData)
head(inputData)
library(Matrix)
library(lme4)
options("scipen"=100, "digits"=4)
# lets force R to treat the variables as factors 
inputData$Genotype <- as.factor(Genotype)
inputData$Block <- as.factor(Block)
Check <- as.factor(Check)
Check.Gen <- as.factor(Check.Gen)
inputData$LRPL <- as.numeric(LRPL)
str(inputData)
library(lmerTest)
str(inputData);  # Display the structure of an arbitrary R Object compactly
options(lmerControl=list(check.nobs.vs.rankZ = "warning",
                         check.nobs.vs.nlev = "warning",
                         check.nobs.vs.nRE = "warning",
                         check.nlev.gtreq.5 = "warning",
                         check.nlev.gtr.1 = "warning"))
#Random model to calculate BLUPS and variances considering checks fixed
LRPLvarcomp <- lmer(data=inputData, LRPL ~ Check.Gen + (1|Block) + (1|Genotype:Check), REML=T)
modelSum <- summary(LRPLvarcomp)
print (modelSum)
anovaout<-anova(LRPLvarcomp)
print(anovaout)
sink("summaryLRPL.csv")
summary (LRPLvarcomp)
LRPLMODEL <- lmer(data=inputData, LRPL ~ (1|Genotype) + (1|Block) + (1|Genotype%in%Block), REML=T)
LRPLBLUP=ranef(LRPLMODEL)
str(LRPLBLUP)
LRPLBLUPGenotype = LRPLBLUP$Genotype
#write.csv(LRPLBLUP, file = "LRPLBLUPS.csv")
BLUPS = LRPLBLUPGenotype[,1]
#hist(RILBLUP, col="brown")
lmean=tapply(LRPL, Genotype, na.rm =T, mean)
plot(BLUPS, lmean, main="LRPLBLUPS",col ="blue")
sink()


