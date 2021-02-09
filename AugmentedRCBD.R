# Script was used to compare mean and visualize distribution of Augmented design for RSA triats. The codes were modified from The R package augmentedRCBD (Aravind et al., 2020) 
setwd("C:/Users/~")
data <- read.csv("~.csv",header=T)
attach(data)
library(Matrix)
library(augmentedRCBD)

options("scipen"=100, "digits"=4)
# lets force R to treat the variables as factors rather than integers
data$Block <- as.factor(data$Block)
data$Genotype <- as.factor(data$Genotype)
#Results for variable y1 (checks inferred)
#Analyzing Augmented Design with replicated checks
out1 <- augmentedRCBD(data$Block, data$Genotype, data$RV, method.comp = "lsd",
                      alpha = 0.05, group = TRUE, console = TRUE)
#Analyzing Frequency Distributions
freq1 <- freqdist.augmentedRCBD(out1, xlab = "Root Volume",
                                highlight.check = FALSE)
class(freq1)
plot(freq1)


