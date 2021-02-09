# The script was used to learn Bayesian Network and compute cross validated genomic predictive ability (GPA)####
########### of 18 agronomic and 13 RSA traits.##################################################################
##############The script were modified from Bnlearn R package (Scutari et al, 2014)############################
library(lme4)
library(parallel)
library(DAAG)
library(tidyverse)
library(BGLR)
library(dplyr)
library(bnlearn)
library(visNetwork)
library(Rgraphviz)
setwd("C:/Users/~")
tildata <- read.table("~.txt",
                      header=T,
                      sep="\t",
                      stringsAsFactors=F,
                      strip.white=T,
                      na.strings='NA',
                      quote="")
#tildata=scale(tildata1,scale=TRUE,center=TRUE) 
#write.csv(tildata,file="~.csv")
#to standardize and normalize data
#tildata<-scale(tildata[2:33])
#tildata <- subset(tildata, select = -c(PercentChalk, GrainLengthWidthRatio,AvgYldPerTill, SdWtMinusPrimPan, SeedWtPrimPan, SeedCntPrimPan))
tildata <- subset(tildata,select= -c(CORT))

# Organize traits into tiers according to temportal order

tier1 <- c("LRL","CRL","LRT","CRSA","LRSA",TSA")
tier2<- c("RV","CR","LRSAP","LRPL","LRD","TRD","LRAL")
tier3<- c("LFNV","LFWDV","LFLV","GRLFEV","GREV","TLFAV")
tier4<- c("TAR","TNR","GRER","GRVR")
tier5<- c("HDT","CHLF","LAIR")
tier6<- c("PHT","TAT","TNT","GRRT")
traittiers <- c(list(tier1),list(tier2),list(tier3),list(tier4),list(tier5),list(tier6))
ids <- names(tildata)[1]
traits <- as.vector(names(tildata)[2:32])
genes <- names(tildata)[33:ncol(tildata)]

#making sure genes are numeric (0 or 2 and 1 for het)##########
tildata[,genes] <- lapply(tildata[,genes] , as.character)
tildata[,genes] <- lapply(tildata[,genes] , as.numeric)

#removing incomplete rows############
partial <- tildata[!complete.cases(tildata), ]
tils <- tildata[complete.cases(tildata), ]

#fitting model###########
fit.the.model = function(data, alpha) {
  cpc = vector(length(traits), mode = "list")
  names(cpc) = traits
    # find the parents of each trait (may be genes or other traits).
  for (t in seq_along(traits)) {
    
    # BLUP away the family structure.
    # No need to do this for RILs, BILs, etc.
    #m = lmer(as.formula(paste(traits[t], "~ (1|Rep)")), data = data)
    #data[!is.na(data[, traits[t]]), traits[t]] =
    #data[, traits[t]] - ranef(m)[[1]][as.character(data$Rep), 1]
    
    # Learn the Markov blanket for each trait; feature selection
    # Turn off debug for less output
    cpc[[t]] = learn.nbr(data[, c(traits, genes)], node = traits[t], debug = TRUE,
                         method = "si.hiton.pc", test = "cor", alpha = alpha)
  }#FOR
  
  # merge the relevant variables to use for learning.
  nodes = unique(c(traits, unlist(cpc)))
  # yield has no children, and genes cannot depend on traits.
  resultvar = "GW"
  traitlist=c(traittiers,list(resultvar))
  # Combine marker nodes and trait nodes
  # Markers can depend on each other (epistasis?  Feature not a bug.)
  modellist = c(list(nodes[!(nodes %in% traits)]),traitlist)
  blacklist = tiers2blacklist(modellist)
  
  # build the Bayesian network.
  bn = hc(data[, nodes], blacklist = blacklist)
  
  return(bn)
  
}#FIT.THE.MODEL

# k-fold cross-validation
xval.the.model = function(data, k = 10, cluster, alpha, ridge) {
  n = nrow(data)
  predcor = numeric(length(traits))
  names(predcor) = traits
  postcor = numeric(length(traits))
  names(postcor) = traits
  #added
  #observed = numeric(length(traits))
  #names(observed) = traits
  
  # shuffle the data to get unbiased splits.
  kcv = split(sample(n), seq_len(k))
  # store the length of each test set.
  kcv.length = sapply(kcv, length)

  # use parLapply to use clusters 
  #predicted = parLapply(kcv, cl = cluster, function(test) {
  predicted = lapply(kcv, function(test) {
    # create a matrix to store the predicted values.
    pred = matrix(0, nrow = length(test), ncol = length(traits))
    colnames(pred) = traits
    # create a matrix to store posterior estimates.
    post = matrix(0, nrow = length(test), ncol = length(traits))
    colnames(post) = traits
    
    cat("* beginning cross-validation fold.\n")
    
    # split training and test.
    dtraining = data[-test, ]
    dtest = data[test, ]

    # fit the model on the training data.
    model = fit.the.model(dtraining, alpha = alpha)
    fitted = bn.fit(model, dtraining[, nodes(model)])
    # maybe re-fit with ridge regression.
    if (ridge) {
      
      library(penalized)
      
      for (no in nodes(fitted)) {
        
        node.parents = parents(fitted, no)
        
        if (length(node.parents) < 3)
          next
        
        opt.lambda = optL2(response = dtraining[, no],
                           penalized = dtraining[, node.parents],
                           model = "linear", trace = FALSE,
                           minlambda2 = 10e-5, maxlambda = 500)$lambda
        fitted[[no]] = penalized(response = dtraining[, no],
                                 penalized = dtraining[, node.parents],
                                 model = "linear", trace = FALSE,
                                 lambda1 = 0, lambda2 = opt.lambda)
        
      }#FOR
      
    }#THEN
    
    # subset the test data.
    dtest = dtest[, nodes(model)]
    
    cat("  > model has", length(nodes(model)), "nodes.\n")
    
    # predict each trait in turn, given all the parents.
    for (t in traits) {
      pred[, t] = predict(fitted, node = t, data = dtest[, nodes(model)])
    }
    
    for (i in seq(nrow(dtest))) {
      #This will crash if no markers are included.  Adjusting alpha is one fix
      post[i, traits] = colMeans(cpdist(fitted, nodes = traits, evidence = as.list(dtest[i, names(dtest) %in% genes]), method = "lw", n = 1000))
      #post[i, traits] = colMeans(cpdist(fitted, nodes = traits, evidence = as.list(dtest[i, names(dtest)]), method = "lw", n = 1000))
    }
    return(list(model = fitted, pred = pred, post = post))
  })
  
  # merge all the predicted values.
  posterior = do.call(rbind, lapply(predicted, `[[`, "post"))
  causal = do.call(rbind, lapply(predicted, `[[`, "pred"))
  
  cat("* overall cross-validated correlations:\n")
  for (t in traits) {
    
    predcor[t] = cor(causal[, t], data[unlist(kcv), t])
    cat("  > PREDCOR(", t, "):", predcor[t], "\n")
    postcor[t] = cor(posterior[, t], data[unlist(kcv), t])
    cat("  > POSTCOR(", t, "):", postcor[t], "\n")
    
  }#FOR
  
  return(list(predicted = causal, posterior = posterior,
              #observed = data[unlist(kcv), t], #something seems wrong.  no t defined
              observed = as.matrix(data[unlist(kcv), traits]), #maybe this is right?
              predcor = predcor, postcor = postcor,
              models = lapply(predicted, `[[`, "model")))
  
}#XVAL.THE.MODEL



#cl = makeCluster(1)
#invisible(clusterEvalQ(cl, library(bnlearn)))
#invisible(clusterEvalQ(cl, library(lme4)))
#clusterExport(cl = cl, c("traits", "genes", "ids", "fit.the.model", "is.string.vector"))

pr001 = vector(10, mode = "list")

for (j in 1:10) {
  pr001[[j]] = xval.the.model(tils, cluster = cl, alpha = .001, ridge = TRUE)
}

#stopCluster(cl)

pred.summary = sapply(pr001, `[[`, "predcor")
print(rowMeans(pred.summary))
write.csv(pred.summary,file="Predsummary.csv")
post.summary = sapply(pr001, `[[`, "postcor")
print(rowMeans(post.summary))
write.csv(post.summary,file="PostSummary.csv")
arclist = list()
for (i in seq_along(pr001)) {
  
  run = pr001[[i]]$models
  for (j in seq_along(run))
    arclist[[length(arclist) + 1]] = arcs(run[[j]])
}#FOR


nodes = unique(unlist(arclist))

strength = custom.strength(arclist, nodes = nodes)
summary(pr001)
averaged = averaged.network(strength)
relevant.nodes = nodes(averaged)[sapply(nodes, degree, object = averaged) > 0]
averaged2 = subgraph(averaged, relevant.nodes)
strength2 = strength[(strength$from %in% relevant.nodes) & (strength$to %in% relevant.nodes), ]
tiff("C:/Users/~.tiff", height = 14, width = 16, units = 'cm', compression = "lzw", res = 1000)
gR = strength.plot(averaged2, strength2, shape = "rectangle", layout = "dot")
#gR = strength.plot(averaged2, strength2, shape = "rectangle", layout = "circo", threshold=.8)
#gR = strength.plot(averaged2, strength2, shape = "rectangle", layout = "fdp", threshold=.8)
#gR = strength.plot(averaged2, strength2, shape = "ellipse", layout = "dot")
# gR = strength.plot(averaged2, strength2, shape = "ellipse", layout = "fdp", threshold=.8)
#gR = strength.plot(averaged2, strength2, shape = "ellipse", layout = "fdp",threshold=.8)
#gR = strength.plot(averaged2, strength2, shape = "ellipse", layout = "fdp")
# add colors
nodeRenderInfo(gR)$fill = "lightblue"
nodeRenderInfo(gR)$fill = "lightblue"
nodeRenderInfo(gR)$col = "darkblue"
nodeRenderInfo(gR)$fill[traits] = "green"
nodeRenderInfo(gR)$col[traits] = "red"
a = arcs(subgraph(averaged, traits))
a = as.character(interaction(a[, "from"], a[, "to"], sep = "~"))
edgeRenderInfo(gR)$col = "lightgray"
edgeRenderInfo(gR)$col[a] = "darkgreen"
renderGraph(gR)
dev.off()

