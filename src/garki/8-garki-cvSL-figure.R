#-------------------------------
# 9-garki-cvSL-figure.R
#
# plot cross-validated risk
# from the super learner
# and consituent algorithms
#-------------------------------



#-------------------------------
# preamble
#-------------------------------
rm(list=ls())
library(RColorBrewer)


#-------------------------------
# load the results
#-------------------------------
load("~/dropbox/articles/antibody-curves/results/raw/garki-cvSL.RData")


#-------------------------------
# get summary estimates
#-------------------------------

# summary estimates for intervention round 5 (no covariates)
sum.tfit <- summary(cv.tfit)


# summary estimates for multivariate example (round 3-5)
sum.mfit <- summary(cv.mfit)


# plot the estimates
plot.slab_cvSL <- function(cvSL,...) {
  
  # get table of estimates
  res <- summary(cvSL)$Table
  
  # remove the discrete super learner
  res$Algorithm <- as.character(res$Algorithm)
  res <- subset(res,Algorithm!="Discrete SL")
  
  # for some of the commonly used algorithms, clean up the names
  res$Algorithm[res$Algorithm=="SL.mean_All"] <- "Simple Mean"
  res$Algorithm[res$Algorithm=="SL.glm_All"] <- "GLM"
  res$Algorithm[res$Algorithm=="SL.bayesglm_All"] <- "Bayes GLM"
  res$Algorithm[res$Algorithm=="SL.gam_All"] <- "GAM"
  res$Algorithm[res$Algorithm=="SL.loess_All"] <- "LOESS"
  res$Algorithm[res$Algorithm=="SL.polymars_All"] <- "MARS"
  res$Algorithm[res$Algorithm=="SL.glmnet_All"] <- "Lasso (glmnet)"
  res$Algorithm[res$Algorithm=="SL.Yman2016_All"] <- "Yman2016"
  res$Algorithm[grep("SL.randomForest",res$Algorithm)] <- "Random Forest"
  
  # sort by CV risk
  res <- res[order(res$Ave,decreasing=TRUE),]
  
  # plot
  op <- par(mar=c(10,4,4,2)+0.1)
  Midpts <- barplot(res$Ave,horiz=TRUE,col=NA,border=NA,names.arg=res$Algorithm,las=1,xlim=c(min(res$Ave-2*res$se),max(res$Ave+2*res$se)))
  arrows(x0=(res$Ave-2*res$se),x1=(res$Ave+2*res$se),y0=Midpts,angle=90,length=0.1,code=3)
  points(res$Ave,Midpts,pch=21,bg="white")
  
  par(op)
  
}


