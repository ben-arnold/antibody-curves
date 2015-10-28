#---------------------------------------
# SLAb-cvRF.R
#
# ben arnold (benarnold@berkeley.edu)
#
# select optimal tree depth for
# randomForest() using cross-validation
# by tuning the nodesize parameter
#
# if you don't do this, RF will usually
# grow the trees that are too deep 
# (e.g., to nodes with 1 obs)
# and this will overfit the data
#
# In the context of Age-antibody curves,
# it will make extremely jagged curves
# that are clear overfits. This adds 
# the one-step optimization of that parameter
# that should probably be included in the
# default implementation of the
# randomForest package, but isn't
#
# this function is designed to be 
# called by SLAb.curve or SLAb.tmle
# if they include SL.randomForest in the
# ensemble library
# it leverages SuperLearner() to estimate
# cross-validated Risks
#---------------------------------------


SLAb.cvRF <- function(Y,X,id,SLlib) {
  # arguments
  # Y     : outcome
  # X     : matrix of features that predict Y
  # id    : id variable to identify independent observations for V-splits
  # SLlib : SuperLearner library
  #
  # returns an updated SuperLearner library

  cat("\nThe ensemble library includes SL.randomForest.")
  cat("\nThe default R implementation of randomForest tends to overfit the data")
  cat("\nby growing trees that are too deep.")
  cat("\nWe can prevent overfitting by properly tuning the node size.")
  cat("\nSelecting the optimal node size (tree depth)")
  cat("\nfrom 5,10,15,...,40 using V-fold cross-validation.")
  cat("\nThis could take a few minutes, depending on the size of your dataset...\n")
  nodesizes <- seq(5,40,by=5)
  create.SL.randomForest <- function(tune = list(nodesize = nodesizes)) {
    for(mm in seq(length(tune$nodesize))) { 
      eval(parse(file = "", text = paste("SL.randomForest.ns", tune$nodesize[mm], "<- function(...,nodesize = ", tune$nodesize[mm], ") SL.randomForest(..., nodesize = nodesize)", sep = "")), envir = .GlobalEnv)
    }
    invisible(TRUE)
  }
  create.SL.randomForest()
  
  # estimate cross-validated Risks for different node sizes
  cvRisks <- rep(NA,length(nodesizes))
  for(nn in seq(length(nodesizes))) {
    fit <- SuperLearner(Y=Y,X=X,id=id,SL.library=paste("SL.randomForest.ns",nodesizes[nn],sep=""))
    cvRisks[nn] <- fit$cvRisk
  }
  # identify the lowest risk
  ns.opt <- nodesizes[order(cvRisks)][1]
  cvr.tab <- cbind(nodesizes,cvRisks)
  colnames(cvr.tab) <- c("nodesize","CVRisk")
  
  # update the library
  SLlib2 <- gsub("SL.randomForest",paste("SL.randomForest.ns",ns.opt,sep=""),SLlib)
  
  cat("\n-----------------------------------")
  cat("\nThe optimal node size is: ",ns.opt,"\n")
  print(cvr.tab)
  cat("-----------------------------------\n")
  
  # return results
  return(list(SLlib=SLlib2,ns.opt=ns.opt,cvRisks=cvr.tab))
  
}








