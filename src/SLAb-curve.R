#-------------------------------
# SLAb-curve.R
#
# ben arnold (benarnold@berkeley.edu)
#
# general wrapper function
# to estimate a age-dependent 
# antibody response curves 
# using machine learning
# prediction with the SuperLearner
# ensemble machine learning algorithm
#
#-------------------------------

#---------------------------------
# functions from 
# https://github.com/ecpolley/SuperLearnerExtra/tree/master/SL
# to create new SL GAM
# algorithms with different parameters
# (higher degrees of 3) 
#---------------------------------
create.SL.gam <- function(deg.gam = c(3, 4)) {
  for(mm in seq(length(deg.gam))){
    eval(parse(text = paste('SL.gam.', deg.gam[mm], '<- function(..., deg.gam = ', deg.gam[mm], ') SL.gam(..., deg.gam = deg.gam)', sep = '')), envir = .GlobalEnv)
  }
  invisible(TRUE)
}

create.SL.gam()

#-------------------------------
# Wrapper function:
# fit age-specific antibody curves
# by country using SuperLearner
# E(Y_x,a) = E(Y | X=x, A=a)
#-------------------------------

SLAb.curve <-function(Y,Age,W=NULL,id,SLlib= c("SL.mean","SL.glm","SL.bayesglm","SL.loess", "SL.gam","SL.gam.3","SL.glmnet","SL.randomForest")) {
	# wrapper function for the SuperLearner algorithm, subset to a specific
	# dataset or subgroup (X=x) of interest
	#
	# generates predicted means of Y at each age (A=a) in the dataset
	# Y   : outcome_i (log10 antibody)
	# Age : age_i
  # W    : matrix of covariates (optional)
	# id  : ID variable for individuals in the dataset (for repeated observations)
	# SLlib : SuperLearner algorithm library (if different from the default, above)
	# 
	# the function returns a data frame with id, Y, Age, W, and pY (SL predictions of marginally averaged Y at A=a)
	
	require(SuperLearner)
  
	# if W is null, create a row of 1s so that the SL.glmnet algorithm will run
	# this will make the SL.glm library throw warnings (due to a rank-deficient model). this is not a problem
  # create a warning handler to muffle that specific warning
  muffw <- function(w) if( any( grepl( "prediction from a rank-deficient fit may be misleading", w) ) ) invokeRestart( "muffleWarning" )
  nullW <- FALSE
	if (is.null(W)) {
	  Wdesign <- rep(1,length(Y))
	  nullW <- TRUE
	} else{
	  # convert W into a design matrix (SuperLearner does not acommodate factor variables)
	  Wdesign <- design.matrix(W)
	}
  
	# restrict dataset to non-missing observations
	fitd <- data.frame(id,Y,Age,Wdesign)
	fitd <- fitd[complete.cases(fitd),]
	X <- subset(fitd,select=-c(1:2) )
	
  # If randomForest is included in the library, 
  # select optimal node size (tree depth) using cross-validated risk
  # and then update the ensemble library to include the optimal node size
  if (grep("SL.randomForest",SLlib)>0) {
    cvRF <- SLAb.cvRF(Y=fitd$Y,X=X,id=fitd$id,SLlib=SLlib)
    SLlib <- cvRF$SLlib
  }
	
	# Fit SuperLearner
	SL.fit <- withCallingHandlers( SuperLearner(Y=fitd$Y,X=X,id=fitd$id,SL.library=SLlib), warning = muffw)
	print(SL.fit)
	
	# obtain marginally averaged values of Y at A=a:  E_W[E(Y|A=a, X=x, W)] 
	# with X=x implied by the subset of data used to fit the function
	# this parsing below is just to speed up computation in the case of
	# no adjustment covariates W=NULL
	# if W is not null, then need to do marginal averaging at each age
	# and merge that back to the analysis data frame, which is slower
	if(nullW==TRUE) {
	  res <- fitd[,1:3]
	  res$pY <- withCallingHandlers(predict(SL.fit)$pred, warning = muffw)
	} else{
	  As <- unique(X$Age)
	  pY <- rep(NA,length(As))
	  for(i in 1:length(As)) {
	    X$Age <- As[i]
	    pYs <- withCallingHandlers( predict(SL.fit,newdata=X)$pred, warning = muffw)
	    pY[i] <- mean(pYs)
	  }
	  res <- merge(fitd,data.frame(Age=As,pY=pY),by="Age",all.x=T,all.y=T)
	}
	
  res <- res[order(res$Age),]
	return(res)
}



# --------------------------------------
# automatic transform of a covariate
# data.frame with factors into one
# with indicator variables (and an 
# ommitted category)
# --------------------------------------
design.matrix <- function(W) {
  # W : data frame of covariates that might include factors
  if(class(W)!="matrix" & class(W)!="data.frame"){
    #cat("\n-----------------------------------------\nThe design matrix you supplied is not a matrix or a data.frame\nAssuming that it is a single variable\n-----------------------------------------\n")
    W <- data.frame(W)
  }
  ncolW <- ncol(W)
  flist <- numeric()
  for(i in 1:ncolW) {
    if(class(W[,i])!="factor"){
      next
    } else {
      flist <- c(flist,i)
      # strip out extra levels
      W[,i] <- factor(W[,i])
      # create a design matrix, remove the first level
      mm <- model.matrix(~-1+W[,i])
      mW <- mm[,-c(1)]
      # format the names of the indicator variables
      # and add them to the design matrix
      levs <- gsub(" ","",levels(W[,i]) )[-c(1)]
      if(length(levs)<2) mW <- matrix(mW,ncol=1)
      colnames(mW) <- paste(names(W)[i],levs,sep="")
      W <- data.frame(W,mW)
    }
  }
  # now drop the factors that have been replaced by indicators
  W <- W[,-c(flist)]
  return(W)
}


