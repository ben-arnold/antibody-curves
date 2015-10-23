#-------------------------------
# SLAb-curve.R
#
# ben arnold (benarnold@berkeley.edu)
#
# general wrapper function
# to estimate a machine learning
# prediction using the SuperLearner
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

SLAb.curve <-function(Y,Age,id,SLlib= c("SL.mean","SL.glm","SL.bayesglm","SL.loess", "SL.gam","SL.gam.3","SL.glmnet")) {
	# wrapper function for the SuperLearner algorithm, subset to a specific
	# dataset or subgroup (X=x) of interest
	#
	# generates predicted means of Y at each age (A=a) in the dataset
	# Y   : outcome_i (log10 antibody)
	# Age : age_i
	# id  : ID variable for individuals in the dataset (for repeated observations)
	# SLlib : SuperLearner algorithm library (if different from the default, above)
	# 
	# the function returns a data frame with Y, Age, id, and pY (SL predictions)
	
	require(SuperLearner)
	
	
	# note that the W matrix includes a row of 1s to 
	#   get the SL.glmnet algorithm to run w/o additional covariates
	# this will make the SL.glm module throw warnings 
	#   (due to a rank-deficient model). this is not a problem
	
	# restrict dataset to non-missing observations
	fitd <- data.frame(Y,Age,id)
	fitd <- fitd[complete.cases(fitd),]
	
	# Fit SuperLearner
	SL.fit <-SuperLearner(Y=fitd$Y,X=data.frame(fitd$Age,rep(1,nrow(fitd))),SL.library=SLlib,id=fitd$id)
	print(SL.fit)
	
	# return predicted values of Y at A=a:  E(Y|A=a, X=x) with X=x implied by
	# the subset of data used to fit the function
	fitd$pY <- predict(SL.fit)$pred
	return(fitd[order(fitd$Age),])
}





