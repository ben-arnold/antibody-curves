#-------------------------------
# 0-SLAb-base-functions.R
#
# ben arnold (benarnold@berkeley.edu)
#
# base functions for the analysis
# of age-specific antibody 
# response curves
#-------------------------------



#-------------------------------
# Wrapper function:
# fit age-specific antibody curves
# by country using SuperLearner
# E(Y_x,a) = E(Y | X=x, A=a)
#-------------------------------

SLAb.curve <-function(Y,Age,id) {
	# wrapper function for the SuperLearner algorithm, subset to a specific
	# dataset or subgroup (X=x) of interest
	#
	# generates predicted means of Y at each age (A=a) in the dataset
	# Y   : outcome_i (log10 antibody)
	# Age : age_i
	# id  : ID variable for individuals in the dataset (for repeated observations)
	SLlib <- c("SL.mean","SL.glm","SL.bayesglm","SL.loess","SL.gam","SL.glmnet")
	# note that the W matrix includes a row of 1s to get the SL.glmnet algorithm to run w/o additional covariates
	# this will make the SL.glm module throw warnings (due to a rank-deficient model). this is not a problem
	
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

#-------------------------------
# wrapper function to call TMLE 
# for antibody response data
# calculating either the mean
# within group X (binary 1/0) or the
# difference between groups 
# if diff=FALSE, estimates psi = E(Y_x)
# if diff=TRUE, estimates  psi = E[ E(Y_1) - E(Y_0) ]
#-------------------------------
SLAb.tmle <- function(Y,Age,X=NULL,W=NULL,id,diff=FALSE) {
	# Y    : outcome
	# Age  : age
	# X    : comparison group (must be binary, 0/1, for tmle()). 
	#        can be null if diff=FALSE
	# W    : matrix of covariates (optional)
	# id   : ID variable, unique for each individual (to identify repeated obs)
	# diff : logical. calculate the difference between treatment groups?
	#           if FALSE (default) it returns the mean
	
	SLlib <- c("SL.mean","SL.glm","SL.bayesglm","SL.loess","SL.gam","SL.glmnet")
	
	# if X is null, create a row of 1s
	if (is.null(X)) {
		X <- rep(1,length(Y))
	}
	# if W is null, create a row of 1s so that the SL.glmnet algorithm will run
	if (is.null(W)) {
		W <- rep(1,length(Y))
	}

	# extract objects to make the calculations easier
	# subset to non-missing data for TMLE fit
	fitd <- data.frame(Y,X,Age,W,id)
	fitd <- fitd[complete.cases(fitd),]
	# estimate either the difference (A=fitd$X) or the adjusted mean (A=NULL)
	if (diff==TRUE) {
		tmle.fit <- tmle(Y=fitd$Y,
			A=fitd$X,
			W=data.frame(fitd$Age,fitd$W),
			id=fitd$id,
			Q.SL.library=SLlib
		)
		print(tmle.fit)
		psi  <- tmle.fit$estimates$ATE$psi
		se   <- sqrt(tmle.fit$estimates$ATE$var.psi)
		ci   <- tmle.fit$estimates$ATE$CI
		p    <- tmle.fit$estimates$ATE$pvalue
	} else {
		tmle.fit <- tmle(Y=fitd$Y,
			A=NULL,
			W=data.frame(fitd$Age,fitd$W),
			id=fitd$id,
			Q.SL.library=SLlib
		)
		print(tmle.fit)
		psi  <- tmle.fit$estimates$EY1$psi
		se   <- sqrt(tmle.fit$estimates$EY1$var.psi)
		ci   <- tmle.fit$estimates$EY1$CI
		p    <- tmle.fit$estimates$EY1$pvalue
	}
	# return estimate, SE, 95% CI, and P-value
	list(psi=psi,se=se,lb=ci[1],ub=ci[2],p=p)
}






