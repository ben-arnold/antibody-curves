






#--------------------------------------
# SL functions to model age-dependent
# antibody curves using the antibody
# acquisition model from Yman 2016
#--------------------------------------

# model function
SL.yman2016 <- function(Y,X,newX=NULL) {
  
  # objective function
  LL <- function(theta, X, Y) {
    # theta : parameters for optimization, length 3, including:
      # alpha : antibody acquisition rate
      # r     : antibody loss rate
      # sigma : sd of the log-normal distribution
    # X     : age of individual i
    # Y     : log antibody level of individual i
    
    alpha <- theta[1]
    r     <- theta[2]
    sigma <- theta[3]
    
    # get age-specific predicted geomean antibody level given alpha, r
    Aa <- (alpha/r)*(1-exp(-r*X))
    
    # likelihood function for log Ab levels
    L <- dnorm(Y,mean=log(Aa),sd=sigma)
    
    # negative log likelihood function (to minimize with optim)
    -sum(log(L))
  }
  
  # close out if X is more than 1 dimension
  X <- data.frame(X)
  if(ncol(X)>1) {
    warning("For the SL.yman2016 model, X can only include one variable (age)\n")
  }
  
  # ML fit given X and Y
  mlfit <- optim(c(0.1,0.01,1),fn=LL,X=X[,1],Y=Y)
  fit <- list(object = mlfit)
  class(fit) <- "SL.yman2016"
  
  # predicted values
  if(is.null(newX)) newX <- X[,1]
  pred <- predict(fit,newdata=newX)
 
  # return results
  return(list(pred=pred,fit=fit))
}


# prediction method for SL.yman2016
predict.SL.yman2016 <- function(object,newdata,...) 
  {
  pred <- (object$object$par[1]/object$object$par[2])*(1-exp(-object$object$par[2]*newdata[,1]))
  pred
}

#--------------------------------------
# testing below
#--------------------------------------
# # simulate some data
# set.seed(1123)
# age <- runif(1000,min=1,max=20)
# alpha <- 0.5
# r <- 0.07
# sigma <- log(1.7)
# abtrue <- (alpha/r)*(1-exp(-r*age))
# abobs <- exp( log(abtrue) +rnorm(length(age),mean=0,sd=sigma) )
# summary(abobs)
# 
# # ML fit
# yman2016fit <- SL.yman2016(Y=log(abobs),X=age)
# 
# 
# # compare against truth
# cbind(c(alpha,r,sigma),yman2016fit$fit$object$par)
# 
# plot(age,log(abobs),cex=0.25,col="gray40",pch=16)
# lines(age[order(age)],log(abtrue[order(age)]),lwd=3,lty=2)
# lines(age[order(age)],log(yman2016fit$pred[order(age)]),lwd=1)



