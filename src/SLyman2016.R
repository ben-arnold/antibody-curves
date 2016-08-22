




rm(list=ls())

#--------------------------------------
# SL functions to model age-dependent
# antibody curves using the antibody
# acquisition model from Yman 2016
#--------------------------------------
# SL model function
SL.yman2016 <- function(Y,X,newX=NULL,...) {
  
  # objective function
  LL <- function(theta, Age, logAb) {
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
    Aa <- (alpha/r)*(1-exp(-r*Age))
    
    # likelihood function for log Ab levels
    L <- dnorm(logAb,mean=log(Aa),sd=sigma)
    
    # negative log likelihood function (to minimize with optim)
    -sum(log(L))
  }
  
  if(is.null(newX)) {
    newX <- X
  }
  
  # throw a warning if X has more than 1 variable
  if(!is.null(ncol(X)) ) {
    #warning("For the SL.yman2016 model, X can only include one variable (age)\nAssuming that the first column in X is Age and discarding other information.")
    X <- X[,1]
    newX <- newX[,1]
  }
  
  # ML fit of L(theta | X, Y)
  mlfit <- optim(c(0.1,0.01,1),fn=LL,Age=X,logAb=Y)
  fit <- list(object = mlfit)
  class(fit) <- "SL.yman2016"
  
  # predicted values
  pred <- predict(fit,newdata=newX)
  
  # return results
  out <- list(pred=pred, fit=fit)
  return(out)
}


# prediction method for SL.yman2016
predict.SL.yman2016 <- function(object,newdata,...) {
  pred <-  log( (object$object$par[1]/object$object$par[2])*(1-exp(-object$object$par[2]*newdata)) )
  pred
}

#--------------------------------------
# testing below
#--------------------------------------
# simulate some data
set.seed(1123)
age <- runif(1000,min=1,max=20)
alpha <- 0.5
r <- 0.07
sigma <- log(1.7)
abtrue <- (alpha/r)*(1-exp(-r*age))
abobs <- exp( log(abtrue) +rnorm(length(age),mean=0,sd=sigma) )
summary(abobs)

# ML fit
yman2016fit <- SL.yman2016(Y=log(abobs),X=age)

# compare against truth
cbind(c(alpha,r,sigma),yman2016fit$fit$object$par)
plot(age,log(abobs),cex=0.25,col="gray40",pch=16)
lines(age[order(age)],log(abtrue[order(age)]),lwd=3,lty=2)
lines(age[order(age)],yman2016fit$pred[order(age)],lwd=1)

# loessfit <- SL.loess(Y=log(abobs),X=data.frame(age),family=gaussian(),obsWeights=rep(1,length(abobs)),newX=data.frame(age))
# lines(age[order(age)],loessfit$pred[order(age)],lwd=1,col="blue")

#--------------------------------------
# super learner fit
#--------------------------------------
require(SuperLearner)

# works fine with others
SL.library <- c("SL.glm","SL.loess","SL.gam")
SLfit <- SuperLearner(Y=log(abobs),X=data.frame(age=age),SL.library=SL.library)
SLfit

# but not with the SL.yman2016
SL.library <- c("SL.glm","SL.loess","SL.gam","SL.yman2016")
SLfit2 <- SuperLearner(Y=log(abobs),X=data.frame(age=age),SL.library=SL.library)
SLfit2

# ... debugging: there is something about this call :
#
# Z[unlist(validRows, use.names = FALSE), ] <- do.call("rbind", 
# lapply(validRows, FUN = .crossValFUN, Y = Y, dataX = X, id = id, 
#        obsWeights = obsWeights, library = library, kScreen = kScreen, 
#        k = k, p = p, libraryNames = libraryNames))
#
# that converts the Z matrix of results (1000 x 4) into a list of length 4000
# it then throws an error in the next line because Z no longer has dimension
# errorsInCVLibrary <- apply(Z, 2, function(x) any(is.na(x)))

# testAlg <- try(do.call("SL.yman2016", 
#             list(Y = log(abobs), X = data.frame(age=age), family = gaussian(), 
#                  id = 1:length(abobs), obsWeights = rep(1,length(abobs))) ) )
# 
# out <- matrix(NA,nrow=1000,ncol=4)
# 
# out[,4] <- testAlg$pred

