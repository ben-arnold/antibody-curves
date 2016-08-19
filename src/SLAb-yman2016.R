






#--------------------------------------
# estimate antibody acquisition rate using
# model in Yman 2016
#--------------------------------------

# objective function
yman2016LL <- function(theta, age, ab) {
  # theta : parameters for optimization, length 3, including:
    # alpha : antibody acquisition rate
    # r     : antibody loss rate
    # sigma : sd of the log-normal distribution
  # age   : age of individual i
  # ab    : log antibody level of individual i

  alpha <- theta[1]
  r     <- theta[2]
  sigma <- theta[3]

  # get age-specific predicted geomean antibody level given alpha, r
  Aa <- (alpha/r)*(1-exp(-r*age))

  # likelihood function for log Ab levels
  L <- dnorm(ab,mean=log(Aa),sd=sigma)

  # negative log likelihood function (to minimize with optim)
  -sum(log(L))
}

# model function
yman2016 <- function(age,Y) {
  mlfit <- optim(c(0.1,0.01,1),fn=yman2016LL,age=age,ab=Y)
  par <- mlfit$par
  pred <- (par[1]/par[2])*(1-exp(-par[2]*age))
  cat("\n-----------------------------------\n","Model parameters:\n","Antibody acquisition rate (alpha):",par[1],"\n Antibody loss rate (r):",par[2],"\n SD on log scale:",exp(par[3]), "\n-----------------------------------\n")
  return(list(par=par,pred=pred,age=age,Y=Y))

}



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
fit <- yman2016(age,log(abobs))

# compare against truth
cbind(c(alpha,r,sigma),fit$par)

# predicted values from ML fit
alphahat <- fit$par[1]
rhat <- fit$par[2]
pAb <- (alphahat/rhat)*(1-exp(-rhat*age))

plot(age,log(abobs),cex=0.25,col="gray40",pch=16)
lines(age[order(age)],log(abtrue[order(age)]),lwd=3,lty=2)
lines(age[order(age)],log(pAb[order(age)]),lwd=1)



