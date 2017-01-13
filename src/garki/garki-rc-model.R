

#-------------------------------
# garki-rc-model.R
#
# try to fit seroprevalence
# for Pf IFA using a reversible
# catalytic model
#
#-------------------------------



#-------------------------------
# preamble
#-------------------------------
rm(list=ls())

# source the base functions for
# SL fits of age-antibody curves
# and TMLE estimates of mean differences
source("~/SLAbcurves/src/SLAb-curve.R")
source("~/SLAbcurves/src/SLAb-tmle.R")
source("~/SLAbcurves/src/SLAb-cvRF.R")


#--------------------------------------
# estimate seroconversion rate using
# a reversible catalytic model
#--------------------------------------
Lrcm <- function(theta,data) {
  # theta : vector of length two, including h, r (parameters to be maximized)
  # data  : data frame with 1 row per age group. cols = age / n / k
  h <- rep(theta[1],nrow(data))
  r <- rep(theta[2],nrow(data))
  t <- data[,1]
  n <- data[,2]
  k <- data[,3]
  p <- h/(r+h)*(1-exp(-(r+h)*t))
  # negative log likelihood function (to minimize with optim)
  sum( - (k)*log(p) - (n-k)*log(1-p) )
}

# reversible catalytic model probability function
rcp <- function(h,r,t) h/(r+h)*(1-exp(-(r+h)*t))

# function to get CIs for a and b from the max likelihood fit (2 parameters)
hrci <- function(llobj) {

  I <- solve(llobj$hessian)
  ahat <- llobj$par[1]
  ahatse <- sqrt(I)[1]
  ahatlb <- ahat-1.96*ahatse
  ahatub <- ahat+1.96*ahatse

  bhat <- llobj$par[2]
  bhatse <- sqrt(I)[2]
  bhatlb <- bhat-1.96*bhatse
  bhatub <- bhat+1.96*bhatse

  res <- cbind(c(ahat,ahatse,ahatlb,ahatub),c(bhat,bhatse,bhatlb,bhatub) )
  rownames(res) <- c("est","se","lb","ub")
  colnames(res) <- c("h","r")
  return(res)
}

#-------------------------------
# load the serology dataset
#-------------------------------
d <- read.csv("~/dropbox/articles/antibody-curves/data/garki/final/garki-sero.csv")

d$mdate <- as.Date(d$mdate,"%d %b %Y")


# add 2 control village names
d$vname <- factor(d$vname,levels=c(levels(d$vname),"Nabanawa","Ajura"))
d$vname[d$village==552] <- "Nabanawa"
d$vname[d$village==553] <- "Ajura"
d$vname <- factor(d$vname)

# identify sero-positive individuals (Table 21 of Molineaux 1980, p 175)
d$ifapfpos <- ifelse(d$ifatpftitre>=20,1,0)

# subset to non-missing values and just a few variables
ad <- subset(d,!is.na(ifapfpos),select=c('tr','serosvy','ageyrs','ifapfpos'))

ad$n <- rep(1,nrow(ad))

fitC <- optim(c(0.3,0.01),fn=Lrcm,data=ad[ad$tr=="Control" & ad$serosvy==5,c("ageyrs","n","ifapfpos")],method="L-BFGS-B", lower=c(0.00001,0.0000001),upper=c(100,100),hessian=TRUE)

fitI <- optim(c(4,0.01),fn=Lrcm,data=ad[ad$tr=="Intervention" & ad$serosvy==5 ,c("ageyrs","n","ifapfpos")],method="L-BFGS-B", lower=c(0.00001,0.0000001),upper=c(4,100),hessian=TRUE)

# brainstorming ideas here
# can get wildy different results, depending on the constraints placed on the model. non-unique solutions!!!
# repeat with an upper bound of 4 and of 10.  Show how they get the same predicted fits, which actually don't even fit the data that well. and
# result in fitted parameters that make no sense at all.
# limit comparison to control and intervention villages at survey 5 to make a point

hrC <- hrci(fitC)
hrI <- hrci(fitI)
ts <- seq(1,70,by=0.1)
psC <- rcp(h=hrC[1,1],r=hrC[1,2],t=ts)
psI <- rcp(h=hrI[1,1],r=hrI[1,2],t=ts)

plot(ts,psC,type="l",ylim=c(0,1),bty="l",las=1)
lines(ts,psI,col="blue")




