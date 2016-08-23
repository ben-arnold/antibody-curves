
#-------------------------------
# 8-garki-cvSL.R
#
# Compute the cross-validated
# risk for the super leaner
# and its constituent algorithms
#
# do calculations for survey 5
# for both control and intervention
# groups
#
#-------------------------------



#-------------------------------
# preamble
#-------------------------------
rm(list=ls())
library(SuperLearner)
library(tmle)

# source the base functions for
# SL fits of age-antibody curves
# and TMLE estimates of mean differences
source("~/SLAbcurves/src/SLAb-curve.R")
source("~/SLAbcurves/src/SLAb-tmle.R")
source("~/SLAbcurves/src/SLAb-cvRF.R")
source("~/SLAbcurves/src/SL.Yman2016.R")
source("~/SLAbcurves/src/slab_cvSL.R")


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

# set sex to factor
d$sex <- as.factor(d$sex)

# for age exactly equal to 0, set it equal to 0.001
# to prevent the Yman 2016 model from blowing up
d$ageyrs[d$ageyrs<=0] <- 0.001


#-------------------------------
# Serological survey timing:
# 1-2 : pre-intervention
# 3-5 : intervention period
# 6-8 : post-intervention
#-------------------------------

# subset the data by group and survey round for convenience
# limit to survey round 5 (end of intervention period)
d.c5 <- d[d$tr=="Control" & d$serosvy==5,]
  d.c5 <- subset(d.c5,!is.na(d.c5$ageyrs) & !is.na(d.c5$ifatpftitre))
d.tr5 <- d[d$tr=="Intervention" & d$serosvy==5,]
  d.tr5 <- subset(d.tr5,!is.na(d.tr5$ageyrs) & !is.na(d.tr5$ifatpftitre))

  
#-------------------------------
# fit cross-validated SL
#-------------------------------

SL.library <- c("SL.mean","SL.glm","SL.bayesglm","SL.loess","SL.gam","SL.randomForest","SL.Yman2016")
# "SL.glmnet"


# fit the cross-validated super learner, with just Age as the predictor (in Round 5)
set.seed(23752)
cv.cfit <- slab_cvSL(Y=log10(d.c5$ifatpftitre+1),Age=d.c5$ageyrs,family=gaussian(),V=10,SL.library=SL.library)
cv.tfit <- slab_cvSL(Y=log10(d.tr5$ifatpftitre+1),Age=d.tr5$ageyrs,family=gaussian(),V=10,SL.library=SL.library)


# fit the cross-validated super learner, across all post-intervention survey rounds (3-5)
# accounting for other covariates
Y <- d$ifatpftitre[d$serosvy>=3|d$serosvy<=5]
Age <- d$ageyrs[d$serosvy>=3|d$serosvy<=5]
W <- d[d$serosvy>=3|d$serosvy<=5,c("tr","wetseason","sex","vname")]
cv.mfit <- slab_cvSL(Y=log10(Y+1),Age=Age,W=W,family=gaussian(),V=10,SL.library=SL.library)


#------------------------------
# print results to log
#------------------------------
# round 5, control villages, age only
summary(cv.cfit)

# round 5, intervention villages, age only
summary(cv.tfit)

# rounds 3-5, intervention + control villages, multivariate
summary(cv.mfit)

#-------------------------------
# save down the results
#-------------------------------
save(cv.cfit,cv.tfit,cv.mfit,file="~/dropbox/articles/antibody-curves/results/raw/garki-cvSL.RData")



