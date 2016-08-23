
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
source("~/SLAbcurves/src/plot_slab_cvSL.R")


#-------------------------------
# load data
#-------------------------------
d <- read.csv("~/dropbox/articles/antibody-curves/data/miton/haiti2-malaria-miton-public.csv")

# for 6 negative MSP_1 values, replace as 0
d$msp1 <- d$msp13d7
d$msp1[d$msp1<0] <- 0


  
#-------------------------------
# fit cross-validated SL
#-------------------------------

SL.library <- c("SL.mean","SL.glm","SL.loess","SL.gam","SL.randomForest","SL.Yman2016")
# "SL.glmnet"


# fit the cross-validated super learner, with just Age as the predictor (in Round 5)
set.seed(32423)
cv.mitonfit <- slab_cvSL(Y=log10(d$msp1+1),Age=d$age,family=gaussian(),V=10,SL.library=SL.library)


#------------------------------
# print results to log
#------------------------------
summary(cv.mitonfit)

#-------------------------------
# plot the estimates
#-------------------------------
pdf("~/dropbox/articles/antibody-curves/results/figs/miton-cvSL.pdf")
plot_slab_cvSL(cv.mitonfit,title="P. falciparum, low transmission (age only)")
dev.off()

#-------------------------------
# save down the results
#-------------------------------
save(cv.mitonfit,file="~/dropbox/articles/antibody-curves/results/raw/miton-cvSL.RData")



