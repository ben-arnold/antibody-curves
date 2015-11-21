#-------------------------------
# 1-miton-msp1-analysis
#
# Fit an age-specific antibody
# curve for MSP-1_19 in Miton, Haiti
# and estimate TMLE means by age group
#
# version 1 (21 nov 2015)
#-------------------------------

#-------------------------------
# input files:
#   haiti2-malaria-miton.dta
#
# output files:
#   miton-msp1-analysis.RData
#-------------------------------


#-------------------------------
# preamble
#-------------------------------

rm(list=ls())
library(SuperLearner)
library(tmle)
library(foreign)


# source the base functions for
# SL fits of age-antibody curves
# and TMLE estimates of mean differences
source("~/SLAbcurves/src/SLAb-curve.R")
source("~/SLAbcurves/src/SLAb-tmle.R")
source("~/SLAbcurves/src/SLAb-cvRF.R")

#-------------------------------------------
# load the Miton data
#-------------------------------------------

d <- read.dta("~/dropbox/haiti2/data/final/haiti2-malaria-miton.dta")

d$msp1 <- d$msp13d7

# for 6 negative MSP_1 values, replace as 0
d$msp1[d$msp1<0] <- 0

# create age categories for stratified analyses
# 5 year age categories (1-20 y)
d$agecat <- cut(d$age,breaks=c(0,5,10,15,20),labels=c("1-5","6-10","11-15","16-20"))

#-------------------------------------------
# estimate a marginal Ab curve E(Y_x,a)
#-------------------------------------------
set.seed(25234)
msp1.EYxa <- SLAb.curve(Y=log10(d$msp1+1),Age=d$age,id=d$id)

#-------------------------------------------
# estimate age-adjusted means E(Y_x)
# by 5 year age category
#-------------------------------------------
agegrps <-c("1-5","6-10","11-15","16-20") 
msp1.EYx <- sapply(agegrps, function(x) 
  SLAb.tmle(Y=log10(d$msp1[d$agecat==x]+1),Age=d$age[d$agecat==x],id=d$id[d$agecat==x]) 
)

#-------------------------------------------
# save results
#-------------------------------------------
save.image(file="~/SLAbcurves/results/raw/miton-msp1-analysis.RData")




