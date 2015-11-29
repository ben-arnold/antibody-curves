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
table(d$agecat)
agegrps <-c("1-5","6-10","11-15","16-20") 
msp1.EYx <- sapply(agegrps, function(x) 
  SLAb.tmle(Y=log10(d$msp1[d$agecat==x]+1),Age=d$age[d$agecat==x],id=d$id[d$agecat==x]) 
)

#-------------------------------------------
# randomly downsample the population to
# N=30, 25, 20, 15, 10 per 5-year age group
#------------------------------------------- 
strat.samp <- function(x,n) {
  d1 <- x[x$agecat=="1-5", ][sample(1:nrow(x[x$agecat=="1-5", ]),n),]
  d2 <- x[x$agecat=="6-10",][sample(1:nrow(x[x$agecat=="6-10",]),n),]
  d3 <- x[x$agecat=="11-15",][sample(1:nrow(x[x$agecat=="11-15",]),n),]
  d4 <- x[x$agecat=="16-20",][sample(1:nrow(x[x$agecat=="16-20",]),n),]
  dd <- rbind(d1,d2,d3,d4)
  return(dd)
}

set.seed(327234)
dlt20 <- subset(d,age<=20)
d30 <- strat.samp(dlt20,30)
d25 <- strat.samp(dlt20,25)
d20 <- strat.samp(dlt20,20)
d15 <- strat.samp(dlt20,15)
d10 <- strat.samp(dlt20,10)


#-------------------------------------------
# re-estimate the curves and 
# age-stratified means with downsampled
# data to see if they deteriorate
#-------------------------------------------

set.seed(33141)
msp1.EYxa.30 <- SLAb.curve(Y=log10(d30$msp1+1),Age=d30$age,id=d30$id)
msp1.EYxa.25 <- SLAb.curve(Y=log10(d25$msp1+1),Age=d25$age,id=d25$id)
msp1.EYxa.20 <- SLAb.curve(Y=log10(d20$msp1+1),Age=d20$age,id=d20$id)
msp1.EYxa.15 <- SLAb.curve(Y=log10(d15$msp1+1),Age=d15$age,id=d15$id)
msp1.EYxa.10 <- SLAb.curve(Y=log10(d10$msp1+1),Age=d10$age,id=d10$id)


# use TMLE to estimate group means for each subsample
agegrps <-c("1-5","6-10","11-15","16-20")
agecatmu <- function(data) {
  print(table(data$agecat))
  msp1.EYx <- sapply(agegrps, function(x) 
    SLAb.tmle(Y=log10(data$msp1[data$agecat==x]+1),Age=data$age[data$agecat==x],id=data$id[data$agecat==x]) 
  )
  return(msp1.EYx)
}
msp1.EYx.30 <- agecatmu(d30)
msp1.EYx.25 <- agecatmu(d25)
msp1.EYx.20 <- agecatmu(d20)
msp1.EYx.15 <- agecatmu(d15)
msp1.EYx.10 <- agecatmu(d10)

#-------------------------------------------
# save results
#-------------------------------------------
save.image(file="~/SLAbcurves/results/raw/miton-msp1-analysis.RData")




