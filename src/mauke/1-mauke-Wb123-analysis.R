

#-------------------------------------------
# 1-mauke-Wb123-analysis.R
# Ben Arnold
#
# SuperLearner fits of age-specific antibody
# response on Mauke island before and after
# MDA, and TMLE estimates of geometric mean
# (area under the curve) of antibody responses
#
#
# version 1 (15 oct 2015)
#-------------------------------------------

#-------------------------------------------
# input files:
#   mauke1974.dta
#   mauke1992.dta
#
# output files:
#   mauke-Wb123-analysis.RData
#-------------------------------------------



#-------------------------------------------
# preamble
#-------------------------------------------

rm(list=ls())
library(foreign)
library(tmle)
library(SuperLearner)

# source the base functions for
# SL fits of age-antibody curves
# and TMLE estimates of mean differences
source("~/SLAbcurve/src/0-SLAb-base-functions.R")

#-------------------------------------------
# load the Mauke data from 1974(1975) and 1992
#-------------------------------------------

d75 <- read.dta("~/dropbox/mauke/data/final/mauke1974.dta")
d92 <- read.dta("~/dropbox/mauke/data/final/mauke1992.dta")

# drop 7 children aged 0 (all Wb123+) in 1975 due to maternal antibodies
a75 <- subset(d75,age>0)
a92 <- d92

# add 0.5 years to age to remove bias (on average) due to rounding to year
a75$ager <- a75$age+0.5
a92$ager <- a92$age+0.5

# create age categories for stratified analyses
# a75$agecat <- cut(a75$age,breaks=c(0,2,4,6,8,10),labels=c("1-2","3-4","5-6","7-8","9-10"))
# a92$agecat <- cut(a92$age,breaks=c(0,2,4,6,8,10),labels=c("1-2","3-4","5-6","7-8","9-10"))
a75$agecat <- cut(a75$age,breaks=c(0,5,10,15,20),labels=c("1-5","6-10","11-15","16-20"))
a92$agecat <- cut(a92$age,breaks=c(0,5,10,15,20),labels=c("1-5","6-10","11-15","16-20"))


# create a standard ID variable before appending
a75$id <- as.integer(a75$id74)
a92$id <- ifelse(is.na(a92$id74),a92$id92,a92$id74)

# identify the pre- vs. post-MDA measurements
a75$mda <- 0
a92$mda <- 1

# append data with common variables
common.vars <- c("id","ager","agecat","wb123","mda")
a7592 <- rbind(subset(a75,select=common.vars),subset(a92,select=common.vars))

#--------------------------------------
# All Ages
#--------------------------------------
set.seed(0237234)

# SuperLearner fits of antibody levels
mauke75 <- SLAb.curve(Y=log10(a75$wb123),Age=a75$ager,id=a75$id)
mauke92 <- SLAb.curve(Y=log10(a92$wb123),Age=a92$ager,id=a92$id)

# estimate group means
EYx.mauke75 <- SLAb.tmle(Y=log10(a75$wb123),Age=a75$ager,id=a75$id)
EYx.mauke92 <- SLAb.tmle(Y=log10(a92$wb123),Age=a92$ager,id=a92$id)

# estimate difference in means
diff.mauke  <- SLAb.tmle(Y=log10(a7592$wb123),Age=a7592$ager,id=a7592$id,X=a7592$mda,diff=TRUE)


#--------------------------------------
# Estimate means and differences between
# time points in
# 2 year age bands from ages 1-10
#--------------------------------------
# agegrps <-c("1-2","3-4","5-6","7-8","9-10") 
agegrps <-c("1-5","6-10","11-15","16-20") 
EYx.mauke75kids <- sapply(agegrps, function(x) 
	SLAb.tmle(Y=log10(a75$wb123[a75$agecat==x]),Age=a75$ager[a75$agecat==x],id=a75$id[a75$agecat==x]) 
	)
EYx.mauke92kids <- sapply(agegrps, function(x) 
	SLAb.tmle(Y=log10(a92$wb123[a92$agecat==x]),Age=a92$ager[a92$agecat==x],id=a92$id[a92$agecat==x]) 
	)
diff.maukekids <- sapply(agegrps, function(x) 
	SLAb.tmle(Y=log10(a7592$wb123[a7592$agecat==x]), Age=a7592$ager[a7592$agecat==x], id=a7592$id[a7592$agecat==x], X=a7592$mda[a7592$agecat==x], diff=TRUE) 
	)



#--------------------------------------
# store results for later summary
# and plotting
#--------------------------------------
save.image("~/SLAbcurve/results/raw/mauke-Wb123-analysis.RData")




