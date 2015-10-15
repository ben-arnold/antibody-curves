

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

# create a standard ID variable before appending
a75$id <- as.integer(a75$id74)
a92$id <- ifelse(is.na(a92$id74),a92$id92,a92$id74)

# identify the year
a75$mda <- 0
a92$mda <- 1

# append data with common variables
common.vars <- c("id","age","wb123","mda")
a7592 <- rbind(subset(a75,select=common.vars),subset(a92,select=common.vars))


#--------------------------------------
# SuperLearner fits of antibody levels
#--------------------------------------
set.seed(0237234)
mauke75 <- SLAb.curve(Y=log10(a75$wb123),Age=a75$age,id=a75$id)
mauke92 <- SLAb.curve(Y=log10(a92$wb123),Age=a92$age,id=a92$id)

#--------------------------------------
# estimate group means
#--------------------------------------
EYx.mauke75 <- SLAb.tmle(Y=log10(a75$wb123),Age=a75$age,id=a75$id)
EYx.mauke92 <- SLAb.tmle(Y=log10(a92$wb123),Age=a92$age,id=a92$id)

#--------------------------------------
# estimate difference in means
#--------------------------------------
diff.mauke  <- SLAb.tmle(Y=log10(a7592$wb123),Age=a7592$age,id=a7592$id,X=a7592$mda,diff=TRUE)


#--------------------------------------
# store results for later summary
# and plotting
#--------------------------------------
save.image("~/SLAbcurve/results/raw/mauke-Wb123-analysis.RData")




