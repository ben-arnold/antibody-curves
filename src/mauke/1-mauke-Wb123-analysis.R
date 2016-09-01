

#-------------------------------------------
# 1-mauke-Wb123-analysis.R
# Ben Arnold
#
# SuperLearner fits of age-specific antibody
# response on Mauke island before and after
# MDA, and TMLE estimates of geometric mean
# (area under the curve) of antibody responses
#
#-------------------------------------------

#-------------------------------------------
# input files:
#   mauke1975-public.csv
#   mauke1992-public.csv
#
# output files:
#   mauke-Wb123-analysis.RData
#-------------------------------------------



#-------------------------------------------
# preamble
#-------------------------------------------

rm(list=ls())
library(tmle)
library(SuperLearner)
library(tmleAb)


#-------------------------------------------
# load the Mauke data from 1974(1975) and 1992
#-------------------------------------------

d75 <- read.csv("~/dropbox/articles/antibody-curves/data/mauke/mauke1975-public.csv")
d92 <- read.csv("~/dropbox/articles/antibody-curves/data/mauke/mauke1992-public.csv")

a75 <- d75
a92 <- d92

# add 0.5 years to age to remove bias (on average) due to rounding to year
a75$ager <- a75$age+0.5
a92$ager <- a92$age+0.5

# create age categories for stratified analyses
# 5 year age categories (1-20 y)
a75$agecat <- cut(a75$ager,breaks=c(0,5,10,15,20),labels=c("1-5","6-10","11-15","16-20"))
a92$agecat <- cut(a92$ager,breaks=c(0,5,10,15,20),labels=c("1-5","6-10","11-15","16-20"))

# 2 year age categories (1-10 y)
a75$agecat2 <- cut(a75$ager,breaks=c(0,2,4,6,8,10),labels=c("1-2","3-4","5-6","7-8","9-10"))
a92$agecat2 <- cut(a92$ager,breaks=c(0,2,4,6,8,10),labels=c("1-2","3-4","5-6","7-8","9-10"))

# create a standard ID variable before appending
a75$id <- as.integer(a75$id75)
a92$id <- ifelse(is.na(a92$id75),a92$id92,a92$id75)

# identify the pre- vs. post-MDA measurements
a75$mda <- 0
a92$mda <- 1

# append data with common variables
common.vars <- c("id","ager","agecat","agecat2","wb123","mda")
a7592 <- rbind(subset(a75,select=common.vars),subset(a92,select=common.vars))


#--------------------------------------
# All Ages
#--------------------------------------
# SL library
SL.library <- c("SL.mean","SL.glm","SL.Yman2016","SL.gam","SL.loess")

set.seed(0237234)

# SuperLearner fits of antibody levels
mauke75 <- ab_agecurve(Y=log10(a75$wb123),Age=a75$ager,id=a75$id,SL.library=SL.library)
mauke92 <- ab_agecurve(Y=log10(a92$wb123),Age=a92$ager,id=a92$id,SL.library=SL.library)

# estimate group means
EYx.mauke75 <- ab_tmle(Y=log10(a75$wb123),Age=a75$ager,id=a75$id,SL.library=SL.library)
EYx.mauke92 <- ab_tmle(Y=log10(a92$wb123),Age=a92$ager,id=a92$id,SL.library=SL.library)

# estimate difference in means
diff.mauke  <- ab_tmle(Y=log10(a7592$wb123),Age=a7592$ager,id=a7592$id,X=a7592$mda,SL.library=SL.library,diff=TRUE)


#--------------------------------------
# Estimate means and differences between
# time points in
# 5 year age bands from ages 1-20
#--------------------------------------
agegrps <-c("1-5","6-10","11-15","16-20") 
EYx.mauke75kids <- sapply(agegrps, function(x) 
	ab_tmle(Y=log10(a75$wb123[a75$agecat==x]),Age=a75$ager[a75$agecat==x],id=a75$id[a75$agecat==x],SL.library=SL.library) 
	)
EYx.mauke92kids <- sapply(agegrps, function(x) 
	ab_tmle(Y=log10(a92$wb123[a92$agecat==x]),Age=a92$ager[a92$agecat==x],id=a92$id[a92$agecat==x],SL.library=SL.library) 
	)
diff.maukekids <- sapply(agegrps, function(x) 
	ab_tmle(Y=log10(a7592$wb123[a7592$agecat==x]), Age=a7592$ager[a7592$agecat==x], id=a7592$id[a7592$agecat==x], X=a7592$mda[a7592$agecat==x],SL.library=SL.library, diff=TRUE) 
	)


#--------------------------------------
# store results for later summary
# and plotting
#--------------------------------------
save.image("~/dropbox/articles/antibody-curves/results/raw/mauke-Wb123-analysis.RData")




