

#-------------------------------------------
# 5-mauke-Wb123-long.R
# Ben Arnold
#
# TMLE means of WB123 antibody in Mauke
# among individuals measured at both time
# points, stratified by age
#
#-------------------------------------------

#-------------------------------------------
# input files:
#   mauke1974.dta
#   mauke1992.dta
#
# output files:
#   mauke-Wb123-long.RData
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
source("~/SLAbcurves/src/SLAb-curve.R")
source("~/SLAbcurves/src/SLAb-tmle.R")
source("~/SLAbcurves/src/SLAb-cvRF.R")

#-------------------------------------------
# load the Mauke data from 1974(1975) and 1992
#-------------------------------------------

d75 <- read.dta("~/dropbox/mauke/data/final/mauke1974.dta")
d92 <- read.dta("~/dropbox/mauke/data/final/mauke1992.dta")


# drop 7 children aged 0 (all Wb123+) in 1975 due to maternal antibodies
a75 <- subset(d75,age>0,select=c("id74","age","CAg","wb123"))
  names(a75) <- c("id","age","CAg","wb123")
a75$mda <- 0
a92 <- subset(d92,select=c("id74","age","CAg","wb123"))
  names(a92) <- c("id","age","CAg","wb123")
a92$id <- paste(a92$id)
a92$mda <- 1

# create an appended dataset (long format)
d7592 <- rbind(a75,a92)

# create a merged dataset (wide format)
names(a75) <- c("id","age.75","CAg.75","wb123.75","mda.75")
names(a92) <- c("id","age.92","CAg.92","wb123.92","mda.92")
d <- merge(a75,a92,by="id")

# restrict long format dataset to individuals 
# with measuremments in both 1975 and 1992
a7592 <- merge(d7592,subset(d,select="id"),by="id",all.x=F)

# create an age category in 1975 (not used)
# d$agecat <- cut(d$age.75,breaks=c(1,20,40,80))
# quantile(d$age.75)
# d$agecat <- cut(d$age.75,breaks=c(0,18,33,43,65))

# create antigen status variables in 1975 and 1992
d$ab75 <- as.factor(ifelse(d$CAg.75>32,"Pos","Neg"))
d$ab92 <- as.factor(ifelse(d$CAg.92>32,"Pos","Neg"))
d$Abstatus <- factor(rep("Neg-Neg",nrow(d)),levels=c("Neg-Neg","Pos-Neg","Pos-Pos","Neg-Pos"))
  d$Abstatus[d$ab75=="Pos"&d$ab92=="Neg"] <- "Pos-Neg"
  d$Abstatus[d$ab75=="Pos"&d$ab92=="Pos"] <- "Pos-Pos"
  d$Abstatus[d$ab75=="Neg"&d$ab92=="Pos"] <- "Neg-Pos"
  
# merge antigen status at both time points into the long format data as well
a7592 <- merge(a7592,subset(d,select=c("id","Abstatus")),by="id",all.x=T)

#--------------------------------------
# Estimate means and differences between
# time points for each antigen+/- 
# combo
#--------------------------------------
set.seed(2194543)
Abstatus <- c("Neg-Neg","Pos-Neg","Pos-Pos")
EYx.75.Abstatus <- sapply(Abstatus, function(x) 
	SLAb.tmle(Y=log10(d$wb123.75[d$Abstatus==x]),Age=d$age.75[d$Abstatus==x],id=d$id[d$Abstatus==x]) 
	)
EYx.92.Abstatus <- sapply(Abstatus, function(x) 
	SLAb.tmle(Y=log10(d$wb123.92[d$Abstatus==x]),Age=d$age.92[d$Abstatus==x],id=d$id[d$Abstatus==x]) 
	)
diff.Abstatus <- sapply(Abstatus, function(x) 
	SLAb.tmle(Y=log10(a7592$wb123[a7592$Abstatus==x]), Age=a7592$age[a7592$Abstatus==x], id=a7592$id[a7592$Abstatus==x], X=a7592$mda[a7592$Abstatus==x], diff=TRUE) 
	)


#--------------------------------------
# store results for later summary
# and plotting
#--------------------------------------
save.image("~/SLAbcurves/results/raw/mauke-Wb123-long.RData")




