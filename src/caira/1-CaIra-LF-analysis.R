
#------------------------------------
# 1-CaIra-LF-analysis.R
#
# LF transmission changes in Ca Ira
# looking at the age-specific
# antibody curves
#
# version 1 (12 Nov 2015)
#
#------------------------------------

#------------------------------------
# input files
#	  CaIra.csv
#
# output files
#   caira-LF-analysis.RData
#------------------------------------

#------------------------------------
# preamble
#------------------------------------

rm(list=ls())
library(SuperLearner)
library(tmle)
library(RColorBrewer)
library(scales)


# source the base functions for
# SL fits of age-antibody curves
# and TMLE estimates of mean differences
source("~/SLAbcurves/src/SLAb-curve.R")
source("~/SLAbcurves/src/SLAb-tmle.R")
source("~/SLAbcurves/src/SLAb-cvRF.R")



#------------------------------------
# load data
#------------------------------------

d <- read.csv("~/dropbox/caira/data/final/caira.csv")
d$id <- as.character(d$id)
d$personid <- as.character(d$personid)
d$sampledate <- as.Date(d$sampledate,format="%m/%d/%Y")

# Wb123 recode 2 negative values to be positive
d$wb123 <- ifelse(d$wb123<0,-d$wb123,d$wb123)

# Bm14. recode 4 negative values to be positive
d$bm14 <- ifelse(d$bm14<0,-d$bm14,d$bm14)

# exclude children who are <1 year old
# to remove potential for maternal antibodies
d <- subset(d,age>=1)

# add a half year to remove bias from rounded age measurement
d$age <- d$age+0.5

# note: round 2 was partial round (not much data), so for differences
# compute round 3 versus round 1
d <- subset(d,round!=2)
d$mda <- ifelse(d$round==3,1,0)


#-------------------------------
# Step 1: 
# Age-specific antibody response
# curves to the Wb123, Bm14, and
# Bm33 antigens for each round 
# (note: round 2 was partial round)
#------------------------------------
set.seed(1742)
wb123.fit1 <- SLAb.curve(Y=log10(d$wb123[d$round==1]),Age=d$age[d$round==1],id=d$id[d$round==1])
wb123.fit3 <- SLAb.curve(Y=log10(d$wb123[d$round==3]),Age=d$age[d$round==3],id=d$id[d$round==3])

bm14.fit1 <- SLAb.curve(Y=log10(d$bm14[d$round==1]),Age=d$age[d$round==1],id=d$id[d$round==1])
bm14.fit3 <- SLAb.curve(Y=log10(d$bm14[d$round==3]),Age=d$age[d$round==3],id=d$id[d$round==3])

bm33.fit1 <- SLAb.curve(Y=log10(d$bm33[d$round==1]),Age=d$age[d$round==1],id=d$id[d$round==1])
bm33.fit3 <- SLAb.curve(Y=log10(d$bm33[d$round==3]),Age=d$age[d$round==3],id=d$id[d$round==3])


#-------------------------------
# Step 2: 
# Calculate age-adjusted means
# and difference in means
# for Wb123, Bm14, Bm33
# only compare for ages that
# overlap: 4.5 - 11.5 years
#-------------------------------

ad <- subset(d,age>=4.5 & age<=11.5)

set.seed(24734)
mu.wb123 <- sapply(c(1,3),function(x) SLAb.tmle(
  Y=log10(ad$wb123[ad$round==x]),
  Age=ad$age[ad$round==x],
  id=ad$id[ad$round==x]
  )
)

diff.wb123 <- SLAb.tmle(
  Y=log10(ad$wb123),
  Age=ad$age,
  id=ad$id,
  X=ad$mda,
  diff=TRUE
  )

mu.bm14 <- sapply(c(1,3),function(x) SLAb.tmle(
  Y=log10(ad$bm14[ad$round==x]),
  Age=ad$age[ad$round==x],
  id=ad$id[ad$round==x]
  )
)

diff.bm14 <- SLAb.tmle(
  Y=log10(ad$bm14),
  Age=ad$age,
  id=ad$id,
  X=ad$mda,
  diff=TRUE
)

mu.bm33 <- sapply(c(1,3),function(x) SLAb.tmle(
  Y=log10(ad$bm33[ad$round==x]),
  Age=ad$age[ad$round==x],
  id=ad$id[ad$round==x]
)
)

diff.bm33 <- SLAb.tmle(
  Y=log10(ad$bm33),
  Age=ad$age,
  id=ad$id,
  X=ad$mda,
  diff=TRUE
)

#-------------------------------
# save down the results
#-------------------------------
save.image("~/SLAbcurves/results/raw/caira-LF-analysis.RData")




