
#-------------------------------
# 1-haiti2-USA-enterics-SL-curves.R
#
# fit enterics luminex reponse by
# age and country
#
# version 2 (9 oct 2015)
# cleaned up code
# separated out plotting script
#
# version 2 (19 oct 2015)
# updated pointer to base functions
#
# version 1 (23 sep 2015)
#-------------------------------


#-------------------------------
# input files:
#   haiti2-long.dta
#   milwaukee2011.csv
#
# output files:
#   haiti2-USA-enterics-SL-curves.pdf
#-------------------------------

#-------------------------------
# preamble
#-------------------------------

rm(list=ls())
library(foreign)
library(SuperLearner)
library(tmle)
library(RColorBrewer)
library(scales)

# source the base functions for
# SL fits of age-antibody curves
# and TMLE estimates of mean differences
source("~/SLAbcurves/src/SLAb-curve.R")
source("~/SLAbcurves/src/SLAb-tmle.R")



#-------------------------------
# load the datasets for analysis
#-------------------------------

d.usa <- read.csv("~/dropbox/haiti2/data/final/milwaukee2011.csv")
d.haiti <- read.dta("~/dropbox/haiti2/data/final/haiti2-long.dta")

# add sequential ids for USA
d.usa$id2 <- 1:nrow(d.usa)

#-------------------------------
# recode ages by adding 1/2 year
# in USA children >1 y to avoid
# bias
#-------------------------------
d.usa$age[d.usa$age>=1] <- d.usa$age[d.usa$age>=1]+0.5


#-------------------------------
# Limit the haiti data to children
# <= 5.5 years old, which is the
# maximum age for the children
# in the USA sample
#-------------------------------
d.hai  <- subset(d.haiti,agey<=5.5)


#-------------------------------
# append the Haiti and USA 
# data to calculate differences 
# in means
#-------------------------------
d.usa$id <- d.usa$id2
d.usa$agey <- d.usa$age
d.usa$salb <- d.usa$sallpsb
d.usa$haiti <- 0
d.hai$haiti <- 1
common.vars <- c("haiti","id","agey","cp17","cp23","vsp5","leca","etec","salb","norogi","norogii")
d.all <- rbind(subset(d.hai,select=common.vars),subset(d.usa,select=common.vars)	)


#-------------------------------
# crypto Cp17
#-------------------------------
set.seed(4634523)

# fit age specific antibody curves
haiti.cp17 <- SLAb.curve(Y=log10(d.hai$cp17+1),Age=d.hai$agey,id=d.hai$id)
usa.cp17   <- SLAb.curve(Y=log10(d.usa$cp17+1),Age=d.usa$agey,id=d.usa$id)

# estimate group means
EYx.haiti.cp17 <- SLAb.tmle(Y=log10(d.hai$cp17+1),Age=d.hai$agey,id=d.hai$id)
EYx.usa.cp17   <- SLAb.tmle(Y=log10(d.usa$cp17+1),Age=d.usa$agey,id=d.usa$id)

# estimate difference in means
diff.cp17      <- SLAb.tmle(Y=log10(d.all$cp17+1),Age=d.all$agey,id=d.all$id,X=d.all$haiti,diff=TRUE)


#-------------------------------
# crypto Cp23
#-------------------------------
set.seed(6323424)

# fit age specific antibody curves
haiti.cp23 <- SLAb.curve(Y=log10(d.hai$cp23+1),Age=d.hai$agey,id=d.hai$id)
usa.cp23   <- SLAb.curve(Y=log10(d.usa$cp23+1),Age=d.usa$agey,id=d.usa$id)

# estimate group means
EYx.haiti.cp23 <- SLAb.tmle(Y=log10(d.hai$cp23+1),Age=d.hai$agey,id=d.hai$id)
EYx.usa.cp23   <- SLAb.tmle(Y=log10(d.usa$cp23+1),Age=d.usa$agey,id=d.usa$id)

# estimate difference in means
diff.cp23      <- SLAb.tmle(Y=log10(d.all$cp23+1),Age=d.all$agey,id=d.all$id,X=d.all$haiti,diff=TRUE)


#-------------------------------
# Giardia VSP-5
#-------------------------------
set.seed(54635)

# fit age specific antibody curves
haiti.giar <- SLAb.curve(Y=log10(d.hai$vsp5+1),Age=d.hai$agey,id=d.hai$id)
usa.giar   <- SLAb.curve(Y=log10(d.usa$vsp5+1),Age=d.usa$agey,id=d.usa$id)

# estimate group means
EYx.haiti.giar <- SLAb.tmle(Y=log10(d.hai$vsp5+1),Age=d.hai$agey,id=d.hai$id)
EYx.usa.giar   <- SLAb.tmle(Y=log10(d.usa$vsp5+1),Age=d.usa$agey,id=d.usa$id)

# estimate difference in means
diff.giar      <- SLAb.tmle(Y=log10(d.all$vsp5+1),Age=d.all$agey,id=d.all$id,X=d.all$haiti,diff=TRUE)


#-------------------------------
# E. hist LecA
#-------------------------------
set.seed(23523412)

# fit age specific antibody curves
haiti.leca <- SLAb.curve(Y=log10(d.hai$leca+1),Age=d.hai$agey,id=d.hai$id)
usa.leca   <- SLAb.curve(Y=log10(d.usa$leca+1),Age=d.usa$agey,id=d.usa$id)

# estimate group means
EYx.haiti.leca <- SLAb.tmle(Y=log10(d.hai$leca+1),Age=d.hai$agey,id=d.hai$id)
EYx.usa.leca   <- SLAb.tmle(Y=log10(d.usa$leca+1),Age=d.usa$agey,id=d.usa$id)

# estimate difference in means
diff.leca      <- SLAb.tmle(Y=log10(d.all$leca+1),Age=d.all$agey,id=d.all$id,X=d.all$haiti,diff=TRUE)


#-------------------------------
# ETEC
#-------------------------------
set.seed(2967234)

# fit age specific antibody curves
haiti.etec <- SLAb.curve(Y=log10(d.hai$etec+1),Age=d.hai$agey,id=d.hai$id)
usa.etec   <- SLAb.curve(Y=log10(d.usa$etec+1),Age=d.usa$agey,id=d.usa$id)

# estimate group means
EYx.haiti.etec <- SLAb.tmle(Y=log10(d.hai$etec+1),Age=d.hai$agey,id=d.hai$id)
EYx.usa.etec   <- SLAb.tmle(Y=log10(d.usa$etec+1),Age=d.usa$agey,id=d.usa$id)

# estimate difference in means
diff.etec      <- SLAb.tmle(Y=log10(d.all$etec+1),Age=d.all$agey,id=d.all$id,X=d.all$haiti,diff=TRUE)


#-------------------------------
# Salmonella LPS B
#-------------------------------
set.seed(987234)

# fit age specific antibody curves
haiti.salb <- SLAb.curve(Y=log10(d.hai$salb+1),Age=d.hai$agey,id=d.hai$id)
usa.salb   <- SLAb.curve(Y=log10(d.usa$salb+1),Age=d.usa$agey,id=d.usa$id)

# estimate group means
EYx.haiti.salb <- SLAb.tmle(Y=log10(d.hai$salb+1),Age=d.hai$agey,id=d.hai$id)
EYx.usa.salb   <- SLAb.tmle(Y=log10(d.usa$salb+1),Age=d.usa$agey,id=d.usa$id)

# estimate difference in means
diff.salb      <- SLAb.tmle(Y=log10(d.all$salb+1),Age=d.all$agey,id=d.all$id,X=d.all$haiti,diff=TRUE)

#-------------------------------
# noro Group I
#-------------------------------
set.seed(5234)
# fit age specific antibody curves
haiti.norogi <- SLAb.curve(Y=log10(d.hai$norogi+1),Age=d.hai$agey,id=d.hai$id)
usa.norogi   <- SLAb.curve(Y=log10(d.usa$norogi+1),Age=d.usa$agey,id=d.usa$id)

# estimate group means
EYx.haiti.norogi <- SLAb.tmle(Y=log10(d.hai$norogi+1),Age=d.hai$agey,id=d.hai$id)
EYx.usa.norogi   <- SLAb.tmle(Y=log10(d.usa$norogi+1),Age=d.usa$agey,id=d.usa$id)

# estimate difference in means
diff.norogi      <- SLAb.tmle(Y=log10(d.all$norogi+1),Age=d.all$agey,id=d.all$id,X=d.all$haiti,diff=TRUE)


#-------------------------------
# noro Group II
#-------------------------------
set.seed(34534)
# fit age specific antibody curves
haiti.norogii <- SLAb.curve(Y=log10(d.hai$norogii+1),Age=d.hai$agey,id=d.hai$id)
usa.norogii   <- SLAb.curve(Y=log10(d.usa$norogii+1),Age=d.usa$agey,id=d.usa$id)

# estimate group means
EYx.haiti.norogii <- SLAb.tmle(Y=log10(d.hai$norogii+1),Age=d.hai$agey,id=d.hai$id)
EYx.usa.norogii   <- SLAb.tmle(Y=log10(d.usa$norogii+1),Age=d.usa$agey,id=d.usa$id)

# estimate difference in means
diff.norogii      <- SLAb.tmle(Y=log10(d.all$norogii+1),Age=d.all$agey,id=d.all$id,X=d.all$haiti,diff=TRUE)

	

#-------------------------------
# save down results for later
# plotting
#-------------------------------
save.image(file="~/SLAbcurves/results/raw/haiti2-usa-enterics-analysis.RData")








