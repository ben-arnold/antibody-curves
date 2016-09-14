
#-------------------------------
# 1-haiti2-USA-enterics-SL-curves.R
#
# fit enterics luminex reponse by
# age and country
#
#-------------------------------


#-------------------------------
# input files:
#   haiti-enterics-public.csv
#   usa-enterics-public.csv
#
# output files:
#   haiti2-usa-enterics-analysis.RData
#-------------------------------

#-------------------------------
# preamble
#-------------------------------

rm(list=ls())
library(SuperLearner)
library(tmle)
library(tmleAb)


#-------------------------------
# load the datasets for analysis
#-------------------------------

# d.usa <- read.csv("~/dropbox/articles/antibody-curves/data/enterics/usa-enterics-public.csv")
# d.haiti <- read.csv("~/dropbox/articles/antibody-curves/data/enterics/haiti-enterics-public.csv")
data(usa_enterics)
data(haiti_enterics)
d.usa <- usa_enterics
d.haiti <- haiti_enterics


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
d.usa$agey <- d.usa$age
d.usa$haiti <- 0
d.hai$haiti <- 1
common.vars <- c("haiti","id","agey","cp17","cp23","vsp5","leca","etec","salb","norogi","norogii")
d.all <- rbind(subset(d.hai,select=common.vars),subset(d.usa,select=common.vars)	)


# SL library
SL.library <- c("SL.mean","SL.glm","SL.Yman2016","SL.gam","SL.loess")


#-------------------------------
# crypto Cp17
#-------------------------------
set.seed(4634523)

# fit age specific antibody curves
haiti.cp17 <- agecurveAb(Y=log10(d.hai$cp17+1),Age=d.hai$agey,id=d.hai$id,SL.library=SL.library,gamdf=2:6)
usa.cp17   <- agecurveAb(Y=log10(d.usa$cp17+1),Age=d.usa$agey,id=d.usa$id,SL.library=SL.library,gamdf=2:6)

# estimate group means
EYx.haiti.cp17 <- tmleAb(Y=log10(d.hai$cp17+1),W=data.frame(Age=d.hai$agey),id=d.hai$id,SL.library=SL.library,gamdf=2:6)
EYx.usa.cp17   <- tmleAb(Y=log10(d.usa$cp17+1),W=data.frame(Age=d.usa$agey),id=d.usa$id,SL.library=SL.library,gamdf=2:6)

# estimate difference in means
diff.cp17      <- tmleAb(Y=log10(d.all$cp17+1),W=data.frame(Age=d.all$agey),id=d.all$id,X=d.all$haiti,SL.library=SL.library,gamdf=2:6)


#-------------------------------
# crypto Cp23
#-------------------------------
set.seed(6323424)

# fit age specific antibody curves
haiti.cp23 <- agecurveAb(Y=log10(d.hai$cp23+1),Age=d.hai$agey,id=d.hai$id,SL.library=SL.library,gamdf=2:6)
usa.cp23   <- agecurveAb(Y=log10(d.usa$cp23+1),Age=d.usa$agey,id=d.usa$id,SL.library=SL.library,gamdf=2:6)

# estimate group means
EYx.haiti.cp23 <- tmleAb(Y=log10(d.hai$cp23+1),W=data.frame(Age=d.hai$agey),id=d.hai$id,SL.library=SL.library,gamdf=2:6)
EYx.usa.cp23   <- tmleAb(Y=log10(d.usa$cp23+1),W=data.frame(Age=d.usa$agey),id=d.usa$id,SL.library=SL.library,gamdf=2:6)

# estimate difference in means
diff.cp23      <- tmleAb(Y=log10(d.all$cp23+1),W=data.frame(Age=d.all$agey),id=d.all$id,X=d.all$haiti,SL.library=SL.library,gamdf=2:6)


#-------------------------------
# Giardia VSP-5
#-------------------------------
set.seed(54635)

# fit age specific antibody curves
haiti.giar <- agecurveAb(Y=log10(d.hai$vsp5+1),Age=d.hai$agey,id=d.hai$id,SL.library=SL.library,gamdf=2:6)
usa.giar   <- agecurveAb(Y=log10(d.usa$vsp5+1),Age=d.usa$agey,id=d.usa$id,SL.library=SL.library,gamdf=2:6)

# estimate group means
EYx.haiti.giar <- tmleAb(Y=log10(d.hai$vsp5+1),W=data.frame(Age=d.hai$agey),id=d.hai$id,SL.library=SL.library,gamdf=2:6)
EYx.usa.giar   <- tmleAb(Y=log10(d.usa$vsp5+1),W=data.frame(Age=d.usa$agey),id=d.usa$id,SL.library=SL.library,gamdf=2:6)

# estimate difference in means
diff.giar      <- tmleAb(Y=log10(d.all$vsp5+1),W=data.frame(Age=d.all$agey),id=d.all$id,X=d.all$haiti,SL.library=SL.library,gamdf=2:6)


#-------------------------------
# E. hist LecA
#-------------------------------
set.seed(23523412)

# fit age specific antibody curves
haiti.leca <- agecurveAb(Y=log10(d.hai$leca+1),Age=d.hai$agey,id=d.hai$id,SL.library=SL.library,gamdf=2:6)
usa.leca   <- agecurveAb(Y=log10(d.usa$leca+1),Age=d.usa$agey,id=d.usa$id,SL.library=SL.library,gamdf=2:6)

# estimate group means
EYx.haiti.leca <- tmleAb(Y=log10(d.hai$leca+1),W=data.frame(Age=d.hai$agey),id=d.hai$id,SL.library=SL.library,gamdf=2:6)
EYx.usa.leca   <- tmleAb(Y=log10(d.usa$leca+1),W=data.frame(Age=d.usa$agey),id=d.usa$id,SL.library=SL.library,gamdf=2:6)

# estimate difference in means
diff.leca      <- tmleAb(Y=log10(d.all$leca+1),W=data.frame(Age=d.all$agey),id=d.all$id,X=d.all$haiti,SL.library=SL.library,gamdf=2:6)


#-------------------------------
# ETEC
#-------------------------------
set.seed(2967234)

# fit age specific antibody curves
haiti.etec <- agecurveAb(Y=log10(d.hai$etec+1),Age=d.hai$agey,id=d.hai$id,SL.library=SL.library,gamdf=2:6)
usa.etec   <- agecurveAb(Y=log10(d.usa$etec+1),Age=d.usa$agey,id=d.usa$id,SL.library=SL.library,gamdf=2:6)

# estimate group means
EYx.haiti.etec <- tmleAb(Y=log10(d.hai$etec+1),W=data.frame(Age=d.hai$agey),id=d.hai$id,SL.library=SL.library,gamdf=2:6)
EYx.usa.etec   <- tmleAb(Y=log10(d.usa$etec+1),W=data.frame(Age=d.usa$agey),id=d.usa$id,SL.library=SL.library,gamdf=2:6)

# estimate difference in means
diff.etec      <- tmleAb(Y=log10(d.all$etec+1),W=data.frame(Age=d.all$agey),id=d.all$id,X=d.all$haiti,SL.library=SL.library,gamdf=2:6)


#-------------------------------
# Salmonella LPS B
#-------------------------------
set.seed(987234)

# fit age specific antibody curves
haiti.salb <- agecurveAb(Y=log10(d.hai$salb+1),Age=d.hai$agey,id=d.hai$id,SL.library=SL.library,gamdf=2:6)
usa.salb   <- agecurveAb(Y=log10(d.usa$salb+1),Age=d.usa$agey,id=d.usa$id,SL.library=SL.library,gamdf=2:6)

# estimate group means
EYx.haiti.salb <- tmleAb(Y=log10(d.hai$salb+1),W=data.frame(Age=d.hai$agey),id=d.hai$id,SL.library=SL.library,gamdf=2:6)
EYx.usa.salb   <- tmleAb(Y=log10(d.usa$salb+1),W=data.frame(Age=d.usa$agey),id=d.usa$id,SL.library=SL.library,gamdf=2:6)

# estimate difference in means
diff.salb      <- tmleAb(Y=log10(d.all$salb+1),W=data.frame(Age=d.all$agey),id=d.all$id,X=d.all$haiti,SL.library=SL.library,gamdf=2:6)

#-------------------------------
# noro Group I
#-------------------------------
set.seed(5234)
# fit age specific antibody curves
haiti.norogi <- agecurveAb(Y=log10(d.hai$norogi+1),Age=d.hai$agey,id=d.hai$id,SL.library=SL.library,gamdf=2:6)
usa.norogi   <- agecurveAb(Y=log10(d.usa$norogi+1),Age=d.usa$agey,id=d.usa$id,SL.library=SL.library,gamdf=2:6)

# estimate group means
EYx.haiti.norogi <- tmleAb(Y=log10(d.hai$norogi+1),W=data.frame(Age=d.hai$agey),id=d.hai$id,SL.library=SL.library,gamdf=2:6)
EYx.usa.norogi   <- tmleAb(Y=log10(d.usa$norogi+1),W=data.frame(Age=d.usa$agey),id=d.usa$id,SL.library=SL.library,gamdf=2:6)

# estimate difference in means
diff.norogi      <- tmleAb(Y=log10(d.all$norogi+1),W=data.frame(Age=d.all$agey),id=d.all$id,X=d.all$haiti,SL.library=SL.library,gamdf=2:6)


#-------------------------------
# noro Group II
#-------------------------------
set.seed(34534)
# fit age specific antibody curves
haiti.norogii <- agecurveAb(Y=log10(d.hai$norogii+1),Age=d.hai$agey,id=d.hai$id,SL.library=SL.library,gamdf=2:6)
usa.norogii   <- agecurveAb(Y=log10(d.usa$norogii+1),Age=d.usa$agey,id=d.usa$id,SL.library=SL.library,gamdf=2:6)

# estimate group means
EYx.haiti.norogii <- tmleAb(Y=log10(d.hai$norogii+1),W=data.frame(Age=d.hai$agey),id=d.hai$id,SL.library=SL.library,gamdf=2:6)
EYx.usa.norogii   <- tmleAb(Y=log10(d.usa$norogii+1),W=data.frame(Age=d.usa$agey),id=d.usa$id,SL.library=SL.library,gamdf=2:6)

# estimate difference in means
diff.norogii      <- tmleAb(Y=log10(d.all$norogii+1),W=data.frame(Age=d.all$agey),id=d.all$id,X=d.all$haiti,SL.library=SL.library,gamdf=2:6)

	

#-------------------------------
# save down results for later
# plotting
#-------------------------------
save.image(file="~/dropbox/articles/antibody-curves/results/raw/haiti2-usa-enterics-analysis.RData")








