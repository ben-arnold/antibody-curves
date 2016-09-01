
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

d.usa <- read.csv("~/dropbox/articles/antibody-curves/data/enterics/usa-enterics-public.csv")
d.haiti <- read.csv("~/dropbox/articles/antibody-curves/data/enterics/haiti-enterics-public.csv")


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
haiti.cp17 <- ab_agecurve(Y=log10(d.hai$cp17+1),Age=d.hai$agey,id=d.hai$id,SL.library=SL.library)
usa.cp17   <- ab_agecurve(Y=log10(d.usa$cp17+1),Age=d.usa$agey,id=d.usa$id,SL.library=SL.library)

# estimate group means
EYx.haiti.cp17 <- ab_tmle(Y=log10(d.hai$cp17+1),Age=d.hai$agey,id=d.hai$id,SL.library=SL.library)
EYx.usa.cp17   <- ab_tmle(Y=log10(d.usa$cp17+1),Age=d.usa$agey,id=d.usa$id,SL.library=SL.library)

# estimate difference in means
diff.cp17      <- ab_tmle(Y=log10(d.all$cp17+1),Age=d.all$agey,id=d.all$id,X=d.all$haiti,SL.library=SL.library,diff=TRUE)


#-------------------------------
# crypto Cp23
#-------------------------------
set.seed(6323424)

# fit age specific antibody curves
haiti.cp23 <- ab_agecurve(Y=log10(d.hai$cp23+1),Age=d.hai$agey,id=d.hai$id,SL.library=SL.library)
usa.cp23   <- ab_agecurve(Y=log10(d.usa$cp23+1),Age=d.usa$agey,id=d.usa$id,SL.library=SL.library)

# estimate group means
EYx.haiti.cp23 <- ab_tmle(Y=log10(d.hai$cp23+1),Age=d.hai$agey,id=d.hai$id,SL.library=SL.library)
EYx.usa.cp23   <- ab_tmle(Y=log10(d.usa$cp23+1),Age=d.usa$agey,id=d.usa$id,SL.library=SL.library)

# estimate difference in means
diff.cp23      <- ab_tmle(Y=log10(d.all$cp23+1),Age=d.all$agey,id=d.all$id,X=d.all$haiti,SL.library=SL.library,diff=TRUE)


#-------------------------------
# Giardia VSP-5
#-------------------------------
set.seed(54635)

# fit age specific antibody curves
haiti.giar <- ab_agecurve(Y=log10(d.hai$vsp5+1),Age=d.hai$agey,id=d.hai$id,SL.library=SL.library)
usa.giar   <- ab_agecurve(Y=log10(d.usa$vsp5+1),Age=d.usa$agey,id=d.usa$id,SL.library=SL.library)

# estimate group means
EYx.haiti.giar <- ab_tmle(Y=log10(d.hai$vsp5+1),Age=d.hai$agey,id=d.hai$id,SL.library=SL.library)
EYx.usa.giar   <- ab_tmle(Y=log10(d.usa$vsp5+1),Age=d.usa$agey,id=d.usa$id,SL.library=SL.library)

# estimate difference in means
diff.giar      <- ab_tmle(Y=log10(d.all$vsp5+1),Age=d.all$agey,id=d.all$id,X=d.all$haiti,SL.library=SL.library,diff=TRUE)


#-------------------------------
# E. hist LecA
#-------------------------------
set.seed(23523412)

# fit age specific antibody curves
haiti.leca <- ab_agecurve(Y=log10(d.hai$leca+1),Age=d.hai$agey,id=d.hai$id,SL.library=SL.library)
usa.leca   <- ab_agecurve(Y=log10(d.usa$leca+1),Age=d.usa$agey,id=d.usa$id,SL.library=SL.library)

# estimate group means
EYx.haiti.leca <- ab_tmle(Y=log10(d.hai$leca+1),Age=d.hai$agey,id=d.hai$id,SL.library=SL.library)
EYx.usa.leca   <- ab_tmle(Y=log10(d.usa$leca+1),Age=d.usa$agey,id=d.usa$id,SL.library=SL.library)

# estimate difference in means
diff.leca      <- ab_tmle(Y=log10(d.all$leca+1),Age=d.all$agey,id=d.all$id,X=d.all$haiti,SL.library=SL.library,diff=TRUE)


#-------------------------------
# ETEC
#-------------------------------
set.seed(2967234)

# fit age specific antibody curves
haiti.etec <- ab_agecurve(Y=log10(d.hai$etec+1),Age=d.hai$agey,id=d.hai$id,SL.library=SL.library)
usa.etec   <- ab_agecurve(Y=log10(d.usa$etec+1),Age=d.usa$agey,id=d.usa$id,SL.library=SL.library)

# estimate group means
EYx.haiti.etec <- ab_tmle(Y=log10(d.hai$etec+1),Age=d.hai$agey,id=d.hai$id,SL.library=SL.library)
EYx.usa.etec   <- ab_tmle(Y=log10(d.usa$etec+1),Age=d.usa$agey,id=d.usa$id,SL.library=SL.library)

# estimate difference in means
diff.etec      <- ab_tmle(Y=log10(d.all$etec+1),Age=d.all$agey,id=d.all$id,X=d.all$haiti,SL.library=SL.library,diff=TRUE)


#-------------------------------
# Salmonella LPS B
#-------------------------------
set.seed(987234)

# fit age specific antibody curves
haiti.salb <- ab_agecurve(Y=log10(d.hai$salb+1),Age=d.hai$agey,id=d.hai$id,SL.library=SL.library)
usa.salb   <- ab_agecurve(Y=log10(d.usa$salb+1),Age=d.usa$agey,id=d.usa$id,SL.library=SL.library)

# estimate group means
EYx.haiti.salb <- ab_tmle(Y=log10(d.hai$salb+1),Age=d.hai$agey,id=d.hai$id,SL.library=SL.library)
EYx.usa.salb   <- ab_tmle(Y=log10(d.usa$salb+1),Age=d.usa$agey,id=d.usa$id,SL.library=SL.library)

# estimate difference in means
diff.salb      <- ab_tmle(Y=log10(d.all$salb+1),Age=d.all$agey,id=d.all$id,X=d.all$haiti,SL.library=SL.library,diff=TRUE)

#-------------------------------
# noro Group I
#-------------------------------
set.seed(5234)
# fit age specific antibody curves
haiti.norogi <- ab_agecurve(Y=log10(d.hai$norogi+1),Age=d.hai$agey,id=d.hai$id,SL.library=SL.library)
usa.norogi   <- ab_agecurve(Y=log10(d.usa$norogi+1),Age=d.usa$agey,id=d.usa$id,SL.library=SL.library)

# estimate group means
EYx.haiti.norogi <- ab_tmle(Y=log10(d.hai$norogi+1),Age=d.hai$agey,id=d.hai$id,SL.library=SL.library)
EYx.usa.norogi   <- ab_tmle(Y=log10(d.usa$norogi+1),Age=d.usa$agey,id=d.usa$id,SL.library=SL.library)

# estimate difference in means
diff.norogi      <- ab_tmle(Y=log10(d.all$norogi+1),Age=d.all$agey,id=d.all$id,X=d.all$haiti,SL.library=SL.library,diff=TRUE)


#-------------------------------
# noro Group II
#-------------------------------
set.seed(34534)
# fit age specific antibody curves
haiti.norogii <- ab_agecurve(Y=log10(d.hai$norogii+1),Age=d.hai$agey,id=d.hai$id,SL.library=SL.library)
usa.norogii   <- ab_agecurve(Y=log10(d.usa$norogii+1),Age=d.usa$agey,id=d.usa$id,SL.library=SL.library)

# estimate group means
EYx.haiti.norogii <- ab_tmle(Y=log10(d.hai$norogii+1),Age=d.hai$agey,id=d.hai$id,SL.library=SL.library)
EYx.usa.norogii   <- ab_tmle(Y=log10(d.usa$norogii+1),Age=d.usa$agey,id=d.usa$id,SL.library=SL.library)

# estimate difference in means
diff.norogii      <- ab_tmle(Y=log10(d.all$norogii+1),Age=d.all$agey,id=d.all$id,X=d.all$haiti,SL.library=SL.library,diff=TRUE)

	

#-------------------------------
# save down results for later
# plotting
#-------------------------------
save.image(file="~/dropbox/articles/antibody-curves/results/raw/haiti2-usa-enterics-analysis.RData")








