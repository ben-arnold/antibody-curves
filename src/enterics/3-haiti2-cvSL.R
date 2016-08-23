
#-------------------------------
# 8-haiti2-cvSL.R
#
# Compute the cross-validated
# risk for the super leaner
# and its constituent algorithms
#
# repeat calculations for multiple
# enteric pathogens in the Haiti-2
# cohort
#-------------------------------



#-------------------------------
# preamble
#-------------------------------
rm(list=ls())
library(SuperLearner)
library(tmle)

# source the base functions for
# SL fits of age-antibody curves
# and TMLE estimates of mean differences
source("~/SLAbcurves/src/SLAb-curve.R")
source("~/SLAbcurves/src/SLAb-tmle.R")
source("~/SLAbcurves/src/SLAb-cvRF.R")
source("~/SLAbcurves/src/SL.Yman2016.R")
source("~/SLAbcurves/src/slab_cvSL.R")
source("~/SLAbcurves/src/plot_slab_cvSL.R")


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


#-------------------------------
# fit cross-validated SL
#-------------------------------

SL.library <- c("SL.mean","SL.glm","SL.loess","SL.gam","SL.randomForest","SL.Yman2016")
# "SL.glmnet"


# crypto
set.seed(56234)
CVcp23 <- slab_cvSL(Y=log10(d.hai$cp23+1),Age=d.hai$agey,id=d.hai$id,family=gaussian(),V=10,SL.library=SL.library)

# giardia
CVvsp5 <- slab_cvSL(Y=log10(d.hai$vsp5+1),Age=d.hai$agey,id=d.hai$id,family=gaussian(),V=10,SL.library=SL.library)

# salmonella LPS B
CVsalB <- slab_cvSL(Y=log10(d.hai$salb+1),Age=d.hai$agey,id=d.hai$id,family=gaussian(),V=10,SL.library=SL.library)

# noro GII
CVnorogii <- slab_cvSL(Y=log10(d.hai$norogii+1),Age=d.hai$agey,id=d.hai$id,family=gaussian(),V=10,SL.library=SL.library)

#------------------------------
# print results to log
#------------------------------
# crypto
summary(CVcp23)

# giardia
summary(CVvsp5)

# salmonella LPS B
summary(CVsalB)

# noro GII
summary(CVnorogii)

#-------------------------------
# plot the estimates
#-------------------------------
pdf("~/dropbox/articles/antibody-curves/results/figs/haiti-cvSL-cp23.pdf")
plot_slab_cvSL(CVcp23,title="Cryptosporidium parvum Cp23, Haiti")
dev.off()
pdf("~/dropbox/articles/antibody-curves/results/figs/haiti-cvSL-vsp5.pdf")
plot_slab_cvSL(CVvsp5,title="Giardia intestinalis VSP-5, Haiti")
dev.off()
pdf("~/dropbox/articles/antibody-curves/results/figs/haiti-cvSL-salB.pdf")
plot_slab_cvSL(CVsalB,title="Salmonella sp. LPS Group B, Haiti")
dev.off()
pdf("~/dropbox/articles/antibody-curves/results/figs/haiti-cvSL-ngii.pdf")
plot_slab_cvSL(CVnorogii,title="Norovirus GII.4 NO, Haiti")
dev.off()


#-------------------------------
# save down the results
#-------------------------------
save(cv.mitonfit,file="~/dropbox/articles/antibody-curves/results/raw/miton-cvSL.RData")



