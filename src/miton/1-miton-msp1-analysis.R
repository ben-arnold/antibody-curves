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
#   haiti2-malaria-miton-public.csv
#
# output files:
#   miton-msp1-analysis.RData
#-------------------------------


#-------------------------------
# preamble
#-------------------------------

rm(list=ls())
# library(SuperLearner)
# library(tmle)
library(SLAb)


#-------------------------------------------
# load the Miton data
#-------------------------------------------

d <- read.csv("~/dropbox/articles/antibody-curves/data/miton/haiti2-malaria-miton-public.csv")

d$msp1 <- d$msp13d7

# for 6 negative MSP_1 values, replace as 0.001
d$msp1[d$msp1<0] <- 0

# create age categories for stratified analyses
# 5 year age categories (1-20 y)
d$agecat <- cut(d$age,breaks=c(0,5,10,15,20),labels=c("1-5","6-10","11-15","16-20"))

#-------------------------------------------
# set the library of models / algorithms
#-------------------------------------------
SL.library <- c("SL.mean","SL.glm","SL.loess","SL.gam","SL.randomForest","SL.Yman2016")

#-------------------------------------------
# estimate a marginal Ab curve E(Y_x,a)
#-------------------------------------------
set.seed(25234)
msp1_EYxa <- slab_curve(Y=log10(d$msp1+1),Age=d$age,family="gaussian",SL.library=SL.library)

#-------------------------------------------
# estimate age-adjusted means E(Y_x)
# by 5 year age category
#-------------------------------------------
agegrps <-c("1-5","6-10","11-15","16-20")
msp1_EYx <- sapply(agegrps, function(x)
  slab_tmle(Y=log10(d$msp1[d$agecat==x]+1),Age=data.frame(Age=d$age[d$agecat==x]),SL.library=SL.library,family="gaussian")
)


#-------------------------------------------
# save results
#-------------------------------------------
save.image(file="~/dropbox/articles/antibody-curves/results/raw/miton-msp1-analysis.RData")




