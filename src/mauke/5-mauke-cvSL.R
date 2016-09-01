
#-------------------------------
# 5-mauke-cvSL.R
#
# Compute the cross-validated
# risk for the super leaner
# and its constituent algorithms
#
# in the Mauke 1992 survey
#-------------------------------

#-------------------------------
# preamble
#-------------------------------
rm(list=ls())
library(SuperLearner)
library(tmleAb)
library(r2weight)
library(RColorBrewer)

#-------------------------------
# load data
#-------------------------------
d <- read.csv("~/dropbox/articles/antibody-curves/data/mauke/mauke1992-public.csv")

# add 0.5 years to age to remove bias (on average) due to rounding to year
d$age <- d$age+0.5


#-------------------------------
# fit cross-validated SL
# with age as the only feature
#-------------------------------
SL.library <- c("SL.mean","SL.glm","SL.Yman2016","SL.gam","SL.loess","SL.randomForest","SL.polymars","SL.svm","SL.nnet")

set.seed(32423)
CVmauke <- ab_cvSL(Y=log10(d$wb123),Age=d$age,family=gaussian(),V=10,SL.library=SL.library)

summary(CVmauke)


#-------------------------------
# plot the CV MSE estimates
#-------------------------------
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cols <- cbPalette[c(1:4,6:8)]
cols <- c(cols,brewer.pal(8,"Dark2")[1:3])

pdf("~/dropbox/articles/antibody-curves/results/figs/mauke-cvSL.pdf")
ab_plot_cvSL(CVmauke,col=cols,ylab="10-fold Cross-validated MSE Estimate",title="W. bancrofti, Mauke 1992")
dev.off()

#-------------------------------
# convert CV MSE into R2
# using the r2weight package
# and plot the R2 estimates
#-------------------------------
pdf("~/dropbox/articles/antibody-curves/results/figs/mauke-cvR2.pdf")
ab_plot_cvR2(CVmauke,data.frame(Age=d$age),col=cols,ylab="10-fold Cross-validated R-squared",title="W. bancrofti, Mauke 1992",ylim=c(0,0.6))
dev.off()


#-------------------------------
# save down the results
#-------------------------------
save(CVmauke,file="~/dropbox/articles/antibody-curves/results/raw/mauke-cvSL.RData")



