#-------------------------------
# 9-garki-cvSL-figure.R
#
# plot cross-validated risk
# from the super learner
# and consituent algorithms
#-------------------------------



#-------------------------------
# preamble
#-------------------------------
rm(list=ls())
library(SuperLearner)
# library(ggplot2)
library(RColorBrewer)

source("~/SLAbcurves/src/plot_slab_cvSL.R")

#-------------------------------
# load the results
#-------------------------------
load("~/dropbox/articles/antibody-curves/results/raw/garki-cvSL.RData")


#-------------------------------
# get summary estimates
#-------------------------------

# summary estimates for intervention round 5 (no covariates)
summary(cv.tfit)

# summary estimates for multivariate example (round 3-5)
summary(cv.mfit)

# plot the estimates
pdf("~/dropbox/articles/antibody-curves/results/figs/garki-cvSL-age.pdf")
plot_slab_cvSL(cv.tfit,title="P. falciparum, high transmission (age only)")
dev.off()
pdf("~/dropbox/articles/antibody-curves/results/figs/garki-cvSL-mv.pdf")
plot_slab_cvSL(cv.mfit,title="P. falciparum, high transmission (multivariate)")
dev.off()




