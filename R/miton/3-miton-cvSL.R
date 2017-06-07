
#-------------------------------
# 3-miton-cvSL.R
#
# Compute the cross-validated
# risk for the super leaner
# and its constituent algorithms
#
# in the Miton Haiti survey
# (low transmission P. falciparum MSP-1)
#-------------------------------

#-------------------------------
# preamble
#-------------------------------
rm(list=ls())
library(SuperLearner)
library(tmleAb)


# load general plotting functions
library(ggplot2)
library(r2weight)
library(RColorBrewer)
source("~/slabcurves/src/ab_plot_cvSL.R")
source("~/slabcurves/src/ab_plot_cvR2.R")

#-------------------------------
# load data
#-------------------------------
# d <- read.csv("~/dropbox/articles/antibody-curves/data/miton/haiti2-malaria-miton-public.csv")
d <- miton_malaria

# for 6 negative MSP_1 values, replace as 0
d$msp1 <- d$msp13d7
d$msp1[d$msp1<0] <- 0


#-------------------------------
# fit cross-validated SL
# with age as the only feature
#-------------------------------
SL.library <- c("SL.mean","SL.glm","SL.Yman2016","SL.gam","SL.loess","SL.randomForest","SL.polymars")

set.seed(32423)
CVmiton <- cvSLAb(Y=log10(d$msp1+1),X=data.frame(Age=d$age),family=gaussian(),V=10,SL.library=SL.library)

summary(CVmiton)


#-------------------------------
# plot the CV MSE estimates
#-------------------------------
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cols <- cbPalette[c(1:4,6:8)]
cols <- c(cols,brewer.pal(8,"Dark2")[3])

# pdf("~/dropbox/articles/antibody-curves/results/figs/miton-cvSL.pdf")
# ab_plot_cvSL(CVmiton,col=cols,ylab="10-fold Cross-validated MSE Estimate",title="P. falciparum, low transmission (age only)")
# dev.off()

#-------------------------------
# convert CV MSE into R2
# using the r2weight package
# and plot the R2 estimates
#-------------------------------
pdf("~/dropbox/articles/antibody-curves/results/figs/miton-cvR2.pdf")
ab_plot_cvR2(CVmiton,data.frame(Age=d$age),col=cols,ylab="10-fold Cross-validated R-squared",title="P. falciparum, (Miton, Haiti)",ylim=c(0,0.6))
dev.off()


