
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

# load general plotting functions
library(ggplot2)
library(r2weight)
library(RColorBrewer)
source("~/slabcurves/src/ab_plot_cvSL.R")
source("~/slabcurves/src/ab_plot_cvR2.R")


#-------------------------------
# load data
# restrict to 1975 for illustration
#-------------------------------
# d <- read.csv("~/dropbox/articles/antibody-curves/data/mauke/mauke1992-public.csv")
data("mauke_wb123")
d <- subset(mauke_wb123,year=='1992' & age<=70)

# add 0.5 years to age to remove bias (on average) due to rounding to year
d$age <- d$age+0.5


#-------------------------------
# fit cross-validated SL
# with age as the only feature
# full library
#-------------------------------
SL.library <- c("SL.mean","SL.glm","SL.Yman2016","SL.gam","SL.loess","SL.randomForest","SL.polymars")

set.seed(0237234)
maukeCVfull <- cvSLAb(Y=log10(d$wb123),X=data.frame(Age=d$age),family=gaussian(),V=10,SL.library=SL.library,gamdf=2:6)
summary(maukeCVfull)

#-------------------------------
# fit cross-validated SL
# with age as the only feature
# restricted library
#-------------------------------
SL.library <- c("SL.mean","SL.glm","SL.Yman2016","SL.gam","SL.loess")

set.seed(0237234)
maukeCVres <- cvSLAb(Y=log10(d$wb123),X=data.frame(Age=d$age),family=gaussian(),V=10,SL.library=SL.library,gamdf=2:6)
summary(maukeCVres)


#-------------------------------
# plot the CV MSE estimates
#-------------------------------
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cols <- cbPalette[c(1:4,6:8)]
cols <- c(cols,brewer.pal(8,"Dark2")[3])

# pdf("~/dropbox/articles/antibody-curves/results/figs/mauke-cvSL.pdf")
# ab_plot_cvSL(CVmauke,col=cols,ylab="10-fold Cross-validated MSE Estimate",title="W. bancrofti, Mauke 1992")
# dev.off()

#-------------------------------
# convert CV MSE into R2
# using the r2weight package
# and plot the R2 estimates
#-------------------------------
pdf("~/dropbox/articles/antibody-curves/results/figs/mauke-cvR2.pdf")
ab_plot_cvR2(maukeCVfull,data.frame(Age=d$age),ylab="10-fold Cross-validated R-squared",title=expression(paste(italic('Wuchereria bancrofti')," Wb123, Mauke")),col=cols,ylim=c(0,0.6))
dev.off()

#-------------------------------
# plot the age-antibody curves
# for the full SL and the 
# restricted SL library to
# compare EYxa and EYx
# setting the same seed to
# get exactly the same splits
# for cross-validation
#-------------------------------

# full library
SL.library <- c("SL.mean","SL.glm","SL.Yman2016","SL.gam","SL.loess","SL.randomForest","SL.polymars")
set.seed(23423)
maukeSLfull <- agecurveAb(Y=log10(d$wb123),Age=d$age,id=d$id,SL.library=SL.library,gamdf=2:6)

# restricted library
set.seed(834524)
SL.library <- c("SL.mean","SL.glm","SL.Yman2016","SL.gam","SL.loess")
maukeSLres <- agecurveAb(Y=log10(d$wb123),Age=d$age,id=d$id,SL.library=SL.library,gamdf=2:6)



# plot age-dependent antibody curves and means for two SL libraries
pdf("~/dropbox/articles/antibody-curves/results/figs/mauke-Wb123-cvSLcurves.pdf",width=5,height=5)

ytics <- seq(2,6,by=1)
xtics <- seq(0,70,by=10)

set.seed(623) # jitter ages (but not actual fits) to better display the data
op <- par(mar=c(4,5,4,0)+0.1,xpd=T)
plot(jitter(maukeSLfull$Age),maukeSLfull$Y,col=alpha("black",alpha=0.4),pch=16,cex=0.3,
     ylim=c(2,6),ylab="",yaxt="n",
     xlim=c(0,max(xtics)+1),xlab="",xaxt="n",
     bty="n",las=1
)
mtext(expression(paste(italic('W. bancrofti')," Wb123 (light units)")),side=2,line=3)
# mtext("d",adj=1,line=2,at=-3,font=2,cex=1.75)
# mtext("Intervention Period",adj=0,line=3,at=0,cex=1.5)
mtext(expression(paste(italic('Wuchereria bancrofti')," Wb123, Mauke")),side=3,line=2)
mtext("Super learner ensemble prediction",side=3,line=0.5)
mtext("Age, years",side=1,line=2.5)
axis(2,at=2:6,labels=c(
  expression(10^2),
  expression(10^3),
  expression(10^4),
  expression(10^5),
  expression(10^6)
), las=1,cex.axis=1.25
)
# segments(x0=min(xtics),x1=max(xtics),y0=ytics,lty=2,col="gray70")
axis(1,at=xtics,cex.axis=1.5)
lines(maukeSLfull$Age,maukeSLfull$pY,col=cols[1],lwd=1)
lines(maukeSLres$Age,maukeSLres$pY,col=cols[2],lwd=1)

legend(x=70,y=2,xjust=1,yjust=0,legend=c("Full library","Restricted library"),lty=c(1,1), lwd=c(2,2),col=cols[1:2],cex=0.8,bty="n")
par(op)

dev.off()



