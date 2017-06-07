
#-------------------------------
# FigS4-3-haiti2-cvSL.R
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
library(tmleAb)

# load general plotting functions
library(ggplot2)
library(r2weight)
library(RColorBrewer)
source("~/slabcurves/src/ab_plot_cvSL.R")
source("~/slabcurves/src/ab_plot_cvR2.R")



#-------------------------------
# load the dataset for analysis
#-------------------------------
# d.haiti <- read.csv("~/dropbox/articles/antibody-curves/data/enterics/haiti-enterics-public.csv")

d <- haiti_enterics
X <- data.frame(Age=d$agey)

#-------------------------------
# fit cross-validated SL
#-------------------------------

SL.library <- c("SL.mean","SL.glm","SL.Yman2016","SL.gam","SL.loess","SL.randomForest","SL.polymars")

# giardia
set.seed(982375)
CVvsp5 <- cvSLAb(Y=log10(d$vsp5+1),X=X,id=d$id,family=gaussian(),V=10,SL.library=SL.library,gamdf=2:6)

# noro GI
set.seed(2359)
CVnorogi <- cvSLAb(Y=log10(d$norogi+1),X=X,id=d$id,family=gaussian(),V=10,SL.library=SL.library,gamdf=2:6)

#------------------------------
# print results to log
#------------------------------

# giardia
summary(CVvsp5)

# noro GI
summary(CVnorogi)

#-------------------------------
# plot the CV Risk estimates
#-------------------------------
library(RColorBrewer)
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cols <- cbPalette[c(1:4,6:8)]
cols <- c(cols,brewer.pal(8,"Dark2")[3])

# pdf("~/dropbox/articles/antibody-curves/results/figs/haiti-cvSL-vsp5.pdf")
# ab_plot_cvSL(CVvsp5,ylab="10-fold Cross-validated MSE",title="Giardia intestinalis VSP-5, Haiti",col=cols,ylim=c(0.4,1.4))
# dev.off()
# 
# pdf("~/dropbox/articles/antibody-curves/results/figs/haiti-cvSL-ngii.pdf")
# ab_plot_cvSL(CVnorogi,ylab="10-fold Cross-validated MSE",title="Norovirus GII.4 NO, Haiti",col=cols,ylim=c(0.4,1.4))
# dev.off()


#-------------------------------
# plot the CV R2 estimates
#-------------------------------

pdf("~/dropbox/articles/antibody-curves/results/figs/haiti-cvR2-vsp5.pdf")
ab_plot_cvR2(CVvsp5,X=X,ylab="10-fold Cross-validated R-squared",title=expression(paste(italic('Giardia intestinalis')," VSP-5, Haiti")),col=cols,ylim=c(0,0.6))
dev.off()

pdf("~/dropbox/articles/antibody-curves/results/figs/haiti-cvR2-norogi.pdf")
ab_plot_cvR2(CVnorogi,X=X,ylab="10-fold Cross-validated R-squared",title="Norovirus GI.4, Haiti",col=cols,ylim=c(0,0.6))
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

#-------------------------------
# Giardia VSP-5
#-------------------------------
# full library
SL.library <- c("SL.mean","SL.glm","SL.Yman2016","SL.gam","SL.loess","SL.randomForest","SL.polymars")
set.seed(543221)
SLvsp5full <- agecurveAb(Y=log10(d$vsp5+1),Age=d$agey,id=d$id,SL.library=SL.library,gamdf=2:6)

# restricted library
SL.library <- c("SL.mean","SL.glm","SL.Yman2016","SL.gam","SL.loess")
set.seed(543221)
SLvsp5res <- agecurveAb(Y=log10(d$vsp5+1),Age=d$agey,id=d$id,SL.library=SL.library,gamdf=2:6)


# plot age-dependent antibody curves and means for two SL libraries
pdf("~/dropbox/articles/antibody-curves/results/figs/haiti2-cvSLcurves-vsp5.pdf",width=5,height=5)

ytics <- seq(0,5,by=1)
xtics <- seq(0,12,by=2)

# EYax
op <- par(mar=c(4,5,4,0)+0.1,xpd=T)
plot(SLvsp5full$Age,SLvsp5full$Y,col=alpha("black",alpha=0.4),pch=16,cex=0.3,
     ylim=c(0,5),ylab="",yaxt="n",
     xlim=c(0,max(xtics)+1),xlab="",xaxt="n",
     bty="n",las=1
)
mtext("Luminex response (MFI-background)",side=2,line=3)
# mtext("d",adj=1,line=2,at=-2,font=2,cex=1.75)
# mtext("Intervention Period",adj=0,line=3,at=0,cex=1.5)
mtext(expression(paste(italic('Giardia intestinalis')," VSP-5, Haiti")),side=3,line=2)
mtext("Super learner ensemble prediction",side=3,line=0.5)
mtext("Age, years",side=1,line=2.5)
axis(2,at=0:5,labels=c(
  expression(10^0),
  expression(10^1),
  expression(10^2),
  expression(10^3),
  expression(10^4),
  expression(10^5)
), las=1,cex.axis=1.25
)
# segments(x0=min(xtics),x1=max(xtics),y0=ytics,lty=2,col="gray70")
axis(1,at=xtics,cex.axis=1.5)
lines(SLvsp5full$Age,SLvsp5full$pY,col=cols[1],lwd=1)
lines(SLvsp5res$Age,SLvsp5res$pY,col=cols[2],lwd=1)

legend(x=12,y=0,xjust=1,yjust=0,legend=c("Full library","Restricted library"),lty=c(1,1), lwd=c(2,2),col=cols[1:2],cex=0.8,bty="n")
par(op)


dev.off()


#-------------------------------
# Noro GI.4
#-------------------------------
# full library
SL.library <- c("SL.mean","SL.glm","SL.Yman2016","SL.gam","SL.loess","SL.randomForest","SL.polymars")
set.seed(543221)
SLnorofull <- agecurveAb(Y=log10(d$norogi+1),Age=d$agey,id=d$id,SL.library=SL.library,gamdf=2:6)

# restricted library
SL.library <- c("SL.mean","SL.glm","SL.Yman2016","SL.gam","SL.loess")
set.seed(543221)
SLnorores <- agecurveAb(Y=log10(d$norogi+1),Age=d$agey,id=d$id,SL.library=SL.library,gamdf=2:6)


# plot age-dependent antibody curves and means for two SL libraries
pdf("~/dropbox/articles/antibody-curves/results/figs/haiti2-cvSLcurves-norogi.pdf",width=5,height=5)

ytics <- seq(0,5,by=1)
xtics <- seq(0,12,by=2)

# EYax
op <- par(mar=c(4,5,4,0)+0.1,xpd=T)
plot(SLnorofull$Age,SLnorofull$Y,col=alpha("black",alpha=0.4),pch=16,cex=0.3,
     ylim=c(0,5),ylab="",yaxt="n",
     xlim=c(0,max(xtics)+1),xlab="",xaxt="n",
     bty="n",las=1
)
mtext("Luminex response (MFI-background)",side=2,line=3)
# mtext("d",adj=1,line=2,at=-2,font=2,cex=1.75)
# mtext("Intervention Period",adj=0,line=3,at=0,cex=1.5)
mtext("Norovirus GI.4, Haiti",side=3,line=2)
mtext("Super learner ensemble prediction",side=3,line=0.5)
mtext("Age, years",side=1,line=2.5)
axis(2,at=0:5,labels=c(
  expression(10^0),
  expression(10^1),
  expression(10^2),
  expression(10^3),
  expression(10^4),
  expression(10^5)
), las=1,cex.axis=1.25
)
# segments(x0=min(xtics),x1=max(xtics),y0=ytics,lty=2,col="gray70")
axis(1,at=xtics,cex.axis=1.5)
lines(SLnorofull$Age,SLnorofull$pY,col=cols[1],lwd=1)
lines(SLnorores$Age,SLnorores$pY,col=cols[2],lwd=1)

legend(x=12,y=0,xjust=1,yjust=0,legend=c("Full library","Restricted library"),lty=c(1,1), lwd=c(2,2),col=cols[1:2],cex=0.8,bty="n")
par(op)


dev.off()





