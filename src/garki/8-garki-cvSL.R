
#-------------------------------
# 8-garki-cvSL.R
#
# Compute the cross-validated
# risk for the super leaner
# and its constituent algorithms
#
# do calculations for surveys
# 3-5 (intervention period)
# in the intervention villages
#
#-------------------------------

#-------------------------------
# preamble
#-------------------------------
rm(list=ls())
library(SuperLearner)
library(tmleAb)
library(r2weight)
library(RColorBrewer)
library(scales)
library(ggplot2)

#-------------------------------
# load the serology dataset
#-------------------------------
# d <- read.csv("~/dropbox/articles/antibody-curves/data/garki/final/garki-sero.csv")

d <- garki_sero

d$mdate <- as.Date(d$mdate,"%d %b %Y")


# add 2 control village names
d$vname <- factor(d$vname,levels=c(levels(d$vname),"Nabanawa","Ajura"))
d$vname[d$village==552] <- "Nabanawa"
d$vname[d$village==553] <- "Ajura"
d$vname <- factor(d$vname)

# set sex to factor
d$sex <- as.factor(d$sex)

# for age exactly equal to 0, set it equal to 0.001
# to prevent the Yman 2016 model from blowing up
# (model is undefined at age=0)
d$ageyrs[d$ageyrs<=0] <- 0.001

# subset to ages 0-20
d <- subset(d,ageyrs<=20)


#-------------------------------
# Serological survey timing:
# 1-2 : pre-intervention
# 3-5 : intervention period
# 6-8 : post-intervention
#-------------------------------

# subset the data by group and intervention period for convenience
# restrict to non-missing obs wrt Pf titres and age
dc <- d[d$tr=="Control" & (d$serosvy>=3 & d$serosvy<=5),]
  dc <- subset(dc,!is.na(dc$ageyrs) & !is.na(dc$ifatpftitre))
dt <- d[d$tr=="Intervention" & (d$serosvy>=3 & d$serosvy<=5),]
  dt <- subset(dt,!is.na(dt$ageyrs) & !is.na(dt$ifatpftitre))

  
#-------------------------------
# Control Villages
#-------------------------------

# fit cross-validated SL
# with a full library
SL.library <- c("SL.mean","SL.glm","SL.Yman2016","SL.gam","SL.loess","SL.randomForest","SL.polymars")
  
# fit the cross-validated super learner
# with just Age as the predictor
set.seed(62522)
CVcfit <- ab_cvSL(Y=log10(dc$ifatpftitre+1),Age=dc$ageyrs,id=dc$id,family=gaussian(),V=10,SL.library=SL.library)
  
# fit the cross-validated super learner
# accounting for other covariates
Wc <- dc[,c("sex","wetseason","vname")]
set.seed(62522)
CVcfitm <- ab_cvSL(Y=log10(dc$ifatpftitre+1),Age=dc$ageyrs,id=dc$id,W=Wc,family=gaussian(),V=10,SL.library=SL.library)
  
# fit cross-validated SL
# with a restricted library
# used in the primary analysis
SL.library <- c("SL.mean","SL.glm","SL.Yman2016","SL.gam","SL.loess")
  
# fit the cross-validated super learner
# with just Age as the predictor
set.seed(62522)
CVcfitr <- ab_cvSL(Y=log10(dc$ifatpftitre+1),Age=dc$ageyrs,id=dc$id,family=gaussian(),V=10,SL.library=SL.library)
  
# fit the cross-validated super learner
# accounting for other covariates
set.seed(62522)
CVcfitrm <- ab_cvSL(Y=log10(dc$ifatpftitre+1),Age=dc$ageyrs,id=dc$id,W=Wc,family=gaussian(),V=10,SL.library=SL.library)

#-------------------------------
# Intervention Villages
#-------------------------------

# fit cross-validated SL
# with a full library
SL.library <- c("SL.mean","SL.glm","SL.Yman2016","SL.gam","SL.loess","SL.randomForest","SL.polymars")

# fit the cross-validated super learner
# with just Age as the predictor
set.seed(62522)
CVtfit <- ab_cvSL(Y=log10(dt$ifatpftitre+1),Age=dt$ageyrs,id=dt$id,family=gaussian(),V=10,SL.library=SL.library)

# fit the cross-validated super learner
# accounting for other covariates
Wt <- dt[,c("sex","wetseason","vname")]
set.seed(62522)
CVtfitm <- ab_cvSL(Y=log10(dt$ifatpftitre+1),Age=dt$ageyrs,id=dt$id,W=Wt,family=gaussian(),V=10,SL.library=SL.library)

# fit cross-validated SL
# with a restricted library
# used in the primary analysis
SL.library <- c("SL.mean","SL.glm","SL.Yman2016","SL.gam","SL.loess")

# fit the cross-validated super learner
# with just Age as the predictor
set.seed(62522)
CVtfitr <- ab_cvSL(Y=log10(dt$ifatpftitre+1),Age=dt$ageyrs,id=dt$id,family=gaussian(),V=10,SL.library=SL.library)

# fit the cross-validated super learner
# accounting for other covariates
set.seed(62522)
CVtfitrm <- ab_cvSL(Y=log10(dt$ifatpftitre+1),Age=dt$ageyrs,id=dt$id,W=Wt,family=gaussian(),V=10,SL.library=SL.library)

#------------------------------
# print results to log
#------------------------------

# control villages, age only
summary(CVcfit)

# control villages, multivariate
summary(CVcfitm)


# restricted libraries

# control villages, age only
summary(CVcfitr)

# control villages, multivariate
summary(CVcfitrm)


# intervention villages, age only
summary(CVtfit)

# intervention villages, multivariate
summary(CVtfitm)


# restricted libraries

# intervention villages, age only
summary(CVtfitr)

# intervention villages, multivariate
summary(CVtfitrm)



#-------------------------------
# plot the CV Risk estimates
#-------------------------------
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cols <- cbPalette[c(1:4,6:8)]
cols <- c(cols,brewer.pal(8,"Dark2")[3])

# pdf("~/dropbox/articles/antibody-curves/results/figs/garki-cvSL-age.pdf")
# ab_plot_cvSL(CVtfit,ylab="10-fold Cross-validated MSE",title="P. falciparum, high transmission",col=cols)
# dev.off()
# 
# pdf("~/dropbox/articles/antibody-curves/results/figs/garki-cvSL-mv.pdf")
# ab_plot_cvSL(CVtfitm,ylab="10-fold Cross-validated MSE",title="P. falciparum, Garki Project (multivariate)",col=cols)
# dev.off()

#-------------------------------
# plot the CV R2 estimates
#-------------------------------
pWc <- design_matrix(Wc)
Xc <- data.frame(Age=dc$ageyrs)
Xcm <- data.frame(Age=dc$ageyrs,pWc)

pdf("~/dropbox/articles/antibody-curves/results/figs/garki-cvR2-age-control.pdf")
ab_plot_cvR2(CVcfit,X=Xc,ylab="10-fold Cross-validated R-squared",title=expression(paste(italic('P. falciparum'),", Garki Project Control Villages")),col=cols,ylim=c(0,0.6))
dev.off()

pdf("~/dropbox/articles/antibody-curves/results/figs/garki-cvR2-mv-control.pdf")
ab_plot_cvR2(CVcfitm,X=Xcm,ylab="10-fold Cross-validated R-squared",title=expression(paste(italic('P. falciparum'),", Garki Project Control Villages (multivariate)")),col=cols,ylim=c(0,0.6))
dev.off()

pWt <- design_matrix(Wt)
Xt <- data.frame(Age=dt$ageyrs)
Xtm <- data.frame(Age=dt$ageyrs,pWt)

pdf("~/dropbox/articles/antibody-curves/results/figs/garki-cvR2-age-int.pdf")
ab_plot_cvR2(CVtfit,X=Xt,ylab="10-fold Cross-validated R-squared",title=expression(paste(italic('P. falciparum'),", Garki Project Intervention Villages")),col=cols,ylim=c(0,0.6))
dev.off()

pdf("~/dropbox/articles/antibody-curves/results/figs/garki-cvR2-int-mv.pdf")
ab_plot_cvR2(CVtfitm,X=Xtm,ylab="10-fold Cross-validated R-squared",title=expression(paste(italic('P. falciparum'),", Garki Project Intervention Villages (multivariate)")),col=cols,ylim=c(0,0.6))
dev.off()

#-------------------------------
# plot the CV R2 estimates
# restricted library
# intervention only (not used)
#-------------------------------

pdf("~/dropbox/articles/antibody-curves/results/figs/garki-cvR2-ager-int.pdf")
ab_plot_cvR2(CVtfitr,X=Xt,ylab="10-fold Cross-validated R-squared",title="P. falciparum, Garki Project",col=cols[1:6],ylim=c(0,0.6))
dev.off()

pdf("~/dropbox/articles/antibody-curves/results/figs/garki-cvR2-mvr-int.pdf")
ab_plot_cvR2(CVtfitrm,X=Xtm,ylab="10-fold Cross-validated R-squared",title="P. falciparum, Garki Project (multivariate)",col=cols[1:6],ylim=c(0,0.6))
dev.off()


#-------------------------------
# plot the age-antibody curves
# for the full SL and the 
# restricted SL library to
# compare EYxa
# Control
#-------------------------------

# full library
SL.library <- c("SL.mean","SL.glm","SL.Yman2016","SL.gam","SL.loess","SL.randomForest","SL.polymars")
set.seed(62522)
SLcfull <- ab_agecurve(Y=log10(dc$ifatpftitre+1),Age=dc$ageyrs,id=dc$id,SL.library=SL.library)

# restricted library
SL.library <- c("SL.mean","SL.glm","SL.Yman2016","SL.gam","SL.loess")
set.seed(62522)
SLcres <- ab_agecurve(Y=log10(dc$ifatpftitre+1),Age=dc$ageyrs,id=dc$id,SL.library=SL.library)


# plot age-dependent antibody curves and means for two SL libraries
pdf("~/dropbox/articles/antibody-curves/results/figs/garki-cvSLcurves-control.pdf",width=5,height=5)

ytics <- seq(0,4,by=1)
xtics <- seq(0,20,by=5)

op <- par(mar=c(4,5,4,0)+0.1,xpd=T)
plot(SLcfull$Age,SLcfull$Y,col=alpha("black",alpha=0.4),pch=16,cex=0.2,
     ylim=c(0,4),ylab="",yaxt="n",
     xlim=c(0,max(xtics)+1),xlab="",xaxt="n",
     bty="n",las=1
)
mtext(expression(paste(italic('P. falciparum')," IFA antibody titre")),side=2,line=3)
# mtext("b",adj=1,line=2,at=-3,font=2,cex=1.75)
# mtext("Intervention Period",adj=0,line=3,at=0,cex=1.5)
mtext(expression(paste(italic('P. falciparum'),", Garki Project Control Villages")),side=3,line=2)
mtext("Super learner ensemble prediction",side=3,line=0.5)
mtext("Age, years",side=1,line=2.5)
axis(2,at=0:4,labels=c(
  expression(10^0),
  expression(10^1),
  expression(10^2),
  expression(10^3),
  expression(10^4)
), las=1,cex.axis=1.25
)
# segments(x0=min(xtics),x1=max(xtics),y0=ytics,lty=2,col="gray70")
axis(1,at=xtics,cex.axis=1.5)
lines(SLcfull$pYframe$Age,SLcfull$pYframe$pY,col=cols[1],lwd=1)
lines(SLcres$pYframe$Age,SLcres$pYframe$pY,col=cols[2],lwd=1)

legend(x=20,y=0,xjust=1,yjust=0,legend=c("Full library","Restricted library"),lty=c(1,1), lwd=c(2,2),col=cols[1:2],cex=0.8,bty="n")
par(op)

dev.off()


#-------------------------------
# plot the age-antibody curves
# for the full SL and the 
# restricted SL library to
# compare EYxa
# Intervention
#-------------------------------

# full library
SL.library <- c("SL.mean","SL.glm","SL.Yman2016","SL.gam","SL.loess","SL.randomForest","SL.polymars")
set.seed(62522)
SLfull <- ab_agecurve(Y=log10(dt$ifatpftitre+1),Age=dt$ageyrs,id=dt$id,SL.library=SL.library)

# restricted library
SL.library <- c("SL.mean","SL.glm","SL.Yman2016","SL.gam","SL.loess")
set.seed(62522)
SLres <- ab_agecurve(Y=log10(dt$ifatpftitre+1),Age=dt$ageyrs,id=dt$id,SL.library=SL.library)


# plot age-dependent antibody curves and means for two SL libraries
pdf("~/dropbox/articles/antibody-curves/results/figs/garki-cvSLcurves-int.pdf",width=5,height=5)

ytics <- seq(0,4,by=1)
xtics <- seq(0,20,by=5)

op <- par(mar=c(4,5,4,0)+0.1,xpd=T)
plot(SLfull$Age,SLfull$Y,col=alpha("black",alpha=0.4),pch=16,cex=0.2,
     ylim=c(0,4),ylab="",yaxt="n",
     xlim=c(0,max(xtics)+1),xlab="",xaxt="n",
     bty="n",las=1
)
mtext(expression(paste(italic('P. falciparum')," IFA antibody titre")),side=2,line=3)
# mtext("d",adj=1,line=2,at=-3,font=2,cex=1.75)
# mtext("Intervention Period",adj=0,line=3,at=0,cex=1.5)
mtext(expression(paste(italic('P. falciparum'),", Garki Project Intervention Villages")),side=3,line=2)
mtext("Super learner ensemble prediction",side=3,line=0.5)
mtext("Age, years",side=1,line=2.5)
axis(2,at=0:4,labels=c(
  expression(10^0),
  expression(10^1),
  expression(10^2),
  expression(10^3),
  expression(10^4)
), las=1,cex.axis=1.25
)
# segments(x0=min(xtics),x1=max(xtics),y0=ytics,lty=2,col="gray70")
axis(1,at=xtics,cex.axis=1.5)
lines(SLfull$pYframe$Age,SLfull$pYframe$pY,col=cols[1],lwd=1)
lines(SLres$pYframe$Age,SLres$pYframe$pY,col=cols[2],lwd=1)

legend(x=20,y=0,xjust=1,yjust=0,legend=c("Full library","Restricted library"),lty=c(1,1), lwd=c(2,2),col=cols[1:2],cex=0.8,bty="n")
par(op)

dev.off()


