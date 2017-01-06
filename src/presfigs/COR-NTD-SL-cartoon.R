#-----------------------------------
# COR-NTD-SL-cartoon.R
#
# build up an ensemble prediction
# showing each model's fit
#-----------------------------------



#-----------------------------------
# preamble
#-----------------------------------
rm(list=ls())
library(SuperLearner)
library(tmle)
library(tmleAb)
library(r2weight)
library(scales)

source("~/slabcurves/src/ab_plot_cvR2.R")

#-------------------------------------------
# load the Hait2 longitudinal data
#-------------------------------------------
library(foreign)
d <- read.dta("~/dropbox/haiti2/data/final/haiti2-long.dta")

#-------------------------------------------
# fit age-antibody curves
#-------------------------------------------
SL.libfull <- c("SL.mean","SL.glm","SL.Yman2016","SL.loess","SL.gam","SL.polymars","SL.randomForest")
SL.libres <- c("SL.mean","SL.glm","SL.gam","SL.loess")

set.seed(3234723)
bm14mean <- agecurveAb(Y=log10(d$bm14),Age=d$agey,id=d$id,SL.library="SL.mean")
set.seed(3234723)
bm14glm <- agecurveAb(Y=log10(d$bm14),Age=d$agey,id=d$id,SL.library="SL.glm")
set.seed(3234723)
bm14yman <- agecurveAb(Y=log10(d$bm14),Age=d$agey,id=d$id,SL.library="SL.Yman2016")
set.seed(3234723)
bm14lo <- agecurveAb(Y=log10(d$bm14),Age=d$agey,id=d$id,SL.library="SL.loess")
set.seed(3234723)
bm14gam <- agecurveAb(Y=log10(d$bm14),Age=d$agey,id=d$id,SL.library="SL.gam")
set.seed(3234723)
bm14mars <- agecurveAb(Y=log10(d$bm14),Age=d$agey,id=d$id,SL.library="SL.polymars")
set.seed(3234723)
bm14rf <- agecurveAb(Y=log10(d$bm14),Age=d$agey,id=d$id,SL.library="SL.randomForest")
set.seed(3234723)
bm14sl <- agecurveAb(Y=log10(d$bm14),Age=d$agey,id=d$id,SL.library=SL.libfull)


#-------------------------------------------
# compute cross-validated risk and R2
#-------------------------------------------
Xa <- data.frame(Age=d$agey)
set.seed(3234723)
bm14cvsl <- cvSLAb(Y=log10(d$bm14),X=Xa,id=d$id,SL.library=SL.libfull,V=10)
bm14cvr2 <- ab_plot_cvR2(cvsl=bm14cvsl,X=Xa)

#-------------------------------------------
# Bm14 curves
#-------------------------------------------

### general plot, to call repeatedly below
bm14shell <- function(x){
  op <- par(mar=c(5,6,2,0)+0.1)
  xtics <- seq(0,12,by=2)
  
  plot(x$Age, x$Y,type="n",
       ylab="",yaxt="n",
       ylim=c(0,5),
       xlab="",xlim=range(xtics),xaxt="n",
       bty="n",las=1
  )
  
  # header
  mtext(expression(paste(italic("W. bancrofti")," Bm14 Luminex Response")),cex=1.5,side=3,line=0)
  
  # axes
  mtext("Age, years",side=1,line=3,cex=1.5)
  #mtext("Luminex Response (MFI-Background)",side=2,line=3.5,cex=1.25)
  axis(1,at=xtics,cex.axis=1.5)
  axis(2,at=0:5,labels=c(
    # expression(10^-1),
    expression(10^0),
    expression(10^1),
    expression(10^2),
    expression(10^3),
    expression(10^4),
    expression(10^5)
  ), las=1,cex.axis=1.5
  )
  
  # plot data (points and fitted lines)
  points(x$Age, x$Y,col=alpha('black',alpha=0.3), pch=19,cex=0.5)
}


cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cols <- cbPalette

# brighter color blind palette:  https://personal.sron.nl/~pault/ 
cblack <- "#000004FF"
cblue <- "#3366AA"
cteal <- "#11AA99"
cgreen <- "#66AA55"
cchartr <- "#CCCC55"
cmagent <- "#992288"
cred <- "#EE3333"
corange <- "#EEA722"
cyellow <- "#FFEE33"
cgrey <- "#777777"
cols=c(cblue,cteal,cgreen,cchartr,corange,cred,cmagent,cblack)

pdf("~/dropbox/articles/antibody-curves/presentations/cor-ntd-2016/figs/slcartoon-0.pdf")
bm14shell(bm14mean)
dev.off()

pdf("~/dropbox/articles/antibody-curves/presentations/cor-ntd-2016/figs/slcartoon-1.pdf")
bm14shell(bm14mean)
lines(bm14mean$Age,bm14mean$pY,col=cols[1],lwd=2)
dev.off()

pdf("~/dropbox/articles/antibody-curves/presentations/cor-ntd-2016/figs/slcartoon-2.pdf")
bm14shell(bm14mean)
lines(bm14mean$Age,bm14mean$pY,col=cols[1],lwd=2)
lines(bm14glm$Age,bm14glm$pY,col=cols[2],lwd=2)
dev.off()


pdf("~/dropbox/articles/antibody-curves/presentations/cor-ntd-2016/figs/slcartoon-3.pdf")
bm14shell(bm14mean)
lines(bm14mean$Age,bm14mean$pY,col=cols[1],lwd=2)
lines(bm14glm$Age,bm14glm$pY,col=cols[2],lwd=2)
lines(bm14yman$Age,bm14yman$pY,col=cols[3],lwd=2)
dev.off()

pdf("~/dropbox/articles/antibody-curves/presentations/cor-ntd-2016/figs/slcartoon-yman.pdf")
bm14shell(bm14mean)
lines(bm14yman$Age,bm14yman$pY,col=cols[3],lwd=3)
dev.off()



pdf("~/dropbox/articles/antibody-curves/presentations/cor-ntd-2016/figs/slcartoon-4.pdf")
bm14shell(bm14mean)
lines(bm14mean$Age,bm14mean$pY,col=cols[1],lwd=2)
lines(bm14glm$Age,bm14glm$pY,col=cols[2],lwd=2)
lines(bm14yman$Age,bm14yman$pY,col=cols[3],lwd=2)
lines(bm14lo$Age,bm14lo$pY,col=cols[4],lwd=2)
dev.off()

pdf("~/dropbox/articles/antibody-curves/presentations/cor-ntd-2016/figs/slcartoon-5.pdf")
bm14shell(bm14mean)
lines(bm14mean$Age,bm14mean$pY,col=cols[1],lwd=2)
lines(bm14glm$Age,bm14glm$pY,col=cols[2],lwd=2)
lines(bm14yman$Age,bm14yman$pY,col=cols[3],lwd=2)
lines(bm14lo$Age,bm14lo$pY,col=cols[4],lwd=2)
lines(bm14gam$Age,bm14gam$pY,col=cols[5],lwd=2)
dev.off()

pdf("~/dropbox/articles/antibody-curves/presentations/cor-ntd-2016/figs/slcartoon-6.pdf")
bm14shell(bm14mean)
lines(bm14mean$Age,bm14mean$pY,col=cols[1],lwd=2)
lines(bm14glm$Age,bm14glm$pY,col=cols[2],lwd=2)
lines(bm14yman$Age,bm14yman$pY,col=cols[3],lwd=2)
lines(bm14lo$Age,bm14lo$pY,col=cols[4],lwd=2)
lines(bm14gam$Age,bm14gam$pY,col=cols[5],lwd=2)
lines(bm14mars$Age,bm14mars$pY,col=cols[6],lwd=2)
dev.off()

pdf("~/dropbox/articles/antibody-curves/presentations/cor-ntd-2016/figs/slcartoon-7.pdf")
bm14shell(bm14mean)
lines(bm14mean$Age,bm14mean$pY,col=cols[1],lwd=2)
lines(bm14glm$Age,bm14glm$pY,col=cols[2],lwd=2)
lines(bm14yman$Age,bm14yman$pY,col=cols[3],lwd=2)
lines(bm14lo$Age,bm14lo$pY,col=cols[4],lwd=2)
lines(bm14gam$Age,bm14gam$pY,col=cols[5],lwd=2)
lines(bm14mars$Age,bm14mars$pY,col=cols[6],lwd=2)
lines(bm14rf$Age,bm14rf$pY,col=cols[7],lwd=2)
dev.off()


pdf("~/dropbox/articles/antibody-curves/presentations/cor-ntd-2016/figs/slcartoon-rf.pdf")
bm14shell(bm14mean)
lines(bm14rf$Age,bm14rf$pY,col=cols[7],lwd=3)
dev.off()

pdf("~/dropbox/articles/antibody-curves/presentations/cor-ntd-2016/figs/slcartoon-8.pdf")
bm14shell(bm14mean)
lines(bm14mean$Age,bm14mean$pY,col=cols[1],lwd=2)
lines(bm14glm$Age,bm14glm$pY,col=cols[2],lwd=2)
lines(bm14yman$Age,bm14yman$pY,col=cols[3],lwd=2)
lines(bm14lo$Age,bm14lo$pY,col=cols[4],lwd=2)
lines(bm14gam$Age,bm14gam$pY,col=cols[5],lwd=2)
lines(bm14mars$Age,bm14mars$pY,col=cols[6],lwd=2)
lines(bm14rf$Age,bm14rf$pY,col=cols[7],lwd=2)
lines(bm14sl$Age,bm14sl$pY,col=cols[8],lwd=3)
dev.off()

pdf("~/dropbox/articles/antibody-curves/presentations/cor-ntd-2016/figs/slcartoon-9.pdf")
bm14shell(bm14mean)
lines(bm14sl$Age,bm14sl$pY,col=cols[8],lwd=2)
dev.off()


