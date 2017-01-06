#-----------------------------------
# COR-NTD-noro-intro.R
#
# plot norovirus GII.4 NO antibody
# data from Haiti and the USA
#-----------------------------------

#-----------------------------------
# preamble
#-----------------------------------
rm(list=ls())
library(SuperLearner)
library(tmle)
library(tmleAb)
library(beeswarm)
library(scales)


#-------------------------------------------
# Load Norovirus data
#-------------------------------------------
data(usa_enterics)
data(haiti_enterics)
d.usa <- usa_enterics
d.haiti <- haiti_enterics

# recode ages by adding 1/2 year
# in USA children >1 y to avoid bias
d.usa$age[d.usa$age>=1] <- d.usa$age[d.usa$age>=1]+0.5

# Limit the haiti data to children
# <= 5.5 years old, which is the
# maximum age for the children
# in the USA sample
d.hai  <- subset(d.haiti,agey<=5.5)

# add 1 to any norovirus MFI values that are zero
d.hai$norogii[d.hai$norogii<=0] <- 1

# append the Haiti and USA 
# data to calculate differences 
# in means and for beeswarm plot
d.usa$agey <- d.usa$age
d.usa$haiti <- 0
d.hai$haiti <- 1
common.vars <- c("haiti","id","agey","cp17","cp23","vsp5","leca","etec","salb","norogi","norogii")
d.all <- rbind(subset(d.hai,select=common.vars),subset(d.usa,select=common.vars)	)
d.all$country <- as.factor(ifelse(d.all$haiti==1,"Haiti","USA"))
d.all$country <- factor(d.all$country,levels=c("USA","Haiti"))


#-------------------------------
# beeswarm plot - Haiti
#-------------------------------

cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cols <- cbPalette[c(7,6)]

# brighter color blind palette:  https://personal.sron.nl/~pault/ 
# colblack <- "#000004FF"
# colblue <- "#3366AA"
# colteal <- "#11AA99"
# colgreen <- "#66AA55"
# colchartr <- "#CCCC55"
# colmagent <- "#992288"
# colred <- "#EE3333"
# colorange <- "#EEA722"
# colyellow <- "#FFEE33"
# colgrey <- "#777777"
# cols=c(colred,colblue)

pdf("~/dropbox/articles/antibody-curves/presentations/cor-ntd-2016/figs/haiti2-norogii-beeswarm.pdf",width=4,height=5)
ytics <- 0:5
ns <- table(d.all$country)
op <- par(mar=c(1,6,3,0)+0.1)
midpts <- beeswarm(log10(d.all$norogii[d.all$country=="Haiti"]),
                   #corral='random',corralWidth = 0.9,
                   labels=NA,
                   col=alpha(cols[2],alpha=0.8),
                   pch=16,cex=0.5,
                   ylab="", yaxt="n",ylim=range(ytics),
                   xlab="",
                   bty="n"
)
mtext("Norovirus Group II, Haiti",side=3,line=1,cex=1.25)
mtext("Luminex Response (MFI-Background)",side=2,line=3.5,cex=1.25)
# mtext(c("USA","Haiti"),side=3,line=-1,at=c(1,2),col=cols)
# mtext(paste("N =",ns),side=1,line=0,at=c(1,2),col="gray40",cex=0.9)
axis(2,at=0:5,labels=c(
  expression(10^0),
  expression(10^1),
  expression(10^2),
  expression(10^3),
  expression(10^4),
  expression(10^5)
), las=1,cex.axis=1.5
)
dev.off()


#-------------------------------------------
# SL fit of noro data for Haiti and USA
#-------------------------------------------

# SL library used in eLife paper
SL.library <- c("SL.mean","SL.glm","SL.Yman2016","SL.gam","SL.loess")

set.seed(12345)
haiti_sl <- agecurveAb(Y=log10(d.hai$norogii),Age=d.hai$agey,id=d.hai$id,SL.library=SL.library)
set.seed(12345)
usa_sl <- agecurveAb(Y=log10(d.usa$norogii),Age=d.usa$agey,id=d.usa$id,SL.library=SL.library)

#-------------------------------------------
# Haiti points, no lines
#-------------------------------------------
pdf("~/dropbox/articles/antibody-curves/presentations/cor-ntd-2016/figs/haiti2-norogii-age1.pdf",width=6,height=5)

op <- par(mar=c(5,6,3,1)+0.1)
xtics <- seq(0,6,by=1)

plot(d.hai$agey, log10(d.hai$norogii),type="n",
     ylab="",yaxt="n",
     ylim=c(0,5),
     xlab="",xlim=range(xtics),xaxt="n",
     bty="n",las=1
)

# header
mtext("Norovirus Group II, Haiti",side=3,line=1,cex=1.25)

# axes
mtext("Age, years",side=1,line=3,cex=1.5)
mtext("Luminex Response (MFI-Background)",side=2,line=3.5,cex=1)
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

# plot data
points(d.hai$agey, log10(d.hai$norogii),col=alpha(cols[2],alpha=0.7), pch=16,cex=0.6)

dev.off()


#-------------------------------------------
# Haiti points, Haiti line
#-------------------------------------------
pdf("~/dropbox/articles/antibody-curves/presentations/cor-ntd-2016/figs/haiti2-norogii-age2.pdf",width=6,height=5)

op <- par(mar=c(5,6,3,1)+0.1)
xtics <- seq(0,6,by=1)

plot(d.hai$agey, log10(d.hai$norogii),type="n",
     ylab="",yaxt="n",
     ylim=c(0,5),
     xlab="",xlim=range(xtics),xaxt="n",
     bty="n",las=1
)

# header
mtext("Norovirus Group II, Haiti",side=3,line=1,cex=1.25)

# axes
mtext("Age, years",side=1,line=3,cex=1.5)
mtext("Luminex Response (MFI-Background)",side=2,line=3.5,cex=1)
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

# plot data
points(d.hai$agey, log10(d.hai$norogii),col=alpha(cols[2],alpha=0.3), pch=16,cex=0.6)

lines(haiti_sl$Age,haiti_sl$pY,col=cols[2],lwd=2)

dev.off()

#-------------------------------------------
# All points, both lines
#-------------------------------------------
pdf("~/dropbox/articles/antibody-curves/presentations/cor-ntd-2016/figs/haiti2-norogii-age3.pdf",width=6,height=5)

op <- par(mar=c(5,6,3,1)+0.1,xpd=TRUE)
xtics <- seq(0,6,by=1)

plot(d.hai$agey, log10(d.hai$norogii),type="n",
     ylab="",yaxt="n",
     ylim=c(0,5),
     xlab="",xlim=range(xtics),xaxt="n",
     bty="n",las=1
)

# header
mtext("Norovirus Group II",side=3,line=1,cex=1.25)

# axes
mtext("Age, years",side=1,line=3,cex=1.5)
mtext("Luminex Response (MFI-Background)",side=2,line=3.5,cex=1)
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

# plot data
points(d.hai$agey, log10(d.hai$norogii),col=alpha(cols[2],alpha=0.3), pch=16,cex=0.6)
points(jitter(d.usa$agey,3), log10(d.usa$norogii),col=alpha(cols[1],alpha=0.3), pch=16,cex=0.6)

lines(haiti_sl$Age,haiti_sl$pY,col=cols[2],lwd=2)
lines(usa_sl$Age,usa_sl$pY,col=cols[1],lwd=2)
text(5.5,usa_sl$pY[nrow(usa_sl)],"USA",col=cols[1],adj=0,pos=4)
text(5.5,haiti_sl$pY[nrow(haiti_sl)],"Haiti",col=cols[2],adj=0,pos=4)

dev.off()

#-------------------------------------------
# ETEC
#-------------------------------------------
pdf("~/dropbox/articles/antibody-curves/presentations/cor-ntd-2016/figs/haiti2-etec.pdf",width=6,height=5)

op <- par(mar=c(5,6,3,1)+0.1,xpd=TRUE)
xtics <- seq(0,6,by=1)

plot(haiti.etec$Age, haiti.etec$Y,type="n",
     ylab="",yaxt="n",
     ylim=c(0,5),
     xlab="",xlim=range(xtics),xaxt="n",
     bty="n",las=1
)

# header
mtext(expression(paste("ETEC heat labile toxin ",beta," subunit")),side=3,line=1,cex=1.25)

# axes
mtext("Age, years",side=1,line=3,cex=1.5)
mtext("Luminex Response (MFI-Background)",side=2,line=3.5,cex=1)
axis(1,at=xtics,cex.axis=1.5)
axis(2,at=0:5,labels=c(
  expression(10^0),
  expression(10^1),
  expression(10^2),
  expression(10^3),
  expression(10^4),
  expression(10^5)
), las=1,cex.axis=1.5
)

# plot data
points(haiti.etec$Age, haiti.etec$Y,col=alpha(cols[2],alpha=0.3), pch=16,cex=0.6)
points(jitter(usa.etec$Age,3), usa.etec$Y,col=alpha(cols[1],alpha=0.3), pch=16,cex=0.6)

lines(haiti.etec$Age,haiti.etec$pY,col=cols[2],lwd=2)
lines(usa.etec$Age,usa.etec$pY,col=cols[1],lwd=2)
text(5.5,usa.etec$pY[nrow(usa.etec)],"USA",col=cols[1],adj=0,pos=4)
text(5.5,haiti.etec$pY[nrow(haiti.etec)],"Haiti",col=cols[2],adj=0,pos=4)

dev.off()


#-------------------------------------------
# Giardia VSP-5
#-------------------------------------------
pdf("~/dropbox/articles/antibody-curves/presentations/cor-ntd-2016/figs/haiti2-giardia.pdf",width=6,height=5)

op <- par(mar=c(5,6,3,1)+0.1,xpd=TRUE)
xtics <- seq(0,6,by=1)

plot(haiti.giar$Age, haiti.giar$Y,type="n",
     ylab="",yaxt="n",
     ylim=c(0,5),
     xlab="",xlim=range(xtics),xaxt="n",
     bty="n",las=1
)

# header
mtext(expression(paste(italic('Giardia intestinalis'), " VSP-5")),side=3,line=1,cex=1.25)

# axes
mtext("Age, years",side=1,line=3,cex=1.5)
mtext("Luminex Response (MFI-Background)",side=2,line=3.5,cex=1)
axis(1,at=xtics,cex.axis=1.5)
axis(2,at=0:5,labels=c(
  expression(10^0),
  expression(10^1),
  expression(10^2),
  expression(10^3),
  expression(10^4),
  expression(10^5)
), las=1,cex.axis=1.5
)

# plot data
points(haiti.giar$Age, haiti.giar$Y,col=alpha(cols[2],alpha=0.3), pch=16,cex=0.6)
points(jitter(usa.giar$Age,3), usa.giar$Y,col=alpha(cols[1],alpha=0.3), pch=16,cex=0.6)

lines(haiti.giar$Age,haiti.giar$pY,col=cols[2],lwd=2)
lines(usa.giar$Age,usa.giar$pY,col=cols[1],lwd=2)
text(5.5,usa.giar$pY[nrow(usa.giar)],"USA",col=cols[1],adj=0,pos=4)
text(5.5,haiti.giar$pY[nrow(haiti.giar)],"Haiti",col=cols[2],adj=0,pos=4)

dev.off()


#-------------------------------------------
# Crypto - cp23
#-------------------------------------------
pdf("~/dropbox/articles/antibody-curves/presentations/cor-ntd-2016/figs/haiti2-crypto.pdf",width=6,height=5)

op <- par(mar=c(5,6,3,1)+0.1,xpd=TRUE)
xtics <- seq(0,6,by=1)

plot(haiti.cp23$Age, haiti.cp23$Y,type="n",
     ylab="",yaxt="n",
     ylim=c(0,5),
     xlab="",xlim=range(xtics),xaxt="n",
     bty="n",las=1
)

# header
mtext(expression(paste(italic('Cryptosporidium parvum'), " Cp23")),side=3,line=1,cex=1.25)

# axes
mtext("Age, years",side=1,line=3,cex=1.5)
mtext("Luminex Response (MFI-Background)",side=2,line=3.5,cex=1)
axis(1,at=xtics,cex.axis=1.5)
axis(2,at=0:5,labels=c(
  expression(10^0),
  expression(10^1),
  expression(10^2),
  expression(10^3),
  expression(10^4),
  expression(10^5)
), las=1,cex.axis=1.5
)

# plot data
points(haiti.cp23$Age, haiti.cp23$Y,col=alpha(cols[2],alpha=0.3), pch=16,cex=0.6)
points(jitter(usa.cp23$Age,3), usa.cp23$Y,col=alpha(cols[1],alpha=0.3), pch=16,cex=0.6)

lines(haiti.cp23$Age,haiti.cp23$pY,col=cols[2],lwd=2)
lines(usa.cp23$Age,usa.cp23$pY,col=cols[1],lwd=2)
text(5.5,usa.cp23$pY[nrow(usa.cp23)],"USA",col=cols[1],adj=0,pos=4)
text(5.5,haiti.cp23$pY[nrow(haiti.cp23)],"Haiti",col=cols[2],adj=0,pos=4)

dev.off()





#-------------------------------------------
# EXTRA CODE BELOW
#-------------------------------------------
#-------------------------------------------
# density distribution plot
#-------------------------------------------
# pdf("~/dropbox/articles/antibody-curves/presentations/cor-ntd-2016/figs/haiti2-norogii-density.pdf",width=6,height=5)
# du <- density(log10(d.usa$norogii))
# dh <- density(log10(d.hai$norogii))
# xtics <- 0:5
# plot(dh,type="n",
#      main="",
#      ylim=c(0,0.6),
#      xlim=range(xtics),xaxt="n",xlab="",
#      las=1,bty="n"
# )
# polygon(du,col=alpha(cols[1],alpha=1),border=cols[1])
# polygon(dh,col=alpha(cols[2],alpha=0.3),border=cols[2])
# 
# axis(1,at=0:5,labels=c(
#   expression(10^0),
#   expression(10^1),
#   expression(10^2),
#   expression(10^3),
#   expression(10^4),
#   expression(10^5)
# ), las=1,cex.axis=1.5
# )
# mtext("Norovirus GII.4 NO",side=3,line=1,cex=1.25)
# mtext("Luminex Response (MFI-Background)",side=1,line=3,cex=1.25)
# text(1.5,0.45,"USA",col=cols[1])
# text(4.5,0.5,"Haiti",col=cols[2])
# dev.off()



