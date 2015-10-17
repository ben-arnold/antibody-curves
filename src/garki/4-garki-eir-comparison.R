
#-------------------------------
# 4-garki-eir-comparison.R
#
# Compare IFAT curves and 
# age-adjusted mean values
# with EIR estimates from
# the original Garki Project
# (Molineaux 1980, Table 4, p 65)
# for the 3 villages with both
# serological and entomological
# analyses at multiple wet season
# time points
#
# version 2 (16 oct 2015)
# updated script to rely on SLAb-curve.R and SLAb-tmle.R
# base functions
#
# version 1 (24 sep 2015)
#-------------------------------



#-------------------------------
# preamble
#-------------------------------

rm(list=ls())
library(RColorBrewer)
library(scales)
library(SuperLearner)
library(tmle)


# source the base functions for
# SL fits of age-antibody curves
# and TMLE estimates of mean differences
source("~/SLAbcurves/src/SLAb-curve.R")
source("~/SLAbcurves/src/SLAb-tmle.R")



#-------------------------------
# input EIR values from 
# Table 4 (p 65) of Molineaux 1980
# for villages:
# Ajura(553), Rafin-Marke(154), Nasakar(218)
#-------------------------------
eir <- c( 37,25,28,NA,NA,
		  18, 5, 2, 4, 6,
		 129, 0, 4,16,24)
village <- rep(c(553,154,218),c(5,5,5))
vname <- rep(c("Ajura","Rafin Marke","Nasakar"),c(5,5,5))
eirdates <- rep(c("1971-06-21 to 1971-11-07","1972-05-22 to 1972-10-22","1973-06-18 to 1973-11-04","1974-07-29 to 1974-12-15","1975-07-14 to 1975-11-30"),3)
wetseason <- rep(1971:1975,3)
deir <- data.frame(village,vname,wetseason,eirdates,eir)


#-------------------------------
# load the serology dataset
#-------------------------------
d <- read.csv("~/dropbox/garki/data/final/garki-sero.csv")

d$mdate <- as.Date(d$mdate,"%d %b %Y")

#-------------------------------
# subset the dataset to
# the 3 villages with EIR ests
#-------------------------------
ad <- subset(d,village==553|village==154|village==218)
ad$vname <- factor(ad$vname)
ad$vname <- factor(ad$vname,levels=c("Ajura","Nasakar","Rafin Marke"))
ad$vname[ad$village==553] <- "Ajura"

#-------------------------------
# Identify observations in the
# wet seasons
#-------------------------------
ad$wetseason <- NA
ad$wetseason[ad$mdate>="1971-06-21" & ad$mdate<="1971-11-07" ] <- 1971
ad$wetseason[ad$mdate>="1972-05-22" & ad$mdate<="1972-10-22" ] <- 1972
ad$wetseason[ad$mdate>="1973-06-18" & ad$mdate<="1973-11-04" ] <- 1973
ad$wetseason[ad$mdate>="1974-07-29" & ad$mdate<="1974-12-15" ] <- 1974
ad$wetseason[ad$mdate>="1975-07-14" & ad$mdate<="1975-11-30" ] <- 1975

table(ad$vname,ad$wetseason)


#-------------------------------
# IFAT-Pf
# SuperLearner curve fits
#-------------------------------

set.seed(2343242)
ajura.1971.SL <- SLAb.curve(
		Y=log10(ad$ifatpftitre[ad$vname=="Ajura" & ad$wetseason==1971]+1),
		Age=ad$ageyrs[ad$vname=="Ajura" & ad$wetseason==1971],
		id=ad$id[ad$vname=="Ajura" & ad$wetseason==1971]
	)
ajura.1972.SL <- SLAb.curve(
		Y=log10(ad$ifatpftitre[ad$vname=="Ajura" & ad$wetseason==1972]+1),
		Age=ad$ageyrs[ad$vname=="Ajura" & ad$wetseason==1972],
		id=ad$id[ad$vname=="Ajura" & ad$wetseason==1972]
	)
ajura.1973.SL <- SLAb.curve(
		Y=log10(ad$ifatpftitre[ad$vname=="Ajura" & ad$wetseason==1973]+1),
		Age=ad$ageyrs[ad$vname=="Ajura" & ad$wetseason==1973],
		id=ad$id[ad$vname=="Ajura" & ad$wetseason==1973]
	)	

rafin.1971.SL <- SLAb.curve(
		Y=log10(ad$ifatpftitre[ad$vname=="Rafin Marke" & ad$wetseason==1971]+1),
		Age=ad$ageyrs[ad$vname=="Rafin Marke" & ad$wetseason==1971],
		id=ad$id[ad$vname=="Rafin Marke" & ad$wetseason==1971]
	)
rafin.1972.SL <- SLAb.curve(
		Y=log10(ad$ifatpftitre[ad$vname=="Rafin Marke" & ad$wetseason==1972]+1),
		Age=ad$ageyrs[ad$vname=="Rafin Marke" & ad$wetseason==1972],
		id=ad$id[ad$vname=="Rafin Marke" & ad$wetseason==1972]
	)
rafin.1973.SL <- SLAb.curve(
		Y=log10(ad$ifatpftitre[ad$vname=="Rafin Marke" & ad$wetseason==1973]+1),
		Age=ad$ageyrs[ad$vname=="Rafin Marke" & ad$wetseason==1973],
		id=ad$id[ad$vname=="Rafin Marke" & ad$wetseason==1973]
	)
rafin.1974.SL <- SLAb.curve(
		Y=log10(ad$ifatpftitre[ad$vname=="Rafin Marke" & ad$wetseason==1974]+1),
		Age=ad$ageyrs[ad$vname=="Rafin Marke" & ad$wetseason==1974],
		id=ad$id[ad$vname=="Rafin Marke" & ad$wetseason==1974]
	)
rafin.1975.SL <- SLAb.curve(
		Y=log10(ad$ifatpftitre[ad$vname=="Rafin Marke" & ad$wetseason==1975]+1),
		Age=ad$ageyrs[ad$vname=="Rafin Marke" & ad$wetseason==1975],
		id=ad$id[ad$vname=="Rafin Marke" & ad$wetseason==1975]
	)

nasak.1971.SL <- SLAb.curve(
		Y=log10(ad$ifatpftitre[ad$vname=="Nasakar" & ad$wetseason==1971]+1),
		Age=ad$ageyrs[ad$vname=="Nasakar" & ad$wetseason==1971],
		id=ad$id[ad$vname=="Nasakar" & ad$wetseason==1971]
	)
nasak.1972.SL <- SLAb.curve(
		Y=log10(ad$ifatpftitre[ad$vname=="Nasakar" & ad$wetseason==1972]+1),
		Age=ad$ageyrs[ad$vname=="Nasakar" & ad$wetseason==1972],
		id=ad$id[ad$vname=="Nasakar" & ad$wetseason==1972]
	)
nasak.1973.SL <- SLAb.curve(
		Y=log10(ad$ifatpftitre[ad$vname=="Nasakar" & ad$wetseason==1973]+1),
		Age=ad$ageyrs[ad$vname=="Nasakar" & ad$wetseason==1973],
		id=ad$id[ad$vname=="Nasakar" & ad$wetseason==1973]
	)
nasak.1974.SL <- SLAb.curve(
		Y=log10(ad$ifatpftitre[ad$vname=="Nasakar" & ad$wetseason==1974]+1),
		Age=ad$ageyrs[ad$vname=="Nasakar" & ad$wetseason==1974],
		id=ad$id[ad$vname=="Nasakar" & ad$wetseason==1974]
	)
nasak.1975.SL <- SLAb.curve(
		Y=log10(ad$ifatpftitre[ad$vname=="Nasakar" & ad$wetseason==1975]+1),
		Age=ad$ageyrs[ad$vname=="Nasakar" & ad$wetseason==1975],
		id=ad$id[ad$vname=="Nasakar" & ad$wetseason==1975]
	)

#-------------------------------
# TMLE estimates of age-adjusted
# mean antibody titres
#-------------------------------

### LEFT OFF HERE
### NEED TO MODIFY CODE BELOW TO USE SLAb.tmle()
### AND THEN MOVE THE PLOTTING CODE TO A SEPARATE SCRIPT


set.seed(5463452)
ajura.1971.tmle <- tmle.wrap("Ajura",1971,ad)
ajura.1972.tmle <- tmle.wrap("Ajura",1972,ad)
ajura.1973.tmle <- tmle.wrap("Ajura",1973,ad)
ajura.tmle <- list(ajura.1971.tmle,ajura.1972.tmle,ajura.1973.tmle)

rafin.1971.tmle <- tmle.wrap("Rafin Marke",1971,ad)
rafin.1972.tmle <- tmle.wrap("Rafin Marke",1972,ad)
rafin.1973.tmle <- tmle.wrap("Rafin Marke",1973,ad)
rafin.1974.tmle <- tmle.wrap("Rafin Marke",1974,ad)
rafin.1975.tmle <- tmle.wrap("Rafin Marke",1975,ad)
rafin.tmle <- list(rafin.1971.tmle,rafin.1972.tmle,rafin.1973.tmle,rafin.1974.tmle,rafin.1975.tmle)


nasak.1971.tmle <- tmle.wrap("Nasakar",1971,ad)
nasak.1972.tmle <- tmle.wrap("Nasakar",1972,ad)
nasak.1973.tmle <- tmle.wrap("Nasakar",1973,ad)
nasak.1974.tmle <- tmle.wrap("Nasakar",1974,ad)
nasak.1975.tmle <- tmle.wrap("Nasakar",1975,ad)
nasak.tmle <- list(nasak.1971.tmle,nasak.1972.tmle,nasak.1973.tmle,nasak.1974.tmle,nasak.1975.tmle)


#-------------------------------
# summarize the means + 95% CIs
#-------------------------------
ajura.mus <- t(sapply(ajura.tmle,function(x) c(x$mu,x$ci)) )
rafin.mus <- t(sapply(rafin.tmle,function(x) c(x$mu,x$ci)) )
nasak.mus <- t(sapply(nasak.tmle,function(x) c(x$mu,x$ci)) )
rownames(ajura.mus) <- 1971:1973
rownames(rafin.mus) <- rownames(nasak.mus) <- 1971:1975
colnames(ajura.mus) <- colnames(rafin.mus) <- colnames(nasak.mus) <- c("mu","lb","ub")
all.mus <- rbind(ajura.mus,rep(NA,3),rep(NA,3),rafin.mus,nasak.mus)

dmus <- data.frame(
		vname=rep(c("Ajura","Rafin Marke","Nasakar"),c(5,5,5)),
		wetseason=rep(1971:1975,3),
		all.mus
)

#-------------------------------
# merge the means to the EIR data
#-------------------------------
md <- merge(deir,dmus,by=c("vname","wetseason"))
md <- subset(md,!is.na(eir))

# calculate log_10 eir.  add 1 to eir=0
md$log10eir <- log10(md$eir)
md$log10eir[md$eir==0] <- log10(1)

#-------------------------------
# calculate spearman's rank 
# correlation test statistic
#-------------------------------

sp.rho <- cor.test(md$log10eir,md$mu,method="spearman")
sp.rho
rho.text <- substitute(paste("Spearman's ",rho," = ",rho.txt ),list(rho.txt=sprintf("%1.2f",sp.rho$estimate)))

#-------------------------------
# plot the values
#-------------------------------

pdf("~/dropbox/garki/figs/garki-IFATPf-EIR.pdf",width=5,height=5)
op <- par(mar=c(5,4,1,2)+0.1)
cols <- c(brewer.pal(8,"Dark2")[c(8,4)],brewer.pal(8,"Set1")[2])
ytics <- seq(1,4,by=1)
xtics <- seq(0,3,by=1)
plot(1,1,type="n",bty="n",
	xaxt="n",xlab="",xlim=c(0,3),
	yaxt="n",ylab="",ylim=range(ytics),
	las=1
	)
	axis(1,at=xtics,labels=c(1,10,100,1000))
	axis(2,at=ytics,labels=c(
		expression(10^1),
		expression(10^2),
		expression(10^3),
		expression(10^4)),
		las=1
	)
	mtext("Age-adjusted Geometric Mean",side=2,line=3)
	mtext(expression(paste(italic('P. falciparum')," IFAT antibody titre")),side=2,line=2)
	mtext("Entomological Innoculation Rate\n(Cumulative Wet Season Infectious Bites per Person)",side=1,line=3.5)
	text(3,1,rho.text,adj=1,cex=0.8)
	
	# Ajura
	points(md$log10eir[md$vname=="Ajura"],md$mu[md$vname=="Ajura"], pch=16,cex=1.5,col=alpha(cols[1],alpha=0.6))

	# Rafin Marke
	points(md$log10eir[md$vname=="Rafin Marke"],md$mu[md$vname=="Rafin Marke"], pch=16, cex=1.5,col=alpha(cols[2],alpha=0.6))

	# Nasakar
	points(md$log10eir[md$vname=="Nasakar"],md$mu[md$vname=="Nasakar"], pch=16,cex=1.5,col=alpha(cols[3],alpha=0.6))

	# circle pre-intervention measures
	points(md$log10eir[md$wetseason==1971],md$mu[md$wetseason==1971],cex=1.5)
	
	
	# label the villages
	ajura.x <- md$log10eir[md$vname=="Ajura" & md$wetseason==1971]
	ajura.y <- md$mu[md$vname=="Ajura" & md$wetseason==1971]
	segments(x0=ajura.x,y0=ajura.y+0.1,y1=ajura.y+0.3,col="gray40")
	text(ajura.x,ajura.y+0.3,"Ajura (control)",col=cols[1],pos=3,cex=0.7)
	
	rafin.x <- md$log10eir[md$vname=="Rafin Marke" & md$wetseason==1971]
	rafin.y <- md$mu[md$vname=="Rafin Marke" & md$wetseason==1971]
	segments(x0=rafin.x-0.07,x1=rafin.x-0.3,y0=rafin.y+0.07,y1=rafin.y+0.25,col="gray40")
	text(rafin.x-0.3,rafin.y+0.3,"Rafin Marke",col=cols[2],pos=2,cex=0.7)
	
	nasak.x <- md$log10eir[md$vname=="Nasakar" & md$wetseason==1971]
	nasak.y <- md$mu[md$vname=="Nasakar" & md$wetseason==1971]
	segments(x0= nasak.x+0.07,x1= nasak.x+0.3,y0= nasak.y +0.07,y1= nasak.y +0.25,col="gray40") 
	text(nasak.x+0.3, nasak.y +0.3,"Nasakar",col=cols[3],pos=4,cex=0.7)
	
par(op)
dev.off()




