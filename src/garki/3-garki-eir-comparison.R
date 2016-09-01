
#-------------------------------
# 3-garki-eir-comparison.R
#
# Compare IFA titre 
# age-adjusted mean values
# with EIR estimates from
# the original Garki Project
# (Molineaux 1980, Table 4, p 65)
# for the 3 villages with both
# serological and entomological
# analyses at multiple wet season
# time points
#
#-------------------------------



#-------------------------------
# preamble
#-------------------------------

rm(list=ls())
library(RColorBrewer)
library(scales)
library(SuperLearner)
library(tmle)
library(tmleAb)


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
d <- read.csv("~/dropbox/articles/antibody-curves/data/garki/final/garki-sero.csv")

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

# SL library
SL.library <- c("SL.mean","SL.glm","SL.Yman2016","SL.gam","SL.loess")


set.seed(2343242)
ajura.1971.SL <- ab_agecurve(
		Y=log10(ad$ifatpftitre[ad$vname=="Ajura" & ad$wetseason==1971]+1),
		Age=ad$ageyrs[ad$vname=="Ajura" & ad$wetseason==1971],
		id=ad$id[ad$vname=="Ajura" & ad$wetseason==1971],
		SL.library=SL.library
	)
ajura.1972.SL <- ab_agecurve(
		Y=log10(ad$ifatpftitre[ad$vname=="Ajura" & ad$wetseason==1972]+1),
		Age=ad$ageyrs[ad$vname=="Ajura" & ad$wetseason==1972],
		id=ad$id[ad$vname=="Ajura" & ad$wetseason==1972],
		SL.library=SL.library
	)
ajura.1973.SL <- ab_agecurve(
		Y=log10(ad$ifatpftitre[ad$vname=="Ajura" & ad$wetseason==1973]+1),
		Age=ad$ageyrs[ad$vname=="Ajura" & ad$wetseason==1973],
		id=ad$id[ad$vname=="Ajura" & ad$wetseason==1973],
		SL.library=SL.library
	)	

rafin.1971.SL <- ab_agecurve(
		Y=log10(ad$ifatpftitre[ad$vname=="Rafin Marke" & ad$wetseason==1971]+1),
		Age=ad$ageyrs[ad$vname=="Rafin Marke" & ad$wetseason==1971],
		id=ad$id[ad$vname=="Rafin Marke" & ad$wetseason==1971],
		SL.library=SL.library
	)
rafin.1972.SL <- ab_agecurve(
		Y=log10(ad$ifatpftitre[ad$vname=="Rafin Marke" & ad$wetseason==1972]+1),
		Age=ad$ageyrs[ad$vname=="Rafin Marke" & ad$wetseason==1972],
		id=ad$id[ad$vname=="Rafin Marke" & ad$wetseason==1972],
		SL.library=SL.library
	)
rafin.1973.SL <- ab_agecurve(
		Y=log10(ad$ifatpftitre[ad$vname=="Rafin Marke" & ad$wetseason==1973]+1),
		Age=ad$ageyrs[ad$vname=="Rafin Marke" & ad$wetseason==1973],
		id=ad$id[ad$vname=="Rafin Marke" & ad$wetseason==1973],
		SL.library=SL.library
	)
rafin.1974.SL <- ab_agecurve(
		Y=log10(ad$ifatpftitre[ad$vname=="Rafin Marke" & ad$wetseason==1974]+1),
		Age=ad$ageyrs[ad$vname=="Rafin Marke" & ad$wetseason==1974],
		id=ad$id[ad$vname=="Rafin Marke" & ad$wetseason==1974],
		SL.library=SL.library
	)
rafin.1975.SL <- ab_agecurve(
		Y=log10(ad$ifatpftitre[ad$vname=="Rafin Marke" & ad$wetseason==1975]+1),
		Age=ad$ageyrs[ad$vname=="Rafin Marke" & ad$wetseason==1975],
		id=ad$id[ad$vname=="Rafin Marke" & ad$wetseason==1975],
		SL.library=SL.library
	)

nasak.1971.SL <- ab_agecurve(
		Y=log10(ad$ifatpftitre[ad$vname=="Nasakar" & ad$wetseason==1971]+1),
		Age=ad$ageyrs[ad$vname=="Nasakar" & ad$wetseason==1971],
		id=ad$id[ad$vname=="Nasakar" & ad$wetseason==1971],
		SL.library=SL.library
	)
nasak.1972.SL <- ab_agecurve(
		Y=log10(ad$ifatpftitre[ad$vname=="Nasakar" & ad$wetseason==1972]+1),
		Age=ad$ageyrs[ad$vname=="Nasakar" & ad$wetseason==1972],
		id=ad$id[ad$vname=="Nasakar" & ad$wetseason==1972],
		SL.library=SL.library
	)
nasak.1973.SL <- ab_agecurve(
		Y=log10(ad$ifatpftitre[ad$vname=="Nasakar" & ad$wetseason==1973]+1),
		Age=ad$ageyrs[ad$vname=="Nasakar" & ad$wetseason==1973],
		id=ad$id[ad$vname=="Nasakar" & ad$wetseason==1973],
		SL.library=SL.library
	)
nasak.1974.SL <- ab_agecurve(
		Y=log10(ad$ifatpftitre[ad$vname=="Nasakar" & ad$wetseason==1974]+1),
		Age=ad$ageyrs[ad$vname=="Nasakar" & ad$wetseason==1974],
		id=ad$id[ad$vname=="Nasakar" & ad$wetseason==1974],
		SL.library=SL.library
	)
nasak.1975.SL <- ab_agecurve(
		Y=log10(ad$ifatpftitre[ad$vname=="Nasakar" & ad$wetseason==1975]+1),
		Age=ad$ageyrs[ad$vname=="Nasakar" & ad$wetseason==1975],
		id=ad$id[ad$vname=="Nasakar" & ad$wetseason==1975],
		SL.library=SL.library
	)

#-------------------------------
# TMLE estimates of age-adjusted
# mean antibody titres in each
# year (wet season only)
#-------------------------------

set.seed(5463452)
ajura.tmle <- sapply(c(1971,1972,1973),function(x) ab_tmle(
  Y=log10(ad$ifatpftitre[ad$vname=="Ajura" & ad$wetseason==x]+1),
  Age=ad$ageyrs[ad$vname=="Ajura" & ad$wetseason==x],
  id=ad$id[ad$vname=="Ajura" & ad$wetseason==x],
  SL.library=SL.library)
  )

rafin.tmle <- sapply(c(1971,1972,1973,1974,1975),function(x) ab_tmle(
  Y=log10(ad$ifatpftitre[ad$vname=="Rafin Marke" & ad$wetseason==x]+1),
  Age=ad$ageyrs[ad$vname=="Rafin Marke" & ad$wetseason==x],
  id=ad$id[ad$vname=="Rafin Marke" & ad$wetseason==x],
  SL.library=SL.library)
)

nasak.tmle <- sapply(c(1971,1972,1973,1974,1975),function(x) ab_tmle(
  Y=log10(ad$ifatpftitre[ad$vname=="Nasakar" & ad$wetseason==x]+1),
  Age=ad$ageyrs[ad$vname=="Nasakar" & ad$wetseason==x],
  id=ad$id[ad$vname=="Nasakar" & ad$wetseason==x],
  SL.library=SL.library)
)

#-------------------------------
# summarize the means + 95% CIs
#-------------------------------
ajura.mus <- matrix(unlist(ajura.tmle[c(1,3,4),]),nrow=3,ncol=3,byrow=T)
rafin.mus <- matrix(unlist(rafin.tmle[c(1,3,4),]),nrow=5,ncol=3,byrow=T)
nasak.mus <- matrix(unlist(nasak.tmle[c(1,3,4),]),nrow=5,ncol=3,byrow=T)
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

pdf("~/dropbox/articles/antibody-curves/results/figs/garki-IFAPf-EIR.pdf",width=6,height=6)
op <- par(mar=c(5,5,3,2)+0.1)
# cols <- c(brewer.pal(8,"Dark2")[c(8,4)],brewer.pal(8,"Set1")[2])
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cols <- c(cbPalette[c(7,6,1)])

# lo <- layout(mat=matrix(1:3,nrow=1,ncol=3))
# plot(1,type="n",bty="n",xlab="",ylab="",xaxt="n",yaxt="n")
# plot(1,type="n",bty="n",xlab="",ylab="",xaxt="n",yaxt="n")


i.cols <- rep(cbPalette[6],8)
ytics <- seq(1,4,by=1)
xtics <- seq(0,2,by=1)
plot(1,1,type="n",bty="n",
	xaxt="n",xlab="",xlim=c(0,2.2),
	yaxt="n",ylab="",ylim=range(ytics),
	las=1
	)
	axis(1,at=xtics,labels=c(1,10,100),cex.axis=1.25)
	axis(2,at=ytics,labels=c(
		expression(10^1),
		expression(10^2),
		expression(10^3),
		expression(10^4)),
		las=1,cex.axis=1.25
	)
  mtext("d",adj=1,line=0.5,at=-0.4,font=2,cex=2)
	mtext("Wet season village geometric mean",side=2,line=4,cex=1.25)
  mtext(expression(paste(italic('P. falciparum')," IFA antibody titre, ", italic(E(Y[x])) )) ,side=2,line=2.5,cex=1.25)
	mtext("Entomological Inoculation Rate\n(cumulative wet season infectious bites per person)",side=1,line=3.5,cex=1.25)
	text(2,1,rho.text,adj=1,cex=1.25)
	
	# Ajura
	points(md$log10eir[md$vname=="Ajura"],md$mu[md$vname=="Ajura"], pch=16,cex=2,col=cols[1])

	# Rafin Marke
	points(md$log10eir[md$vname=="Rafin Marke"],md$mu[md$vname=="Rafin Marke"], pch=16, cex=2,col=cols[2])

	# Nasakar, exclude 1972 b/c it is out of the scale
	points(md$log10eir[md$vname=="Nasakar" & md$wetseason!=1972],md$mu[md$vname=="Nasakar" & md$wetseason!=1972], pch=16,cex=2,col=cols[3])

	# circle pre-intervention measures
	points(md$log10eir[md$wetseason==1971],md$mu[md$wetseason==1971],cex=2)
	
	# label the villages
	ajura.x <- md$log10eir[md$vname=="Ajura" & md$wetseason==1971]
	ajura.y <- md$mu[md$vname=="Ajura" & md$wetseason==1971]
	segments(x0=ajura.x,y0=ajura.y+0.1,y1=ajura.y+0.2,col="gray40")
	text(ajura.x,ajura.y+0.2,"Ajura (control)",col=cols[1],pos=3,cex=1)
	
	rafin.x <- md$log10eir[md$vname=="Rafin Marke" & md$wetseason==1971]
	rafin.y <- md$mu[md$vname=="Rafin Marke" & md$wetseason==1971]
	segments(x0=rafin.x-0.07,x1=rafin.x-0.3,y0=rafin.y+0.07,y1=rafin.y+0.25,col="gray40")
	text(rafin.x-0.3,rafin.y+0.3,"Rafin Marke",col=cols[2],pos=2,cex=1)
	
	nasak.x <- md$log10eir[md$vname=="Nasakar" & md$wetseason==1971]
	nasak.y <- md$mu[md$vname=="Nasakar" & md$wetseason==1971]
  segments(x0= nasak.x,y0= nasak.y +0.1,y1= nasak.y +0.4,col="gray40") 
  text(nasak.x, nasak.y +0.4,"Nasakar",col=cols[3],pos=3,cex=1)
	
par(op)
dev.off()




#-------------------------------
# save the output
#-------------------------------
rm(d)
save.image("~/dropbox/articles/antibody-curves/results/raw/garki-eir-comparison.RData")



