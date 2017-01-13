
#-------------------------------
# 3b-garki-eir-comparison-seroprev.R
#
# Compare IFA seroprevalence
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


# bright color blind palette:  https://personal.sron.nl/~pault/ 
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
# d <- read.csv("~/dropbox/articles/antibody-curves/data/garki/final/garki-sero.csv")
data("garki_sero")
d <- garki_sero

d$mdate <- as.Date(d$mdate,"%d %b %Y")

# identify sero-positive individuals (Table 21 of Molineaux 1980, p 175)
d$ifapfpos <- ifelse(d$ifatpftitre>=20,1,0)

# subset to ages 0-20
# d <- subset(d,ageyrs<=5)

# drop observations with missing outcome or age information
d <- subset(d,!is.na(d$ifapfpos) & !is.na(d$ageyrs))

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
# TMLE estimates of 
# mean seroprevalence in each
# year (wet season only)
#-------------------------------

# SL library
SL.library <- c("SL.mean","SL.glm","SL.gam","SL.loess")


set.seed(5463452)
# 1971 and 1972 are 100% positive, so TMLE estimation impossible
# c(1971,1972,1973)
ajura.tmle <- sapply(c(1973),function(x) tmleAb(
  Y=ad$ifapfpos[ad$vname=="Ajura" & ad$wetseason==x],
  W=data.frame(Age=ad$ageyrs[ad$vname=="Ajura" & ad$wetseason==x]),
  id=ad$id[ad$vname=="Ajura" & ad$wetseason==x],
  SL.library=SL.library)[c("psi","se","lb","ub","p")]
  )

rafin.tmle <- sapply(c(1971,1972,1973,1974,1975),function(x) tmleAb(
  Y=ad$ifapfpos[ad$vname=="Rafin Marke" & ad$wetseason==x],
  W=data.frame(Age=ad$ageyrs[ad$vname=="Rafin Marke" & ad$wetseason==x]),
  id=ad$id[ad$vname=="Rafin Marke" & ad$wetseason==x],
  SL.library=SL.library)[c("psi","se","lb","ub","p")]
)

# 1971 is 100% positive, so TMLE estimation impossible
# c(1971,1972,1973,1974,1975)
nasak.tmle <- sapply(c(1972,1973,1974,1975),function(x) tmleAb(
  Y=ad$ifapfpos[ad$vname=="Nasakar" & ad$wetseason==x],
  W=data.frame(Age=ad$ageyrs[ad$vname=="Nasakar" & ad$wetseason==x]),
  id=ad$id[ad$vname=="Nasakar" & ad$wetseason==x],
  SL.library=SL.library)[c("psi","se","lb","ub","p")]
)

#-------------------------------
# add back in the 100% means
#-------------------------------
ajura.tmle <- cbind(c(1,0,1,1,0),c(1,0,1,1,0),ajura.tmle)
nasak.tmle <- cbind(c(1,0,1,1,0),nasak.tmle)
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

pdf("~/dropbox/articles/antibody-curves/results/figs/garki-IFAPf-EIR-seroprev.pdf",width=6,height=6)
op <- par(mar=c(5,6,3,2)+0.1,xpd=TRUE)
# cols <- c(brewer.pal(8,"Dark2")[c(8,4)],brewer.pal(8,"Set1")[2])
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cols <- c(cbPalette[c(7,6,1)])

# lo <- layout(mat=matrix(1:3,nrow=1,ncol=3))
# plot(1,type="n",bty="n",xlab="",ylab="",xaxt="n",yaxt="n")
# plot(1,type="n",bty="n",xlab="",ylab="",xaxt="n",yaxt="n")


i.cols <- rep(cbPalette[6],8)
ytics <- seq(0.7,1,by=0.05)
xtics <- seq(0,2,by=1)
plot(1,1,type="n",bty="n",
	xaxt="n",xlab="",xlim=c(0,2.2),
	yaxt="n",ylab="",ylim=range(ytics),
	las=1
	)
	axis(1,at=xtics,labels=c(1,10,100),cex.axis=1.25)
	axis(2,at=ytics,labels=sprintf("%1.0f",ytics*100),las=1)
  mtext("d",adj=1,line=0.5,at=-0.4,font=2,cex=2)
	mtext(expression(paste("Wet season village ",italic('P. falciparum'))),side=2,line=4,cex=1.25)
  mtext(expression(paste("IFA seroprevalence, ", italic(E(Y[x])) )) ,side=2,line=2.5,cex=1.25)
	mtext("Entomological Inoculation Rate\n(cumulative wet season infectious bites per person)",side=1,line=3.5,cex=1.25)
	text(2,0.7,rho.text,adj=1,cex=1.25)
	
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
	segments(x0=ajura.x,y0=ajura.y+0.01,y1=ajura.y+0.02,col="gray40")
	text(ajura.x,ajura.y+0.02,"Ajura (control)",col=cols[1],pos=3,cex=1)
	
	rafin.x <- md$log10eir[md$vname=="Rafin Marke" & md$wetseason==1971]
	rafin.y <- md$mu[md$vname=="Rafin Marke" & md$wetseason==1971]
	segments(x0=rafin.x-0.1,x1=rafin.x-0.2,y0=rafin.y,col="gray40")
	text(rafin.x-0.2,rafin.y,"Rafin Marke",col=cols[2],pos=2,cex=1)
	
	nasak.x <- md$log10eir[md$vname=="Nasakar" & md$wetseason==1971]
	nasak.y <- md$mu[md$vname=="Nasakar" & md$wetseason==1971]
  segments(x0= nasak.x,y0= nasak.y-0.01,y1= nasak.y-0.02,col="gray40") 
  text(nasak.x, nasak.y-0.02,"Nasakar",col=cols[3],pos=1,cex=1)
	
par(op)
dev.off()




#-------------------------------
# save the output
#-------------------------------
rm(d)
save.image("~/dropbox/articles/antibody-curves/results/raw/garki-eir-comparison-seroprev.RData")



