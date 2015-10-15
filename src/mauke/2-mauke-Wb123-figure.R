

#-------------------------------------------
# mauke-LF-auc-super-learner.R
# Ben Arnold
#
# calculate area under the curve for LF
# antibody data. Estimate the curve using
# TMLE
#
# version 2 (16 apr 2015)
# added points to the chart
#
# version 1 (26 mar 2015)
# 
#
#-------------------------------------------

#-------------------------------------------
# input files:
#   mauke1974.dta
#   mauke1992.dta
#
# output files:
#   xxx
#-------------------------------------------



#-------------------------------------------
# preamble
#-------------------------------------------

rm(list=ls())
library(RColorBrewer)
library(foreign)
library(zoo)
library(scales)
library(SuperLearner)

#-------------------------------------------
# load the Mauke data from 1974 and 1992
#-------------------------------------------

d75 <- read.dta("~/dropbox/mauke/data/final/mauke1974.dta")
d92 <- read.dta("~/dropbox/mauke/data/final/mauke1992.dta")

# drop 7 children aged 0 (all Wb123+) due to maternal antibodies
a75 <- subset(d75,age>0)
a75$wb123p <- ifelse(a75$wb123pos=="Positive",1,0)

a92 <- d92
a92$wb123p <- ifelse(a92$wb123pos=="Positive",1,0)

#--------------------------------------
# loess fits of antibody levels
#--------------------------------------
afit75 <- loess(log10(wb123)~age,data=a75)
afit92 <- loess(log10(wb123)~age,data=a92)

#--------------------------------------
# SuperLearner fits of antibody levels
#--------------------------------------
set.seed(0237234)
# SL.glmnet",
SLlib <- c("SL.loess","SL.glm","SL.bayesglm","SL.gam","SL.mean")

SLfit75 <- SuperLearner(Y=log10(a75$wb123),X=data.frame(a75$age),SL.library=SLlib,id=a75$id74)
	SLfit75
SLfit92 <- SuperLearner(Y=log10(a92$wb123),X=data.frame(a92$age),SL.library=SLlib,id=a92$id92)
	SLfit92

#--------------------------------------
# calculate area under the curve
#--------------------------------------

# area under the curve (trapezoidal method)
auc <- function(x,y) {
	# calculate the area under the curve using
	# rectangles and triangles
	# x : x-coordinates of the curve
	# y : y-coordinates of the curve
	#
	require(zoo)
	xs <- diff(x)
	ys <- rollmean(y,2)
	return( sum(xs*ys) )
}

#--------------------------------------
# loess AUC
#--------------------------------------
# get predicted values at 1:70
loess.p75 <- predict(afit75,newdata=data.frame(age=1:70))
loess.p92 <- predict(afit92,newdata=data.frame(age=1:70))

# area under the curve for ages 1 - 70
loess.auc75 <- auc(1:70,loess.p75)
loess.auc92 <- auc(1:70,loess.p92)
cbind(loess.auc75,loess.auc92,loess.auc92/loess.auc75)

# area under the curve for ages 1 - 10
loess.auc75.10 <- auc(1:10,loess.p75[1:10])
loess.auc92.10 <- auc(1:10,loess.p92[1:10])
cbind(loess.auc75.10,loess.auc92.10,loess.auc92.10/loess.auc75.10)

#--------------------------------------
# Super Learner AUC
#--------------------------------------
# get predicted values at 1:70
sl.p75 <- predict(SLfit75,newdata=data.frame(a75.age=1:70))$pred
sl.p92 <- predict(SLfit92,newdata=data.frame(a92.age=1:70))$pred

# area under the curve for ages 1 - 70
sl.auc75 <- auc(1:70,sl.p75)
sl.auc92 <- auc(1:70,sl.p92)
cbind(sl.auc75,sl.auc92,sl.auc92/sl.auc75)

# area under the curve for ages 1 - 10
sl.auc75.10 <- auc(1:10,sl.p75[1:10])
sl.auc92.10 <- auc(1:10,sl.p92[1:10])
cbind(sl.auc75.10,sl.auc92.10,sl.auc92.10/sl.auc75.10)


#--------------------------------------
# Plot Loess and SL fits for comparison
#--------------------------------------
lfit75.x <- afit75$x[order(afit75$x)]
lfit75.f <- afit75$fitted[order(afit75$x)]
lfit92.x <- afit92$x[order(afit92$x)]
lfit92.f <- afit92$fitted[order(afit92$x)]
slfit75.x <- a75$age[order(a75$age)]
slfit75.f <- SLfit75$SL.predict[order(a75$age)]
slfit92.x <- a92$age[order(a92$age)]
slfit92.f <- SLfit92$SL.predict[order(a92$age)]

#--------------------------------------
# 1974 vs 1992 Wb123 antibody levels
#--------------------------------------
pdf("~/dropbox/mauke/figs/mauke-Wb123-Ab-1975v1992-loess-vs-SL-fits.pdf",width=7,height=7)
cols1 <- brewer.pal(8,"Set1")[c(3)]
cols2 <- brewer.pal(11,"Spectral")[c(11)]
brewcols <- c(cols2,cols1)

op <- par(mar=c(4,3,3,1)+0.1)
ylim <- c(2,6)
plot(lfit75.x,lfit75.f,type="n",
	xlab="",xaxt="n",xlim=c(0,70),
	ylab="",yaxt="n",ylim=ylim,
	main="",
	las=1,bty="n"
	)
	axis(1,at=seq(0,70,by=10),cex.axis=1.5)
	axis(2,at=2:6,labels=c(
		expression(10^2),
		expression(10^3),
		expression(10^4),
		expression(10^5),
		expression(10^6)
		), las=1,cex.axis=1.5
	)
	# abline(log(10986),0,lty=2,col="gray60")
	
	# points(jitter(a75$age),log10(a75$wb123),cex=0.45,pch=16,col=brewcols[1])
	# points(jitter(a92$age),log10(a92$wb123),cex=0.45,pch=16,col=brewcols[2])
	
	lines(lfit75.x,lfit75.f,col=brewcols[1],lwd=3)
	lines(lfit92.x,lfit92.f,col=brewcols[2],lwd=3)
	
	lines(slfit75.x,slfit75.f,col=brewcols[1],lwd=1.5,lty=2)
	lines(slfit92.x,slfit92.f,col=brewcols[2],lwd=1.5,lty=2)
	text(30,4,"Heavy solid lines show loess fit\nDashed lines show SuperLearner fit",adj=0)
	
	mtext(substitute(paste("Wb123 (Light Units), ",italic('W. bancrofti'))),side=3,line=0,at=-10,adj=0,cex=1.5)
	# mtext("Wb123 (Light Units)",side=3,line=0,at=-10,adj=0,cex=1.5)
	mtext("Age, years",side=1,line=2.5,cex=1.5)
	
	op <- par(xpd=NA)
	legend(40,6.5,legend=c("1975 (pre MDA)","1992 (post MDA)"),lty=c(1,1),lwd=c(2,2),col=brewcols[1:2],bty="n",cex=1.5,bg="white")
	
par(op)
dev.off()

pdf("~/dropbox/mauke/figs/mauke-Wb123-Ab-1975v1992-SL-fits.pdf",width=7,height=7)
cols1 <- brewer.pal(8,"Set1")[c(3)]
cols2 <- brewer.pal(11,"Spectral")[c(11)]
brewcols <- c(cols2,cols1)

op <- par(mar=c(4,3,3,1)+0.1)
ylim <- c(2,6)
plot(lfit75.x,lfit75.f,type="n",
	xlab="",xaxt="n",xlim=c(0,70),
	ylab="",yaxt="n",ylim=ylim,
	main="",
	las=1,bty="n"
	)
	axis(1,at=seq(0,70,by=10),cex.axis=1.5)
	axis(2,at=2:6,labels=c(
		expression(10^2),
		expression(10^3),
		expression(10^4),
		expression(10^5),
		expression(10^6)
		), las=1,cex.axis=1.5
	)
	# abline(log(10986),0,lty=2,col="gray60")
	
	points(jitter(a75$age),log10(a75$wb123),cex=0.5,pch=16,col=alpha(brewcols[1],0.5))
	points(jitter(a92$age),log10(a92$wb123),cex=0.5,pch=16,col=alpha(brewcols[2],0.5))
	
	# lines(lfit75.x,lfit75.f,col=brewcols[1],lwd=3)
	# lines(lfit92.x,lfit92.f,col=brewcols[2],lwd=3)
	
	lines(slfit75.x,slfit75.f,col=brewcols[1],lwd=3)
	lines(slfit92.x,slfit92.f,col=brewcols[2],lwd=3)
	
	mtext(substitute(paste("Wb123 (Light Units), ",italic('W. bancrofti'))),side=3,line=0,at=-10,adj=0,cex=1.5)
	# mtext("Wb123 (Light Units)",side=3,line=0,at=-10,adj=0,cex=1.5)
	mtext("Age, years",side=1,line=2.5,cex=1.5)
	
	op <- par(xpd=NA)
	legend(40,6.5,legend=c("1975 (pre MDA)","1992 (post MDA)"),lty=c(1,1),lwd=c(2,2),col=brewcols[1:2],bty="n",cex=1.5,bg="white")
	
par(op)
dev.off()




# Not used yet

# # #--------------------------------------
# # # bootstrap the ratio for inference
# # #--------------------------------------
# # set.seed(5872935)

# # iter <- 1000
# # aucs <- rep(NA,iter)
# # aucs.10 <- rep(NA,iter)
# # for (bb in 1:iter) {
	# # b75 <- a75[sample(1:nrow(a75),nrow(a75),replace=TRUE),]
	# # b92 <- a92[sample(1:nrow(a92),nrow(a92),replace=TRUE),]
	# # bfit75 <- loess(log10(wb123)~age,data=b75)
	# # bfit92 <- loess(log10(wb123)~age,data=b92)
	# # bp75 <- predict(bfit75,newdata=data.frame(age=1:70))
	# # bp92 <- predict(bfit92,newdata=data.frame(age=1:70))
	# # aucs[bb] <- auc(1:70,bp92)/auc(1:70,bp75)
	# # aucs.10[bb] <- auc(1:10,bp92[1:10])/auc(1:10,bp75[1:10])
# # }

# # # calculate percentile CIs
# # auc.ci <- quantile(aucs,probs=c(0.025,0.975),na.rm=T)
# # auc.10.ci <- quantile(aucs.10,probs=c(0.025,0.975),na.rm=T)

# # # print results
# # c(auc75,auc92,auc92/auc75,auc.ci)
# # c(auc75.10,auc92.10,auc92.10/auc75.10,auc.10.ci)











