
#-------------------------------
# 3-garki-village-EY-by-round
#
# Calculate age-adjusted mean
# IFAT antibody titres by
# village and survey round
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


#-------------------------------
# load the serology dataset
#-------------------------------
d <- read.csv("~/dropbox/garki/data/final/garki-sero.csv")

d$mdate <- as.Date(d$mdate,"%d %b %Y")


# add 2 control village names
d$vname <- factor(d$vname,levels=c(levels(d$vname),"Nabanawa","Ajura"))
d$vname[d$village==552] <- "Nabanawa"
d$vname[d$village==553] <- "Ajura"
d$vname <- factor(d$vname)

#-------------------------------
# Serological survey timing:
# 1-2 : pre-intervention
# 3-5 : intervention period
# 6-8 : post-intervention
#-------------------------------

#-------------------------------
# wrapper function to call TMLE 
# for each village and survey round
#-------------------------------
tmle.wrap <- function(serosvy,data) {
	# serosvy : integer, survey round
	# data    : data frame used for estimation (must include vars extracted below)
	
	SLlib <- c("SL.mean","SL.glm","SL.bayesglm","SL.loess","SL.gam","SL.glmnet")
	# note that the W matrix includes a row of 1s to get the SL.glmnet algorithm to run
	# extract objects to make the calculations easier
	# subset to non-missing data for SL fit
	ifatpf <- log10(data$ifatpftitre+1)
	vil <- data$village
	svy <- data$serosvy
	age <- data$ageyrs
	id  <- data$id
	fitd <- data.frame(ifatpf,vil,svy,age,id)
	fitd <- fitd[complete.cases(fitd),]
	
	mu.fit <- tmle(Y=fitd$ifatpf[fitd$svy==serosvy],
			A=NULL,
			W=data.frame(fitd$age[fitd$svy==serosvy], 
				 rep(1,length(fitd$age[fitd$svy==serosvy]))),
			id=fitd$id[fitd$svy==serosvy],
			Q.SL.library=SLlib
	)
	print(mu.fit)
	mu <- mu.fit$estimates$EY1$psi
	se <- sqrt(mu.fit$estimates$EY1$var.psi)
	ci <- mu.fit$estimates$EY1$CI
	list(mu=mu,se=se,lb=ci[1],ub=ci[2])
}

#-------------------------------
# Estimate mean IFAT P. falciparm
# titre by village and survey round
#-------------------------------

set.seed(5463452)

## Control Villages (no measurement in round 6)
# Nabanawa + Ajura
v552 <- sapply(c(1:5,7,8),tmle.wrap,data=d[d$village==552|d$village==553,])
	# add NA column for survey round 6
	v552 <- cbind(v552[,1:5],rep(NA,4),v552[,6:7])



### Spraying (Propoxur) + MDA

# village cluster 5
# Kawari
v153 <- sapply(c(1:8),tmle.wrap,data=d[d$village==153,])
# Rafin Marke
v154 <- sapply(c(1:8),tmle.wrap,data=d[d$village==154,])
# Kukar Maikiva
v155 <- sapply(c(1:8),tmle.wrap,data=d[d$village==155,])

# village cluster 7
# Kargo Kudu
v213 <- sapply(c(1:8),tmle.wrap,data=d[d$village==213,])
# Nasakar
v218 <- sapply(c(1:8),tmle.wrap,data=d[d$village==218,])
# Bakan Sabara
v220 <- sapply(c(1:8),tmle.wrap,data=d[d$village==220,])



#-------------------------------
# plot means over time
#-------------------------------


# plotting schema, repeated for each intervention village
EYplot <- function(x,cols,vname,header=FALSE,footer=FALSE) {
	# x    : intervention village/survey round results, from above
	# cols : colors for control and intervention points
	# vname: village name for printing
	# header: logical. print header text?
	# footer: logical. print footer text?
	
	# set up an empty plot
	ytics <- seq(1,4)
	MidPts <- barplot(1:8,names.arg=NA,border=NA,col=NA,
		ylim=range(ytics),ylab="",yaxt="n",
		las=1,bty="n"
	)
	axis(2,at=1:4,labels=c(
		expression(10^1),
		expression(10^2),
		expression(10^3),
		expression(10^4)
		), las=1,cex.axis=1.25
	)
	mtext(1:8,side=1,line=1,at=MidPts,cex=0.8,col="gray40")
	if(header==TRUE) mtext(c("Pre-Intervention","Intervention Phase","Post-Intervention"),side=3,line=1,at=c(mean(MidPts[1:2]),mean(MidPts[3:5]),mean(MidPts[6:8])) )
	segments(x0=c(mean(MidPts[2:3]), mean(MidPts[5:6])),y0=min(ytics), y1=max(ytics), lty=2,col="gray80",lwd=2 )
	
	mtext(vname,side=2,line=3,adj=1,col=cols[2],las=1)
	
	
	# control 
	segments(x0=MidPts,y0=unlist(v552[3,]),y1=unlist(v552[4,]),lwd=2,col=alpha(cols[1],alpha=0.6))
	points(MidPts,v552[1,],pch=16,cex=1.5,col=alpha(cols[1],alpha=0.6))

	# intervention
	segments(x0=MidPts,y0=unlist(x[3,]),y1=unlist(x[4,]),lwd=2,col=alpha(cols[2],alpha=0.6))
	points(MidPts,x[1,],pch=16,cex=1.5,col=alpha(cols[2],alpha=0.6))
}


pdf("~/dropbox/garki/figs/garki-IFATpf-by-village-svy.pdf",width=7,height=10)
brewcols1 <- brewer.pal(8,"Dark2")
brewcols2 <- brewer.pal(9,"Set1")
cols <- c(brewcols1[c(8,4)],brewcols2[c(1,5,3,2,4)])
lo <- layout(mat=matrix(1:8,nrow=8,ncol=1),heights=c(0.3,rep(1,6),0.3))
# header
op <- par(mar=c(0,10,0,0)+0.1,xpd=TRUE)
MidPts <- barplot(rep(1,8),names.arg=NA,border=NA,col=NA,ylab="",yaxt="n",bty="n")
text(c(mean(MidPts[1:2]),mean(MidPts[3:5]),mean(MidPts[6:8])),rep(0.5,3),c("Pre-Intervention","Intervention Phase","Post-Intervention"),cex=1.5)
mtext(expression(italic('P. falciparum')),side=2,at=0.5,adj=1,las=1,cex=0.9 )
mtext("IFAT antibody titre",side=2,at=0,adj=1,las=1,cex=0.9 )
# figures
op <- par(mar=c(2,10,1,0)+0.1)
EYplot(v154,cols=cols[c(1,2)],vname="Rafin\nMarke")
EYplot(v153,cols=cols[c(1,3)],vname="Kawari")
EYplot(v155,cols=cols[c(1,4)],vname="Kukar\nMaikiva")
EYplot(v213,cols=cols[c(1,5)],vname="Kargo\nKudu")
EYplot(v218,cols=cols[c(1,6)],vname="Nasakar")
EYplot(v220,cols=cols[c(1,7)],vname="Bakan\nSabara")
# footer
op <- par(mar=c(0,10,0,0)+0.1)
plot(1:8,rep(1,8),type="n",bty="n",xlab="",xaxt="n",ylab="",yaxt="n")
text(mean(1:8),1,"Survey Round",cex=1.5)
par(op)
dev.off()



