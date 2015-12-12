#-------------------------------
# 2-garki-main-analysis-figure.R
#
# summarize the results in a figure
#
# version 2 (11 oct 2015)
# drew on common source code for
# estimation
# migrated everything to git
#
# version 1 (26 sep 2015)
#-------------------------------





#-------------------------------
# preamble
#-------------------------------

rm(list=ls())
library(RColorBrewer)
library(scales)

# load the analysis results
load("~/SLAbcurves/results/raw/garki-main-analysis.RData")



#-------------------------------
# Plot antibody response 
# curves for the 3 study phases
#-------------------------------
pdf("~/SLAbcurves/results/figs/garki-antibody-curves-IFATPf.pdf",height=5,width=15)

# general plotting parameters
ytics <- seq(0,4,by=1)
xtics <- seq(0,70,by=10)
# i.cols <- c(brewer.pal(9,"YlGnBu")[7:6],brewer.pal(9,"BuPu")[9:7],brewer.pal(9,"YlOrRd")[6:8])
# c.cols <- brewer.pal(8,"Dark2")[8]
# c.cols <- "black"
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
i.cols <- rep(cbPalette[6],8)
c.cols <- cbPalette[7]

lo <- layout(mat=matrix(1:6,nrow=1,ncol=6),widths=c(1,0.25,1,0.25,1,0.25))

# Pre-intervention age-antibody curves
op <- par(mar=c(4,5,7,1)+0.1,xpd=T)
plot(p.c12$Age,p.c12$pY,type="l",lwd=2,col=c.cols,
	ylim=c(0.5,4),ylab="",yaxt="n",
	xlim=c(0,71),xlab="",xaxt="n",
	bty="n",las=1
	)
	mtext(expression(paste(italic('P. falciparum')," IFAT antibody titre")),side=2,line=3)
	mtext("a",adj=1,line=3,at=-5,font=2,cex=1.75)
	mtext("Pre-Intervention Period",adj=0,line=3,at=0,cex=1.5)
	mtext(expression(paste(italic(E),"(",italic(Y[x][","][a]),")")),side=3,line=0)
	mtext("Age, years",side=1,line=2.5)
	axis(2,at=1:4,labels=c(
		expression(10^1),
		expression(10^2),
		expression(10^3),
		expression(10^4)
		), las=1,cex.axis=1.25
	)
	# segments(x0=min(xtics),x1=max(xtics),y0=ytics,lty=2,col="gray70")
	axis(1,at=xtics,cex.axis=1.5)
	lines(p.tr1$Age,p.tr1$pY,col=i.cols[1])
	lines(p.tr2$Age,p.tr2$pY,col=i.cols[2])
	
	text(70,3.1,"Control",cex=1,adj=1,font=2,col=c.cols)
	text(70,3.6,"Intervention\n(survey 1, 2)",cex=1,adj=1,font=1,col=i.cols[1])
	text(3,3.4,"1",cex=0.75,col=i.cols[1])
	text(10,3.6,"2",cex=0.75,col=i.cols[2])
	
# Pre-intervention  mean estimates
op <- par(mar=c(4,0,7,4)+0.1)
plot(1:3,1:3,type="n",
	xlim=c(0.5,3.5),xaxt="n",xlab="",
	ylim=c(0.5,4),ylab="",yaxt="n",
	las=1,bty="n"
)
# axis(4,at=1:4,labels=c(
	# expression(10^1),
	# expression(10^2),
	# expression(10^3),
	# expression(10^4)
	# ), las=1,cex.axis=1.25
# )
mtext(c("1","2"),side=1,line=1,at=c(1.5,2.5),cex=1,col=i.cols[1:2])
mtext(expression(paste(italic(E),"(",italic(Y[x]),")")),side=3,line=0)
mtext("Survey",side=1,line=2.5)
# control 
	arrows(x0=c(1.5,2.5),y0=unlist(mu.c[3,1:2]), y1=unlist(mu.c[4,1:2]),lwd=1,col=c.cols,length=0.05,angle=90,code=3)
	points(c(1.5,2.5),mu.c[1,1:2], pch=21,cex=1.5, lwd=1,bg=c.cols,col=c.cols)
# intervention
	arrows(x0=c(1.5,2.5),y0=unlist(mu.i[3,1:2]), y1=unlist(mu.i[4,1:2]),lwd=1,col=i.cols[1:2],length=0.05,angle=90,code=3)
	points(c(1.5,2.5),mu.i[1,1:2], pch=21,bg="white",cex=1.5,lwd=1, col=i.cols[1:2])
	

# Active Intervention  age-antibody curves
op <- par(mar=c(4,5,7,1)+0.1)
plot(p.c345$Age,p.c345$pY,type="l",lwd=2,col=c.cols,
	ylim=c(0.5,4),ylab="",yaxt="n",
	xlim=c(0,71),xlab="",xaxt="n",
	bty="n",las=1
	)
	# mtext(expression(paste(italic('P. falciparum')," IFAT antibody titre")),side=2,line=3)
	mtext("b",adj=1,line=3,at=-5,font=2,cex=1.75)
	mtext("Active Intervention Period",adj=0,line=3,at=0,cex=1.5)
	mtext(expression(paste(italic(E),"(",italic(Y[x][","][a]),")")),side=3,line=0)
	mtext("Age, years",side=1,line=2.5)
	axis(2,at=1:4,labels=c(
		expression(10^1),
		expression(10^2),
		expression(10^3),
		expression(10^4)
		), las=1,cex.axis=1.25
	)
	# segments(x0=min(xtics),x1=max(xtics),y0=ytics,lty=2,col="gray70")
	axis(1,at=xtics,cex.axis=1.5)
	lines(p.tr3$Age,p.tr3$pY,col=i.cols[3])
	lines(p.tr4$Age,p.tr4$pY,col=i.cols[4])
	lines(p.tr5$Age,p.tr5$pY,col=i.cols[5])
	
	text(70,3.4,"Control",cex=1,adj=1,font=2,col=c.cols)
	text(70,2.5,"Intervention\n(survey 3, 4, 5)",cex=1,adj=1,font=1,col=i.cols[1])
	text(11,3.0,"3",cex=0.75,col=i.cols[3])
	text(14,2.85,"4",cex=0.75,col=i.cols[4])
	text(19,2.6,"5",cex=0.75,col=i.cols[5],adj=0)
	
# Active intervention  mean estimates
op <- par(mar=c(4,0,7,4)+0.1)
plot(1:3,1:3,type="n",
	xlim=c(0.5,3.5),xaxt="n",xlab="",
	ylim=c(0.5,4),ylab="",yaxt="n",
	las=1,bty="n"
)
# axis(4,at=1:4,labels=c(
	# expression(10^1),
	# expression(10^2),
	# expression(10^3),
	# expression(10^4)
	# ), las=1,cex.axis=1.25
# )
mtext(c(3:5),side=1,line=1,at=1:3,cex=1,col=i.cols[3:5])
mtext(expression(paste(italic(E),"(",italic(Y[x]),")")),side=3,line=0)
mtext("Survey",side=1,line=2.5)
# control 
	arrows(x0=1:3,y0=unlist(mu.c[3,3:5]), y1=unlist(mu.c[4,3:5]),lwd=1,col=c.cols,length=0.05,angle=90,code=3)
	points(1:3,mu.c[1,3:5], pch=21,cex=1.5, lwd=1,bg=c.cols,col=c.cols)
# intervention
	arrows(x0=1:3,y0=unlist(mu.i[3,3:5]), y1=unlist(mu.i[4,3:5]),lwd=1,col=i.cols[3:5],length=0.05,angle=90,code=3)
	points(1:3,mu.i[1,3:5], pch=21,bg="white",cex=1.5,lwd=1, col=i.cols[3:5])
# mark periods with P<0.01 with Bonferroni corrected pvalues
	text(1:3,0.6,"*",cex=2)

# Post-Intervention  age-antibody curves
op <- par(mar=c(4,5,7,1)+0.1)
plot(p.c78$Age,p.c78$pY,type="l",lwd=2,col=c.cols,
	ylim=c(0.5,4),ylab="",yaxt="n",
	xlim=c(0,71),xlab="",xaxt="n",
	bty="n",las=1
	)
	# mtext(expression(paste(italic('P. falciparum')," IFAT antibody titre")),side=2,line=3)
	mtext("c",adj=1,line=3,at=-5,font=2,cex=1.75)
	mtext("Post Intervention Period",adj=0,line=3,at=0,cex=1.5)
	mtext(expression(paste(italic(E),"(",italic(Y[x][","][a]),")")),side=3,line=0)
	mtext("Age, years",side=1,line=2.5)
	axis(2,at=1:4,labels=c(
		expression(10^1),
		expression(10^2),
		expression(10^3),
		expression(10^4)
		), las=1,cex.axis=1.25
	)
	# segments(x0=min(xtics),x1=max(xtics),y0=ytics,lty=2,col="gray70")
	axis(1,at=xtics,cex.axis=1.5)
	lines(p.tr6$Age,p.tr6$pY,col=i.cols[6])
	lines(p.tr7$Age,p.tr7$pY,col=i.cols[7])
	lines(p.tr8$Age,p.tr8$pY,col=i.cols[8])
	
	text(70,3.4,"Control",cex=1,adj=1,font=2,col=c.cols)
	text(70,2.5,"Intervention\n(survey 6, 7, 8)",cex=1,adj=1,font=1,col=i.cols[1])
	text(20,2.7,"6",cex=0.75,col=i.cols[6])
	text(12,3.1,"7",cex=0.75,col=i.cols[7])
	text(15,2.9,"8",cex=0.75,col=i.cols[8],adj=0)

# Post-intervention  mean estimates
op <- par(mar=c(4,0,7,4)+0.1)
plot(1:3,1:3,type="n",
	xlim=c(0.5,3.5),xaxt="n",xlab="",
	ylim=c(0.5,4),ylab="",yaxt="n",
	las=1,bty="n"
)
# axis(4,at=1:4,labels=c(
	# expression(10^1),
	# expression(10^2),
	# expression(10^3),
	# expression(10^4)
	# ), las=1,cex.axis=1.25
# )
mtext(c(6:8),side=1,line=1,at=1:3,cex=1,col=i.cols[6:8])
mtext(expression(paste(italic(E),"(",italic(Y[x]),")")),side=3,line=0)
mtext("Survey",side=1,line=2.5)
# control 
	arrows(x0=2:3,y0=unlist(mu.c[3,6:7]), y1=unlist(mu.c[4,6:7]),lwd=1,col=c.cols,length=0.05,angle=90,code=3)
	points(2:3,mu.c[1,6:7], pch=21,cex=1.5, lwd=1,bg=c.cols,col=c.cols)
# intervention
	arrows(x0=1:3,y0=unlist(mu.i[3,6:8]), y1=unlist(mu.i[4,6:8]),lwd=1,col=i.cols[6:8],length=0.05,angle=90,code=3)
	points(1:3,mu.i[1,6:8], pch=21,bg="white",cex=1.5,lwd=1, col=i.cols[6:8])
# mark periods with P<0.01 with Bonferroni corrected pvalues
	text(2:3,0.6,"*",cex=2)

par(op)
dev.off()


#-------------------------------
# plot means over time
# not used 
#-------------------------------
# # # P-values for differences in means, with Bonferroni correction for 7 tests
# sprintf("%1.4f",unlist(diff.psi[5,])*7)

# # Number of observations by survey round and treatment group
# trNs <- table(d$serosvy[!is.na(d$ifatpftitre)],d$tr[!is.na(d$ifatpftitre)])


# pdf("~/SLAbcurve/results/figs/garki-IFATpf-by-svy.pdf",width=10,height=4.5)
# op <- par(mar=c(7,6,4,5)+0.1,xpd=TRUE)
# # use the same colors for intervention estimates as for the other plots
# # (actually not used)
# i.cols <- c(brewer.pal(9,"YlGnBu")[7:6],brewer.pal(9,"BuPu")[8:6],brewer.pal(9,"YlOrRd")[6:8])
# i.cols <- i.cols[1]
# c.cols <- brewer.pal(8,"Dark2")[8]

# # set up an empty plot
# ytics <- seq(1,4)
# MidPts <- barplot(1:8,names.arg=NA,border=NA,col=NA,
	# ylim=range(ytics),ylab="",yaxt="n",
	# las=1,bty="n"
# )
# axis(2,at=1:4,labels=c(
	# expression(10^1),
	# expression(10^2),
	# expression(10^3),
	# expression(10^4)
	# ), las=1,cex.axis=1.25
# )
# mtext(1:8,side=1,line=1,at=MidPts,cex=1)
# mtext(c("Pre-Intervention Period","Active Intervention Period","Post-Intervention Period"),side=3,line=1,at=c(mean(MidPts[1:2]),mean(MidPts[3:5]),mean(MidPts[6:8])) )
# segments(x0=c(mean(MidPts[2:3]), mean(MidPts[5:6])),y0=min(ytics), y1=max(ytics), lty=2,col="gray80",lwd=2 )
# mtext("Survey Round",side=1,line=2.5,at=MidPts[4])
# mtext(expression(paste(italic('P. falciparum')," IFAT Antibody Titre")),side=2,line=4)
# mtext("Age-adjusted Geometric Mean",side=2,line=3)
	
# # control 
# segments(x0=MidPts[-c(6)],y0=unlist(mu.c[3,]),y1=unlist(mu.c[4,]),lwd=2,col=alpha(c.cols,alpha=0.6))
	# points(MidPts[-c(6)],mu.c[1,],pch=21,cex=1.25,lwd=2,bg=alpha(c.cols,alpha=0.6),col=alpha(c.cols,alpha=0.6))

	# # intervention
	# segments(x0=MidPts,y0=unlist(mu.i[3,]), y1=unlist(mu.i[4,]),lwd=2,col=alpha(i.cols,alpha=0.6))
	# points(MidPts,mu.i[1,],pch=21,bg="white",cex=1.25,lwd=2,col=alpha(i.cols,alpha=0.6))

	# # labels
	# mtext("Control",side=4,line=-1,at=mu.c[1,7],las=1,col=c.cols[1])
	# mtext("Intervention",side=4,line=-1,at=mu.i[1,8],las=1,col=i.cols[1])
	
	# # mark periods with P<0.01 with Bonferroni corrected pvalues
	# text(MidPts[c(3:5,7:8)],1,"*",cex=1.25)
	
	# # make a table of Ns at the bottom
	# mtext(format(trNs[,1],big.mark=","),side=1,line=4,at=MidPts,cex=0.8,col=c.cols)
	# mtext(format(trNs[,2],big.mark=","),side=1,line=5,at=MidPts,cex=0.8,col=i.cols)
	# mtext("N Control",side=1,line=4,at=MidPts[8]+0.5,adj=0,col=c.cols[1])
	# mtext("N Intervention",side=1,line=5,at=MidPts[8]+0.5,adj=0,col=i.cols[1])

# par(op)
# dev.off()

