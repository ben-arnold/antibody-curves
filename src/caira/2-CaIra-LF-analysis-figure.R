
#--------------------------------------
# caira-LF-analysis.RData
#
# plot SL and TMLE estimates
#
# version 1 (12 nov 2015)
#--------------------------------------

#--------------------------------------
# preamble
#--------------------------------------
rm(list=ls())
library(RColorBrewer)
library(scales)


#--------------------------------------
# load the Ca Ira LF analysis output
#--------------------------------------
load("~/SLAbcurves/results/raw/caira-LF-analysis.RData")


#-------------------------------
# Plot schema that will be
# repeated for each antibody in
# the composite figure
#-------------------------------
SLAb.plotEYax <- function(SLfit0,SLfit1,main,letter,xlabel=FALSE,ylabel=FALSE) {
	# plot the age-adjusted antibody curve 
	
	# SLfit0 : SL fit object for group 0 (pre MDA) returned from SLAb.curve()
	# SLfit1 : SL fit object for group 1 (post MDA) returned from SLAb.curve()
	# main : plot title
	# letter: letter for multi-panel plots (e.g., "A")
	# xlabel: logical. print a label for the X-axis?
	# ylabel: logical. print a label for the Y-axis
	
	# plotting parameters and empty plot
	op <- par(mar=c(5,6,4,0)+0.1)
	xtics <- seq(0,14,by=1)
	cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
	cols <- cbPalette[c(7,6)]
	
	
	plot(SLfit0$Age, SLfit0$Y,type="n",
		ylab="",yaxt="n",
		ylim=c(0,5),
		xlab="",xlim=range(xtics),xaxt="n",
		bty="n",las=1
		)
	
	# header
	mtext(letter,side=3,line=1.25,font=2,at=-1.7,cex=1.75)
	mtext(main,cex=1.25,line=1.5)
	mtext(expression(paste(italic(E),"(",italic(Y[x][","][a]),")")),side=3,line=-0.5)
	
	# axes
	if (xlabel==TRUE) mtext("Age, years",side=1,line=3,cex=1.5)
	if (ylabel==TRUE) mtext("Luminex Response (MFI-Background)",side=2,line=3.5,cex=1.25)
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
  # jitter age data slightly
	points(jitter(SLfit0$Age,1.5), SLfit0$Y,col=alpha(cols[1],alpha=0.6), pch=16,cex=0.7)
	points(jitter(SLfit1$Age,1.5), SLfit1$Y,col=alpha(cols[2],alpha=0.6), pch=16,cex=0.7)
	lines(SLfit0$Age,SLfit0$pY,col=cols[1],lwd=2)
	lines(SLfit1$Age,SLfit1$pY,col=cols[2],lwd=2)
	
	# survey labels
# 	mtext("USA",side=4,line=1,adj=1,at=SLfit0$pY[length(SLfit0$pY)], col=cols[1],cex=1.25,las=1)
# 	mtext("Haiti",side=4,line=1,adj=1,at=SLfit1$pY[length(SLfit1$pY)], col=cols[2],cex=1.25,las=2)
# 	par(op)
	
}

# skinny E(Yx) plot to add next to the age-specific antibody curves
SLAb.plotEYx <- function(EY0,EY1,Ediff) {
	# plot the adjusted means
	
	# EY0    : Adjusted mean Y for group 0 (pre MDA), returned from SLAb.tmle()
	# EY1    : Adjusted mean Y for group 1 (post MDA), returned from SLAb.tmle()
	# Ediff  : Adjusted mean difference for pre-post MDA, returned from SLAb.tmle()
	
	# plot parameters and empty plot
	op <- par(mar=c(5,0,4,0)+0.1)
	cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
	cols <- cbPalette[c(7,6)]
	plot(1,1,type="n",
		xlim=c(0,1),xaxt="n",xlab="",
		ylim=c(0,5),ylab="",yaxt="n",
		las=1,bty="n"
	)
	
	# Y-axis
	# axis(4,at=0:5,labels=c(
		# # expression(10^-1),
		# expression(10^0),
		# expression(10^1),
		# expression(10^2),
		# expression(10^3),
		# expression(10^4),
		# expression(10^5)
		# ), las=1,cex.axis=1.5
	# )
	mtext(expression(paste(italic(E),"(",italic(Y[x]),")")),side=1,line=1)
	mtext(expression(paste(italic(E),"(",italic(Y[x]),")")),side=3,line=-0.5)
	
	# plot data
	arrows(x0=0.5,y0=EY0$lb, y1=EY0$ub,lwd=1,col=cols[1],length=0.05,angle=90,code=3)
	points(0.5,EY0$psi, pch=16,cex=1.75, lwd=1,bg="white",col=cols[1])
	
	arrows(x0=0.5,y0=EY1$lb, y1=EY1$ub,lwd=1,col=cols[2],length=0.05,angle=90,code=3)
	points(0.5,EY1$psi, pch=21,cex=1.75, lwd=1,bg="white",col=cols[2])
	
  
  # label points
  mtext("pre-\nMDA",side=2,line=-2,adj=1,las=1,col=cols[1],at=EY0$psi,cex=0.8)
	mtext("post-\nMDA",side=2,line=-2,adj=1,las=1,col=cols[2],at=EY1$psi,cex=0.8)
  
	# if(Ediff$p<0.01) text(0.5,0.1,"*",cex=2)
	
}


set.seed(27313)
pdf("~/SLAbcurves/results/figs/caira-LF-curves.pdf",height=5,width=15)
lo <- layout(mat=matrix(1:6,nrow=1,ncol=6,byrow=TRUE),widths=c(1,0.2,1,0.2,1,0.2))

SLAb.plotEYax(wb123.fit1,wb123.fit3,expression(paste(italic('W. bancrofti'), " Wb123")),"a",xlabel=T,ylabel=T)
SLAb.plotEYx(mu.wb123[,1],mu.wb123[,2],diff.wb123) 

SLAb.plotEYax(bm14.fit1,bm14.fit3,expression(paste(italic('W. bancrofti'), " Bm14")),"b",xlabel=T)
SLAb.plotEYx(mu.bm14[,1],mu.bm14[,2],diff.bm14) 

SLAb.plotEYax(bm33.fit1,bm33.fit3,expression(paste(italic('W. bancrofti'), " Bm33")),"c",xlabel=T)
SLAb.plotEYx(mu.bm33[,1],mu.bm33[,2],diff.bm33) 

dev.off()

	