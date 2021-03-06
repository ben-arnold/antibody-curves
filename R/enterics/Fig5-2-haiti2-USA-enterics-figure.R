
#--------------------------------------
# Fig5-2-haiti2-USA-enterics-figure.R
#
# ben arnold (benarnold@berkeley.edu)
#
# plot SL and TMLE estimates for the
# enteric antibody analyses
#
#--------------------------------------

#--------------------------------------
# preamble
#--------------------------------------
rm(list=ls())
library(RColorBrewer)
library(scales)


#--------------------------------------
# input files:
#   haiti2-usa-enterics-analysis.RData
#
# output files:
#   haiti2-USA-enterics-Ab-curves.pdf
#--------------------------------------


#--------------------------------------
# load the enterics analysis output
#--------------------------------------
load("~/dropbox/articles/antibody-curves/results/raw/haiti2-usa-enterics-analysis.RData")


#-------------------------------
# Plot schema that will be
# repeated for each antibody in
# the composite figure
#-------------------------------
SLAb.plotEYax <- function(SLfit0,SLfit1,main,letter,xlabel=FALSE,ylabel=FALSE) {
	# plot the age-adjusted antibody curve 
	
	# SLfit0 : SL fit object for group 0 (USA) returned from SLAb.curve()
	# SLfit1 : SL fit object for group 1 (Haiti) returned from SLAb.curve()
	# main : plot title
	# letter: letter for multi-panel plots (e.g., "A")
	# xlabel: logical. print a label for the X-axis?
	# ylabel: logical. print a label for the Y-axis
	
	# plotting parameters and empty plot
	op <- par(mar=c(5,6,5,0)+0.1)
	xtics <- seq(0,6,by=1)
	# cols1 <- brewer.pal(8,"Set1")[c(3)]
	# cols2 <- brewer.pal(11,"Spectral")[c(11)]
	# cols <- c(cols1,cols2)
	cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
	cols <- cbPalette[c(6,7)]
	
	
	plot(SLfit0$Age, SLfit0$Y,type="n",
		ylab="",yaxt="n",
		ylim=c(0,5),
		xlab="",xlim=range(xtics),xaxt="n",
		bty="n",las=1
		)
	
	# header
	mtext(letter,side=3,line=2,font=2,at=-0.7,cex=1.75)
	mtext(main,cex=1.5,line=2,at=0,adj=0)
	mtext(expression(paste("Age-dependent mean, ",italic(E),"(",italic(Y[a][","][x]),")")),side=3,line=-0.5,at=0,adj=0,cex=1.25)
	
	# axes
	if (xlabel==TRUE) mtext("Age, years",side=1,line=3.75,cex=1.5)
	if (ylabel==TRUE) mtext("Luminex response (MFI-background)",side=2,line=4,cex=1.25)
	axis(1,at=xtics,cex.axis=2)
	axis(2,at=0:5,labels=c(
		# expression(10^-1),
		expression(10^0),
		expression(10^1),
		expression(10^2),
		expression(10^3),
		expression(10^4),
		expression(10^5)
		), las=1,cex.axis=2
	)
	
	# plot data (points and fitted lines)
	points(SLfit0$Age, SLfit0$Y,col=alpha(cols[1],alpha=0.6), pch=16,cex=0.7)
	points(SLfit1$Age, SLfit1$Y,col=alpha(cols[2],alpha=0.6), pch=16,cex=0.7)
	lines(SLfit0$Age,SLfit0$pY,col=cols[1],lwd=2)
	lines(SLfit1$Age,SLfit1$pY,col=cols[2],lwd=2)
	
	# country labels
	mtext("USA",side=4,line=1,adj=1,at=SLfit0$pY[length(SLfit0$pY)], col=cols[1],cex=1.25,las=1)
	mtext("Haiti",side=4,line=1,adj=1,at=SLfit1$pY[length(SLfit1$pY)], col=cols[2],cex=1.25,las=2)
	par(op)
	
}

# skinny E(Yx) plot to add next to the age-specific antibody curves
SLAb.plotEYx <- function(EY0,EY1,Ediff) {
	# plot the adjusted means
	
	# EY0    : Adjusted mean Y for group 0 (USA), returned from SLAb.tmle()
	# EY1    : Adjusted mean Y for group 1 (Haiti), returned from SLAb.tmle()
	# Ediff  : Adjusted mean difference for Haiti-USA, returned from SLAb.tmle()
	
	# plot parameters and empty plot
	op <- par(mar=c(5,0,4,0.5)+0.1)
	# cols1 <- brewer.pal(8,"Set1")[c(3)]
	# cols2 <- brewer.pal(11,"Spectral")[c(11)]
	# cols <- c(cols1,cols2)
	cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
	cols <- cbPalette[c(6,7)]
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
	mtext(expression(paste(italic(E),"(",italic(Y[x]),")")),side=1,line=1,cex=1.25)
	# mtext(expression(paste("Mean, ",italic(E),"(",italic(Y[x]),")")),side=3,line=-0.5,cex=1.25)
	
	# plot data
	arrows(x0=0.5,y0=EY0$lb, y1=EY0$ub,lwd=1,col=cols[1],length=0.05,angle=90,code=3)
	points(0.5,EY0$psi, pch=21,cex=1.75, lwd=1,bg="white",col=cols[1])
	
	arrows(x0=0.5,y0=EY1$lb, y1=EY1$ub,lwd=1,col=cols[2],length=0.05,angle=90,code=3)
	points(0.5,EY1$psi, pch=16,cex=1.75, lwd=1,bg="white",col=cols[2])
	
	# if(Ediff$p<0.01) text(0.5,0.1,"*",cex=2)
	
}



pdf("~/dropbox/articles/antibody-curves/results/figs/haiti2-USA-enterics-Ab-curves.pdf",height=10,width=20)
lo <- layout(mat=matrix(1:16,nrow=2,ncol=8,byrow=TRUE),widths=rep(c(1,0.2),4))
# lo <- layout(mat=matrix(1:2,nrow=1,ncol=2),widths=c(1,0.2))

SLAb.plotEYax(usa.cp17,haiti.cp17,expression(paste(italic('Cryptosporidium parvum'), " Cp17")),"a",ylabel=T)
SLAb.plotEYx(EYx.usa.cp17,EYx.haiti.cp17,diff.cp17) 

SLAb.plotEYax(usa.cp23,haiti.cp23,expression(paste(italic('Cryptosporidium parvum'), " Cp23")),"b")
SLAb.plotEYx(EYx.usa.cp23,EYx.haiti.cp23,diff.cp23) 

SLAb.plotEYax(usa.giar,haiti.giar,expression(paste(italic('Giardia intestinalis'), " VSP-5")),"c")
SLAb.plotEYx(EYx.usa.giar,EYx.haiti.giar,diff.giar) 

SLAb.plotEYax(usa.leca,haiti.leca,expression(paste(italic('Entamoeba histolytica'), " LecA")),"d")
SLAb.plotEYx(EYx.usa.leca,EYx.haiti.leca,diff.leca) 

SLAb.plotEYax(usa.etec,haiti.etec,expression(paste("ETEC heat labile toxin ",beta," subunit")),"e",ylabel=T,xlabel=T)
SLAb.plotEYx(EYx.usa.etec,EYx.haiti.etec,diff.etec) 

SLAb.plotEYax(usa.salb,haiti.salb,expression(paste(italic('Salmonella sp.'), " LPS Group B")),"f",xlabel=T)
SLAb.plotEYx(EYx.usa.salb,EYx.haiti.salb,diff.salb) 

SLAb.plotEYax(usa.norogi,haiti.norogi,"Norovirus GI.4","g",xlabel=T)
SLAb.plotEYx(EYx.usa.norogi,EYx.haiti.norogi,diff.norogi) 

SLAb.plotEYax(usa.norogii,haiti.norogii,"Norovirus GII.4 NO","h",xlabel=T)
SLAb.plotEYx(EYx.usa.norogii,EYx.haiti.norogii,diff.norogii) 

dev.off()

	