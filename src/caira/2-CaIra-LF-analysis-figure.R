
#--------------------------------------
# 3-caira-LF-analysis-figure.R
#
# plot SL and TMLE estimates
#
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
# reshape the long form data
# to wide format, merging by
# individual ID, to plot
# individual antibody trajectories
# note: this plot is limited
# to the 73 children who
# were between 4.5 and 11.5 years
# old at measurement round 1
# (for consistency with the EY_x
# estimates, which are limited
# to the age range of overlap
# between measurement rounds)
#-------------------------------
d1 <- subset(ad,round==1,select=c("id","age","wb123","bm14","bm33"))
d3 <- subset(ad,round==3,select=c("id","age","wb123","bm14","bm33"))
dl <- merge(d1,d3,by="id",all.x=TRUE)
dl <- subset(dl,!is.na(age.y))

# summarize the number of children
# for whom antibody levels declined between
# rounds. N=65 children with 2 measures in age range of overlap
table(dl$wb123.y<dl$wb123.x)
table(dl$bm14.y<dl$bm14.x)
table(dl$bm33.y<dl$bm33.x)

#-------------------------------
# antibody curve 
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
	
}

#-------------------------------
# skinny E(Yx) plot schema 
# to add next to age-specific antibody curves
#-------------------------------
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


SLAb.plotLong <- function(Ab1,Ab2,mus,main,letter,ylabel=FALSE) {
  # plot individual level trajectories between measurements
  
  # Ab1 : log10 antibody level for each individual at measurement 1
  # Ab1 : log10 antibody level for each individual at measurement 2
  # main : plot title
  # letter: letter for multi-panel plots (e.g., "a")
  # ylabel: logical. print a label for the Y-axis
  
  op <- par(mar=c(5,6,4,0)+0.1)
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  cols <- cbPalette[c(7,6)]

plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",xlim=c(0,2),ylim=c(0,5))
  xs <- c(0.5,1.5)
  # axes
  mtext(c("2011\npre-MDA","2013\npost-MDA"),side=1,line=2,at=xs,col=cols)
  if (ylabel==TRUE) mtext("Luminex Response (MFI-Background)",side=2,line=3.5,cex=1.25)
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

  # header
  mtext(letter,side=3,line=1.25,font=2,at=-0.2,cex=1.75)
  mtext(main,cex=1.25,line=1.5)

  # plot individual trajectories
  segments(x0=xs[1],x1=xs[2],y0=Ab1[Ab2<Ab1],y1=Ab2[Ab2<Ab1],col=cols[1])
  segments(x0=xs[1],x1=xs[2],y0=Ab1[Ab2>=Ab1],y1=Ab2[Ab2>=Ab1],col=cols[2])
  points(rep(xs[1],length(Ab1)),Ab1,col=cols[1],pch=16)
  points(rep(xs[2],length(Ab2)),Ab2,col=cols[2],pch=16)
  # plot geometric means
  arrows(x0=c(0.25,1.75),y0=unlist(mus[3,]),y1=unlist(mus[4,]),col=cols,lwd=2,length=0.05,angle=90,code=3)
  points(c(0.25,1.75),mus[1,],col=cols,cex=2,pch=21,bg=c(cols[1],"white"))
  par(op)

}


#-------------------------------
# make the figure, using the
# general schema above
#-------------------------------

set.seed(27313)
pdf("~/SLAbcurves/results/figs/caira-LF-curves.pdf",height=15,width=11)
lo <- layout(mat=matrix(1:9,nrow=3,ncol=3,byrow=TRUE),widths=c(1,0.2,0.7))

SLAb.plotEYax(wb123.fit1,wb123.fit3,expression(paste(italic('W. bancrofti'), " Wb123")),"a",xlabel=F,ylabel=T)
SLAb.plotEYx(mu.wb123[,1],mu.wb123[,2],diff.wb123)
SLAb.plotLong(Ab1=log10(dl$wb123.x+1),Ab2=log10(dl$wb123.y+1),mus=mu.wb123,expression(paste(italic('W. bancrofti'), " Wb123")),"d")

SLAb.plotEYax(bm14.fit1,bm14.fit3,expression(paste(italic('W. bancrofti'), " Bm14")),"b",xlabel=F,ylabel=T)
SLAb.plotEYx(mu.bm14[,1],mu.bm14[,2],diff.bm14)
SLAb.plotLong(Ab1=log10(dl$bm14.x+1),Ab2=log10(dl$bm14.y+1),mus=mu.bm14,expression(paste(italic('W. bancrofti'), " Bm14")),"e")


SLAb.plotEYax(bm33.fit1,bm33.fit3,expression(paste(italic('W. bancrofti'), " Bm33")),"c",xlabel=T,ylabel=T)
SLAb.plotEYx(mu.bm33[,1],mu.bm33[,2],diff.bm33) 
SLAb.plotLong(Ab1=log10(dl$bm33.x+1),Ab2=log10(dl$bm33.y+1),mus=mu.bm33,expression(paste(italic('W. bancrofti'), " Bm33")),"f")

dev.off()



