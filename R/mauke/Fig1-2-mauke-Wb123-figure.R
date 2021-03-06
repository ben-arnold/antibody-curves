

#-------------------------------------------
# Fig1-2-mauke-Wb123-figure.R
# Ben Arnold
#
# Plot age-specific antibody curves from 
# Mauke and age-stratified means
#
#-------------------------------------------

#-------------------------------------------
# input files:
#   mauke-Wb123-analysis.RData
#
# output files:
#   mauke-Wb123-analysis.pdf
#-------------------------------------------



#-------------------------------------------
# preamble
#-------------------------------------------
rm(list=ls())
library(RColorBrewer)
library(scales)

#-------------------------------------------
# load the Mauke analysis results
#-------------------------------------------

load("~/dropbox/articles/antibody-curves/results/raw/mauke-Wb123-analysis.RData")


#-------------------------------------------
# plot
#-------------------------------------------
# cross-tab of groups for sample sizes
table(d$agecat,d$mda)

pdf("~/dropbox/articles/antibody-curves/results/figs/mauke-Wb123-analysis.pdf",width=12,height=6)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cols <- cbPalette[c(7,6)]

lo <- layout(mat=matrix(1:2,nrow=1,ncol=2,byrow=TRUE),widths=c(1,0.75))


# Panel A age-specific antibody curves E(Y_x,a)
op <- par(mar=c(5,5,3,4)+0.1)
ytics <- 2:6
plot(mauke75$Age,mauke75$pY,type="n",
	xlab="",xaxt="n",xlim=c(0,70),
	ylab="",yaxt="n",ylim=range(ytics),
	main="",
	las=1,bty="n"
	)
	axis(1,at=seq(0,70,by=10),cex.axis=1.5)
	axis(2,at=ytics,labels=c(
		expression(10^2),
		expression(10^3),
		expression(10^4),
		expression(10^5),
		expression(10^6)
		), las=1,cex.axis=1.5
	)
	# abline(log(10986),0,lty=2,col="gray60")
	set.seed(623) # jitter ages (but not actual fits) to better display the data
	points(jitter(mauke75$Age),mauke75$Y,cex=0.45,pch=16,col=alpha(cols[1],alpha=0.6))
	points(jitter(mauke92$Age),mauke92$Y,cex=0.45,pch=16,col=alpha(cols[2],alpha=0.6))
	
	lines(mauke75$Age,mauke75$pY,col=cols[1],lwd=2)
	lines(mauke92$Age,mauke92$pY,col=cols[2],lwd=2)
	
	# Axis labels
	mtext(expression(paste(italic('W. bancrofti'),"  Wb123 (light units)")),side=2,line=3,cex=1.5)
	mtext("Age, years",side=1,line=3,cex=1.5)
	mtext("a",line=1,at=-10,adj=0,font=2,cex=2)
	mtext(expression(paste("Age-dependent mean, ",italic(E),"(",italic(Y[a][","][x]),")")),line=1,cex=1.5)
	
	# Group labels
	mtext("1975",side=4,line=0.5,adj=0,at=5.1,col=cols[1],cex=1.25,las=1)
	mtext("1992",side=4,line=0.5,adj=0,at=4.75,col=cols[2],cex=1.25,las=1)
	mtext("(5y Post-MDA)",side=4,line=0.5,adj=0,at=4.5,col=cols[2],cex=1,las=1)
	
par(op)

# skinny plot of E(Y_x) for all ages (Not used)
	# op <- par(mar=c(5,0,3,0)+0.1)
	# plot(1,1,type="n",
		# xlim=c(0,1),xaxt="n",xlab="",
		# ylim=range(ytics),ylab="",yaxt="n",
		# las=1,bty="n"
	# )
	# mtext(expression(paste(italic(E),"(",italic(Y[x]),")")),side=1,line=1)
	# mtext(expression(paste(italic(E),"(",italic(Y[x]),")")),side=3,line=0)
	
	# # plot data
	# segments(x0=0.5,y0=EYx.mauke75$lb, y1=EYx.mauke75$ub,lwd=1,col=cols[1])
	# points(0.5,EYx.mauke75$psi, pch=21,cex=1.75, lwd=1,bg="white",col=cols[1])
	
	# segments(x0=0.5,y0=EYx.mauke92$lb, y1=EYx.mauke92$ub,lwd=1,col=cols[2])
	# points(0.5,EYx.mauke92$psi, pch=21,cex=1.75, lwd=1,bg="white",col=cols[2])
	
# par(op)


# dev.off()

# Panel B age category E(Y_x) plot
op <- par(mar=c(5,5,3,0)+0.1)

plot(1:4,1:4,type="n",
     ylab="",yaxt="n",ylim=range(ytics),
     xlab="",xaxt="n",xlim=c(0.5,4.5),
     bty="n"
     )
axis(side=2,at=ytics,labels=c(
  expression(10^2),
  expression(10^3),
  expression(10^4),
  expression(10^5),
  expression(10^6)
  ), las=1,cex.axis=1.5
)
# X labels and line segments
mtext(levels(d$agecat),side=1,line=1,at=1:4,cex=1.5)
mtext("Age category, years",side=1,line=3,cex=1.5)

# Y label
mtext("b",line=1,at=-0.3,adj=0,font=2,cex=2)
mtext(expression(paste("Mean, ",italic(E),"(",italic(Y[x]),"), stratified by child age")),line=1,cex=1.5)
# mtext(c("1975","1992"),at=c(1,2),col=cols[1:2],side=3,line=-0.5)

# add in geometric means
arrows(x0=c(1:4), y0=unlist(EYx_mauke75kids[3,]), y1=unlist(EYx_mauke75kids[4,]), col=cols[1],lwd=2,length=0.05,angle=90,code=3)
points(c(1:4),unlist(EYx_mauke75kids[1,]),pch=16,cex=1.75,bg="white",col=cols[1],lwd=2)

arrows(x0=c(1:4), y0=unlist(EYx_mauke92kids[3,]), y1=unlist(EYx_mauke92kids[4,]), col=cols[2],lwd=2,length=0.05,angle=90,code=3)
points(c(1:4),unlist(EYx_mauke92kids[1,]),pch=21,cex=1.75,bg="white",col=cols[2],lwd=2)

# label data series
text(4,EYx_mauke75kids[1,4],"1975",col=cols[1],pos=4)
text(4,EYx_mauke92kids[1,4],"1992",col=cols[2],pos=4)

par(op)
dev.off()




