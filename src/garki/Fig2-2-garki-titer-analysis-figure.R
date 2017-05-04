#-------------------------------
# Fig2-2-garki-titer-analysis-figure.R
#
# summarize the garki antibody curve
# results in a figure

#-------------------------------


#-------------------------------
# preamble
#-------------------------------

rm(list=ls())
library(RColorBrewer)
library(scales)

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
# load the analysis results
#-------------------------------
load("~/dropbox/articles/antibody-curves/results/raw/garki-main-analysis.RData")


#-------------------------------
# Plot antibody response 
# curves for the 3 study phases
#-------------------------------
pdf("~/dropbox/articles/antibody-curves/results/figs/garki-antibody-curves-IFAPf.pdf",height=5,width=15)

# general plotting parameters
ytics <- seq(0,4,by=1)
xtics <- seq(0,20,by=5)

c.cols=c(cblack)
i.cols=c(cmagent,cblue,cteal,cgreen,cchartr,corange,cred,cmagent)

lo <- layout(mat=matrix(1:6,nrow=1,ncol=6,byrow=T),widths=c(1,0.25,1,0.25,1,0.25))

# Panel A: Pre-intervention age-antibody curves
op <- par(mar=c(4,5,7,1)+0.1,xpd=T)
plot(p.c12$Age,p.c12$pY,type="l",lwd=2,col=c.cols,
	ylim=c(0,4),ylab="",yaxt="n",
	xlim=c(0,max(xtics)+1),xlab="",xaxt="n",
	bty="n",las=1
	)
	mtext(expression(paste(italic('P. falciparum')," IFA antibody titer")),side=2,line=3,cex=1.25)
	# mtext("a",adj=1,line=3,at=-2,font=2,cex=1.75)
	mtext("Pre-Intervention Period",adj=0,line=3,at=0,cex=1.5)
	mtext(expression(paste("Age-dependent mean, ",italic(E),"(",italic(Y[a][","][x]),")")),side=3,line=0,adj=0,at=0)
	mtext("Age, years",side=1,line=2.5,cex=1.25)
	axis(2,at=0:4,labels=c(
	  expression(10^0),
	  expression(10^1),
	  expression(10^2),
	  expression(10^3),
	  expression(10^4)
	), las=1,cex.axis=1.75
	)
	# segments(x0=min(xtics),x1=max(xtics),y0=ytics,lty=2,col="gray70")
	axis(1,at=xtics,cex.axis=1.5)
	lines(p.tr1$Age,p.tr1$pY,col=i.cols[1])
	lines(p.tr2$Age,p.tr2$pY,col=i.cols[2])
	
	text(20,3.1,"Control",cex=1.25,adj=1,font=2,col=c.cols)
	text(4,2.6,"Intervention\n(survey 1, 2)",cex=1.25,adj=0,font=1,col=i.cols[2])
	text(1,3.3,"1",cex=1.25,col=i.cols[1])
	text(2,2.5,"2",cex=1.25,col=i.cols[2])
	
# Pre-intervention  mean estimates
op <- par(mar=c(4,0,7,4)+0.1)
plot(1:3,1:3,type="n",
	xlim=c(0.5,3.5),xaxt="n",xlab="",
	ylim=c(0,4),ylab="",yaxt="n",
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

# vertical guides
segments(x0=c(1.5,2.5),y0=0,y1=4,lty=2,col="gray80")

# control 
	arrows(x0=c(1.5,2.5),y0=unlist(mu_c[3,1:2]), y1=unlist(mu_c[4,1:2]),lwd=1,col=c.cols,length=0.05,angle=90,code=3)
	points(c(1.5,2.5),mu_c[1,1:2], pch=21,cex=1.5, lwd=1,bg=c.cols,col=c.cols)
# intervention
	arrows(x0=c(1.5,2.5),y0=unlist(mu_i[3,1:2]), y1=unlist(mu_i[4,1:2]),lwd=1,col=i.cols[1:2],length=0.05,angle=90,code=3)
	points(c(1.5,2.5),mu_i[1,1:2], pch=21,bg="white",cex=1.5,lwd=1, col=i.cols[1:2])
	

# Panel B: Active Intervention  age-antibody curves
op <- par(mar=c(4,5,7,1)+0.1)
plot(p.c345$Age,p.c345$pY,type="l",lwd=2,col=c.cols,
	ylim=c(0,4),ylab="",yaxt="n",
	xlim=c(0,max(xtics)+1),xlab="",xaxt="n",
	bty="n",las=1
	)
	# mtext(expression(paste(italic('P. falciparum')," IFA antibody titer")),side=2,line=3)
	# mtext("b",adj=1,line=3,at=-2,font=2,cex=1.75)
	mtext("Active Intervention Period",adj=0,line=3,at=0,cex=1.5)
	mtext(expression(paste("Age-dependent mean, ",italic(E),"(",italic(Y[a][","][x]),")")),side=3,line=0,adj=0,at=0)
	mtext("Age, years",side=1,line=2.5,cex=1.25)
	axis(2,at=0:4,labels=c(
	  expression(10^0),
	  expression(10^1),
	  expression(10^2),
	  expression(10^3),
	  expression(10^4)
	), las=1,cex.axis=1.75
	)
	# segments(x0=min(xtics),x1=max(xtics),y0=ytics,lty=2,col="gray70")
	axis(1,at=xtics,cex.axis=1.5)
	lines(p.tr3$Age,p.tr3$pY,col=i.cols[3])
	lines(p.tr4$Age,p.tr4$pY,col=i.cols[4])
	lines(p.tr5$Age,p.tr5$pY,col=i.cols[5])
	
	text(20,3.6,"Control",cex=1.25,adj=1,font=2,col=c.cols)
	text(20,2.3,"Intervention\n(survey 3, 4, 5)",cex=1.25,adj=1,font=1,col=i.cols[3])
	text(7,2.85,"3",cex=1.25,col=i.cols[3])
	text(9,2.6,"4",cex=1.25,col=i.cols[4])
	text(10.5,2.3,"5",cex=1.25,col=i.cols[5],adj=0)
	
# Active intervention  mean estimates
op <- par(mar=c(4,0,7,4)+0.1)
plot(1:3,1:3,type="n",
	xlim=c(0.5,3.5),xaxt="n",xlab="",
	ylim=c(0,4),ylab="",yaxt="n",
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

# vertical guides
segments(x0=c(1:3),y0=0,y1=4,lty=2,col="gray80")

# control 
	arrows(x0=1:3,y0=unlist(mu_c[3,3:5]), y1=unlist(mu_c[4,3:5]),lwd=1,col=c.cols,length=0.05,angle=90,code=3)
	points(1:3,mu_c[1,3:5], pch=21,cex=1.5, lwd=1,bg=c.cols,col=c.cols)
# intervention
	arrows(x0=1:3,y0=unlist(mu_i[3,3:5]), y1=unlist(mu_i[4,3:5]),lwd=1,col=i.cols[3:5],length=0.05,angle=90,code=3)
	points(1:3,mu_i[1,3:5], pch=21,bg="white",cex=1.5,lwd=1, col=i.cols[3:5])
# mark periods with P<0.01 with Bonferroni corrected pvalues
	# text(1:3,0.1,"*",cex=2)

# Panel C: Post-Intervention  age-antibody curves
op <- par(mar=c(4,5,7,1)+0.1)
plot(p.c78$Age,p.c78$pY,type="l",lwd=2,col=c.cols,
	ylim=c(0,4),ylab="",yaxt="n",
	xlim=c(0,max(xtics)+1),xlab="",xaxt="n",
	bty="n",las=1
	)
	# mtext(expression(paste(italic('P. falciparum')," IFA antibody titer")),side=2,line=3)
	# mtext("c",adj=1,line=3,at=-2,font=2,cex=1.75)
	mtext("Post Intervention Period",adj=0,line=3,at=0,cex=1.5)
	mtext(expression(paste("Age-dependent mean, ",italic(E),"(",italic(Y[a][","][x]),")")),side=3,line=0,adj=0,at=0)
	mtext("Age, years",side=1,line=2.5,cex=1.25)
	axis(2,at=0:4,labels=c(
	  expression(10^0),
		expression(10^1),
		expression(10^2),
		expression(10^3),
		expression(10^4)
		), las=1,cex.axis=1.75
	)
	# segments(x0=min(xtics),x1=max(xtics),y0=ytics,lty=2,col="gray70")
	axis(1,at=xtics,cex.axis=1.5)
	lines(p.tr6$Age,p.tr6$pY,col=i.cols[6])
	lines(p.tr7$Age,p.tr7$pY,col=i.cols[7])
	lines(p.tr8$Age,p.tr8$pY,col=i.cols[8])
	
	text(20,3.7,"Control",cex=1.25,adj=1,font=2,col=c.cols)
	text(20,2.3,"Intervention\n(survey 6, 7, 8)",cex=1.25,adj=1,font=1,col=i.cols[8])
	text(5,1.2,"6",cex=1.25,col=i.cols[6])
	text(4,2.0,"7",cex=1.25,col=i.cols[7])
	text(3,2.8,"8",cex=1.25,col=i.cols[8],adj=0)

# Post-intervention  mean estimates
op <- par(mar=c(4,0,7,4)+0.1)
plot(1:3,1:3,type="n",
	xlim=c(0.5,3.5),xaxt="n",xlab="",
	ylim=c(0,4),ylab="",yaxt="n",
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

# vertical guides
segments(x0=c(1:3),y0=0,y1=4,lty=2,col="gray80")

# control 
	arrows(x0=2:3,y0=unlist(mu_c[3,6:7]), y1=unlist(mu_c[4,6:7]),lwd=1,col=c.cols,length=0.05,angle=90,code=3)
	points(2:3,mu_c[1,6:7], pch=21,cex=1.5, lwd=1,bg=c.cols,col=c.cols)
# intervention
	arrows(x0=1:3,y0=unlist(mu_i[3,6:8]), y1=unlist(mu_i[4,6:8]),lwd=1,col=i.cols[6:8],length=0.05,angle=90,code=3)
	points(1:3,mu_i[1,6:8], pch=21,bg="white",cex=1.5,lwd=1, col=i.cols[6:8])
# mark periods with P<0.01 with Bonferroni corrected pvalues
	# text(2:3,0.1,"*",cex=2)

par(op)


dev.off()



