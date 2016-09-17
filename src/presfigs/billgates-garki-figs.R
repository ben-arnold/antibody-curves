
#-------------------------------
# garki-figs.R
#
# make formatted figures for
# the garki progject results
# in the Bill Gates handout
#-------------------------------



#-------------------------------
# preamble
#-------------------------------

rm(list=ls())
library(RColorBrewer)
library(scales)


#-------------------------------
# Plot antibody response 
# curves for the 3 study phases
#-------------------------------
load("~/dropbox/articles/antibody-curves/results/raw/garki-main-analysis.RData")

pdf("~/dropbox/k01/presentations/160503-billgates/figs/garki-curves.pdf",height=5,width=5)

# general plotting parameters
ytics <- seq(0,4,by=1)
xtics <- seq(0,70,by=10)
# i.cols <- c(brewer.pal(9,"YlGnBu")[7:6],brewer.pal(9,"BuPu")[9:7],brewer.pal(9,"YlOrRd")[6:8])
# c.cols <- brewer.pal(8,"Dark2")[8]
# c.cols <- "black"
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
i.cols <- rep(cbPalette[6],8)
c.cols <- cbPalette[7]

lo <- layout(mat=matrix(1:2,nrow=1,ncol=2,byrow=T),widths=c(1,0.15))

# Panel B: Active Intervention  age-antibody curves
op <- par(mar=c(4,5,2,1)+0.1)
plot(p.c345$Age,p.c345$pY,type="l",lwd=2,col=c.cols,
     ylim=c(0.5,4),ylab="",yaxt="n",
     xlim=c(0,71),xlab="",xaxt="n",
     bty="n",las=1
)
mtext(expression(paste(italic('P. falciparum')," IFA antibody titre")),side=2,line=3,cex=1.25)
mtext("a",adj=1,line=0.5,at=-10,font=2,cex=1.75)
# mtext("Active Intervention Period",adj=0,line=3,at=0,cex=1.5)
mtext(expression(paste(italic(E),"(",italic(Y[a][","][x]),")")),side=3,line=0,cex=1.25)
mtext("Age, years",side=1,line=2.5)
axis(2,at=1:4,labels=c(
  expression(10^1),
  expression(10^2),
  expression(10^3),
  expression(10^4)
), las=1,cex.axis=1.25
)
# segments(x0=min(xtics),x1=max(xtics),y0=ytics,lty=2,col="gray70")
axis(1,at=xtics,cex.axis=1)
lines(p.tr3$Age,p.tr3$pY,col=i.cols[3])
lines(p.tr4$Age,p.tr4$pY,col=i.cols[4])
lines(p.tr5$Age,p.tr5$pY,col=i.cols[5])

text(70,3.4,"Control",cex=1,adj=1,font=2,col=c.cols)
text(70,2.5,"Intervention\n(survey 3, 4, 5)",cex=1,adj=1,font=1,col=i.cols[1])
text(11,3.0,"3",cex=0.75,col=i.cols[3])
text(14,2.85,"4",cex=0.75,col=i.cols[4])
text(19,2.6,"5",cex=0.75,col=i.cols[5],adj=0)

# Active intervention  mean estimates
op <- par(mar=c(4,0,2,1)+0.1)
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
mtext(expression(paste(italic(E),"(",italic(Y[x]),")")),side=3,line=0,cex=1.25)
mtext("Survey",side=1,line=2.5)
# control 
arrows(x0=1:3,y0=unlist(mu.c[3,3:5]), y1=unlist(mu.c[4,3:5]),lwd=1,col=c.cols,length=0.05,angle=90,code=3)
points(1:3,mu.c[1,3:5], pch=21,cex=1.5, lwd=1,bg=c.cols,col=c.cols)
# intervention
arrows(x0=1:3,y0=unlist(mu.i[3,3:5]), y1=unlist(mu.i[4,3:5]),lwd=1,col=i.cols[3:5],length=0.05,angle=90,code=3)
points(1:3,mu.i[1,3:5], pch=21,bg="white",cex=1.5,lwd=1, col=i.cols[3:5])
# mark periods with P<0.01 with Bonferroni corrected pvalues
# text(1:3,0.6,"*",cex=2)

par(op)
dev.off()





#-------------------------------
# EIR plot
#-------------------------------
load("~/dropbox/articles/antibody-curves/results/raw/garki-eir-comparison.RData")


pdf("~/dropbox/k01/presentations/160503-billgates/figs/garki-EIR.pdf",width=6,height=6)
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
  mtext("b",adj=1,line=0.5,at=-0.4,font=2,cex=2)
	mtext("Wet season village geometric mean",side=2,line=4,cex=1.25)
  mtext(expression(paste(italic('P. falciparum')," IFA antibody titre, ", italic(E(Y[x])) )) ,side=2,line=2.5,cex=1.25)
	mtext("Entomological Inoculation Rate\n(cumulative wet season infectious bites per person)",side=1,line=3.5,cex=1.25)
	text(2,1,rho.text,adj=1,cex=1.25)
	
	# Ajura
	points(md$log10eir[md$vname=="Ajura"],md$mu[md$vname=="Ajura"], pch=16,cex=2,col=cols[2])

	# Rafin Marke
	points(md$log10eir[md$vname=="Rafin Marke"],md$mu[md$vname=="Rafin Marke"], pch=16, cex=2,col=cols[2])

	# Nasakar, exclude 1972 b/c it is out of the scale
	points(md$log10eir[md$vname=="Nasakar" & md$wetseason!=1972],md$mu[md$vname=="Nasakar" & md$wetseason!=1972], pch=16,cex=2,col=cols[2])

	# circle pre-intervention measures
	# points(md$log10eir[md$wetseason==1971],md$mu[md$wetseason==1971],cex=2)
	
	# label the villages
# 	ajura.x <- md$log10eir[md$vname=="Ajura" & md$wetseason==1971]
# 	ajura.y <- md$mu[md$vname=="Ajura" & md$wetseason==1971]
# 	segments(x0=ajura.x,y0=ajura.y+0.1,y1=ajura.y+0.2,col="gray40")
# 	text(ajura.x,ajura.y+0.2,"Ajura (control)",col=cols[1],pos=3,cex=1)
# 	
# 	rafin.x <- md$log10eir[md$vname=="Rafin Marke" & md$wetseason==1971]
# 	rafin.y <- md$mu[md$vname=="Rafin Marke" & md$wetseason==1971]
# 	segments(x0=rafin.x-0.07,x1=rafin.x-0.3,y0=rafin.y+0.07,y1=rafin.y+0.25,col="gray40")
# 	text(rafin.x-0.3,rafin.y+0.3,"Rafin Marke",col=cols[2],pos=2,cex=1)
# 	
# 	nasak.x <- md$log10eir[md$vname=="Nasakar" & md$wetseason==1971]
# 	nasak.y <- md$mu[md$vname=="Nasakar" & md$wetseason==1971]
#   segments(x0= nasak.x,y0= nasak.y +0.1,y1= nasak.y +0.4,col="gray40") 
#   text(nasak.x, nasak.y +0.4,"Nasakar",col=cols[3],pos=3,cex=1)
	
par(op)
dev.off()



