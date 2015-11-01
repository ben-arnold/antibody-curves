

#-------------------------------------------
# 3-mauke-Wb123-supp-figure
# Ben Arnold
#
# Mauke beeswarm plot of age-stratified
# means in 2 year age bands
#
# version 1 (19 oct 2015)
# 
#
#-------------------------------------------

#-------------------------------------------
# input files:
#   mauke-Wb123-analysis.RData
#
# output files:
#   mauke-Wb123-analysis-2y-beeswarm.pdf
#-------------------------------------------



#-------------------------------------------
# preamble
#-------------------------------------------
rm(list=ls())
library(RColorBrewer)
library(scales)
library(beeswarm)

#-------------------------------------------
# load the Mauke analysis results
#-------------------------------------------

load("~/SLAbcurves/results/raw/mauke-Wb123-analysis.RData")


#-------------------------------------------
# plot E(Y_x) in 2 year age bands
# for a supplemental figure
#-------------------------------------------

# cross-tab of groups for sample sizes
Ntab <- table(a7592$agecat2,a7592$mda)
Ntab 

pdf("~/SLAbcurves/results/figs/mauke-Wb123-analysis-2y.pdf",width=7,height=6)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cols <- cbPalette[c(7,6)]

op <- par(mar=c(7,5,3,1)+0.1)
ytics <- 2:6

plot(1:5,1:5,type="n",
     ylab="",yaxt="n",ylim=range(ytics),
     xlab="",xaxt="n",xlim=c(0.5,5.5),
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
# X labels 
mtext(levels(a7592$agecat2),side=1,line=3.5,at=1:5,cex=1.5)
mtext("Age Category, Years",side=1,line=5.5,cex=1.5)

# Y label
mtext(expression(paste(italic('Wuchereria bancrofti')," Wb123 (Light Units)")),side=2,line=3,cex=1.25)

# header
mtext(expression(paste(italic(E),"(",italic(Y[x]),") stratified by child age")),line=1,cex=1.5)

# add in geometric means	
arrows(x0=1:5, y0=unlist(EYx.mauke75kids2y[3,]), y1=unlist(EYx.mauke75kids2y[4,]), col=cols[1],lwd=2,length=0.07,angle=90,code=3)
points(1:5,unlist(EYx.mauke75kids2y[1,]),pch=16,cex=1.75,bg="white",col=cols[1],lwd=2)

arrows(x0=1:5, y0=unlist(EYx.mauke92kids2y[3,]), y1=unlist(EYx.mauke92kids2y[4,]), col=cols[2],lwd=2,length=0.07,angle=90,code=3)
points(1:5,unlist(EYx.mauke92kids2y[1,]),pch=21,cex=1.75,bg="white",col=cols[2],lwd=2)

# add in Ns
mtext("N obs, 1975",side=1,line=1,at=0.5,adj=1,col=cols[1])
  mtext(Ntab[,1],side=1,line=1,at=1:5,col=cols[1])
mtext("N obs, 1992",side=1,line=2,at=0.5,adj=1,col=cols[2])
  mtext(Ntab[,2],side=1,line=2,at=1:5,col=cols[2])

# label data series
text(5,EYx.mauke75kids2y[1,5],"1975",col=cols[1],pos=4)
text(5,EYx.mauke92kids2y[1,5],"1992",col=cols[2],pos=4)

# add Bonferroni-corrected p-values (multiplied by 5 for 5 tests)
ptext <- ifelse(unlist(diff.maukekids2y[5,])*5<0.001,"p < 0.001",paste("p =",sprintf("%1.3f",unlist(diff.maukekids2y[5,])*5)))
mtext(ptext,side=1,line=-2,at=1:5,cex=0.75)

par(op)
dev.off()




