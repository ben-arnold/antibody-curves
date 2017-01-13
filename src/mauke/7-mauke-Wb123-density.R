

#-------------------------------------------
# 7-mauke-Wb123-density.R
#
# estimate the kernel density of the Wb123
# MFI values in 1975 and 1992
#
#-------------------------------------------

#-------------------------------------------
# input files:
#   mauke_wb123.RData
#
# output files:
#   mauke-Wb123-density.pdf
#-------------------------------------------



#-------------------------------------------
# preamble
#-------------------------------------------

rm(list=ls())
library(tmleAb)
library(scales)

#-------------------------------------------
# load the Mauke data from 1974(1975) and 1992
#-------------------------------------------

data("mauke_wb123")
d  <- mauke_wb123

# subset to each year as well, for convenience
d75 <- subset(d,year=="1975")
d92 <- subset(d,year=="1992")



#-------------------------------------------
# density distribution plot
#-------------------------------------------
pdf("~/dropbox/articles/antibody-curves/results/figs/mauke-Wb123-density.pdf",width=6,height=5)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cols <- cbPalette[c(7,6)]


den75 <- density(log10(d75$wb123))
den92 <- density(log10(d92$wb123))
xtics <- 1:7
plot(den75,type="n",
     main="",
     ylim=c(0,1.2),
     xlim=range(xtics),xaxt="n",xlab="",
     las=1,bty="n"
)
polygon(den75,col=alpha(cols[1],alpha=0.3),border=cols[1])
polygon(den92,col=alpha(cols[2],alpha=0.3),border=cols[2])

axis(1,at=1:6,labels=c(
  expression(10^1),
  expression(10^2),
  expression(10^3),
  expression(10^4),
  expression(10^5),
  expression(10^6)
), las=1,cex.axis=1
)

# labels
mtext("Wb123",side=3,line=1,cex=1.25)
mtext("LIPS Light Units",side=1,line=2.5,cex=1)
text(6.1,0.8,"1975",col=cols[1],cex=0.75)
text(3,0.35,"1992",col=cols[2],cex=0.75)
segments(x0=log10(10968),y0=0,y1=1.2,lty=2)
text(log10(10968),1.2,"Seropositive",pos=4,cex=0.75)

dev.off()






