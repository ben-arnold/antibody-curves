
#-------------------------------
# 7-garki-wet-dry-analysis-figure.R
#
# Plot antibody curves by
# wet and dry seasons in non-
# intervention conditions
#
#-------------------------------



#-------------------------------
# preamble
#-------------------------------
rm(list=ls())
library(RColorBrewer)
library(scales)

#-------------------------------
# load analysis output
#-------------------------------
load("~/dropbox/articles/antibody-curves/results/raw/garki-wet-dry-analysis.RData")

#-------------------------------
# grab the number of children
# in each age category / survey
# to include in the figure
#-------------------------------
Ns <- table(d$agecat[d$tr=="Control" & !is.na(d$ifatpftitre)],d$serosvy[d$tr=="Control" & !is.na(d$ifatpftitre)])
Ns

#-------------------------------
# NOTE:
# only plotting summaries for
# the lowest 2 age groups
# to streamline the presentation
# since main comparison is the
# pattern in  ages 0-4 vs older
# also, the Ns in the 10-14 age
# category are very small (<25),
# so not very well-estimated
# compared to younger ages, which
# have between 84-119 obs
#-------------------------------

#-------------------------------
# grab p-values for test of
# differences between rounds
# correct the P-values for
# 8 comparisons (4 comparisons
# in each of 2 age groups)
#-------------------------------
p12 <- unlist(wet.diff.12[5,1:2])*8
p23 <- unlist(wet.diff.23[5,1:2])*8
p34 <- unlist(wet.diff.34[5,1:2])*8
p45 <- unlist(wet.diff.45[5,1:2])*8
pvalues <- rbind(p12,p23,p34,p45)

pform <- function(p) {
  ifelse(p<0.001,"p < 0.001",paste("p =",sprintf("%1.3f",p)))
}

pformat <- pform(pvalues)


#-------------------------------
# plot antibody curves and means
#-------------------------------

pdf("~/dropbox/articles/antibody-curves/results/figs/garki-wet-dry-analysis.pdf",width=15,height=5)

# general plotting parameters
ytics <- seq(0,4,by=1)
xtics <- seq(0,10,by=2)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cols <- cbPalette[2:8]
# subset and re-order the colors to correspond to the figure
cols <- cols[c(2,1,3,6,5)]

lo <- layout(mat=matrix(1:3,nrow=1,ncol=3),widths=c(0.6,0.6,1))

# Antibody curves, E(Y_x,a) by survey round (season)
op <- par(mar=c(4,5,7,1)+0.1,xpd=F)
plot(p.c1$pYframe$Age,p.c1$pYframe$pY,type="n",lwd=2,
     ylim=c(0.5,4),ylab="",yaxt="n",
     xlim=range(xtics),xlab="",xaxt="n",
     bty="n",las=1
)
mtext(expression(paste(italic('P. falciparum')," IFA antibody titre")),side=2,line=3)
mtext("a",adj=1,line=4,at=-1.5,font=2,cex=1.75)
mtext("Control villages\n1971-72 wet and dry seasons",adj=0,line=3,at=0,cex=1)
mtext(expression(paste(italic(E),"(",italic(Y[a][","][x]),")")),side=3,line=0)
mtext("Age, years",side=1,line=2.5)
axis(2,at=1:4,labels=c(
  expression(10^1),
  expression(10^2),
  expression(10^3),
  expression(10^4)
  ), las=1,cex.axis=1.25
)
axis(1,at=xtics,cex.axis=1.5)
lines(p.c1$pYframe$Age,p.c1$pYframe$pY,col=cols[1])
lines(p.c2$pYframe$Age,p.c2$pYframe$pY,col=cols[2])

text(1,3.4,"1971 wet",cex=1,adj=0,font=1,col=cols[1])
text(3.5,2.9,"1972 dry season",cex=1,adj=0,font=1,col=cols[2])


plot(p.c3$pYframe$Age,p.c3$pYframe$pY,type="n",lwd=2,
     ylim=c(0.5,4),ylab="",yaxt="n",
     xlim=range(xtics),xlab="",xaxt="n",
     bty="n",las=1
)
# mtext(expression(paste(italic('P. falciparum')," IFA antibody titre")),side=2,line=3)
# mtext("e",adj=1,line=3,at=-1,font=2,cex=1.75)
mtext("Control villages\n1972-73 wet and dry seasons",adj=0,line=3,at=0,cex=1)
mtext(expression(paste(italic(E),"(",italic(Y[a][","][x]),")")),side=3,line=0)
mtext("Age, years",side=1,line=2.5)
axis(2,at=1:4,labels=c(
  expression(10^1),
  expression(10^2),
  expression(10^3),
  expression(10^4)
  ), las=1,cex.axis=1.25
)
axis(1,at=xtics,cex.axis=1.5)
lines(p.c3$pYframe$Age,p.c3$pYframe$pY,col=cols[3])
lines(p.c4$pYframe$Age,p.c4$pYframe$pY,col=cols[4])
lines(p.c5$pYframe$Age,p.c5$pYframe$pY,col=cols[5])

segments(x0=1.3,y0=2.9,y1=3.2,col="gray50",lty=2)
text(1.3,3.3,"1972 wet",cex=1,adj=1,font=1,col=cols[3])
text(4.5,2.9,"1973 dry season",cex=1,adj=0,font=1,col=cols[4])
text(8.5,3.7,"1973 wet",cex=1,adj=0,font=1,col=cols[5])


# Summary mean estimates, E(Y_x) by survey round (season), stratified by child age category
xsep <- c(-0.3,-0.15,0,0.15,0.3)
xs.0to4   <- 1+xsep
xs.5to9   <- 2+xsep
xs <- c(xs.0to4,xs.5to9)
plot(1:2,1:2,type="n",
     xlim=c(0.5,2.5),xaxt="n",xlab="",
     ylim=c(0.5,4),ylab="",yaxt="n",
     las=1,bty="n"
)
axis(2,at=1:4,labels=c(
  expression(10^1),
  expression(10^2),
  expression(10^3),
  expression(10^4)
  ), las=1,cex.axis=1.25
)
# panel title info
mtext("b",adj=1,line=4,at=0.3,font=2,cex=1.75)
mtext("Control villages",side=3,line=4.5,adj=0)
mtext(expression(paste(italic(E),"(",italic(Y[x]),"), stratified by child age and survey")),side=3,line=2.5,adj=0)

# age group labels
mtext(c("<5","5 to 10"),side=1,line=1,at=1:2,cex=1.25)
mtext("Age category, years",side=1,line=2.75,cex=1.25)

# children ages 0-4
arrows(x0=xs.0to4,y0=unlist(mu.0to4[3,1:5]), y1=unlist(mu.0to4[4,1:5]),lwd=1,col=cols,length=0.05,angle=90,code=3)
points(xs.0to4,mu.0to4[1,1:5], pch=21,cex=1.5, lwd=1,bg=cols,col=cols)

# children ages 5-9
arrows(x0=xs.5to9,y0=unlist(mu.5to9[3,1:5]), y1=unlist(mu.5to9[4,1:5]),lwd=1,col=cols,length=0.05,angle=90,code=3)
points(xs.5to9,mu.5to9[1,1:5], pch=21,cex=1.5, lwd=1,bg=cols,col=cols)

# header labels
mtext(c("wet","dry","wet","dry","wet"),side=3,line=1,col=cols,at=xs)
mtext(c("1971","1972","1972","1973","1973"),side=3,line=-0.5,col=cols,at=xs)

# footer labels
mtext("Survey",side=1,line=-2.5,at=0.5,adj=1)
mtext(c(1:5),side=1,line=-2.5,col=cols,at=xs)
mtext("N obs",side=1,line=-1,at=0.5,adj=1)
mtext(c(sprintf("%1.0f",Ns[1,1:5]),sprintf("%1.0f",Ns[2,1:5])),side=1,line=-1,at=xs,col="gray40")


# label significant differences
plabel.y <- 3.8
text(mean(xs.0to4[1:2]),y=plabel.y,labels=pformat[1,1],pos=3,cex=0.7,col="gray40")
segments(x0=xs.0to4[1],x1=xs.0to4[2],y0=plabel.y,col="gray40")
segments(x0=xs.0to4[1],y0=plabel.y,y1=plabel.y-0.1,col="gray40")
segments(x0=xs.0to4[2],y0=plabel.y,y1=plabel.y-0.1,col="gray40")

text(mean(xs.0to4[3:4]),y=plabel.y,labels=pformat[3,1],pos=3,cex=0.7,col="gray40")
segments(x0=xs.0to4[3],x1=xs.0to4[4]-0.01,y0=plabel.y,col="gray40")
segments(x0=xs.0to4[3],y0=plabel.y,y1=plabel.y-0.1,col="gray40")
segments(x0=xs.0to4[4]-0.01,y0=plabel.y,y1=plabel.y-0.1,col="gray40")

text(mean(xs.0to4[4:5]),y=plabel.y,labels=pformat[4,1],pos=3,cex=0.7,col="gray40")
segments(x0=xs.0to4[4]+0.01,x1=xs.0to4[5],y0=plabel.y,col="gray40")
segments(x0=xs.0to4[4]+0.01,y0=plabel.y,y1=plabel.y-0.1,col="gray40")
segments(x0=xs.0to4[5],y0=plabel.y,y1=plabel.y-0.1,col="gray40")

text(mean(xs.5to9[4:5]),y=plabel.y,labels=pformat[4,2],pos=3,cex=0.7,col="gray40")
segments(x0=xs.5to9[4]+0.01,x1=xs.5to9[5],y0=plabel.y,col="gray40")
segments(x0=xs.5to9[4]+0.01,y0=plabel.y,y1=plabel.y-0.1,col="gray40")
segments(x0=xs.5to9[5],y0=plabel.y,y1=plabel.y-0.1,col="gray40")

par(op)
dev.off()

