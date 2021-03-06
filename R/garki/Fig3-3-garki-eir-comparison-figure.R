#-------------------------------
# Fig3-3-garki-eir-comparison-figure.R
#
# Compare IFA titer 
# age-adjusted mean values
# with EIR estimates from
# the original Garki Project
# (Molineaux 1980, Table 4, p 65)
# for the 3 villages with both
# serological and entomological
# analyses at multiple wet season
# time points
#
#-------------------------------



#-------------------------------
# preamble
#-------------------------------

rm(list=ls())
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
# load the EIR / quantitative
# analysis results
#-------------------------------
load("~/dropbox/articles/antibody-curves/results/raw/garki-eir-comparison.RData")

# rename a few objects
qrho <- sp.rho
qmd <- md

#-------------------------------
# load the EIR / seroprev
# analysis results
#-------------------------------
load("~/dropbox/articles/antibody-curves/results/raw/garki-eir-comparison-seroprev.RData")

# rename a few objects
prho <- sp.rho
pmd <- md

#-------------------------------
# calculate spearman's rho
# between means and seroprev
#-------------------------------

pqrho <- cor.test(qmd$mu,pmd$mu,method="spearman")
pqrho

#-------------------------------
# format the rho estimates
#-------------------------------
q_rho.text <- substitute(paste(rho," = ",rho.txt ),list(rho.txt=sprintf("%1.2f",qrho$estimate)))
p_rho.text <- substitute(paste(rho," = ",rho.txt ),list(rho.txt=sprintf("%1.2f",prho$estimate)))
pq_rho.text <- substitute(paste(rho," = ",rho.txt ),list(rho.txt=sprintf("%1.2f",pqrho$estimate)))

#-------------------------------
# plot the values
#-------------------------------

pdf("~/dropbox/articles/antibody-curves/results/figs/garki-IFAPf-EIR.pdf",width=10,height=10)

lo <- layout(mat=matrix(1:4,nrow=2,ncol=2,byrow=TRUE))

op <- par(mar=c(6,6,3,1)+0.1,xpd=TRUE)

cols <- c(cgrey,cteal,corange)
ytics <- seq(1,4,by=1)
xtics <- seq(0,2,by=1)
plot(1,1,type="n",bty="n",
     xaxt="n",xlab="",xlim=c(0,2.2),
     yaxt="n",ylab="",ylim=range(ytics),
     las=1
)
axis(1,at=xtics,labels=c(1,10,100),cex.axis=1.5)
axis(2,at=ytics,labels=c(
  expression(10^1),
  expression(10^2),
  expression(10^3),
  expression(10^4)),
  las=1,cex.axis=1.5
)
mtext("a",adj=1,line=0.5,at=-0.4,font=2,cex=2)
mtext(expression(paste("Wet season village ",italic('P. falciparum'))),side=2,line=4,cex=1.25)
mtext(expression(paste("mean IFA antibody titer, ", italic(E(Y[x])) )) ,side=2,line=2.5,cex=1.25)
# mtext("Entomological Inoculation Rate\n(cumulative wet season infectious bites per person)",side=1,line=3.5,cex=1.25)

# add correlation
text(2,1.25,q_rho.text,adj=1,cex=1.5)

# Ajura
points(qmd$log10eir[qmd$vname=="Ajura"],qmd$mu[qmd$vname=="Ajura"], pch=16,cex=2,col=cols[1])

# Rafin Marke
points(qmd$log10eir[qmd$vname=="Rafin Marke"],qmd$mu[qmd$vname=="Rafin Marke"], pch=16, cex=2,col=cols[2])

# Nasakar, exclude 1972 b/c it is out of the scale
points(qmd$log10eir[qmd$vname=="Nasakar" & qmd$wetseason!=1972],qmd$mu[qmd$vname=="Nasakar" & qmd$wetseason!=1972], pch=16,cex=2,col=cols[3])



# circle pre-intervention measures
points(qmd$log10eir[qmd$wetseason==1971],qmd$mu[qmd$wetseason==1971],cex=2)

# label the villages
ajura.x <- qmd$log10eir[qmd$vname=="Ajura" & qmd$wetseason==1971]
ajura.y <- qmd$mu[qmd$vname=="Ajura" & qmd$wetseason==1971]
segments(x0=ajura.x,y0=ajura.y+0.1,y1=ajura.y+0.2,col="gray40")
text(ajura.x,ajura.y+0.2,"Ajura (control)",col=cols[1],pos=3,cex=1.25)

rafin.x <- qmd$log10eir[qmd$vname=="Rafin Marke" & qmd$wetseason==1971]
rafin.y <- qmd$mu[qmd$vname=="Rafin Marke" & qmd$wetseason==1971]
segments(x0=rafin.x-0.07,x1=rafin.x-0.3,y0=rafin.y+0.07,y1=rafin.y+0.25,col="gray40")
text(rafin.x-0.3,rafin.y+0.3,"Rafin Marke",col=cols[2],pos=2,cex=1.25)

nasak.x <- qmd$log10eir[qmd$vname=="Nasakar" & qmd$wetseason==1971]
nasak.y <- qmd$mu[qmd$vname=="Nasakar" & qmd$wetseason==1971]
segments(x0= nasak.x,y0= nasak.y +0.1,y1= nasak.y +0.4,col="gray40")
text(nasak.x, nasak.y +0.4,"Nasakar",col=cols[3],pos=3,cex=1.25)

#-------------------------------
# quantitative vs. seroprevalence panel
#-------------------------------

ytics <- seq(1,4,by=1)
xtics <- seq(0.7,1,by=0.05)
plot(1,1,type="n",bty="n",
     xaxt="n",xlab="",xlim=range(xtics),
     yaxt="n",ylab="",ylim=range(ytics),
     las=1
)
axis(1,at=xtics,labels=sprintf("%1.0f",xtics*100),cex.axis=1.5)
axis(2,at=ytics,labels=c(
  expression(10^1),
  expression(10^2),
  expression(10^3),
  expression(10^4)),
  las=1,cex.axis=1.5
)
mtext("b",adj=1,line=0.5,at=0.65,font=2,cex=2)
# mtext("Wet season village geometric mean",side=2,line=4,cex=1.25)
# mtext(expression(paste(italic('P. falciparum')," IFA antibody titer, ", italic(E(Y[x])) )) ,side=2,line=2.5,cex=1.25)
mtext(expression(paste(italic('P. falciparum')," IFA seroprevalence (%)")),side=1,line=3.5,cex=1.25)

# add correlation
text(1,1.25,pq_rho.text,adj=1,cex=1.5)

# Ajura
points(pmd$mu[pmd$vname=="Ajura"],qmd$mu[qmd$vname=="Ajura"], pch=16,cex=2,col=cols[1])

# Rafin Marke
points(pmd$mu[pmd$vname=="Rafin Marke"],qmd$mu[qmd$vname=="Rafin Marke"], pch=16, cex=2,col=cols[2])

# Nasakar, exclude 1972 b/c it is out of the scale
# points(pmd$mu[pmd$vname=="Nasakar" & pmd$wetseason!=1972],qmd$mu[qmd$vname=="Nasakar" & qmd$wetseason!=1972], pch=16,cex=2,col=cols[3])
points(pmd$mu[pmd$vname=="Nasakar"],qmd$mu[qmd$vname=="Nasakar"], pch=16,cex=2,col=cols[3])

# circle pre-intervention measures
# points(pmd$mu[pmd$wetseason==1971],qmd$mu[qmd$wetseason==1971],cex=2)

# label the villages
# ajura.x <- qmd$log10eir[qmd$vname=="Ajura" & qmd$wetseason==1971]
# ajura.y <- qmd$mu[qmd$vname=="Ajura" & qmd$wetseason==1971]
# segments(x0=ajura.x,y0=ajura.y+0.1,y1=ajura.y+0.2,col="gray40")
# text(ajura.x,ajura.y+0.2,"Ajura (control)",col=cols[1],pos=3,cex=1)
# 
# rafin.x <- qmd$log10eir[qmd$vname=="Rafin Marke" & qmd$wetseason==1971]
# rafin.y <- qmd$mu[qmd$vname=="Rafin Marke" & qmd$wetseason==1971]
# segments(x0=rafin.x-0.07,x1=rafin.x-0.3,y0=rafin.y+0.07,y1=rafin.y+0.25,col="gray40")
# text(rafin.x-0.3,rafin.y+0.3,"Rafin Marke",col=cols[2],pos=2,cex=1)
# 
# nasak.x <- qmd$log10eir[qmd$vname=="Nasakar" & qmd$wetseason==1971]
# nasak.y <- qmd$mu[qmd$vname=="Nasakar" & qmd$wetseason==1971]
# segments(x0= nasak.x,y0= nasak.y +0.1,y1= nasak.y +0.4,col="gray40") 
# text(nasak.x, nasak.y +0.4,"Nasakar",col=cols[3],pos=3,cex=1)



#-------------------------------
# seroprevalence panel
#-------------------------------
ytics <- seq(0.7,1,by=0.05)
xtics <- seq(0,2,by=1)
plot(1,1,type="n",bty="n",
     xaxt="n",xlab="",xlim=c(0,2.2),
     yaxt="n",ylab="",ylim=range(ytics),
     las=1
)
axis(1,at=xtics,labels=c(1,10,100),cex.axis=1.5)
axis(2,at=ytics,labels=sprintf("%1.0f",ytics*100),las=1,cex.axis=1.5)
mtext("c",adj=1,line=0.5,at=-0.4,font=2,cex=2)
mtext(expression(paste("Wet season village ",italic('P. falciparum'))),side=2,line=4,cex=1.25)
mtext("IFA seroprevalence (%)",side=2,line=2.75,cex=1.25)
mtext("Entomological Inoculation Rate\n(cumulative wet season infectious bites per person)",side=1,line=4,cex=1.25)

# add correlation
text(2,0.72,p_rho.text,adj=1,cex=1.5)

# Ajura
points(pmd$log10eir[pmd$vname=="Ajura"],pmd$mu[pmd$vname=="Ajura"], pch=16,cex=2,col=cols[1])

# Rafin Marke
points(pmd$log10eir[pmd$vname=="Rafin Marke"],pmd$mu[pmd$vname=="Rafin Marke"], pch=16, cex=2,col=cols[2])

# Nasakar, exclude 1972 b/c it is out of the scale
points(pmd$log10eir[pmd$vname=="Nasakar" & pmd$wetseason!=1972],pmd$mu[pmd$vname=="Nasakar" & pmd$wetseason!=1972], pch=16,cex=2,col=cols[3])

# circle pre-intervention measures
# points(pmd$log10eir[pmd$wetseason==1971],pmd$mu[pmd$wetseason==1971],cex=2)

# # label the villages
# ajura.x <- pmd$log10eir[pmd$vname=="Ajura" & pmd$wetseason==1971]
# ajura.y <- pmd$mu[pmd$vname=="Ajura" & pmd$wetseason==1971]
# segments(x0=ajura.x,y0=ajura.y+0.01,y1=ajura.y+0.02,col="gray40")
# text(ajura.x,ajura.y+0.02,"Ajura (control)",col=cols[1],pos=3,cex=1)
# 
# rafin.x <- pmd$log10eir[pmd$vname=="Rafin Marke" & pmd$wetseason==1971]
# rafin.y <- pmd$mu[pmd$vname=="Rafin Marke" & pmd$wetseason==1971]
# segments(x0=rafin.x-0.1,x1=rafin.x-0.2,y0=rafin.y,col="gray40")
# text(rafin.x-0.2,rafin.y,"Rafin Marke",col=cols[2],pos=2,cex=1)
# 
# nasak.x <- pmd$log10eir[pmd$vname=="Nasakar" & pmd$wetseason==1971]
# nasak.y <- pmd$mu[pmd$vname=="Nasakar" & pmd$wetseason==1971]
# segments(x0= nasak.x,y0= nasak.y-0.01,y1= nasak.y-0.02,col="gray40")
# text(nasak.x, nasak.y-0.02,"Nasakar",col=cols[3],pos=1,cex=1)

par(op)
lo <- layout(mat=matrix(1))
dev.off()


