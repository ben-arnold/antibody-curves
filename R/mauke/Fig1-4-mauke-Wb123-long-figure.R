

#-------------------------------------------
# Fig1-4-mauke-Wb123-long-figure.R
# Ben Arnold
#
# plot Wb123 antibody levels for individuals
# measured at both time points
#
#-------------------------------------------

#-------------------------------------------
# input files:
#   mauke-Wb123-long.RData
#
# output files:
#   mauke-Wb123-long.pdf
#-------------------------------------------



#-------------------------------------------
# preamble
#-------------------------------------------

rm(list=ls())
library(RColorBrewer)
library(scales)

#-------------------------------------------
# load the Mauke longitudinal analysis output
#-------------------------------------------
load("~/dropbox/articles/antibody-curves/results/raw/mauke-Wb123-long.RData")


#-------------------------------
# Plot schema that will be
# repeated for each Ab group in
# the composite figure
#-------------------------------

SLAb.plotLong <- function(Ab1,Ab2,mu1,mu2,diff,labels=c("Ab-","Ab-"),letter="",ylabel=FALSE) {
  # plot individual level trajectories between measurements
  
  # Ab1 : log10 antibody level for each individual at measurement 1
  # Ab2 : log10 antibody level for each individual at measurement 2
  # mu1 : geometric mean and 95% CIs at measurement 1
  # mu2 : geometric mean and 95% CIs at measurement 2
  # diff: difference between geometric means at measurement 1 and measurement 2
  # labels: Antibody status group labels, length 2. e.g., c("Ab-","Ab+")
  # letter: letter for multi-panel plots (e.g., "a")
  # ylabel: logical. print a label for the Y-axis
  
  op <- par(mar=c(3,6,4,0)+0.1)
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  cols <- cbPalette[c(7,6,2,3)]
  mus <- cbind(mu1,mu2)
  
  plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",xlim=c(0,2),ylim=c(2,6))
  xs <- c(0.5,1.5)
  # axes
  mtext(labels,side=3,line=0.5,at=xs,col=cols[1:2])
  mtext(c("1975","1992"),side=1,line=0,at=xs,col=cols[1:2])
  if (ylabel==TRUE) mtext(expression(paste(italic('W. bancrofti'), " Wb123 (light units)")),side=2,line=3.5,cex=1.25)
  axis(2,at=2:6,labels=c(
    # expression(10^-1),
    # expression(10^0),
    # expression(10^1),
    expression(10^2),
    expression(10^3),
    expression(10^4),
    expression(10^5),
    expression(10^6)
  ), las=1,cex.axis=1.5
  )
  
  # header
  mtext(letter,side=3,line=1.75,font=2,at=-0.3,cex=1.75)
  
  # plot individual trajectories
  segments(x0=xs[1],x1=xs[2],y0=Ab1[Ab2>=Ab1],y1=Ab2[Ab2>=Ab1],col=cols[4])
  segments(x0=xs[1],x1=xs[2],y0=Ab1[Ab2<Ab1],y1=Ab2[Ab2<Ab1],col=cols[3])
  points(rep(xs[1],length(Ab1)),Ab1,col=cols[1],pch=16)
  points(rep(xs[2],length(Ab2)),Ab2,col=cols[2],pch=16)
  # plot geometric means
  arrows(x0=c(0.25,1.75),y0=unlist(mus[3,]),y1=unlist(mus[4,]),col=cols,lwd=1,length=0.05,angle=90,code=3)
  points(c(0.25,1.75),mus[1,],col=cols,cex=1.5,pch=21,bg=c(cols[1],"white"))
  # plot P-values
#   if(diff[5]*3<0.0001) {
#     ptext <- "p < 0.0001"
#   } else {
#     ptext <- ifelse(diff[5]*3>1,"p = 1.000",paste("p =",sprintf("%1.3f",diff[5]*3)) )
#   }
#   mtext(ptext,side=1,line=-1)
  par(op)
  
}

#-------------------------------
# make the plot
#-------------------------------

pdf("~/dropbox/articles/antibody-curves/results/figs/mauke-Wb123-long.pdf",width=10,height=4)
lo <- layout(mat=matrix(1:3,nrow=1,ncol=3,byrow=TRUE))

SLAb.plotLong(Ab1=log10(d$wb123.75[d$Abstatus=="Pos-Neg"]),Ab2=log10(d$wb123.92[d$Abstatus=="Pos-Neg"]),mu1=unlist(EYx_75_Abstatus[,2]),mu2=unlist(EYx_92_Abstatus[,2]),diff=unlist(diff_Abstatus[,2]),labels=c("Ag+","Ag-"),letter="c",ylabel=TRUE)

SLAb.plotLong(Ab1=log10(d$wb123.75[d$Abstatus=="Neg-Neg"]),Ab2=log10(d$wb123.92[d$Abstatus=="Neg-Neg"]),mu1=unlist(EYx_75_Abstatus[,1]),mu2=unlist(EYx_92_Abstatus[,1]),diff=unlist(diff_Abstatus[,1]),labels=c("Ag-","Ag-"),letter="",ylabel=FALSE)

SLAb.plotLong(Ab1=log10(d$wb123.75[d$Abstatus=="Pos-Pos"]),Ab2=log10(d$wb123.92[d$Abstatus=="Pos-Pos"]),mu1=unlist(EYx_75_Abstatus[,3]),mu2=unlist(EYx_92_Abstatus[,3]),diff=unlist(diff_Abstatus[,3]),labels=c("Ag+","Ag+"),letter="",ylabel=FALSE)

dev.off()







