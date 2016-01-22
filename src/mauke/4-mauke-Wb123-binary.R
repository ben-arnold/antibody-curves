

#-------------------------------------------
# 4-mauke-Wb123-binary.R
#
# estimate age-specific antibody curves
# after reducing the quantitative Wb123
# antibody response to a binary pos/neg
# measurement. the analysis parallels the 
# main analysis using the quantitative response
# (see 1-mauke-Wb123-analysis.R)
#
# Also, compare the results with
# a standard parametric model used for such data --
# the "reversible catalytic model"
#
#-------------------------------------------

#-------------------------------------------
# input files:
#   mauke1974.dta
#   mauke1992.dta
#
# output files:
#   mauke-Wb123-binary.RData
#   mauke-Wb123-binary.pdf
#-------------------------------------------



#-------------------------------------------
# preamble
#-------------------------------------------

rm(list=ls())
library(RColorBrewer)
library(foreign)
library(tmle)
library(SuperLearner)

# source the base functions for
# SL fits of age-antibody curves
# and TMLE estimates of mean differences
source("~/SLAbcurves/src/SLAb-curve.R")
source("~/SLAbcurves/src/SLAb-tmle.R")
source("~/SLAbcurves/src/SLAb-cvRF.R")

#-------------------------------------------
# load the Mauke data from 1974(1975) and 1992
#-------------------------------------------

d75 <- read.dta("~/dropbox/mauke/data/final/mauke1974.dta")
d92 <- read.dta("~/dropbox/mauke/data/final/mauke1992.dta")

# drop 7 children aged 0 (all Wb123+) in 1975 due to maternal antibodies
a75 <- subset(d75,age>0)
a92 <- d92

# add 0.5 years to age to remove bias (on average) due to rounding to year
a75$ager <- a75$age+0.5
a92$ager <- a92$age+0.5

# create age categories for stratified analyses
# 5 year age categories (1-20 y)
a75$agecat <- cut(a75$age,breaks=c(0,5,10,15,20),labels=c("1-5","6-10","11-15","16-20"))
a92$agecat <- cut(a92$age,breaks=c(0,5,10,15,20),labels=c("1-5","6-10","11-15","16-20"))

# 2 year age categories (1-10 y)
a75$agecat2 <- cut(a75$age,breaks=c(0,2,4,6,8,10),labels=c("1-2","3-4","5-6","7-8","9-10"))
a92$agecat2 <- cut(a92$age,breaks=c(0,2,4,6,8,10),labels=c("1-2","3-4","5-6","7-8","9-10"))

# create a standard ID variable before appending
a75$id <- as.integer(a75$id74)
a92$id <- ifelse(is.na(a92$id74),a92$id92,a92$id74)

# identify the pre- vs. post-MDA measurements
a75$mda <- 0
a92$mda <- 1

# recode the dichotomous Wb123 variable to 0/1
a75$wb123p <- ifelse(a75$wb123pos=="Positive",1,0)
a92$wb123p <- ifelse(a92$wb123pos=="Positive",1,0)

# append data with common variables
common.vars <- c("id","ager","agecat","agecat2","wb123","wb123p","mda")
a7592 <- rbind(subset(a75,select=common.vars),subset(a92,select=common.vars))

# add a variable "n" <- 1 for fitting the reversible catalytic model, below
a7592$n <- 1


#--------------------------------------
# All Ages
#--------------------------------------
set.seed(36234)

# SuperLearner fits of antibody levels
mauke75 <- SLAb.curve(Y=a75$wb123p,Age=a75$ager,id=a75$id)
mauke92 <- SLAb.curve(Y=a92$wb123p,Age=a92$ager,id=a92$id)

# estimate group means
EYx.mauke75 <- SLAb.tmle(Y=a75$wb123p,Age=a75$ager,id=a75$id)
EYx.mauke92 <- SLAb.tmle(Y=a92$wb123p,Age=a92$ager,id=a92$id)

# estimate difference in means
diff.mauke  <- SLAb.tmle(Y=a7592$wb123p,Age=a7592$ager,id=a7592$id,X=a7592$mda,diff=TRUE)


#--------------------------------------
# Estimate means and differences between
# time points in
# 5 year age bands from ages 1-20
#--------------------------------------
agegrps <-c("1-5","6-10","11-15","16-20") 
EYx.mauke75kids <- sapply(agegrps, function(x) 
  SLAb.tmle(Y=a75$wb123p[a75$agecat==x],Age=a75$ager[a75$agecat==x],id=a75$id[a75$agecat==x]) 
)
EYx.mauke92kids <- sapply(agegrps, function(x) 
  SLAb.tmle(Y=a92$wb123p[a92$agecat==x],Age=a92$ager[a92$agecat==x],id=a92$id[a92$agecat==x]) 
)
diff.maukekids <- sapply(agegrps, function(x) 
  SLAb.tmle(Y=a7592$wb123p[a7592$agecat==x], Age=a7592$ager[a7592$agecat==x], id=a7592$id[a7592$agecat==x], X=a7592$mda[a7592$agecat==x], diff=TRUE) 
)


#--------------------------------------
# Estimate means and differences between
# time points in
# 2 year age bands from ages 1-10
#-------------------------------------- 
agegrps2 <-c("1-2","3-4","5-6","7-8","9-10")
EYx.mauke75kids2y <- sapply(agegrps2, function(x) 
  SLAb.tmle(Y=a75$wb123p[a75$agecat2==x],Age=a75$ager[a75$agecat2==x],id=a75$id[a75$agecat2==x]) 
)
EYx.mauke92kids2y <- sapply(agegrps2, function(x) 
  SLAb.tmle(Y=a92$wb123p[a92$agecat2==x],Age=a92$ager[a92$agecat2==x],id=a92$id[a92$agecat2==x]) 
)
diff.maukekids2y <- sapply(agegrps2, function(x) 
  SLAb.tmle(Y=a7592$wb123p[a7592$agecat2==x], Age=a7592$ager[a7592$agecat2==x], id=a7592$id[a7592$agecat2==x], X=a7592$mda[a7592$agecat2==x], diff=TRUE) 
)

#-------------------------------------- 
# save results
#-------------------------------------- 
save.image("~/SLAbcurves/results/raw/mauke-Wb123-binary.RData")

#--------------------------------------
# make figure
#--------------------------------------

pdf("~/SLAbcurves/results/figs/mauke-Wb123-binary.pdf",width=12,height=12)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cols <- cbPalette[c(7,6)]

lo <- layout(mat=matrix(1:4,nrow=2,ncol=2,byrow=TRUE))


# Panel A age-specific curve E(Y_x,a), 1975
op <- par(mar=c(5,5,3,4)+0.1)
ytics <- seq(0,1,by=0.1)
plot(mauke75$Age,mauke75$pY,type="n",
     xlab="",xaxt="n",xlim=c(0,70),
     ylab="",yaxt="n",ylim=range(ytics),
     main="",
     las=1,bty="n"
)
axis(1,at=seq(0,70,by=10),cex.axis=1.5)
axis(2,at=ytics,labels=seq(0,100,by=10), las=1,cex.axis=1.5)

lines(mauke75$Age,mauke75$pY,col=cols[1],lwd=2)
set.seed(623) # jitter ages (but not actual fits) to better display the data in the rug
rug.low <- ifelse(mauke75$Y[mauke75$Age<=70.5]==0,0,mauke75$Y[mauke75$Age<=70.5]+0.02)
rug.hgh <- ifelse(mauke75$Y[mauke75$Age<=70.5]==1,1,mauke75$Y[mauke75$Age<=70.5]+0.02)
segments(x0=jitter(mauke75$Age[mauke75$Age<=70.5],1.5),y0=rug.low,y1=rug.hgh,col=cols[1])

# Axis labels
mtext(expression(paste(italic('W. bancrofti')," Wb123 Seroprevalence (%)")),side=2,line=3,cex=1.25)
# mtext("Age, years",side=1,line=3,cex=1.5)
mtext("a",line=1,at=-10,adj=0,font=2,cex=2)
mtext(expression(paste(italic(E),"(",italic(Y[a][","][x]),") in 1975 (pre-MDA)")),line=1,cex=1.5)

par(op)


# Panel B age-specific curve E(Y_x,a), 1992
op <- par(mar=c(5,5,3,0)+0.1)
ytics <- seq(0,1,by=0.1)
plot(mauke92$Age,mauke92$pY,type="n",
     xlab="",xaxt="n",xlim=c(0,70),
     ylab="",yaxt="n",ylim=range(ytics),
     main="",
     las=1,bty="n"
)
axis(1,at=seq(0,70,by=10),cex.axis=1.5)
axis(2,at=ytics,labels=seq(0,100,by=10), las=1,cex.axis=1.5)

lines(mauke92$Age,mauke92$pY,col=cols[2],lwd=2)
set.seed(623) # jitter ages (but not actual fits) to better display the data in the rug
rug.low <- ifelse(mauke92$Y[mauke92$Age<=70.5]==0,0,mauke92$Y[mauke92$Age<=70.5]+0.02)
rug.hgh <- ifelse(mauke92$Y[mauke92$Age<=70.5]==1,1,mauke92$Y[mauke92$Age<=70.5]+0.02)
segments(x0=jitter(mauke92$Age[mauke92$Age<=70.5],1.5),y0=rug.low,y1=rug.hgh,col=cols[2])

# Axis labels
# mtext(expression(paste(italic('W. bancrofti')," Wb123 Seroprevalence (%)")),side=2,line=3,cex=1.25)
mtext("Age, years",side=1,line=3,cex=1.5)
mtext("b",line=1,at=-10,adj=0,font=2,cex=2)
mtext(expression(paste(italic(E),"(",italic(Y[a][","][x]),") in 1992 (post-MDA)")),line=1,cex=1.5)

par(op)

# Panel C, overlay
op <- par(mar=c(5,5,3,4)+0.1)
ytics <- seq(0,1,by=0.1)
plot(mauke92$Age,mauke92$pY,type="n",
     xlab="",xaxt="n",xlim=c(0,70),
     ylab="",yaxt="n",ylim=range(ytics),
     main="",
     las=1,bty="n"
)
axis(1,at=seq(0,70,by=10),cex.axis=1.5)
axis(2,at=ytics,labels=seq(0,100,by=10), las=1,cex.axis=1.5)

lines(mauke75$Age,mauke75$pY,col=cols[1],lwd=2)
lines(mauke92$Age,mauke92$pY,col=cols[2],lwd=2)

# Axis labels
mtext(expression(paste(italic('W. bancrofti')," Wb123 Seroprevalence (%)")),side=2,line=3,cex=1.25)
mtext("Age, years",side=1,line=3,cex=1.5)
mtext("c",line=1,at=-10,adj=0,font=2,cex=2)
mtext(expression(paste(italic(E),"(",italic(Y[a][","][x]),")")),line=1,cex=1.5)

# Group labels
mtext("1975",side=4,line=0.5,adj=0,at=0.92,col=cols[1],cex=1.25,las=1)
mtext("1992",side=4,line=0.5,adj=0,at=0.82,col=cols[2],cex=1.25,las=1)

par(op)


# Panel D age category E(Y_x) plot
op <- par(mar=c(5,5,3,0)+0.1)

plot(1:4,1:4,type="n",
     ylab="",yaxt="n",ylim=range(ytics),
     xlab="",xaxt="n",xlim=c(0.5,4.5),
     bty="n"
)
axis(2,at=ytics,labels=seq(0,100,by=10), las=1,cex.axis=1.5)
# labels and line segments
mtext(levels(a7592$agecat),side=1,line=1,at=1:4,cex=1.5)
mtext("Age Category, Years",side=1,line=3,cex=1.5)

# Y label
mtext("d",line=1,at=-0.3,adj=0,font=2,cex=2)
mtext(expression(paste(italic(E),"(",italic(Y[x]),") stratified by child age")),line=1,cex=1.5)
# mtext(c("1975","1992"),at=c(1,2),col=cols[1:2],side=3,line=-0.5)

# add in seroprevalence estimates
arrows(x0=c(1:4), y0=unlist(EYx.mauke75kids[3,]), y1=unlist(EYx.mauke75kids[4,]), col=cols[1],lwd=2,length=0.05,angle=90,code=3)
points(c(1:4),unlist(EYx.mauke75kids[1,]),pch=16,cex=1.75,bg="white",col=cols[1],lwd=2)

arrows(x0=c(1:4), y0=unlist(EYx.mauke92kids[3,]), y1=unlist(EYx.mauke92kids[4,]), col=cols[2],lwd=2,length=0.05,angle=90,code=3)
points(c(1:4),unlist(EYx.mauke92kids[1,]),pch=21,cex=1.75,bg="white",col=cols[2],lwd=2)

# label data series
text(4,EYx.mauke75kids[1,4],"1975",col=cols[1],pos=4)
text(4,EYx.mauke92kids[1,4],"1992",col=cols[2],pos=4)

par(op)


dev.off()





