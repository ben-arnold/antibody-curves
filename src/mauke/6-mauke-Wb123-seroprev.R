

#-------------------------------------------
# 6-mauke-Wb123-seroprev.R
#
# estimate age-specific antibody curves
# after reducing the quantitative Wb123
# antibody response to a binary pos/neg
# measurement. the analysis parallels the 
# main analysis using the quantitative response
# (see 1-mauke-Wb123-analysis.R)
#
#-------------------------------------------

#-------------------------------------------
# input files:
#   mauke1975-public.csv
#   mauke1992-public.csv
#
# output files:
#   mauke-Wb123-seroprev.RData
#   mauke-Wb123-seroprev.pdf
#-------------------------------------------



#-------------------------------------------
# preamble
#-------------------------------------------

rm(list=ls())
library(tmle)
library(SuperLearner)
library(tmleAb)

#-------------------------------------------
# load the Mauke data from 1974(1975) and 1992
#-------------------------------------------

data("mauke_wb123")
d  <- mauke_wb123

# add 0.5 years to age to remove bias (on average) due to rounding to year
d$ager <- d$age+0.5

# create age categories for stratified analyses
# 5 year age categories (1-20 y)
d$agecat <- cut(d$ager,breaks=c(0,5,10,15,20),labels=c("1-5","6-10","11-15","16-20"))

# identify the pre- vs. post-MDA measurements
d$mda <- ifelse(d$year=="1992",1,0)

# make an unique individual id variable
d$id <- ifelse(!is.na(d$id75),d$id75,d$id92)

# subset to common variables
common.vars <- c("id","ager","agecat","wb123","wb123pos","year","mda")
d <- subset(d,select=common.vars)

# recode the dichotomous Wb123 variable to 0/1
d$wb123p <- ifelse(d$wb123pos=="Positive",1,0)

# add a variable "n" <- 1 for fitting the reversible catalytic model, below
d$n <- 1


# subset to each year as well, for convenience
d75 <- subset(d,year=="1975")
d92 <- subset(d,year=="1992")


#--------------------------------------
# All Ages
#--------------------------------------
SL.library <- c("SL.mean","SL.glm","SL.gam","SL.loess")



# SuperLearner fits of antibody levels
set.seed(12345)
mauke75 <- agecurveAb(Y=d75$wb123p,Age=d75$ager,id=d75$id,SL.library=SL.library,gamdf=2:3)
set.seed(12345)
mauke92 <- agecurveAb(Y=d92$wb123p,Age=d92$ager,id=d92$id,SL.library=SL.library,gamdf=2:3)

# estimate group means
EYx.mauke75 <- tmleAb(Y=d75$wb123p,W=d75$ager,id=d75$id,SL.library=SL.library)
EYx.mauke92 <- tmleAb(Y=d92$wb123p,W=d92$ager,id=d92$id,SL.library=SL.library)

# estimate difference in means
set.seed(12345)
diff.mauke  <- tmleAb(Y=d$wb123p,W=d$ager,id=d$id,X=d$mda,SL.library=SL.library)


#--------------------------------------
# Estimate means and differences between
# time points in
# 5 year age bands from ages 1-20
#--------------------------------------
agegrps <-c("1-5","6-10","11-15","16-20") 
EYx.mauke75kids <- sapply(agegrps, function(x) 
  tmleAb(Y=d75$wb123p[d75$agecat==x],W=d75$ager[d75$agecat==x],id=d75$id[d75$agecat==x], SL.library=SL.library) 
)
EYx.mauke92kids <- sapply(agegrps, function(x) 
  tmleAb(Y=d92$wb123p[d92$agecat==x],W=d92$ager[d92$agecat==x],id=d92$id[d92$agecat==x], SL.library=SL.library) 
)
diff.maukekids <- sapply(agegrps, function(x) 
  tmleAb(Y=d$wb123p[d$agecat==x], W=d$ager[d$agecat==x], id=d$id[d$agecat==x], X=d$mda[d$agecat==x], SL.library=SL.library) 
)



#-------------------------------------- 
# save results
#-------------------------------------- 
save.image("~/dropbox/articles/antibody-curves/results/raw/mauke-Wb123-seroprev.RData")

#--------------------------------------
# make figure
#--------------------------------------

pdf("~/dropbox/articles/antibody-curves/results/figs/mauke-Wb123-seroprev.pdf",width=12,height=6)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cols <- cbPalette[c(7,6)]

lo <- layout(mat=matrix(1:2,nrow=1,ncol=2,byrow=TRUE))


# Panel a, age dependent curves
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

# set.seed(623) # jitter ages (but not actual fits) to better display the data in the rug
# rug.low <- ifelse(mauke75$Y[mauke75$Age<=70.5]==0,0,mauke75$Y[mauke75$Age<=70.5]+0.02)
# rug.hgh <- ifelse(mauke75$Y[mauke75$Age<=70.5]==1,1,mauke75$Y[mauke75$Age<=70.5]+0.02)
# segments(x0=jitter(mauke75$Age[mauke75$Age<=70.5],1.5),y0=rug.low,y1=rug.hgh,col=cols[1])
# 
# set.seed(623) # jitter ages (but not actual fits) to better display the data in the rug
# rug.low <- ifelse(mauke92$Y[mauke92$Age<=70.5]==0,0,mauke92$Y[mauke92$Age<=70.5]+0.02)
# rug.hgh <- ifelse(mauke92$Y[mauke92$Age<=70.5]==1,1,mauke92$Y[mauke92$Age<=70.5]+0.02)
# segments(x0=jitter(mauke92$Age[mauke92$Age<=70.5],1.5),y0=rug.low,y1=rug.hgh,col=cols[2])

# Axis labels
mtext(expression(paste(italic('W. bancrofti')," Wb123 seroprevalence (%)")),side=2,line=3,cex=1.25)
mtext("Age, years",side=1,line=3,cex=1.5)
mtext("a",line=1,at=-10,adj=0,font=2,cex=2)
mtext(expression(paste("Seroprevalence by age, ", italic(E),"(",italic(Y[a][","][x]),")")),line=1,cex=1.5)

# Group labels
mtext("1975",side=4,line=0.5,adj=0,at=0.92,col=cols[1],cex=1.25,las=1)
mtext("1992",side=4,line=0.5,adj=0,at=0.82,col=cols[2],cex=1.25,las=1)

par(op)


# Panel B age category E(Y_x) plot
op <- par(mar=c(5,5,3,0)+0.1)

plot(1:4,1:4,type="n",
     ylab="",yaxt="n",ylim=range(ytics),
     xlab="",xaxt="n",xlim=c(0.5,4.5),
     bty="n"
)
axis(2,at=ytics,labels=seq(0,100,by=10), las=1,cex.axis=1.5)
# labels and line segments
mtext(levels(d$agecat),side=1,line=1,at=1:4,cex=1.5)
mtext("Age category, years",side=1,line=3,cex=1.5)

# Y label
mtext("b",line=1,at=-0.07,adj=0,font=2,cex=2)
mtext(expression(paste("Seroprevalence by age category, ",italic(E),"(",italic(Y[x]),")")),line=1,cex=1.5)
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





