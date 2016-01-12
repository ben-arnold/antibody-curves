
#-------------------------------
# 6-garki-wet-dry-analysis.R
#
# Calculate age-adjusted mean
# IFA Pf antibody titres by
# treatment group and wet/dry season
# in the pre-intervention period
# (i.e., stratified by svy 1, 2)
#
#-------------------------------



#-------------------------------
# preamble
#-------------------------------
rm(list=ls())
library(SuperLearner)
library(tmle)

# source the base functions for
# SL fits of age-antibody curves
# and TMLE estimates of mean differences
source("~/SLAbcurves/src/SLAb-curve.R")
source("~/SLAbcurves/src/SLAb-tmle.R")
source("~/SLAbcurves/src/SLAb-cvRF.R")

#-------------------------------
# load the serology dataset
#-------------------------------
d <- read.csv("~/dropbox/garki/data/final/garki-sero.csv")

d$mdate <- as.Date(d$mdate,"%d %b %Y")


# add 2 control village names
d$vname <- factor(d$vname,levels=c(levels(d$vname),"Nabanawa","Ajura"))
d$vname[d$village==552] <- "Nabanawa"
d$vname[d$village==553] <- "Ajura"
d$vname <- factor(d$vname)

# factor for child age group categories
d$agecat <- cut(d$ageyrs,c(-1,5,10,15),labels=c("0-4","5-9","10-14"))

# create a numeric 0/1 treatment variable
d$tr01 <- ifelse(d$tr=="Control",0,1)

# create a numeric 0/1 wet season variable
d$wet <- ifelse(d$wetseason=="Wet",1,0)


#-------------------------------
# subset the data to the pre-
# intervention period (svy 1-2)
#-------------------------------
# subset the data by group and survey round for convenience
d.c1  <- d[d$tr=="Control" & (d$serosvy==1),]
d.c2  <- d[d$tr=="Control" & (d$serosvy==2),]
d.c3  <- d[d$tr=="Control" & (d$serosvy==3),]
d.c4  <- d[d$tr=="Control" & (d$serosvy==4),]
d.c5  <- d[d$tr=="Control" & (d$serosvy==5),]
d.c7  <- d[d$tr=="Control" & (d$serosvy==7),]
d.c8  <- d[d$tr=="Control" & (d$serosvy==8),]

d.tr1 <- d[d$tr=="Intervention" & d$serosvy==1,]
d.tr2 <- d[d$tr=="Intervention" & d$serosvy==2,]

# all villages combined (pre-intervention)
d.1 <- d[d$serosvy==1,]
d.2 <- d[d$serosvy==2,]

#-------------------------------
# Step 1: 
# Age-specific antibody response
# curves
#-------------------------------

set.seed(65452)

# Pre-intervention period fitted curves
# stratified by season, averaged over control + intervention
p.1 <-SLAb.curve(Y=log10(d.1$ifatpftitre+1),Age=d.1$ageyrs,id=d.1$id)
p.2 <-SLAb.curve(Y=log10(d.2$ifatpftitre+1),Age=d.2$ageyrs,id=d.2$id)

# surveys 1-2 stratified by intervention assignment and season
# doing this extra leg work to confirm that there
# is no difference between assignment
p.c1 <- SLAb.curve(Y=log10(d.c1$ifatpftitre+1),Age=d.c1$ageyrs,id=d.c1$id)
p.c2 <- SLAb.curve(Y=log10(d.c2$ifatpftitre+1),Age=d.c2$ageyrs,id=d.c2$id)
p.tr1 <- SLAb.curve(Y=log10(d.tr1$ifatpftitre+1),Age=d.tr1$ageyrs,id=d.tr1$id)
p.tr2 <- SLAb.curve(Y=log10(d.tr2$ifatpftitre+1),Age=d.tr2$ageyrs,id=d.tr2$id)

# control village fitted curves for surveys 3,4,5,7,8 (wet, dry, wet, wet, dry)
p.c3 <- SLAb.curve(Y=log10(d.c3$ifatpftitre+1),Age=d.c3$ageyrs,id=d.c3$id)
p.c4 <- SLAb.curve(Y=log10(d.c4$ifatpftitre+1),Age=d.c4$ageyrs,id=d.c4$id)
p.c5 <- SLAb.curve(Y=log10(d.c5$ifatpftitre+1),Age=d.c5$ageyrs,id=d.c5$id)
p.c7 <- SLAb.curve(Y=log10(d.c7$ifatpftitre+1),Age=d.c7$ageyrs,id=d.c7$id)
p.c8 <- SLAb.curve(Y=log10(d.c8$ifatpftitre+1),Age=d.c8$ageyrs,id=d.c8$id)




# temp code to look at curves
# lo <- layout(mat=matrix(1:3,nrow=3,ncol=1))
# cols <- rainbow(n=7,v=0.7)
# plot(p.1$Age,p.1$pY,type="l",bty="n",ylim=c(1,4),xlim=c(0,15),col=cols[1])
# lines(p.2$Age,p.2$pY,col=cols[2])
# 
# plot(p.c3$Age,p.c3$pY,type="l",bty="n",ylim=c(1,4),xlim=c(0,15),col=cols[3])
# lines(p.c4$Age,p.c4$pY,col=cols[4])
# lines(p.c5$Age,p.c5$pY,col=cols[5])
# 
# plot(p.c7$Age,p.c7$pY,type="l",bty="n",ylim=c(1,4),xlim=c(0,15),col=cols[6])
# lines(p.c8$Age,p.c8$pY,col=cols[7])



#-------------------------------
# Step 2: 
# Calculate age-adjusted means
# and difference in means
# by wet/dry season survey
#-------------------------------
set.seed(8969754)

## All villages, rounds 1-2
mu12.0to4 <- sapply(c(1:2),function(x) SLAb.tmle(
	Y=log10(d$ifatpftitre[d$agecat=="0-4" & d$serosvy==x]+1),
	Age=d$ageyrs[d$agecat=="0-4" & d$serosvy==x],
	id=d$id[d$agecat=="0-4" & d$serosvy==x]
	)
)
mu12.5to9 <- sapply(c(1:2),function(x) SLAb.tmle(
  Y=log10(d$ifatpftitre[d$agecat=="5-9" & d$serosvy==x]+1),
  Age=d$ageyrs[d$agecat=="5-9" & d$serosvy==x],
  id=d$id[d$agecat=="5-9" & d$serosvy==x]
  )
)
mu12.10to14 <- sapply(c(1:2),function(x) SLAb.tmle(
  Y=log10(d$ifatpftitre[d$agecat=="10-14" & d$serosvy==x]+1),
  Age=d$ageyrs[d$agecat=="10-14" & d$serosvy==x],
  id=d$id[d$agecat=="10-14" & d$serosvy==x]
  )
)


# Control villages only, round 3-5,7,8
mu.0to4 <- sapply(c(1:5,7,8),function(x) SLAb.tmle(
  Y=log10(d$ifatpftitre[d$agecat=="0-4" & d$tr=="Control" & d$serosvy==x]+1),
  Age=d$ageyrs[d$agecat=="0-4" & d$tr=="Control" & d$serosvy==x],
  id=d$id[d$agecat=="0-4" & d$tr=="Control" & d$serosvy==x]
  )
)
mu.5to9 <- sapply(c(1:5,7,8),function(x) SLAb.tmle(
  Y=log10(d$ifatpftitre[d$agecat=="5-9" & d$tr=="Control" & d$serosvy==x]+1),
  Age=d$ageyrs[d$agecat=="5-9" & d$tr=="Control" & d$serosvy==x],
  id=d$id[d$agecat=="5-9" & d$tr=="Control" & d$serosvy==x]
  )
)
mu.10to14 <- sapply(c(1:5,7,8),function(x) SLAb.tmle(
  Y=log10(d$ifatpftitre[d$agecat=="10-14" & d$serosvy==x]+1),
  Age=d$ageyrs[d$agecat=="10-14" & d$serosvy==x],
  id=d$id[d$agecat=="10-14" & d$serosvy==x]
  )
)

#-------------------------------
# Estimate difference between
# control and intervention villages in
# IFA P. falciparum
# titre by survey round
# to ensure it is reasonable to
# pool them in rounds 1-2
#-------------------------------

set.seed(8234234)
ages <- levels(d$agecat)
tc.diff.psi1 <- sapply(ages,function(x) SLAb.tmle(
	Y=log10(d$ifatpftitre[d$serosvy==1 & d$agecat==x]+1),
	Age=d$ageyrs[d$serosvy==1 & d$agecat==x],
	id=d$id[d$serosvy==1 & d$agecat==x],
	X=d$tr01[d$serosvy==1 & d$agecat==x],
	diff=TRUE
	)
)
tc.diff.psi2 <- sapply(ages,function(x) SLAb.tmle(
  Y=log10(d$ifatpftitre[d$serosvy==2 & d$agecat==x]+1),
  Age=d$ageyrs[d$serosvy==2 & d$agecat==x],
  id=d$id[d$serosvy==2 & d$agecat==x],
  X=d$tr01[d$serosvy==2 & d$agecat==x],
  diff=TRUE
  )
)

# P-values for differences in means, with Bonferroni correction for 6 tests
sprintf("%1.4f",unlist(tc.diff.psi1[5,])*6)
sprintf("%1.4f",unlist(tc.diff.psi2[5,])*6)

#-------------------------------
# Estimate difference between
# successive wet and dry survey 
# rounds in the 2 control villages
#-------------------------------

# small wrapper for the TMLE function
# to streamline between-round comparisons
# for wet v. dry, stratified by child age
tmle.wrap <- function(d) {
  ages <- levels(d$agecat)
  diffs <- sapply(ages,function(x) SLAb.tmle(
    Y=log10(d$ifatpftitre[d$agecat==x]+1),
    Age=d$ageyrs[d$agecat==x],
    id=d$id[d$agecat==x],
    X=d$wet[d$agecat==x],
    diff=TRUE
    )
  )
  return(diffs)
}
d.12 <- rbind(d.c1,d.c2)
d.23 <- rbind(d.c2,d.c3)
d.34 <- rbind(d.c3,d.c4)
d.45 <- rbind(d.c4,d.c5)


set.seed(3413)
wet.diff.12 <- tmle.wrap(d.12)
wet.diff.23 <- tmle.wrap(d.23)
wet.diff.34 <- tmle.wrap(d.34)
wet.diff.45 <- tmle.wrap(d.45)
# note: 5 v 7 and 7 v 8 are not tested b/c all wet season


#-------------------------------
# save down the results
#-------------------------------
save.image("~/SLAbcurves/results/raw/garki-wet-dry-analysis.RData")



