
#-------------------------------
# Fig4-1-garki-wet-dry-analysis.R
#
# Calculate age-adjusted mean
# IFA Pf antibody titres by
# treatment group and wet/dry season
# in the control villages
#
#-------------------------------



#-------------------------------
# preamble
#-------------------------------
rm(list=ls())
library(SuperLearner)
library(tmle)
library(tmleAb)


#-------------------------------
# load the serology dataset
#-------------------------------
# d <- read.csv("~/dropbox/articles/antibody-curves/data/garki/final/garki-sero.csv")

d <- garki_sero

d$mdate <- as.Date(d$mdate,"%d %b %Y")

# subset to observations with non-missing Ab + age measures
d <- subset(d,!is.na(d$ifatpftitre) & !is.na(d$ageyrs))

# add 2 control village names
d$vname <- factor(d$vname,levels=c(levels(d$vname),"Nabanawa","Ajura"))
d$vname[d$village==552] <- "Nabanawa"
d$vname[d$village==553] <- "Ajura"
d$vname <- factor(d$vname)

# set sex to factor
d$sex <- as.factor(d$sex)

# for age exactly equal to 0, set it equal to 0.001
# to prevent the Yman 2016 model from blowing up
# (model is undefined at age=0)
d$ageyrs[d$ageyrs<=0] <- 0.001

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

# model / algorithm library
SL.library <- c("SL.mean","SL.glm","SL.Yman2016","SL.gam","SL.loess")

set.seed(35234)

# Pre-intervention period fitted curves
# stratified by season, averaged over control + intervention
p.1 <-agecurveAb(Y=log10(d.1$ifatpftitre+1),Age=d.1$ageyrs,id=d.1$id,SL.library=SL.library)
p.2 <-agecurveAb(Y=log10(d.2$ifatpftitre+1),Age=d.2$ageyrs,id=d.2$id,SL.library=SL.library)

# surveys 1-2 stratified by intervention assignment and season
# doing this extra leg work to confirm that there
# is no difference between assignment
p.c1 <- agecurveAb(Y=log10(d.c1$ifatpftitre+1),Age=d.c1$ageyrs,id=d.c1$id,SL.library=SL.library)
p.c2 <- agecurveAb(Y=log10(d.c2$ifatpftitre+1),Age=d.c2$ageyrs,id=d.c2$id,SL.library=SL.library)
p.tr1 <- agecurveAb(Y=log10(d.tr1$ifatpftitre+1),Age=d.tr1$ageyrs,id=d.tr1$id,SL.library=SL.library)
p.tr2 <- agecurveAb(Y=log10(d.tr2$ifatpftitre+1),Age=d.tr2$ageyrs,id=d.tr2$id,SL.library=SL.library)

# control village fitted curves for surveys 3,4,5,7,8 (wet, dry, wet, wet, dry)
p.c3 <- agecurveAb(Y=log10(d.c3$ifatpftitre+1),Age=d.c3$ageyrs,id=d.c3$id,SL.library=SL.library)
p.c4 <- agecurveAb(Y=log10(d.c4$ifatpftitre+1),Age=d.c4$ageyrs,id=d.c4$id,SL.library=SL.library)
p.c5 <- agecurveAb(Y=log10(d.c5$ifatpftitre+1),Age=d.c5$ageyrs,id=d.c5$id,SL.library=SL.library)
p.c7 <- agecurveAb(Y=log10(d.c7$ifatpftitre+1),Age=d.c7$ageyrs,id=d.c7$id,SL.library=SL.library)
p.c8 <- agecurveAb(Y=log10(d.c8$ifatpftitre+1),Age=d.c8$ageyrs,id=d.c8$id,SL.library=SL.library)


#-------------------------------
# Step 2:
# Calculate age-adjusted means
# and difference in means
# by wet/dry season survey
#-------------------------------
set.seed(8969754)

## All villages, rounds 1-2
mu12.0to4 <- sapply(c(1:2),function(x) tmleAb(
	Y=log10(d$ifatpftitre[d$agecat=="0-4" & d$serosvy==x]+1),
	W=data.frame(Age=d$ageyrs[d$agecat=="0-4" & d$serosvy==x]),
	id=d$id[d$agecat=="0-4" & d$serosvy==x],
	SL.library=SL.library
	)[c("psi","se","lb","ub","p")]
)
mu12.5to9 <- sapply(c(1:2),function(x) tmleAb(
  Y=log10(d$ifatpftitre[d$agecat=="5-9" & d$serosvy==x]+1),
  W=data.frame(Age=d$ageyrs[d$agecat=="5-9" & d$serosvy==x]),
  id=d$id[d$agecat=="5-9" & d$serosvy==x],
  SL.library=SL.library
  )[c("psi","se","lb","ub","p")]
)
mu12.10to14 <- sapply(c(1:2),function(x) tmleAb(
  Y=log10(d$ifatpftitre[d$agecat=="10-14" & d$serosvy==x]+1),
  W=data.frame(Age=d$ageyrs[d$agecat=="10-14" & d$serosvy==x]),
  id=d$id[d$agecat=="10-14" & d$serosvy==x],
  SL.library=SL.library
  )[c("psi","se","lb","ub","p")]
)


# Control villages only, round 3-5,7,8
mu.0to4 <- sapply(c(1:5,7,8),function(x) tmleAb(
  Y=log10(d$ifatpftitre[d$agecat=="0-4" & d$tr=="Control" & d$serosvy==x]+1),
  W=data.frame(Age=d$ageyrs[d$agecat=="0-4" & d$tr=="Control" & d$serosvy==x]),
  id=d$id[d$agecat=="0-4" & d$tr=="Control" & d$serosvy==x],
  SL.library=SL.library
  )[c("psi","se","lb","ub","p")]
)
mu.5to9 <- sapply(c(1:5,7,8),function(x) tmleAb(
  Y=log10(d$ifatpftitre[d$agecat=="5-9" & d$tr=="Control" & d$serosvy==x]+1),
  W=data.frame(Age=d$ageyrs[d$agecat=="5-9" & d$tr=="Control" & d$serosvy==x]),
  id=d$id[d$agecat=="5-9" & d$tr=="Control" & d$serosvy==x],
  SL.library=SL.library
  )[c("psi","se","lb","ub","p")]
)
mu.10to14 <- sapply(c(1:5,7,8),function(x) tmleAb(
  Y=log10(d$ifatpftitre[d$agecat=="10-14" & d$serosvy==x]+1),
  W=data.frame(Age=d$ageyrs[d$agecat=="10-14" & d$serosvy==x]),
  id=d$id[d$agecat=="10-14" & d$serosvy==x],
  SL.library=SL.library
  )[c("psi","se","lb","ub","p")]
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
tc.diff.psi1 <- sapply(ages,function(x) tmleAb(
	Y=log10(d$ifatpftitre[d$serosvy==1 & d$agecat==x]+1),
	X=d$tr01[d$serosvy==1 & d$agecat==x],
	W=data.frame(Age=d$ageyrs[d$serosvy==1 & d$agecat==x]),
	id=d$id[d$serosvy==1 & d$agecat==x],
	SL.library=SL.library
	)
)
tc.diff.psi2 <- sapply(ages,function(x) tmleAb(
  Y=log10(d$ifatpftitre[d$serosvy==2 & d$agecat==x]+1),
  X=d$tr01[d$serosvy==2 & d$agecat==x],
  W=data.frame(Age=d$ageyrs[d$serosvy==2 & d$agecat==x]),
  id=d$id[d$serosvy==2 & d$agecat==x],
  SL.library=SL.library
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
  diffs <- sapply(ages,function(x) tmleAb(
    Y=log10(d$ifatpftitre[d$agecat==x]+1),
    X=d$wet[d$agecat==x],
    W=data.frame(Age=d$ageyrs[d$agecat==x]),
    id=d$id[d$agecat==x],
    SL.library=SL.library
    )[c("psi","se","lb","ub","p")]
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
save.image("~/dropbox/articles/antibody-curves/results/raw/garki-wet-dry-analysis.RData")



