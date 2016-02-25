
#-------------------------------
# 2-garki-main-analysis.R
#
# Calculate age-adjusted mean
# IFA antibody titres by
# intervention group and survey round
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
d <- read.csv("~/dropbox/articles/antibody-curves/data/garki/final/garki-sero.csv")

d$mdate <- as.Date(d$mdate,"%d %b %Y")


# add 2 control village names
d$vname <- factor(d$vname,levels=c(levels(d$vname),"Nabanawa","Ajura"))
d$vname[d$village==552] <- "Nabanawa"
d$vname[d$village==553] <- "Ajura"
d$vname <- factor(d$vname)

#-------------------------------
# Serological survey timing:
# 1-2 : pre-intervention
# 3-5 : intervention period
# 6-8 : post-intervention
#-------------------------------

# subset the data by group and survey round for convenience
# combine control measures by phase of the study
d.c12  <- d[d$tr=="Control" & (d$serosvy==1|d$serosvy==2),]
d.c345 <- d[d$tr=="Control" & (d$serosvy==3|d$serosvy==4|d$serosvy==5),]
d.c78  <- d[d$tr=="Control" & (d$serosvy==7|d$serosvy==8),]

d.tr1 <- d[d$tr=="Intervention" & d$serosvy==1,]
d.tr2 <- d[d$tr=="Intervention" & d$serosvy==2,]
d.tr3 <- d[d$tr=="Intervention" & d$serosvy==3,]
d.tr4 <- d[d$tr=="Intervention" & d$serosvy==4,]
d.tr5 <- d[d$tr=="Intervention" & d$serosvy==5,]
d.tr6 <- d[d$tr=="Intervention" & d$serosvy==6,]
d.tr7 <- d[d$tr=="Intervention" & d$serosvy==7,]
d.tr8 <- d[d$tr=="Intervention" & d$serosvy==8,]

#-------------------------------
# Step 1: 
# Age-specific antibody response
# curves
#-------------------------------

set.seed(23752)

# Pre-intervention period fitted curves
p.c12 <- SLAb.curve(Y=log10(d.c12$ifatpftitre+1),Age=d.c12$ageyrs,id=d.c12$id)
p.tr1 <- SLAb.curve(Y=log10(d.tr1$ifatpftitre+1),Age=d.tr1$ageyrs,id=d.tr1$id)
p.tr2 <- SLAb.curve(Y=log10(d.tr2$ifatpftitre+1),Age=d.tr2$ageyrs,id=d.tr2$id)

# Intervention phase fitted curves
p.c345 <- SLAb.curve(Y=log10(d.c345$ifatpftitre+1),Age=d.c345$ageyrs,id=d.c345$id)
p.tr3 <- SLAb.curve(Y=log10(d.tr3$ifatpftitre+1),Age=d.tr3$ageyrs,id=d.tr3$id)
p.tr4 <- SLAb.curve(Y=log10(d.tr4$ifatpftitre+1),Age=d.tr4$ageyrs,id=d.tr4$id)
p.tr5 <- SLAb.curve(Y=log10(d.tr5$ifatpftitre+1),Age=d.tr5$ageyrs,id=d.tr5$id)

# post intervention period fitted curves
# note: no control measurement in round 6
p.c78 <- SLAb.curve(Y=log10(d.c78$ifatpftitre+1),Age=d.c78$ageyrs,id=d.c78$id)
p.tr6 <- SLAb.curve(Y=log10(d.tr6$ifatpftitre+1),Age=d.tr6$ageyrs,id=d.tr6$id)
p.tr7 <- SLAb.curve(Y=log10(d.tr7$ifatpftitre+1),Age=d.tr7$ageyrs,id=d.tr7$id)
p.tr8 <- SLAb.curve(Y=log10(d.tr8$ifatpftitre+1),Age=d.tr8$ageyrs,id=d.tr8$id)


#-------------------------------
# Step 2: 
# Calculate age-adjusted means
# and difference in means
# by treatment group and 
# survey round
#-------------------------------

# create a numeric 0/1 treatment variable
d$tr01 <- ifelse(d$tr=="Control",0,1)


set.seed(5463452)

## Control Villages 
# Nabanawa + Ajura
# (no measurement in survey round 6)
mu.c <- sapply(c(1:5,7:8),function(x) SLAb.tmle(
	Y=log10(d$ifatpftitre[d$tr=="Control" & d$serosvy==x]+1),
	Age=d$ageyrs[d$tr=="Control" & d$serosvy==x],
	id=d$id[d$tr=="Control" & d$serosvy==x]
	)
)

### Intervention Villages - Spraying (Propoxur) + MDA
mu.i <- sapply(1:8,function(x) SLAb.tmle(
	Y=log10(d$ifatpftitre[d$tr=="Intervention" & d$serosvy==x]+1),
	Age=d$ageyrs[d$tr=="Intervention" & d$serosvy==x],
	id=d$id[d$tr=="Intervention" & d$serosvy==x]
	)
)


#-------------------------------
# Estimate difference between
# control and intervention villages in
# IFAT P. falciparm
# titre by survey round
#-------------------------------

set.seed(79287234)
diff.psi <- sapply(c(1:5,7:8),function(x) SLAb.tmle(
	Y=log10(d$ifatpftitre[d$serosvy==x]+1),
	Age=d$ageyrs[d$serosvy==x],
	id=d$id[d$serosvy==x],
	X=d$tr01[d$serosvy==x],
	diff=TRUE
	)
)


# P-values for differences in means, with Bonferroni correction for 7 tests
sprintf("%1.4f",unlist(diff.psi[5,])*7)

#-------------------------------
# save down the results
#-------------------------------
save.image("~/dropbox/articles/antibody-curves/results/raw/garki-main-analysis.RData")



