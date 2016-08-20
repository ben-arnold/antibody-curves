
#-------------------------------
# 8-garki-cvSL.R
#
# Compute the cross-validated
# risk for the super leaner
# and its constituent algorithms
#
# do calculations for survey 5
# for both control and intervention
# groups
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
source("~/SLAbcurves/src/SLyman2016.R")
source("~/SLAbcurves/src/slab_cvSL.R")


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
# limit to survey round 5 (end of intervention period)
d.c5 <- d[d$tr=="Control" & d$serosvy==5,]
  d.c5 <- subset(d.c5,!is.na(d.c5$ageyrs) & !is.na(d.c5$ifatpftitre))
d.tr5 <- d[d$tr=="Intervention" & d$serosvy==5,]
  d.tr5 <- subset(d.tr5,!is.na(d.tr5$ageyrs) & !is.na(d.tr5$ifatpftitre))

  
#-------------------------------
# try out this Yman 2016 model
#-------------------------------

ymanfit <- SL.yman2016(Y=log10(d.c5$ifatpftitre+1),X=d.c5$ageyrs)
ymanfit2 <- SL.yman2016(Y=log10(d.tr5$ifatpftitre+1),X=d.tr5$ageyrs)


plot(d.c5$ageyrs[order(d.c5$ageyrs)],ymanfit$pred[order(d.c5$ageyrs)],type='l',ylim=c(0,4))
lines(d.tr5$ageyrs[order(d.tr5$ageyrs)],ymanfit2$pred[order(d.tr5$ageyrs)],lty=2)
  
#-------------------------------
# fit cross-validated SL
#-------------------------------

SL.library <- c("SL.mean","SL.glm","SL.bayesglm","SL.loess", "SL.gam","SL.glmnet","SL.randomForest","SL.yman2016")

SL.library <- c("SL.glm","SL.loess","SL.yman2016")

cfit <- SuperLearner(Y=log10(d.c5$ifatpftitre+1),X=data.frame(Age=d.c5$ageyrs,Age2=d.c5$ageyrs^2),SL.library=SL.library)

tfit <- SuperLearner(Y=log10(d.tr5$ifatpftitre+1),X=data.frame(Age=d.tr5$ageyrs),SL.library=SL.library)


lines(d.c5$ageyrs[order(d.c5$ageyrs)],cfit$SL.predict[order(d.c5$ageyrs)],col="blue")
lines(d.tr5$ageyrs[order(d.tr5$ageyrs)],tfit$SL.predict[order(d.tr5$ageyrs)],col="blue")


# fit the super learner 
set.seed(23752)
cfit <- slab_cvSL(Y=log10(d.c5$ifatpftitre+1),Age=d.c5$ageyrs,id=d.c5$id,SL.library=SL.library)






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



