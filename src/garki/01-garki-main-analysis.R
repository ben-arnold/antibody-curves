
#-------------------------------
# 1-garki-main-analysis.R
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
library(tmleAb)

#-------------------------------
# load the serology dataset
#-------------------------------
# d <- read.csv("~/dropbox/articles/antibody-curves/data/garki/final/garki-sero.csv")
d <- garki_sero

d$mdate <- as.Date(d$mdate,"%d %b %Y")


# add 2 control village names
d$vname <- factor(d$vname,levels=c(levels(d$vname),"Nabanawa","Ajura"))
d$vname[d$village==552] <- "Nabanawa"
d$vname[d$village==553] <- "Ajura"
d$vname <- factor(d$vname)

# for age exactly equal to 0, set it equal to 0.001
# to prevent the Yman 2016 model from blowing up
# (model is undefined at age=0)
d$ageyrs[d$ageyrs<=0] <- 0.001

# subset to observations with non-missing Ab + age measures
d <- subset(d,!is.na(d$ifatpftitre) & !is.na(d$ageyrs))

# subset to ages 0-20
d <- subset(d,ageyrs<=20)

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

# SL library
SL.library <- c("SL.mean","SL.glm","SL.Yman2016","SL.gam","SL.loess")


# Pre-intervention period fitted curves
set.seed(23752)
p.c12 <- agecurveAb(Y=log10(d.c12$ifatpftitre+1),Age=d.c12$ageyrs,id=d.c12$id,SL.library=SL.library)
p.tr1 <- agecurveAb(Y=log10(d.tr1$ifatpftitre+1),Age=d.tr1$ageyrs,id=d.tr1$id,SL.library=SL.library)
p.tr2 <- agecurveAb(Y=log10(d.tr2$ifatpftitre+1),Age=d.tr2$ageyrs,id=d.tr2$id,SL.library=SL.library)

# Intervention phase fitted curves
set.seed(2436)
p.c345 <- agecurveAb(Y=log10(d.c345$ifatpftitre+1),Age=d.c345$ageyrs,id=d.c345$id,SL.library=SL.library)
p.tr3 <- agecurveAb(Y=log10(d.tr3$ifatpftitre+1),Age=d.tr3$ageyrs,id=d.tr3$id,SL.library=SL.library)
p.tr4 <- agecurveAb(Y=log10(d.tr4$ifatpftitre+1),Age=d.tr4$ageyrs,id=d.tr4$id,SL.library=SL.library)
p.tr5 <- agecurveAb(Y=log10(d.tr5$ifatpftitre+1),Age=d.tr5$ageyrs,id=d.tr5$id,SL.library=SL.library)

# post intervention period fitted curves
# note: no control measurement in round 6
set.seed(45234)
p.c78 <- agecurveAb(Y=log10(d.c78$ifatpftitre+1),Age=d.c78$ageyrs,id=d.c78$id,SL.library=SL.library)
p.tr6 <- agecurveAb(Y=log10(d.tr6$ifatpftitre+1),Age=d.tr6$ageyrs,id=d.tr6$id,SL.library=SL.library)
p.tr7 <- agecurveAb(Y=log10(d.tr7$ifatpftitre+1),Age=d.tr7$ageyrs,id=d.tr7$id,SL.library=SL.library)
p.tr8 <- agecurveAb(Y=log10(d.tr8$ifatpftitre+1),Age=d.tr8$ageyrs,id=d.tr8$id,SL.library=SL.library)


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
mu_c <- sapply(c(1:5,7:8),function(x) tmleAb(Y=log10(d$ifatpftitre[d$tr=="Control" & d$serosvy==x]+1),
                                             W=data.frame(Age=d$ageyrs[d$tr=="Control" & d$serosvy==x]),
                                             id=d$id[d$tr=="Control" & d$serosvy==x],
                                             SL.library=SL.library)[c("psi","se","lb","ub","p")] )


### Intervention Villages - Spraying (Propoxur) + MDA
mu_i <- sapply(c(1:8),function(x) tmleAb(Y=log10(d$ifatpftitre[d$tr=="Intervention" & d$serosvy==x]+1),
                                             W=data.frame(Age=d$ageyrs[d$tr=="Intervention" & d$serosvy==x]),
                                             id=d$id[d$tr=="Intervention" & d$serosvy==x],
                                             SL.library=SL.library)[c("psi","se","lb","ub","p")] )

#-------------------------------
# Estimate difference between
# control and intervention villages in
# IFA P. falciparm
# titre by survey round
#-------------------------------

set.seed(79287234)
diff_psi <- sapply(c(1:5,7:8),function(x) tmleAb(Y=log10(d$ifatpftitre[d$serosvy==x]+1),
                                                 X=d$tr01[d$serosvy==x],
                                                 W=data.frame(Age=d$ageyrs[d$serosvy==x]),
                                                 id=d$id[d$serosvy==x],
                                                 SL.library=SL.library)[c("psi","se","lb","ub","p")])

# P-values for differences in means, with Bonferroni correction for 7 tests
sprintf("%1.4f",unlist(diff_psi[5,])*7)

#-------------------------------
# save down the results
#-------------------------------
save.image("~/dropbox/articles/antibody-curves/results/raw/garki-main-analysis.RData")



