

#-------------------------------------------
# Fig1-1-mauke-Wb123-analysis.R
# Ben Arnold
#
# SuperLearner fits of age-specific antibody
# response on Mauke island before and after
# MDA, and TMLE estimates of geometric mean
# (area under the curve) of antibody responses
#
#-------------------------------------------

#-------------------------------------------
# input files:
#   mauke1975-public.csv
#   mauke1992-public.csv
#
# output files:
#   mauke-Wb123-analysis.RData
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

# d75 <- read.csv("~/dropbox/articles/antibody-curves/data/mauke/mauke1975-public.csv")
# d92 <- read.csv("~/dropbox/articles/antibody-curves/data/mauke/mauke1992-public.csv")
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
common.vars <- c("id","ager","agecat","wb123","year","mda")
d <- subset(d,select=common.vars)


# subset to each year as well, for convenience
d75 <- subset(d,year=="1975")
d92 <- subset(d,year=="1992")


#--------------------------------------
# All Ages
#--------------------------------------
# SL library
SL.library <- c("SL.mean","SL.glm","SL.Yman2016","SL.gam","SL.loess")

set.seed(0237234)

# SuperLearner fits of antibody levels
mauke75 <- agecurveAb(Y=log10(d75$wb123),Age=d75$ager,id=d75$id,SL.library=SL.library)
mauke92 <- agecurveAb(Y=log10(d92$wb123),Age=d92$ager,id=d92$id,SL.library=SL.library)

# estimate group means
EYx.mauke75 <- tmleAb(Y=log10(d75$wb123),W=data.frame(Age=d75$ager),id=d75$id,SL.library=SL.library)
EYx.mauke92 <- tmleAb(Y=log10(d92$wb123),W=data.frame(Age=d92$ager),id=d92$id,SL.library=SL.library)

# estimate difference in means
diff.mauke  <- tmleAb(Y=log10(d$wb123),X=d$mda,W=data.frame(Age=d$ager),id=d$id,SL.library=SL.library)


#--------------------------------------
# Estimate means and differences between
# time points in
# 5 year age bands from ages 1-20
#--------------------------------------

EYx75_1 <- tmleAb(Y=log10(d75$wb123[d75$agecat=="1-5"]),
                  W=data.frame(Age=d75$ager[d75$agecat=="1-5"]),
                  id=d75$id[d75$agecat=="1-5"],
                  SL.library=SL.library) 
EYx75_2 <- tmleAb(Y=log10(d75$wb123[d75$agecat=="6-10"]),
                  W=data.frame(Age=d75$ager[d75$agecat=="6-10"]),
                  id=d75$id[d75$agecat=="6-10"],
                  SL.library=SL.library) 
EYx75_3 <- tmleAb(Y=log10(d75$wb123[d75$agecat=="11-15"]),
                  W=data.frame(Age=d75$ager[d75$agecat=="11-15"]),
                  id=d75$id[d75$agecat=="11-15"],
                  SL.library=SL.library) 
EYx75_4 <- tmleAb(Y=log10(d75$wb123[d75$agecat=="16-20"]),
                  W=data.frame(Age=d75$ager[d75$agecat=="16-20"]),
                  id=d75$id[d75$agecat=="16-20"],
                  SL.library=SL.library) 

EYx92_1 <- tmleAb(Y=log10(d92$wb123[d92$agecat=="1-5"]),
                  W=data.frame(Age=d92$ager[d92$agecat=="1-5"]),
                  id=d92$id[d92$agecat=="1-5"],
                  SL.library=SL.library) 
EYx92_2 <- tmleAb(Y=log10(d92$wb123[d92$agecat=="6-10"]),
                  W=data.frame(Age=d92$ager[d92$agecat=="6-10"]),
                  id=d92$id[d92$agecat=="6-10"],
                  SL.library=SL.library) 
EYx92_3 <- tmleAb(Y=log10(d92$wb123[d92$agecat=="11-15"]),
                  W=data.frame(Age=d92$ager[d92$agecat=="11-15"]),
                  id=d92$id[d92$agecat=="11-15"],
                  SL.library=SL.library) 
EYx92_4 <- tmleAb(Y=log10(d92$wb123[d92$agecat=="16-20"]),
                  W=data.frame(Age=d92$ager[d92$agecat=="16-20"]),
                  id=d92$id[d92$agecat=="16-20"],
                  SL.library=SL.library) 

diff_1 <- tmleAb(Y=log10(d$wb123[d$agecat=="1-5"]),W=data.frame(Age=d$ager[d$agecat=="1-5"]),
                 id=d$id[d$agecat=="1-5"],
                 SL.library=SL.library) 
diff_2 <- tmleAb(Y=log10(d$wb123[d$agecat=="6-10"]),
                 W=data.frame(Age=d$ager[d$agecat=="6-10"]),
                 id=d$id[d$agecat=="6-10"],
                 SL.library=SL.library) 
diff_3 <- tmleAb(Y=log10(d$wb123[d$agecat=="11-15"]),
                 W=data.frame(Age=d$ager[d$agecat=="11-15"]),
                 id=d$id[d$agecat=="11-15"],
                 SL.library=SL.library) 
diff_4 <- tmleAb(Y=log10(d$wb123[d$agecat=="16-20"]),
                 W=data.frame(Age=d$ager[d$agecat=="16-20"]),
                 id=d$id[d$agecat=="16-20"],
                 SL.library=SL.library) 

# condense results
getpsi <- function(x) {
  res <- c(x$psi,x$se,x$lb,x$ub)
  names(res) <- c("psi","se","lb","ub")
  return(res)
}
EYx75 <- list(EYx75_1,EYx75_2,EYx75_3,EYx75_4)
EYx92 <- list(EYx92_1,EYx92_2,EYx92_3,EYx92_4)
diffs <- list(diff_1,diff_2,diff_3,diff_4)
EYx_mauke75kids <- sapply(EYx75, getpsi)
EYx_mauke92kids <- sapply(EYx92, getpsi)
diff_maukekids <- sapply(diffs, getpsi)
colnames(EYx_mauke75kids) <- colnames(EYx_mauke92kids) <- colnames(diff_maukekids) <- c("1-5","6-10","11-15","16-20")

#--------------------------------------
# store results for later summary
# and plotting
#--------------------------------------
save.image("~/dropbox/articles/antibody-curves/results/raw/mauke-Wb123-analysis.RData")




