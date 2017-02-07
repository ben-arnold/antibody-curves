
#-------------------------------
# FigS2-3-garki-village-EY-by-round-seroprev.R
#
# Calculate age-specific antibody
# curves and age-adjusted mean
# Pf IFA seroprevalence by
# village and survey round
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
d <- garki_sero

d$mdate <- as.Date(d$mdate,"%d %b %Y")



# add 2 control village names
d$vname <- factor(d$vname,levels=c(levels(d$vname),"Nabanawa","Ajura"))
d$vname[d$village==552] <- "Nabanawa"
d$vname[d$village==553] <- "Ajura"
d$vname <- factor(d$vname)

# identify sero-positive individuals (Table 21 of Molineaux 1980, p 175)
d$ifapfpos <- ifelse(d$ifatpftitre>=20,1,0)

# subset to ages 0-20
d <- subset(d,ageyrs<=20)

# drop observations without outcome or age measurements
d <- subset(d,!is.na(d$ifapfpos) & !is.na(d$ageyrs))

#-------------------------------
# Estimate antibody curves
# by village and survey round
# group 2 control villages into
# the same curve for visual
# presentation
#-------------------------------


# SL library
SL.library <- c("SL.mean","SL.glm","SL.gam","SL.loess")

set.seed(5463452)


# small wrapper function to call
# agecurveAb() by village
# and by survey round
SLAb.wrap <- function(svy,vil,d) {
  # svy : survey round
  # vil : village number
  # d   : dataset
  return( agecurveAb(Y=d$ifapfpos[d$village==vil & d$serosvy==svy],Age=d$ageyrs[d$village==vil & d$serosvy==svy],id=d$id[d$village==vil & d$serosvy==svy],SL.library=SL.library) )
}

# Control Villages 
# Nabanawa + Ajura
# (not using wrapper fn above because need to combine across svy rounds)
c12.EYxa <- agecurveAb(
  Y=d$ifapfpos[d$tr=="Control" & d$serosvy>=1 & d$serosvy<=2],
  Age=d$ageyrs[d$tr=="Control" & d$serosvy>=1 & d$serosvy<=2],
  id=d$id[d$tr=="Control" & d$serosvy>=1 & d$serosvy<=2],
  SL.library=SL.library
  )
c345.EYxa <- agecurveAb(
  Y=d$ifapfpos[d$tr=="Control" & d$serosvy>=3 & d$serosvy<=5],
  Age=d$ageyrs[d$tr=="Control" & d$serosvy>=3 & d$serosvy<=5],
  id=d$id[d$tr=="Control" & d$serosvy>=3 & d$serosvy<=5],
  SL.library=SL.library
)
c678.EYxa <- agecurveAb(
  Y=d$ifapfpos[d$tr=="Control" & d$serosvy>=6 & d$serosvy<=8],
  Age=d$ageyrs[d$tr=="Control" & d$serosvy>=6 & d$serosvy<=8],
  id=d$id[d$tr=="Control" & d$serosvy>=6 & d$serosvy<=8],
  SL.library=SL.library
)

# village cluster 5
# Kawari
v153.EYxa <- sapply(1:8,SLAb.wrap,vil=153,d=d,simplify=FALSE)

# Rafin Marke
v154.EYxa <- sapply(1:8,SLAb.wrap,vil=154,d=d,simplify=FALSE)

# Kukar Maikiva
v155.EYxa <- sapply(1:8,SLAb.wrap,vil=155,d=d,simplify=FALSE)

# village cluster 7
# Kargo Kudu
v213.EYxa <- sapply(1:8,SLAb.wrap,vil=213,d=d,simplify=FALSE)

# Nasakar
v218.EYxa <- sapply(1:8,SLAb.wrap,vil=218,d=d,simplify=FALSE)

# Bakan Sabara
v220.EYxa <- sapply(1:8,SLAb.wrap,vil=220,d=d,simplify=FALSE)



#-------------------------------
# Estimate mean IFA P. falciparm
# seroprevalence by village and survey round
#-------------------------------


# since there are many, many cases
# of villages with 100% seropositive
# for a given survey round, it is
# impossible to estimate a TMLE.
# instead,  estimate a mean
# with exact binomial confidence interval

seroprev_est <- function(x) {
  # x : a binary indicator of seropositive(1) vs. seronegative(0)
  tabx <- table(x)
  if(length(tabx)<2) {
    tabx <- c(0,tabx)
  }
  estx <- binom.test(x=tabx[2],n=sum(tabx))
  res <- c(estx$estimate,estx$conf.int)
  names(res) <- c('mean','min95','max95')
  return(res)
}


#-------------------------------
# Serological survey timing:
# 1-2 : pre-intervention
# 3-5 : intervention period
# 6-8 : post-intervention
#-------------------------------

## Control Villages (no measurement in round 6)
# Nabanawa + Ajura

v552 <-  sapply(c(1:5,7:8),function(x) seroprev_est(d$ifapfpos[d$tr=="Control" & d$serosvy==x]) )
	# add NA column for survey round 6
	v552 <- cbind(v552[,1:5],rep(NA,3),v552[,6:7])



### Spraying (Propoxur) + MDA

# village cluster 5
# Kawari
v153 <- sapply(c(1:8),function(x) seroprev_est(d$ifapfpos[d$village==153 & d$serosvy==x]) )

# Rafin Marke
v154 <- sapply(c(1:8),function(x) seroprev_est(d$ifapfpos[d$village==154 & d$serosvy==x]) )

# Kukar Maikiva
v155 <- sapply(c(1:8),function(x) seroprev_est(d$ifapfpos[d$village==155 & d$serosvy==x]) )


# village cluster 7
# Kargo Kudu
v213 <- sapply(c(1:8),function(x) seroprev_est(d$ifapfpos[d$village==213 & d$serosvy==x]) )

# Nasakar
v218 <- sapply(c(1:8),function(x) seroprev_est(d$ifapfpos[d$village==218 & d$serosvy==x]) )

# Bakan Sabara
v220 <- sapply(c(1:8),function(x) seroprev_est(d$ifapfpos[d$village==220 & d$serosvy==x]) )


#-------------------------------
# save the analysis output
#-------------------------------
rm(d)
save.image("~/dropbox/articles/antibody-curves/results/raw/garki-village-EY-by-round-seroprev.RData")


