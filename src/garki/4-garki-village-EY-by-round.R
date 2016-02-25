
#-------------------------------
# 3-garki-village-EY-by-round
#
# Calculate age-specific antibody
# curves and age-adjusted mean
# Pf IFA antibody titres by
# village and survey round
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
# Estimate antibody curves
# by village and survey round
# group 2 control villages into
# the same curve for visual
# presentation
#-------------------------------
set.seed(5463452)


# small wrapper function to call
# SLAb.curve() by village
# and by survey round
SLAb.wrap <- function(svy,vil,d) {
  # svy : survey round
  # vil : village number
  # d   : dataset
  return( SLAb.curve(Y=log10(d$ifatpftitre[d$village==vil & d$serosvy==svy]+1),Age=d$ageyrs[d$village==vil & d$serosvy==svy],id=d$id[d$village==vil & d$serosvy==svy]) )
}

# Control Villages 
# Nabanawa + Ajura
# (not using wrapper fn above because need to combine across svy rounds)
c12.EYxa <- SLAb.curve(
  Y=log10(d$ifatpftitre[d$tr=="Control" & d$serosvy>=1 & d$serosvy<=2]+1),
  Age=d$ageyrs[d$tr=="Control" & d$serosvy>=1 & d$serosvy<=2],
  id=d$id[d$tr=="Control" & d$serosvy>=1 & d$serosvy<=2]
  )
c345.EYxa <- SLAb.curve(
  Y=log10(d$ifatpftitre[d$tr=="Control" & d$serosvy>=3 & d$serosvy<=5]+1),
  Age=d$ageyrs[d$tr=="Control" & d$serosvy>=3 & d$serosvy<=5],
  id=d$id[d$tr=="Control" & d$serosvy>=3 & d$serosvy<=5]
)
c678.EYxa <- SLAb.curve(
  Y=log10(d$ifatpftitre[d$tr=="Control" & d$serosvy>=6 & d$serosvy<=8]+1),
  Age=d$ageyrs[d$tr=="Control" & d$serosvy>=6 & d$serosvy<=8],
  id=d$id[d$tr=="Control" & d$serosvy>=6 & d$serosvy<=8]
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
# titre by village and survey round
#-------------------------------

#-------------------------------
# Serological survey timing:
# 1-2 : pre-intervention
# 3-5 : intervention period
# 6-8 : post-intervention
#-------------------------------


## Control Villages (no measurement in round 6)
# Nabanawa + Ajura

v552 <-  sapply(c(1:5,7:8),function(x) SLAb.tmle(
  Y=log10(d$ifatpftitre[d$tr=="Control" & d$serosvy==x]+1),
  Age=d$ageyrs[d$tr=="Control" & d$serosvy==x],
  id=d$id[d$tr=="Control" & d$serosvy==x]
  )
)
	# add NA column for survey round 6
	v552 <- cbind(v552[,1:5],rep(NA,5),v552[,6:7])



### Spraying (Propoxur) + MDA

# village cluster 5
# Kawari
v153 <- sapply(c(1:8),function(x) SLAb.tmle(
    Y=log10(d$ifatpftitre[d$village==153 & d$serosvy==x]+1),
    Age=d$ageyrs[d$village==153 & d$serosvy==x],
    id=d$id[d$village==153 & d$serosvy==x]
    )
  )
# Rafin Marke
v154 <- sapply(c(1:8),function(x) SLAb.tmle(
  Y=log10(d$ifatpftitre[d$village==154 & d$serosvy==x]+1),
  Age=d$ageyrs[d$village==154 & d$serosvy==x],
  id=d$id[d$village==154 & d$serosvy==x]
  )
)
# Kukar Maikiva
v155 <- sapply(c(1:8),function(x) SLAb.tmle(
  Y=log10(d$ifatpftitre[d$village==155 & d$serosvy==x]+1),
  Age=d$ageyrs[d$village==155 & d$serosvy==x],
  id=d$id[d$village==155 & d$serosvy==x]
  )
)

# village cluster 7
# Kargo Kudu
v213 <- sapply(c(1:8),function(x) SLAb.tmle(
  Y=log10(d$ifatpftitre[d$village==213 & d$serosvy==x]+1),
  Age=d$ageyrs[d$village==213 & d$serosvy==x],
  id=d$id[d$village==213 & d$serosvy==x]
  )
)
# Nasakar
v218 <- sapply(c(1:8),function(x) SLAb.tmle(
  Y=log10(d$ifatpftitre[d$village==218 & d$serosvy==x]+1),
  Age=d$ageyrs[d$village==218 & d$serosvy==x],
  id=d$id[d$village==218 & d$serosvy==x]
  )
)
# Bakan Sabara
v220 <- sapply(c(1:8),function(x) SLAb.tmle(
  Y=log10(d$ifatpftitre[d$village==220 & d$serosvy==x]+1),
  Age=d$ageyrs[d$village==220 & d$serosvy==x],
  id=d$id[d$village==220 & d$serosvy==x]
  )
)

#-------------------------------
# save the analysis output
#-------------------------------
rm(d)
save.image("~/dropbox/articles/antibody-curves/results/raw/garki-village-EY-by-round.RData")


