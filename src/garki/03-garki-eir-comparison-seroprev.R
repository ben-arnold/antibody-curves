
#-------------------------------
# 3b-garki-eir-comparison-seroprev.R
#
# Compare IFA seroprevalence
# with EIR estimates from
# the original Garki Project
# (Molineaux 1980, Table 4, p 65)
# for the 3 villages with both
# serological and entomological
# analyses at multiple wet season
# time points
#
#-------------------------------



#-------------------------------
# preamble
#-------------------------------

rm(list=ls())
library(RColorBrewer)
library(scales)
library(SuperLearner)
library(tmle)
library(tmleAb)


# bright color blind palette:  https://personal.sron.nl/~pault/ 
cblack <- "#000004FF"
cblue <- "#3366AA"
cteal <- "#11AA99"
cgreen <- "#66AA55"
cchartr <- "#CCCC55"
cmagent <- "#992288"
cred <- "#EE3333"
corange <- "#EEA722"
cyellow <- "#FFEE33"
cgrey <- "#777777"


#-------------------------------
# input EIR values from 
# Table 4 (p 65) of Molineaux 1980
# for villages:
# Ajura(553), Rafin-Marke(154), Nasakar(218)
#-------------------------------
eir <- c( 37,25,28,NA,NA,
		  18, 5, 2, 4, 6,
		 129, 0, 4,16,24)
village <- rep(c(553,154,218),c(5,5,5))
vname <- rep(c("Ajura","Rafin Marke","Nasakar"),c(5,5,5))
eirdates <- rep(c("1971-06-21 to 1971-11-07","1972-05-22 to 1972-10-22","1973-06-18 to 1973-11-04","1974-07-29 to 1974-12-15","1975-07-14 to 1975-11-30"),3)
wetseason <- rep(1971:1975,3)
deir <- data.frame(village,vname,wetseason,eirdates,eir)


#-------------------------------
# load the serology dataset
#-------------------------------
# d <- read.csv("~/dropbox/articles/antibody-curves/data/garki/final/garki-sero.csv")
data("garki_sero")
d <- garki_sero

d$mdate <- as.Date(d$mdate,"%d %b %Y")

# identify sero-positive individuals (Table 21 of Molineaux 1980, p 175)
d$ifapfpos <- ifelse(d$ifatpftitre>=20,1,0)

# subset to ages 0-20
# d <- subset(d,ageyrs<=5)

# drop observations with missing outcome or age information
d <- subset(d,!is.na(d$ifapfpos) & !is.na(d$ageyrs))

#-------------------------------
# subset the dataset to
# the 3 villages with EIR ests
#-------------------------------
ad <- subset(d,village==553|village==154|village==218)
ad$vname <- factor(ad$vname)
ad$vname <- factor(ad$vname,levels=c("Ajura","Nasakar","Rafin Marke"))
ad$vname[ad$village==553] <- "Ajura"

#-------------------------------
# Identify observations in the
# wet seasons
#-------------------------------
ad$wetseason <- NA
ad$wetseason[ad$mdate>="1971-06-21" & ad$mdate<="1971-11-07" ] <- 1971
ad$wetseason[ad$mdate>="1972-05-22" & ad$mdate<="1972-10-22" ] <- 1972
ad$wetseason[ad$mdate>="1973-06-18" & ad$mdate<="1973-11-04" ] <- 1973
ad$wetseason[ad$mdate>="1974-07-29" & ad$mdate<="1974-12-15" ] <- 1974
ad$wetseason[ad$mdate>="1975-07-14" & ad$mdate<="1975-11-30" ] <- 1975

table(ad$vname,ad$wetseason)


#-------------------------------
# TMLE estimates of 
# mean seroprevalence in each
# year (wet season only)
#-------------------------------

# SL library
SL.library <- c("SL.mean","SL.glm","SL.gam","SL.loess")


set.seed(5463452)
# 1971 and 1972 are 100% positive, so TMLE estimation impossible
# c(1971,1972,1973)
ajura.tmle <- sapply(c(1973),function(x) tmleAb(
  Y=ad$ifapfpos[ad$vname=="Ajura" & ad$wetseason==x],
  W=data.frame(Age=ad$ageyrs[ad$vname=="Ajura" & ad$wetseason==x]),
  id=ad$id[ad$vname=="Ajura" & ad$wetseason==x],
  SL.library=SL.library)[c("psi","se","lb","ub","p")]
  )

rafin.tmle <- sapply(c(1971,1972,1973,1974,1975),function(x) tmleAb(
  Y=ad$ifapfpos[ad$vname=="Rafin Marke" & ad$wetseason==x],
  W=data.frame(Age=ad$ageyrs[ad$vname=="Rafin Marke" & ad$wetseason==x]),
  id=ad$id[ad$vname=="Rafin Marke" & ad$wetseason==x],
  SL.library=SL.library)[c("psi","se","lb","ub","p")]
)

# 1971 is 100% positive, so TMLE estimation impossible
# c(1971,1972,1973,1974,1975)
nasak.tmle <- sapply(c(1972,1973,1974,1975),function(x) tmleAb(
  Y=ad$ifapfpos[ad$vname=="Nasakar" & ad$wetseason==x],
  W=data.frame(Age=ad$ageyrs[ad$vname=="Nasakar" & ad$wetseason==x]),
  id=ad$id[ad$vname=="Nasakar" & ad$wetseason==x],
  SL.library=SL.library)[c("psi","se","lb","ub","p")]
)

#-------------------------------
# add back in the 100% means
#-------------------------------
ajura.tmle <- cbind(c(1,0,1,1,0),c(1,0,1,1,0),ajura.tmle)
nasak.tmle <- cbind(c(1,0,1,1,0),nasak.tmle)
#-------------------------------
# summarize the means + 95% CIs
#-------------------------------
ajura.mus <- matrix(unlist(ajura.tmle[c(1,3,4),]),nrow=3,ncol=3,byrow=T)
rafin.mus <- matrix(unlist(rafin.tmle[c(1,3,4),]),nrow=5,ncol=3,byrow=T)
nasak.mus <- matrix(unlist(nasak.tmle[c(1,3,4),]),nrow=5,ncol=3,byrow=T)
colnames(ajura.mus) <- colnames(rafin.mus) <- colnames(nasak.mus) <- c("mu","lb","ub")
all.mus <- rbind(ajura.mus,rep(NA,3),rep(NA,3),rafin.mus,nasak.mus)
dmus <- data.frame(
		vname=rep(c("Ajura","Rafin Marke","Nasakar"),c(5,5,5)),
		wetseason=rep(1971:1975,3),
		all.mus
)

#-------------------------------
# merge the means to the EIR data
#-------------------------------
md <- merge(deir,dmus,by=c("vname","wetseason"))
md <- subset(md,!is.na(eir))

# calculate log_10 eir.  add 1 to eir=0
md$log10eir <- log10(md$eir)
md$log10eir[md$eir==0] <- log10(1)

#-------------------------------
# calculate spearman's rank 
# correlation test statistic
#-------------------------------

sp.rho <- cor.test(md$log10eir,md$mu,method="spearman")
sp.rho


#-------------------------------
# save the output
#-------------------------------
rm(d)
save.image("~/dropbox/articles/antibody-curves/results/raw/garki-eir-comparison-seroprev.RData")



