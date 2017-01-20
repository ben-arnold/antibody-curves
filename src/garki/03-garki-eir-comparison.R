
#-------------------------------
# 3-garki-eir-comparison.R
#
# Compare IFA titre 
# age-adjusted mean values
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

# drop observations with missing outcome or age information
d <- subset(d,!is.na(d$ifatpftitre) & !is.na(d$ageyrs))

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
# IFAT-Pf
# SuperLearner curve fits
#-------------------------------

# SL library
SL.library <- c("SL.mean","SL.glm","SL.Yman2016","SL.gam","SL.loess")


set.seed(2343242)
ajura.1971.SL <- agecurveAb(
		Y=log10(ad$ifatpftitre[ad$vname=="Ajura" & ad$wetseason==1971]+1),
		Age=ad$ageyrs[ad$vname=="Ajura" & ad$wetseason==1971],
		id=ad$id[ad$vname=="Ajura" & ad$wetseason==1971],
		SL.library=SL.library
	)
ajura.1972.SL <- agecurveAb(
		Y=log10(ad$ifatpftitre[ad$vname=="Ajura" & ad$wetseason==1972]+1),
		Age=ad$ageyrs[ad$vname=="Ajura" & ad$wetseason==1972],
		id=ad$id[ad$vname=="Ajura" & ad$wetseason==1972],
		SL.library=SL.library
	)
ajura.1973.SL <- agecurveAb(
		Y=log10(ad$ifatpftitre[ad$vname=="Ajura" & ad$wetseason==1973]+1),
		Age=ad$ageyrs[ad$vname=="Ajura" & ad$wetseason==1973],
		id=ad$id[ad$vname=="Ajura" & ad$wetseason==1973],
		SL.library=SL.library
	)	

rafin.1971.SL <- agecurveAb(
		Y=log10(ad$ifatpftitre[ad$vname=="Rafin Marke" & ad$wetseason==1971]+1),
		Age=ad$ageyrs[ad$vname=="Rafin Marke" & ad$wetseason==1971],
		id=ad$id[ad$vname=="Rafin Marke" & ad$wetseason==1971],
		SL.library=SL.library
	)
rafin.1972.SL <- agecurveAb(
		Y=log10(ad$ifatpftitre[ad$vname=="Rafin Marke" & ad$wetseason==1972]+1),
		Age=ad$ageyrs[ad$vname=="Rafin Marke" & ad$wetseason==1972],
		id=ad$id[ad$vname=="Rafin Marke" & ad$wetseason==1972],
		SL.library=SL.library
	)
rafin.1973.SL <- agecurveAb(
		Y=log10(ad$ifatpftitre[ad$vname=="Rafin Marke" & ad$wetseason==1973]+1),
		Age=ad$ageyrs[ad$vname=="Rafin Marke" & ad$wetseason==1973],
		id=ad$id[ad$vname=="Rafin Marke" & ad$wetseason==1973],
		SL.library=SL.library
	)
rafin.1974.SL <- agecurveAb(
		Y=log10(ad$ifatpftitre[ad$vname=="Rafin Marke" & ad$wetseason==1974]+1),
		Age=ad$ageyrs[ad$vname=="Rafin Marke" & ad$wetseason==1974],
		id=ad$id[ad$vname=="Rafin Marke" & ad$wetseason==1974],
		SL.library=SL.library
	)
rafin.1975.SL <- agecurveAb(
		Y=log10(ad$ifatpftitre[ad$vname=="Rafin Marke" & ad$wetseason==1975]+1),
		Age=ad$ageyrs[ad$vname=="Rafin Marke" & ad$wetseason==1975],
		id=ad$id[ad$vname=="Rafin Marke" & ad$wetseason==1975],
		SL.library=SL.library
	)

nasak.1971.SL <- agecurveAb(
		Y=log10(ad$ifatpftitre[ad$vname=="Nasakar" & ad$wetseason==1971]+1),
		Age=ad$ageyrs[ad$vname=="Nasakar" & ad$wetseason==1971],
		id=ad$id[ad$vname=="Nasakar" & ad$wetseason==1971],
		SL.library=SL.library
	)
nasak.1972.SL <- agecurveAb(
		Y=log10(ad$ifatpftitre[ad$vname=="Nasakar" & ad$wetseason==1972]+1),
		Age=ad$ageyrs[ad$vname=="Nasakar" & ad$wetseason==1972],
		id=ad$id[ad$vname=="Nasakar" & ad$wetseason==1972],
		SL.library=SL.library
	)
nasak.1973.SL <- agecurveAb(
		Y=log10(ad$ifatpftitre[ad$vname=="Nasakar" & ad$wetseason==1973]+1),
		Age=ad$ageyrs[ad$vname=="Nasakar" & ad$wetseason==1973],
		id=ad$id[ad$vname=="Nasakar" & ad$wetseason==1973],
		SL.library=SL.library
	)
nasak.1974.SL <- agecurveAb(
		Y=log10(ad$ifatpftitre[ad$vname=="Nasakar" & ad$wetseason==1974]+1),
		Age=ad$ageyrs[ad$vname=="Nasakar" & ad$wetseason==1974],
		id=ad$id[ad$vname=="Nasakar" & ad$wetseason==1974],
		SL.library=SL.library
	)
nasak.1975.SL <- agecurveAb(
		Y=log10(ad$ifatpftitre[ad$vname=="Nasakar" & ad$wetseason==1975]+1),
		Age=ad$ageyrs[ad$vname=="Nasakar" & ad$wetseason==1975],
		id=ad$id[ad$vname=="Nasakar" & ad$wetseason==1975],
		SL.library=SL.library
	)

#-------------------------------
# TMLE estimates of age-adjusted
# mean antibody titres in each
# year (wet season only)
#-------------------------------

set.seed(5463452)
ajura.tmle <- sapply(c(1971,1972,1973),function(x) tmleAb(
  Y=log10(ad$ifatpftitre[ad$vname=="Ajura" & ad$wetseason==x]+1),
  W=data.frame(Age=ad$ageyrs[ad$vname=="Ajura" & ad$wetseason==x]),
  id=ad$id[ad$vname=="Ajura" & ad$wetseason==x],
  SL.library=SL.library)[c("psi","se","lb","ub","p")]
  )

rafin.tmle <- sapply(c(1971,1972,1973,1974,1975),function(x) tmleAb(
  Y=log10(ad$ifatpftitre[ad$vname=="Rafin Marke" & ad$wetseason==x]+1),
  W=data.frame(Age=ad$ageyrs[ad$vname=="Rafin Marke" & ad$wetseason==x]),
  id=ad$id[ad$vname=="Rafin Marke" & ad$wetseason==x],
  SL.library=SL.library)[c("psi","se","lb","ub","p")]
)

nasak.tmle <- sapply(c(1971,1972,1973,1974,1975),function(x) tmleAb(
  Y=log10(ad$ifatpftitre[ad$vname=="Nasakar" & ad$wetseason==x]+1),
  W=data.frame(Age=ad$ageyrs[ad$vname=="Nasakar" & ad$wetseason==x]),
  id=ad$id[ad$vname=="Nasakar" & ad$wetseason==x],
  SL.library=SL.library)[c("psi","se","lb","ub","p")]
)

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
save.image("~/dropbox/articles/antibody-curves/results/raw/garki-eir-comparison.RData")



