

#-------------------------------------------
# Fig1-3-mauke-Wb123-long.R
# Ben Arnold
#
# TMLE means of WB123 antibody in Mauke
# among individuals measured at both time
# points, stratified by age
#
#-------------------------------------------

#-------------------------------------------
# input files:
#   mauke1975-public.csv
#   mauke1992-public.csv
#
# output files:
#   mauke-Wb123-long.RData
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

data("mauke_wb123")
md <- mauke_wb123

# make an MDA identifier (all obs in 1992 are 5y post MDA)
# subset to relevant vars
md$mda <- ifelse(md$year=='1992',1,0)
md <- subset(md,select=c("id75","age","CAg","wb123","mda"))
names(md) <- c("id","age","CAg","wb123","mda")

# create a wide format dataset
d75 <- subset(md,mda==0)
d92 <- subset(md,mda==1)
names(d75) <- c("id","age.75","CAg.75","wb123.75","mda.75")
names(d92) <- c("id","age.92","CAg.92","wb123.92","mda.92")
d <- merge(d75,d92,by="id")

# restrict long format dataset to individuals 
# with measuremments in both 1975 and 1992
a7592 <- merge(md,subset(d,select="id"),by="id",all.x=F)

# create antigen status variables in 1975 and 1992
d$ab75 <- as.factor(ifelse(d$CAg.75>32,"Pos","Neg"))
d$ab92 <- as.factor(ifelse(d$CAg.92>32,"Pos","Neg"))
d$Abstatus <- factor(rep("Neg-Neg",nrow(d)),levels=c("Neg-Neg","Pos-Neg","Pos-Pos","Neg-Pos"))
  d$Abstatus[d$ab75=="Pos"&d$ab92=="Neg"] <- "Pos-Neg"
  d$Abstatus[d$ab75=="Pos"&d$ab92=="Pos"] <- "Pos-Pos"
  d$Abstatus[d$ab75=="Neg"&d$ab92=="Pos"] <- "Neg-Pos"
  
# merge antigen status at both time points into the long format data as well
a7592 <- merge(a7592,subset(d,select=c("id","Abstatus")),by="id",all.x=T)


#--------------------------------------
# Estimate means and differences between
# time points for each antigen+/- 
# combo
#--------------------------------------
SL.library <- c("SL.mean","SL.glm","SL.Yman2016","SL.gam","SL.loess")
set.seed(2194543)
Abstatus <- c("Neg-Neg","Pos-Neg","Pos-Pos")
EYx_75_Abstatus <- sapply(Abstatus, function(x) 
	tmleAb(Y=log10(d$wb123.75[d$Abstatus==x]),
	       W=data.frame(Age=d$age.75[d$Abstatus==x]),
	       id=d$id[d$Abstatus==x],
	       SL.library = SL.library)[c("psi","se","lb","ub","p")]
)
EYx_92_Abstatus <- sapply(Abstatus, function(x) 
	tmleAb(Y=log10(d$wb123.92[d$Abstatus==x]),
	       W=data.frame(Age=d$age.92[d$Abstatus==x]),
	       id=d$id[d$Abstatus==x],
	       SL.library = SL.library)[c("psi","se","lb","ub","p")]
)
diff_Abstatus <- sapply(Abstatus, function(x) 
	tmleAb(Y=log10(a7592$wb123[a7592$Abstatus==x]), 
	       X=a7592$mda[a7592$Abstatus==x],
	       W=data.frame(Age=a7592$age[a7592$Abstatus==x]), 
	       id=a7592$id[a7592$Abstatus==x],
	       SL.library = SL.library)[c("psi","se","lb","ub","p")] 
)


#--------------------------------------
# store results for later summary
# and plotting
#--------------------------------------
save.image("~/dropbox/articles/antibody-curves/results/raw/mauke-Wb123-long.RData")




