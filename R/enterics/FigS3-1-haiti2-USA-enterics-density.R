

#-------------------------------------------
# FigS3-1-haiti2-USA-enterics-density.R
#
# estimate the kernel density of the enterics
# MFI values in the USA and Haiti
#
# estimate seropositivity cutoffs using
# finite gaussian mixture models
#
#-------------------------------------------

#-------------------------------------------
# input files:
#   usa_enterics.RData
#   haiti_enterics.RData
#
# output files:
#   haiti2-USA-enterics-density.pdf
#-------------------------------------------



#-------------------------------------------
# preamble
#-------------------------------------------

rm(list=ls())
library(tmleAb)
library(scales)
library(mixtools)


#-------------------------------------------
# load the Mauke data from 1974(1975) and 1992
#-------------------------------------------
data("usa_enterics")
data("haiti_enterics")
d.usa <- usa_enterics
d.haiti <- haiti_enterics


# Limit the haiti data to children
# <= 5.5 years old, which is the
# maximum age for the children
# in the USA sample
d.hai  <- subset(d.haiti,agey<=5.5)

# add 1 to any norovirus and giardia VSP-5
# values that are <=0
d.hai$norogii[d.hai$norogii<=0] <- 1
d.hai$vsp5[d.hai$vsp5<=0] <- 1

# append USA and haiti datasets
d.usa$agey <- d.usa$age
d.usa$cntry <- "USA"
d.hai$cntry <- "Haiti"
commonvar <- c("cntry","id","agey","cp17","cp23","vsp5","leca","etec","salb","norogi","norogii")
d <- rbind(d.usa[commonvar],d.hai[commonvar])


#--------------------------------
# use a gaussian mixture model
# to define cutoffs
# combined (USA + Haiti)
# and separately for each population
#--------------------------------

gmmcut <- function(x) {
  mmx <- normalmixEM(x,lambda=0.5,k=2)
  summary(mmx)
  mini <- mmx$mu<=min(mmx$mu)
  cutx <- mmx$mu[mini]+3*mmx$sigma[mini]
  cat("\nCutpoint (mean + 3SD):",cutx,"\n")
  return(cutx)
}

# cp17
cp17cut     <- gmmcut(log10(d$cp17+1))
cp17cut_hai <- gmmcut(log10(d.hai$cp17+1))
cp17cut_usa <- gmmcut(log10(d.usa$cp17+1))
rbind(cp17cut,cp17cut_hai,cp17cut_usa)

# cp23
cp23cut     <- gmmcut(log10(d$cp23+1))
cp23cut_hai <- gmmcut(log10(d.hai$cp23+1))
cp23cut_usa <- gmmcut(log10(d.usa$cp23+1))
rbind(cp23cut,cp23cut_hai,cp23cut_usa)

# vsp5
vsp5cut     <- gmmcut(log10(d$vsp5+1))
vsp5cut_hai <- gmmcut(log10(d.hai$vsp5+1))
vsp5cut_usa <- gmmcut(log10(d.usa$vsp5+1))
rbind(vsp5cut,vsp5cut_hai,vsp5cut_usa)

# leca
lecacut     <- gmmcut(log10(d$leca+1))
lecacut_hai <- gmmcut(log10(d.hai$leca+1))
lecacut_usa <- gmmcut(log10(d.usa$leca+1))
rbind(lecacut,lecacut_hai,lecacut_usa)

# etec
eteccut     <- gmmcut(log10(d$etec+1))
eteccut_hai <- gmmcut(log10(d.hai$etec+1))
eteccut_usa <- gmmcut(log10(d.usa$etec+1))
rbind(eteccut,eteccut_hai,eteccut_usa)

# salb
salbcut     <- gmmcut(log10(d$salb+1))
salbcut_hai <- gmmcut(log10(d.hai$salb+1))
salbcut_usa <- gmmcut(log10(d.usa$salb+1))
rbind(salbcut,salbcut_hai,salbcut_usa)

# norogi
norogicut     <- gmmcut(log10(d$norogi+1))
norogicut_hai <- gmmcut(log10(d.hai$norogi+1))
norogicut_usa <- gmmcut(log10(d.usa$norogi+1))
rbind(norogicut,norogicut_hai,norogicut_usa)


# norogii
norogiicut     <- gmmcut(log10(d$norogii+1))
norogiicut_hai <- gmmcut(log10(d.hai$norogii+1))
norogiicut_usa <- gmmcut(log10(d.usa$norogii+1))
rbind(norogiicut,norogiicut_hai,norogiicut_usa)




#-------------------------------------------
# density distribution plot function
# repeated for each Ab
#-------------------------------------------

uhden <- function(usaAb,haiAb,allcut,usacut,haicut,letter,main){
  
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  cols <- cbPalette[c(6,7)]
  
  denu <- density(log10(usaAb))
  denh <- density(log10(haiAb))
  xtics <- 0:5
  plot(denu,type="n",
       main="",
       ylim=c(-0.02,1),
       xlim=range(xtics),xaxt="n",xlab="",
       las=1,bty="n"
  )
  polygon(denu,col=alpha(cols[1],alpha=0.3),border=cols[1])
  polygon(denh,col=alpha(cols[2],alpha=0.3),border=cols[2])
  
  axis(1,at=0:5,labels=c(
    expression(10^0),
    expression(10^1),
    expression(10^2),
    expression(10^3),
    expression(10^4),
    expression(10^5)
  ), las=1,cex.axis=1
  )
  
  # seropositivity cutpoints
  points(allcut,-0.02,pch=24,cex=2,col="black",bg="white")
  points(usacut,-0.02,pch=24,cex=2,col=cols[1],bg=alpha(cols[1],alpha=0.3))
  points(haicut,-0.02,pch=24,cex=2,,col=cols[2],bg=alpha(cols[2],alpha=0.3))
#   segments(x0=allcut,y0=0,y1=1,col="black",lty=2,lwd=1.5)
#   segments(x0=usacut,y0=0,y1=1,col=cols[1],lty=2,lwd=1.5)
#   segments(x0=haicut,y0=0,y1=1,col=cols[2],lty=2,lwd=1.5)
  
  # labels
  mtext(main,side=3,line=1,cex=1.25)
  mtext(letter,side=3,line=1,adj=0,at=-0.5,cex=1.5,font=2)
  mtext("Luminex Response (MFI-bg)",side=1,line=2.5,cex=1)
  text(0.5,0.8,"USA",col=cols[1],cex=1)
  text(4,0.5,"Haiti",col=cols[2],cex=1)
  
}


pdf("~/dropbox/articles/antibody-curves/results/figs/haiti2-USA-enterics-density.pdf",height=10,width=20)

# op <- par(xpt=T)
lo <- layout(mat=matrix(1:8,nrow=2,ncol=4,byrow=TRUE))

# cp17
uhden(d.usa$cp17,d.hai$cp17,cp17cut,cp17cut_usa,cp17cut_hai,"a",main=expression(paste(italic('Cryptosporidium parvum'), " Cp17")))

# cp23
uhden(d.usa$cp23,d.hai$cp23,cp23cut,cp23cut_usa,cp23cut_hai,"b",main=expression(paste(italic('Cryptosporidium parvum'), " Cp23")))

# vsp5
uhden(d.usa$vsp5,d.hai$vsp5,vsp5cut,vsp5cut_usa,vsp5cut_hai,"c",main=expression(paste(italic('Giardia intestinalis'), " VSP-5")))

# leca
uhden(d.usa$leca,d.hai$leca,lecacut,lecacut_usa,lecacut_hai,"d",main=expression(paste(italic('Entamoeba histolytica'), " LecA")))

# etec
uhden(d.usa$etec,d.hai$etec,eteccut,eteccut_usa,eteccut_hai,"e",main=expression(paste("ETEC heat labile toxin ",beta," subunit")))

# salb
uhden(d.usa$salb,d.hai$salb,salbcut,salbcut_usa,salbcut_hai,"f",main=expression(paste(italic('Salmonella sp.'), " LPS Group B")))

# norogi
uhden(d.usa$norogi,d.hai$norogi,norogicut,norogicut_usa,norogicut_hai,"g",main="Norovirus GI.4")

# norogii
uhden(d.usa$norogii,d.hai$norogii,norogiicut,norogiicut_usa,norogiicut_hai,"h",main="Norovirus GII.4 NO")

# par(op)
dev.off()






