

#-------------------------------------------
# haiti2-enterics-density
#
# estimate the kernel density of the enterics
# MFI values in the USA and Haiti
#
#-------------------------------------------

#-------------------------------------------
# input files:
#   usa_enterics.RData
#   haiti_enterics.RData
#
# output files:
#   xxx
#-------------------------------------------



#-------------------------------------------
# preamble
#-------------------------------------------

rm(list=ls())
library(tmleAb)
library(scales)

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



#-------------------------------------------
# density distribution plot function
# repeated for each Ab
#-------------------------------------------

uhden <- function(usaAb,haiAb,letter,main){
  
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  cols <- cbPalette[c(6,7)]
  
  denu <- density(log10(usaAb))
  denh <- density(log10(haiAb))
  xtics <- 0:5
  plot(denu,type="n",
       main="",
       ylim=c(0,1),
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
  
  # labels
  mtext(main,side=3,line=1,cex=1.25)
  mtext(letter,side=3,line=1,adj=0,at=-0.5,cex=1.5,font=2)
  mtext("Luminex Response (MFI-bg)",side=1,line=2.5,cex=1)
  text(0.5,0.8,"USA",col=cols[1],cex=1)
  text(4,0.5,"Haiti",col=cols[2],cex=1)
  
  
  
}


pdf("~/dropbox/articles/antibody-curves/results/figs/haiti2-USA-enterics-density.pdf",height=10,width=20)

lo <- layout(mat=matrix(1:8,nrow=2,ncol=4,byrow=TRUE))

# cp17
uhden(d.usa$cp17,d.hai$cp17,"a",main=expression(paste(italic('Cryptosporidium parvum'), " Cp17")))

# cp23
uhden(d.usa$cp23,d.hai$cp23,"b",main=expression(paste(italic('Cryptosporidium parvum'), " Cp23")))

# vsp5
uhden(d.usa$vsp5,d.hai$vsp5,"c",main=expression(paste(italic('Giardia intestinalis'), " VSP-5")))

# leca
uhden(d.usa$leca,d.hai$leca,"d",main=expression(paste(italic('Entamoeba histolytica'), " LecA")))

# etec
uhden(d.usa$etec,d.hai$etec,"e",main=expression(paste("ETEC heat labile toxin ",beta," subunit")))

# salb
uhden(d.usa$salb,d.hai$salb,"f",main=expression(paste(italic('Salmonella sp.'), " LPS Group B")))

# norogi
uhden(d.usa$norogi,d.hai$norogi,"g",main="Norovirus GI.4")

# norogii
uhden(d.usa$norogi,d.hai$norogi,"h",main="Norovirus GII.4 NO")


dev.off()






