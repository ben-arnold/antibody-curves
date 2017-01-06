#-------------------------------
# COR-NTD-garki-figs.R
#
# summarize the garki antibody curve
# results in a figure

#-------------------------------


#-------------------------------
# preamble
#-------------------------------

rm(list=ls())
library(RColorBrewer)
library(scales)
library(tmleAb)
library(SuperLearner)
library(tmle)

#-------------------------------
# load the analysis results
#-------------------------------
load("~/dropbox/articles/antibody-curves/results/raw/garki-main-analysis.RData")

# fit a curve for intervention
# surveys 1 and 2 combined
d.tr12 <- d[d$tr=="Intervention" & (d$serosvy==1 | d$serosvy==2),]
set.seed(23752)
p.tr12 <- agecurveAb(Y=log10(d.tr12$ifatpftitre+1),Age=d.tr12$ageyrs,id=d.tr12$id,SL.library=SL.library)

# calculate means and differences for surveys 1+2 combined
set.seed(23752)
muc12 <- tmleAb(Y=log10(d$ifatpftitre[d$tr=="Control" & (d$serosvy==1|d$serosvy==2)]+1),
                          W=data.frame(Age=d$ageyrs[d$tr=="Control" & (d$serosvy==1|d$serosvy==2)]),
                          id=d$id[d$tr=="Control" & (d$serosvy==1|d$serosvy==2)],
                          SL.library=SL.library)
set.seed(23752)
mui12 <- tmleAb(Y=log10(d$ifatpftitre[d$tr=="Intervention" & (d$serosvy==1|d$serosvy==2)]+1),
                W=data.frame(Age=d$ageyrs[d$tr=="Intervention" & (d$serosvy==1|d$serosvy==2)]),
                id=d$id[d$tr=="Intervention" & (d$serosvy==1|d$serosvy==2)],
                SL.library=SL.library)

set.seed(23752)
diff12 <- tmleAb(Y=log10(d$ifatpftitre[d$serosvy==1|d$serosvy==2]+1),
       X=d$tr01[d$serosvy==1|d$serosvy==2],
       W=data.frame(Age=d$ageyrs[d$serosvy==1|d$serosvy==2]),
       id=d$id[d$serosvy==1|d$serosvy==2],
       SL.library=SL.library)

#-------------------------------
# Plot antibody response 
# curves for the 3 study phases
#-------------------------------


# general plotting parameters
ytics <- seq(0,4,by=1)
xtics <- seq(0,20,by=5)
# i.cols <- c(brewer.pal(9,"YlGnBu")[7:6],brewer.pal(9,"BuPu")[9:7],brewer.pal(9,"YlOrRd")[6:8])
# c.cols <- brewer.pal(8,"Dark2")[8]
# c.cols <- "black"
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# i.cols <- cbPalette[c(6,4,8,6,4,8,6,4)]
i.cols <- rep(cbPalette[6],8)
c.cols <- cbPalette[7]

# brighter color blind palette:  https://personal.sron.nl/~pault/ 
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
c.cols=c(cblack)
i.cols=c(cblue,cblue,cblue,cteal,cgreen)





#-------------------------------
# pre-intervention
#-------------------------------
pdf("~/dropbox/articles/antibody-curves/presentations/cor-ntd-2016/figs/garki-svy1.pdf",height=5,width=8)
lo <- layout(mat=matrix(1:2,nrow=1,ncol=2,byrow=T),widths=c(1,0.4))

# Panel A: Pre-intervention age-antibody curves
op <- par(mar=c(4,5,3,1)+0.1,xpd=T)
plot(p.c12$Age,p.c12$pY,type="l",lwd=2,col=c.cols,
	ylim=c(0,4),ylab="",yaxt="n",
	xlim=c(0,max(xtics)+1),xlab="",xaxt="n",
	bty="n",las=1
	)
	mtext(expression(paste(italic('P. falciparum')," IFA antibody titre")),side=2,line=3,cex=1.25)
	#mtext("a",adj=1,line=3,at=-2,font=2,cex=1.75)
	#mtext("Pre-Intervention Period",adj=0,line=3,at=0,cex=1.5)
	#mtext(expression(paste(italic(E),"(",italic(Y[a][","][x]),")")),side=3,line=0)
	mtext("Age, years",side=1,line=2.5,cex=1.5)
	axis(2,at=0:4,labels=c(
	  expression(10^0),
	  expression(10^1),
	  expression(10^2),
	  expression(10^3),
	  expression(10^4)
	), las=1,cex.axis=1.25
	)
	# segments(x0=min(xtics),x1=max(xtics),y0=ytics,lty=2,col="gray70")
	axis(1,at=xtics,cex.axis=1.5)
	lines(p.tr12$Age,p.tr12$pY,col=i.cols[1],lwd=2)
# 	lines(p.tr1$Age,p.tr1$pY,col=i.cols[1])
# 	lines(p.tr2$Age,p.tr2$pY,col=i.cols[2])
	
	text(20,3.1,"Control",cex=1,adj=1,font=2,col=c.cols)
	text(20,3.9,"Intervention",cex=1,adj=1,font=1,col=i.cols[1])
# 	text(3,3.4,"1",cex=0.75,col=i.cols[1])
# 	text(10,3.6,"2",cex=0.75,col=i.cols[2])
	
# Pre-intervention  mean estimates
op <- par(mar=c(4,0,3,4)+0.1)
plot(1:3,1:3,type="n",
	xlim=c(0.5,3.5),xaxt="n",xlab="",
	ylim=c(0,4),ylab="",yaxt="n",
	las=1,bty="n"
)

# mtext(c("1","2"),side=1,line=1,at=c(1.5,2.5),cex=1,col=i.cols[1:2])
mtext(expression(paste(italic(E),"(",italic(Y[x]),")")),side=3,line=0)
# mtext("Survey",side=1,line=2.5)
# control 
	arrows(x0=2,y0=muc12$lb, y1=muc12$ub,lwd=1,col=c.cols,length=0.05,angle=90,code=3)
	points(2,muc12$psi, pch=21,cex=1.5, lwd=1,bg=c.cols,col=c.cols)
# intervention
	arrows(x0=2,y0=mui12$lb, y1=mui12$ub,lwd=1,col=i.cols[1],length=0.05,angle=90,code=3)
	points(2,mui12$psi, pch=21,bg=i.cols[1],cex=1.5,lwd=1, col=i.cols[1])
	
par(op)
dev.off()

#-------------------------------
# intervention, control only svy 3
#-------------------------------
pdf("~/dropbox/articles/antibody-curves/presentations/cor-ntd-2016/figs/garki-svy2.pdf",height=5,width=8)
lo <- layout(mat=matrix(1:2,nrow=1,ncol=2,byrow=T),widths=c(1,0.4))

# age-antibody curves
op <- par(mar=c(4,5,3,1)+0.1)
plot(p.c345$Age,p.c345$pY,type="l",lwd=2,col=c.cols,
     ylim=c(0,4),ylab="",yaxt="n",
     xlim=c(0,max(xtics)+1),xlab="",xaxt="n",
     bty="n",las=1
)
mtext(expression(paste(italic('P. falciparum')," IFA antibody titre")),side=2,line=3,cex=1.25)
# 	mtext("b",adj=1,line=3,at=-2,font=2,cex=1.75)
# 	mtext("Active Intervention Period",adj=0,line=3,at=0,cex=1.5)
# 	mtext(expression(paste(italic(E),"(",italic(Y[a][","][x]),")")),side=3,line=0)
mtext("Age, years",side=1,line=2.5,cex=1.5)
axis(2,at=0:4,labels=c(
  expression(10^0),
  expression(10^1),
  expression(10^2),
  expression(10^3),
  expression(10^4)
), las=1,cex.axis=1.25
)
axis(1,at=xtics,cex.axis=1.5)
# lines(p.tr3$Age,p.tr3$pY,col=i.cols[3],lwd=2)
# lines(p.tr4$Age,p.tr4$pY,col=i.cols[4],lwd=2)
# lines(p.tr5$Age,p.tr5$pY,col=i.cols[5],lwd=2)

text(20,3.6,"Control",cex=1,adj=1,font=2,col=c.cols)
# text(20,2.5,"Intervention\n(survey 3, 4, 5)",cex=1,adj=1,font=1,col=i.cols[1])
# text(7,2.9,"3",cex=0.75,col=i.cols[3])
# text(9,2.65,"4",cex=0.75,col=i.cols[4])
# text(10.5,2.3,"5",cex=0.75,col=i.cols[5],adj=0)

# Active intervention  mean estimates
op <- par(mar=c(4,0,3,4)+0.1)
plot(1:3,1:3,type="n",
     xlim=c(0.5,3.5),xaxt="n",xlab="",
     ylim=c(0,4),ylab="",yaxt="n",
     las=1,bty="n"
)


par(op)
dev.off()


#-------------------------------
# intervention, survey 3
#-------------------------------
pdf("~/dropbox/articles/antibody-curves/presentations/cor-ntd-2016/figs/garki-svy3.pdf",height=5,width=8)
lo <- layout(mat=matrix(1:2,nrow=1,ncol=2,byrow=T),widths=c(1,0.4))

# age-antibody curves
op <- par(mar=c(4,5,3,1)+0.1)
plot(p.c345$Age,p.c345$pY,type="l",lwd=2,col=c.cols,
     ylim=c(0,4),ylab="",yaxt="n",
     xlim=c(0,max(xtics)+1),xlab="",xaxt="n",
     bty="n",las=1
)
mtext(expression(paste(italic('P. falciparum')," IFA antibody titre")),side=2,line=3,cex=1.25)
# 	mtext("b",adj=1,line=3,at=-2,font=2,cex=1.75)
# 	mtext("Active Intervention Period",adj=0,line=3,at=0,cex=1.5)
# 	mtext(expression(paste(italic(E),"(",italic(Y[a][","][x]),")")),side=3,line=0)
mtext("Age, years",side=1,line=2.5,cex=1.5)
axis(2,at=0:4,labels=c(
  expression(10^0),
  expression(10^1),
  expression(10^2),
  expression(10^3),
  expression(10^4)
), las=1,cex.axis=1.25
)
axis(1,at=xtics,cex.axis=1.5)
lines(p.tr3$Age,p.tr3$pY,col=i.cols[3],lwd=2)
# lines(p.tr4$Age,p.tr4$pY,col=i.cols[4],lwd=2)
# lines(p.tr5$Age,p.tr5$pY,col=i.cols[5],lwd=2)

text(20,3.6,"Control",cex=1,adj=1,font=2,col=c.cols)
# text(20,2.5,"Intervention\n(survey 3, 4, 5)",cex=1,adj=1,font=1,col=i.cols[1])
# text(7,2.9,"3",cex=0.75,col=i.cols[3])
# text(9,2.65,"4",cex=0.75,col=i.cols[4])
# text(10.5,2.3,"5",cex=0.75,col=i.cols[5],adj=0)

# Active intervention  mean estimates
op <- par(mar=c(4,0,3,4)+0.1)
plot(1:3,1:3,type="n",
     xlim=c(0.5,3.5),xaxt="n",xlab="",
     ylim=c(0,4),ylab="",yaxt="n",
     las=1,bty="n"
)

mtext(c(20),side=1,line=1,at=1,cex=1.5,col=i.cols[3])
mtext(expression(paste(italic(E),"(",italic(Y[x]),")")),side=3,line=0)
mtext("Weeks",side=1,line=2.5,cex=1.5)
segments(x0=1,y0=rep(0,1),y1=rep(4,1),col="gray80",lty=2)

# control 
arrows(x0=1,y0=unlist(mu_c[3,3]), y1=unlist(mu_c[4,3]),lwd=1,col=c.cols,length=0.05,angle=90,code=3)
points(1,mu_c[1,3], pch=21,cex=1.5, lwd=1,bg=c.cols,col=c.cols)
# intervention
arrows(x0=1,y0=unlist(mu_i[3,3]), y1=unlist(mu_i[4,3]),lwd=1,col=i.cols[3],length=0.05,angle=90,code=3)
points(1,mu_i[1,3], pch=21,bg=i.cols[3],cex=1.5,lwd=1, col=i.cols[3])
par(op)
dev.off()


#-------------------------------
# intervention, survey 3-4
#-------------------------------
pdf("~/dropbox/articles/antibody-curves/presentations/cor-ntd-2016/figs/garki-svy4.pdf",height=5,width=8)
lo <- layout(mat=matrix(1:2,nrow=1,ncol=2,byrow=T),widths=c(1,0.4))

# age-antibody curves
op <- par(mar=c(4,5,3,1)+0.1)
plot(p.c345$Age,p.c345$pY,type="l",lwd=2,col=c.cols,
     ylim=c(0,4),ylab="",yaxt="n",
     xlim=c(0,max(xtics)+1),xlab="",xaxt="n",
     bty="n",las=1
)
mtext(expression(paste(italic('P. falciparum')," IFA antibody titre")),side=2,line=3,cex=1.25)
# 	mtext("b",adj=1,line=3,at=-2,font=2,cex=1.75)
# 	mtext("Active Intervention Period",adj=0,line=3,at=0,cex=1.5)
# 	mtext(expression(paste(italic(E),"(",italic(Y[a][","][x]),")")),side=3,line=0)
mtext("Age, years",side=1,line=2.5,cex=1.5)
axis(2,at=0:4,labels=c(
  expression(10^0),
  expression(10^1),
  expression(10^2),
  expression(10^3),
  expression(10^4)
), las=1,cex.axis=1.25
)
axis(1,at=xtics,cex.axis=1.5)
lines(p.tr3$Age,p.tr3$pY,col=i.cols[3],lwd=2)
lines(p.tr4$Age,p.tr4$pY,col=i.cols[4],lwd=2)
# lines(p.tr5$Age,p.tr5$pY,col=i.cols[5],lwd=2)

text(20,3.6,"Control",cex=1,adj=1,font=2,col=c.cols)

# Active intervention  mean estimates
op <- par(mar=c(4,0,3,4)+0.1)
plot(1:3,1:3,type="n",
     xlim=c(0.5,3.5),xaxt="n",xlab="",
     ylim=c(0,4),ylab="",yaxt="n",
     las=1,bty="n"
)

mtext(c(20,50),side=1,line=1,at=1:2,cex=1.5,col=i.cols[3:4])
mtext(expression(paste(italic(E),"(",italic(Y[x]),")")),side=3,line=0)
mtext("Weeks",side=1,line=2.5,cex=1.5)
segments(x0=1:2,y0=rep(0,2),y1=rep(4,2),col="gray80",lty=2)

# control 
arrows(x0=1:2,y0=unlist(mu_c[3,3:4]), y1=unlist(mu_c[4,3:4]),lwd=1,col=c.cols,length=0.05,angle=90,code=3)
points(1:2,mu_c[1,3:4], pch=21,cex=1.5, lwd=1,bg=c.cols,col=c.cols)
# intervention
arrows(x0=1:2,y0=unlist(mu_i[3,3:4]), y1=unlist(mu_i[4,3:4]),lwd=1,col=i.cols[3:4],length=0.05,angle=90,code=3)
points(1:2,mu_i[1,3:4], pch=21,bg=i.cols[3:4],cex=1.5,lwd=1, col=i.cols[3:4])
par(op)
dev.off()


#-------------------------------
# intervention, survey 3-5
#-------------------------------
pdf("~/dropbox/articles/antibody-curves/presentations/cor-ntd-2016/figs/garki-svy5.pdf",height=5,width=8)
lo <- layout(mat=matrix(1:2,nrow=1,ncol=2,byrow=T),widths=c(1,0.4))

# age-antibody curves
op <- par(mar=c(4,5,3,1)+0.1)
plot(p.c345$Age,p.c345$pY,type="l",lwd=2,col=c.cols,
	ylim=c(0,4),ylab="",yaxt="n",
	xlim=c(0,max(xtics)+1),xlab="",xaxt="n",
	bty="n",las=1
	)
	mtext(expression(paste(italic('P. falciparum')," IFA antibody titre")),side=2,line=3,cex=1.25)
# 	mtext("b",adj=1,line=3,at=-2,font=2,cex=1.75)
# 	mtext("Active Intervention Period",adj=0,line=3,at=0,cex=1.5)
# 	mtext(expression(paste(italic(E),"(",italic(Y[a][","][x]),")")),side=3,line=0)
	mtext("Age, years",side=1,line=2.5,cex=1.5)
	axis(2,at=0:4,labels=c(
	  expression(10^0),
	  expression(10^1),
	  expression(10^2),
	  expression(10^3),
	  expression(10^4)
	), las=1,cex.axis=1.25
	)
	axis(1,at=xtics,cex.axis=1.5)
	lines(p.tr3$Age,p.tr3$pY,col=i.cols[3],lwd=2)
	lines(p.tr4$Age,p.tr4$pY,col=i.cols[4],lwd=2)
	lines(p.tr5$Age,p.tr5$pY,col=i.cols[5],lwd=2)
	
	text(20,3.6,"Control",cex=1,adj=1,font=2,col=c.cols)

	
# Active intervention  mean estimates
op <- par(mar=c(4,0,3,4)+0.1)
plot(1:3,1:3,type="n",
	xlim=c(0.5,3.5),xaxt="n",xlab="",
	ylim=c(0,4),ylab="",yaxt="n",
	las=1,bty="n"
)

mtext(c(20,50,70),side=1,line=1,at=1:3,cex=1.5,col=i.cols[3:5])
mtext(expression(paste(italic(E),"(",italic(Y[x]),")")),side=3,line=0)
mtext("Weeks",side=1,line=2.5,cex=1.5)
segments(x0=1:3,y0=rep(0,3),y1=rep(4,3),col="gray80",lty=2)

# control 
	arrows(x0=1:3,y0=unlist(mu_c[3,3:5]), y1=unlist(mu_c[4,3:5]),lwd=1,col=c.cols,length=0.05,angle=90,code=3)
	points(1:3,mu_c[1,3:5], pch=21,cex=1.5, lwd=1,bg=c.cols,col=c.cols)
# intervention
	arrows(x0=1:3,y0=unlist(mu_i[3,3:5]), y1=unlist(mu_i[4,3:5]),lwd=1,col=i.cols[3:5],length=0.05,angle=90,code=3)
	points(1:3,mu_i[1,3:5], pch=21,bg=i.cols[3:5],cex=1.5,lwd=1, col=i.cols[3:5])
par(op)
dev.off()



