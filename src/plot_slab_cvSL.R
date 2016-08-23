
#--------------------------------
# generalized plotting routine
# for cross-validated super learner
# predictions of antibody data
#--------------------------------

plot_slab_cvSL <- function(x,col='blue',title=NULL,xlim=NULL) {
  
  # summarize the object
  sumx <- summary(x)
  
  # remove the discrete super learner
  sumx$Table <- subset(sumx$Table,Algorithm!="Discrete SL")
  
  # for some of the commonly used algorithms, clean up the names
  sumx$Table$Algorithm <- as.character(sumx$Table$Algorithm)
  sumx$Table$Algorithm[sumx$Table$Algorithm=="SL.mean_All"] <- "Simple Mean"
  sumx$Table$Algorithm[sumx$Table$Algorithm=="SL.glm_All"] <- "GLM"
  sumx$Table$Algorithm[sumx$Table$Algorithm=="SL.bayesglm_All"] <- "Bayes GLM"
  sumx$Table$Algorithm[sumx$Table$Algorithm=="SL.gam_All"] <- "GAM"
  sumx$Table$Algorithm[sumx$Table$Algorithm=="SL.loess_All"] <- "LOESS"
  sumx$Table$Algorithm[sumx$Table$Algorithm=="SL.polymars_All"] <- "MARS"
  sumx$Table$Algorithm[sumx$Table$Algorithm=="SL.glmnet_All"] <- "Lasso (glmnet)"
  sumx$Table$Algorithm[sumx$Table$Algorithm=="SL.Yman2016_All"] <- "Yman2016"
  sumx$Table$Algorithm[grep("SL.randomForest",sumx$Table$Algorithm)] <- "Random Forest"
  sumx$Table$Algorithm <- as.factor(sumx$Table$Algorithm)
  
  # sort algorithm by CV risk
  # always put Super Learner at the top
  sumx$Table$Algorithm <- reorder(sumx$Table$Algorithm, -sumx$Table$Ave)
  if(levels(sumx$Table$Algorithm)[8]!="Super Learner") {
    levels(sumx$Table$Algorithm) <- c(levels(sumx$Table$Algorithm)[-grep("Super Learner",levels(sumx$Table$Algorithm))],"Super Learner")
    
  }
  
  
  # plot
  Mean <- sumx$Table$Ave
  se <- sumx$Table$se
  Lower <- Mean - qnorm(0.975) * se
  Upper <- Mean + qnorm(0.975) * se
  assign("d", data.frame(Y = Mean, X = sumx$Table$Algorithm, 
                         Lower = Lower, Upper = Upper))
  
  # .SLAb.require("ggplot2")
  p <- ggplot(d, aes_string(x = "X", y = "Y", ymin = "Lower",  ymax = "Upper")) + 
    geom_linerange(color=col,alpha=0.5) + geom_point(color=col,size=4,alpha=0.5) + coord_flip() + 
    ylab("V-fold CV Risk Estimate") + xlab("Method") + 
    theme_bw() +
    ggtitle(title)
  
  p
  
  return(p)
  
}