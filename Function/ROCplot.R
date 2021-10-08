# Function for ROC
smooth.binormal <- function(TP,FP, n) {
  df <- data.frame(sp=qnorm(1-FP), se=qnorm(TP))
  df <- df[apply(df, 1, function(x) all(is.finite(x))),]
  # if (dim(df)[1] <= 1) # ROC curve or with only 1 point
  #   stop("ROC curve not smoothable (not enough points).")
  model <- lm(sp~se, df)
  # if(any(is.na(model$coefficients[2])))
  #   stop("ROC curve not smoothable (not enough points).")
  se <- qnorm(seq(0, 1, 1/(n-1)))
  sp <- predict(model, data.frame(se))
  
  return(list(roc_y = pnorm(se),
              roc_x = 1-pnorm(sp)))
}

myPlot<-function(plotDs=plotDs,colorList,smooth=T,printList=seq(2,60,2),textAUC=NULL,main="ROC of 3-year of TCGA",
                 legendLab,legendCol,legendLty){
  
  printListno <- c(1:length(plotDs$times))
  
  for (i in printListno){
    
    FP <- plotDs$FP[,i]
    TP <- plotDs$TP[,i]
    if(smooth){
      smooth_roc <-  smooth.binormal(TP=TP,FP=FP,512)
      roc_x <- smooth_roc$roc_x
      roc_y <- smooth_roc$roc_y
    }else{
      roc_x <- FP
      roc_y <- TP
    }
    
    if(i == 1){
      plot(roc_x,roc_y,
           type="l",xlim=c(0,1),ylim=c(0,1),
           xaxt="n",yaxt="n",
           xlab="",ylab="",
           lty=1,lwd=2.5,col=colorList[i],main = list(main,cex=1.1),cex=1)
    }else{
      par(new=T)
      plot(roc_x,roc_y,
           type="l",xlim=c(0,1),ylim=c(0,1),
           xaxt="n",yaxt="n",
           xlab="",ylab="",
           lty=1,lwd=2.5,col=colorList[i],cex=1)
    }
  }
  
  axis(1,at=seq(0,1,0.2),labels=sprintf("%.1f",seq(0,1,0.2)),cex.axis=1.2)
  mtext("Time-dependent False Positive",1,line=2.2,cex=0.8)
  axis(2,at=seq(0,1,0.2),labels=sprintf("%.1f",seq(0,1,0.2)),cex.axis=1.2)
  mtext("Time-dependent True Positive",2,line=2.2,cex=0.8)
  
  # text(0.7,0.15,textAUC,cex=1.2,col="black")
  abline(0,1,lty=2,lwd=1.9)
  
  legend(x=0.18,y=0.18,seg.len=1.2,legend=legendLab,col=legendCol,lty=legendLty,lwd=2,bty="n",cex=1)
}