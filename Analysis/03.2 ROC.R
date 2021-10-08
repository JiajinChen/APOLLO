# rm(list=ls())

source('ROCplot.R')
model_C <- coxph(Surv(OS,Censor) ~ Clin_4covar,data = TCGA)
model_A <- coxph(Surv(OS,Censor) ~ Clin_4covar+All_4covar,data = TCGA)

# TCGA
lp_TCGA_C <- predict(model_C,type="lp",newdata = TCGA)
ROC_TCGA_C <- timeROC(T=TCGA$OS,
                 delta=TCGA$Censor,marker=lp_TCGA_C,
                 cause=1,weighting="marginal",iid=T,
                 times=c(36,60))
lp_TCGA_A <- predict(model_A,type="lp",newdata = TCGA)
ROC_TCGA_A <- timeROC(T=TCGA$OS,
                 delta=TCGA$Censor,marker=lp_TCGA_A,
                 cause=1,weighting="marginal",iid=T,
                 times=c(36,60))

# CGGA-Seq
lp_CGGASeq_C <- predict(model_C,type="lp",newdata = CGGA1)
ROC_CGGASeq_C <- timeROC(T=CGGA1$OS,
                      delta=CGGA1$Censor,marker=lp_CGGASeq_C,
                      cause=1,weighting="marginal",iid=T,
                      times=c(36,60))
lp_CGGASeq_A <- predict(model_A,type="lp",newdata = CGGA1)
ROC_CGGASeq_A <- timeROC(T=CGGA1$OS,
                      delta=CGGA1$Censor,marker=lp_CGGASeq_A,
                      cause=1,weighting="marginal",iid=T,
                      times=c(36,60))

# CGGA-Array
lp_CGGAArray_C <- predict(model_C,type="lp",newdata = CGGA2)
ROC_CGGAArray_C <- timeROC(T=CGGA2$OS,
                         delta=CGGA2$Censor,marker=lp_CGGAArray_C,
                         cause=1,weighting="marginal",iid=T,
                         times=c(36,60))
lp_CGGAArray_A <- predict(model_A,type="lp",newdata = CGGA2)
ROC_CGGAArray_A <- timeROC(T=CGGA2$OS,
                         delta=CGGA2$Censor,marker=lp_CGGAArray_A,
                         cause=1,weighting="marginal",iid=T,
                         times=c(36,60))

# Rembrandt
lp_Rembrandt_C <- predict(model_C,type="lp",newdata = Rembrandt)
ROC_Rembrandt_C <- timeROC(T=Rembrandt$OS,
                         delta=Rembrandt$Censor,marker=lp_Rembrandt_C,
                         cause=1,weighting="marginal",iid=T,
                         times=c(36,60))
lp_Rembrandt_A <- predict(model_A,type="lp",newdata = Rembrandt)
ROC_Rembrandt_A <- timeROC(T=Rembrandt$OS,
                         delta=Rembrandt$Censor,marker=lp_Rembrandt_A,
                         cause=1,weighting="marginal",iid=T,
                         times=c(36,60))

# GSE61374
lp_GSE61374_C <- predict(model_C,type="lp",newdata = GSE61374)
ROC_GSE61374_C <- timeROC(T=GSE61374$OS,
                           delta=GSE61374$Censor,marker=lp_GSE61374_C,
                           cause=1,weighting="marginal",iid=T,
                           times=c(36,60))
lp_GSE61374_A <- predict(model_A,type="lp",newdata = GSE61374)
ROC_GSE61374_A <- timeROC(T=GSE61374$OS,
                           delta=GSE61374$Censor,marker=lp_GSE61374_A,
                           cause=1,weighting="marginal",iid=T,
                           times=c(36,60))

# GSE16011
lp_GSE16011_C <- predict(model_C,type="lp",newdata = GSE16011)
ROC_GSE16011_C <- timeROC(T=GSE16011$OS,
                          delta=GSE16011$Censor,marker=lp_GSE16011_C,
                          cause=1,weighting="marginal",iid=T,
                          times=c(36,60))
lp_GSE16011_A <- predict(model_A,type="lp",newdata = GSE16011)
ROC_GSE16011_A <- timeROC(T=GSE16011$OS,
                          delta=GSE16011$Censor,marker=lp_GSE16011_A,
                          cause=1,weighting="marginal",iid=T,
                          times=c(36,60))




pdf("Fig_ROC_35_Smooth.pdf",width=9,height = 6)
par(pin = c(9,6),mfrow=c(2,3),mar = rep(3, 4))


# TCGA
plotDs=ROC_TCGA_A
legendLab <- c(paste0('3-year: ',sprintf("%0.3f",round(plotDs$AUC[1],3)),
                      ' (',sprintf("%0.3f",round(confint(ROC_TCGA_A)$CI_AUC[1,1]/100,3)),', ',sprintf("%0.3f",round(confint(ROC_TCGA_A)$CI_AUC[1,2]/100,3)),')'),
               paste0('5-year: ',sprintf("%0.3f",round(plotDs$AUC[2],3)),
                      ' (',sprintf("%0.3f",round(confint(ROC_TCGA_A)$CI_AUC[2,1]/100,3)),', ',sprintf("%0.3f",round(confint(ROC_TCGA_A)$CI_AUC[2,2]/100,3)),')'))
myPlot(plotDs=plotDs,
       colorList=c('darkblue','darkred'),
       smooth=T,
       main="ROC for OS prediction of TCGA",
       legendLab=legendLab,
       legendCol=c('darkblue','darkred'),
       legendLty=c(1,1))

# CGGA-Seq
plotDs=ROC_CGGASeq_A
legendLab <- c(paste0('3-year: ',sprintf("%0.3f",round(plotDs$AUC[1],3)),
                      ' (',sprintf("%0.3f",round(confint(ROC_CGGASeq_A)$CI_AUC[1,1]/100,3)),', ',sprintf("%0.3f",round(confint(ROC_CGGASeq_A)$CI_AUC[1,2]/100,3)),')'),
               paste0('5-year: ',sprintf("%0.3f",round(plotDs$AUC[2],3)),
                      ' (',sprintf("%0.3f",round(confint(ROC_CGGASeq_A)$CI_AUC[2,1]/100,3)),', ',sprintf("%0.3f",round(confint(ROC_CGGASeq_A)$CI_AUC[2,2]/100,3)),')'))
myPlot(plotDs=plotDs,
       colorList=c('darkblue','darkred'),
       smooth=T,
       main="ROC for OS prediction of CGGA1",
       legendLab=legendLab,
       legendCol=c('darkblue','darkred'),
       legendLty=c(1,1))

# CGGA-Array
plotDs=ROC_CGGAArray_A
legendLab <- c(paste0('3-year: ',sprintf("%0.3f",round(plotDs$AUC[1],3)),
                      ' (',sprintf("%0.3f",round(confint(ROC_CGGAArray_A)$CI_AUC[1,1]/100,3)),', ',sprintf("%0.3f",round(confint(ROC_CGGAArray_A)$CI_AUC[1,2]/100,3)),')'),
               paste0('5-year: ',sprintf("%0.3f",round(plotDs$AUC[2],3)),
                      ' (',sprintf("%0.3f",round(confint(ROC_CGGAArray_A)$CI_AUC[2,1]/100,3)),', ',sprintf("%0.3f",round(confint(ROC_CGGAArray_A)$CI_AUC[2,2]/100,3)),')'))
myPlot(plotDs=plotDs,
       colorList=c('darkblue','darkred'),
       smooth=T,
       main="ROC for OS prediction of CGGA2",
       legendLab=legendLab,
       legendCol=c('darkblue','darkred'),
       legendLty=c(1,1))

# Rembrandt
plotDs=ROC_Rembrandt_A
legendLab <- c(paste0('3-year: ',sprintf("%0.3f",round(plotDs$AUC[1],3)),
                      ' (',sprintf("%0.3f",round(confint(ROC_Rembrandt_A)$CI_AUC[1,1]/100,3)),', ',sprintf("%0.3f",round(confint(ROC_Rembrandt_A)$CI_AUC[1,2]/100,3)),')'),
               paste0('5-year: ',sprintf("%0.3f",round(plotDs$AUC[2],3)),
                      ' (',sprintf("%0.3f",round(confint(ROC_Rembrandt_A)$CI_AUC[2,1]/100,3)),', ',sprintf("%0.3f",round(confint(ROC_Rembrandt_A)$CI_AUC[2,2]/100,3)),')'))
myPlot(plotDs=plotDs,
       colorList=c('darkblue','darkred'),
       smooth=T,
       main="ROC for OS prediction of Rembrandt",
       legendLab=legendLab,
       legendCol=c('darkblue','darkred'),
       legendLty=c(1,1))

# GSE61374
plotDs=ROC_GSE61374_A
legendLab <- c(paste0('3-year: ',sprintf("%0.3f",round(plotDs$AUC[1],3)),
                      ' (',sprintf("%0.3f",round(confint(ROC_GSE61374_A)$CI_AUC[1,1]/100,3)),', ',sprintf("%0.3f",round(confint(ROC_GSE61374_A)$CI_AUC[1,2]/100,3)),')'),
               paste0('5-year: ',sprintf("%0.3f",round(plotDs$AUC[2],3)),
                      ' (',sprintf("%0.3f",round(confint(ROC_GSE61374_A)$CI_AUC[2,1]/100,3)),', ',sprintf("%0.3f",round(confint(ROC_GSE61374_A)$CI_AUC[2,2]/100,3)),')'))
myPlot(plotDs=plotDs,
       colorList=c('darkblue','darkred'),
       smooth=T,
       main="ROC for OS prediction of Weller",
       legendLab=legendLab,
       legendCol=c('darkblue','darkred'),
       legendLty=c(1,1))

# GSE16011
plotDs=ROC_GSE16011_A
legendLab <- c(paste0('3-year: ',sprintf("%0.3f",round(plotDs$AUC[1],3)),
                      ' (',sprintf("%0.3f",round(confint(ROC_GSE16011_A)$CI_AUC[1,1]/100,3)),', ',sprintf("%0.3f",round(confint(ROC_GSE16011_A)$CI_AUC[1,2]/100,3)),')'),
               paste0('5-year: ',sprintf("%0.3f",round(plotDs$AUC[2],3)),
                      ' (',sprintf("%0.3f",round(confint(ROC_GSE16011_A)$CI_AUC[2,1]/100,3)),', ',sprintf("%0.3f",round(confint(ROC_GSE16011_A)$CI_AUC[2,2]/100,3)),')'))
myPlot(plotDs=plotDs,
       colorList=c('darkblue','darkred'),
       smooth=T,
       main="ROC for OS prediction of Gravendeel",
       legendLab=legendLab,
       legendCol=c('darkblue','darkred'),
       legendLty=c(1,1))

dev.off()

