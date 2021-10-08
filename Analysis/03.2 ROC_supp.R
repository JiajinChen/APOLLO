# rm(list=ls())
source('ROCplot.R')

model_C <- coxph(Surv(OS,Censor) ~ Clin_4covar,data = TCGA)
model_A <- coxph(Surv(OS,Censor) ~ Clin_4covar+All_4covar,data = TCGA)

# TCGA
lp_TCGA_C <- predict(model_C,type="lp",newdata = TCGA)
ROC_TCGA_C <- timeROC(T=TCGA$OS,
                      delta=TCGA$Censor,marker=lp_TCGA_C,
                      cause=1,weighting="marginal",iid = T,
                      times=c(36,60))
lp_TCGA_A <- predict(model_A,type="lp",newdata = TCGA)
ROC_TCGA_A <- timeROC(T=TCGA$OS,
                      delta=TCGA$Censor,marker=lp_TCGA_A,
                      cause=1,weighting="marginal",iid = T,
                      times=c(36,60))

# CGGA-Seq
lp_CGGASeq_C <- predict(model_C,type="lp",newdata = CGGA1)
ROC_CGGASeq_C <- timeROC(T=CGGA1$OS,
                         delta=CGGA1$Censor,marker=lp_CGGASeq_C,
                         cause=1,weighting="marginal",iid = T,
                         times=c(36,60))
lp_CGGASeq_A <- predict(model_A,type="lp",newdata = CGGA1)
ROC_CGGASeq_A <- timeROC(T=CGGA1$OS,
                         delta=CGGA1$Censor,marker=lp_CGGASeq_A,
                         cause=1,weighting="marginal",iid = T,
                         times=c(36,60))

# CGGA-Array
lp_CGGAArray_C <- predict(model_C,type="lp",newdata = CGGA2)
ROC_CGGAArray_C <- timeROC(T=CGGA2$OS,
                           delta=CGGA2$Censor,marker=lp_CGGAArray_C,
                           cause=1,weighting="marginal",iid = T,
                           times=c(36,60))
lp_CGGAArray_A <- predict(model_A,type="lp",newdata = CGGA2)
ROC_CGGAArray_A <- timeROC(T=CGGA2$OS,
                           delta=CGGA2$Censor,marker=lp_CGGAArray_A,
                           cause=1,weighting="marginal",iid = T,
                           times=c(36,60))

# Rembrandt
lp_Rembrandt_C <- predict(model_C,type="lp",newdata = Rembrandt)
ROC_Rembrandt_C <- timeROC(T=Rembrandt$OS,
                           delta=Rembrandt$Censor,marker=lp_Rembrandt_C,
                           cause=1,weighting="marginal",iid = T,
                           times=c(36,60))
lp_Rembrandt_A <- predict(model_A,type="lp",newdata = Rembrandt)
ROC_Rembrandt_A <- timeROC(T=Rembrandt$OS,
                           delta=Rembrandt$Censor,marker=lp_Rembrandt_A,
                           cause=1,weighting="marginal",iid = T,
                           times=c(36,60))

# GSE61374
lp_GSE61374_C <- predict(model_C,type="lp",newdata = GSE61374)
ROC_GSE61374_C <- timeROC(T=GSE61374$OS,
                          delta=GSE61374$Censor,marker=lp_GSE61374_C,
                          cause=1,weighting="marginal",iid = T,
                          times=c(36,60))
lp_GSE61374_A <- predict(model_A,type="lp",newdata = GSE61374)
ROC_GSE61374_A <- timeROC(T=GSE61374$OS,
                          delta=GSE61374$Censor,marker=lp_GSE61374_A,
                          cause=1,weighting="marginal",iid = T,
                          times=c(36,60))

# GSE16011
lp_GSE16011_C <- predict(model_C,type="lp",newdata = GSE16011)
ROC_GSE16011_C <- timeROC(T=GSE16011$OS,
                          delta=GSE16011$Censor,marker=lp_GSE16011_C,
                          cause=1,weighting="marginal",iid = T,
                          times=c(36,60))
lp_GSE16011_A <- predict(model_A,type="lp",newdata = GSE16011)
ROC_GSE16011_A <- timeROC(T=GSE16011$OS,
                          delta=GSE16011$Censor,marker=lp_GSE16011_A,
                          cause=1,weighting="marginal",iid = T,
                          times=c(36,60))




pdf("Fig_ROC_CA2_Smooth.pdf",width=9.9,height = 13.2)
par(pin = c(9.9,13.2),mfrow=c(4,3),mar = rep(3.2, 4))

# TCGA 3
plotDs=list(ROC_TCGA_C,ROC_TCGA_A);index=1
legendLab <- c(paste0('Covariates: ',sprintf("%0.3f",round(plotDs[[1]]$AUC[index],3)),
                      ' (',sprintf("%0.3f",round(confint(plotDs[[1]])$CI_AUC[index,1]/100,3)),', ',sprintf("%0.3f",round(confint(plotDs[[1]])$CI_AUC[index,2]/100,3)),')'),
               paste0('APOLLO: ',sprintf("%0.3f",round(plotDs[[2]]$AUC[index],3)),
                      ' (',sprintf("%0.3f",round(confint(plotDs[[2]])$CI_AUC[index,1]/100,3)),', ',sprintf("%0.3f",round(confint(plotDs[[2]])$CI_AUC[index,2]/100,3)),')'))
myPlot(plotDs=plotDs,
       colorList=c('darkblue','darkred'),
       index=index,
       smooth=T,
       main="ROC for 36-month OS of TCGA",
       legendLab=legendLab,
       legendCol=c('darkblue','darkred'),
       legendLty=c(1,1))

# CGGA-Seq 3
plotDs=list(ROC_CGGASeq_C,ROC_CGGASeq_A);index=1
legendLab <- c(paste0('Covariates: ',sprintf("%0.3f",round(plotDs[[1]]$AUC[index],3)),
                      ' (',sprintf("%0.3f",round(confint(plotDs[[1]])$CI_AUC[index,1]/100,3)),', ',sprintf("%0.3f",round(confint(plotDs[[1]])$CI_AUC[index,2]/100,3)),')'),
               paste0('APOLLO: ',sprintf("%0.3f",round(plotDs[[2]]$AUC[index],3)),
                      ' (',sprintf("%0.3f",round(confint(plotDs[[2]])$CI_AUC[index,1]/100,3)),', ',sprintf("%0.3f",round(confint(plotDs[[2]])$CI_AUC[index,2]/100,3)),')'))
myPlot(plotDs=plotDs,
       colorList=c('darkblue','darkred'),
       smooth=T,
       main="ROC for 36-month OS of CGGA1",
       legendLab=legendLab,
       legendCol=c('darkblue','darkred'),
       legendLty=c(1,1))

# CGGA-Array 3
plotDs=list(ROC_CGGAArray_C,ROC_CGGAArray_A);index=1
legendLab <- c(paste0('Covariates: ',sprintf("%0.3f",round(plotDs[[1]]$AUC[index],3)),
                      ' (',sprintf("%0.3f",round(confint(plotDs[[1]])$CI_AUC[index,1]/100,3)),', ',sprintf("%0.3f",round(confint(plotDs[[1]])$CI_AUC[index,2]/100,3)),')'),
               paste0('APOLLO: ',sprintf("%0.3f",round(plotDs[[2]]$AUC[index],3)),
                      ' (',sprintf("%0.3f",round(confint(plotDs[[2]])$CI_AUC[index,1]/100,3)),', ',sprintf("%0.3f",round(confint(plotDs[[2]])$CI_AUC[index,2]/100,3)),')'))
myPlot(plotDs=plotDs,
       colorList=c('darkblue','darkred'),
       smooth=T,
       main="ROC for 36-month OS of CGGA2",
       legendLab=legendLab,
       legendCol=c('darkblue','darkred'),
       legendLty=c(1,1))

# Rembrandt 3
plotDs=list(ROC_Rembrandt_C,ROC_Rembrandt_A);index=1
legendLab <- c(paste0('Covariates: ',sprintf("%0.3f",round(plotDs[[1]]$AUC[index],3)),
                      ' (',sprintf("%0.3f",round(confint(plotDs[[1]])$CI_AUC[index,1]/100,3)),', ',sprintf("%0.3f",round(confint(plotDs[[1]])$CI_AUC[index,2]/100,3)),')'),
               paste0('APOLLO: ',sprintf("%0.3f",round(plotDs[[2]]$AUC[index],3)),
                      ' (',sprintf("%0.3f",round(confint(plotDs[[2]])$CI_AUC[index,1]/100,3)),', ',sprintf("%0.3f",round(confint(plotDs[[2]])$CI_AUC[index,2]/100,3)),')'))
myPlot(plotDs=plotDs,
       colorList=c('darkblue','darkred'),
       smooth=T,
       main="ROC for 36-month OS of Rembrandt",
       legendLab=legendLab,
       legendCol=c('darkblue','darkred'),
       legendLty=c(1,1))

# GSE61374 3
plotDs=list(ROC_GSE61374_C,ROC_GSE61374_A);index=1
legendLab <- c(paste0('Covariates: ',sprintf("%0.3f",round(plotDs[[1]]$AUC[index],3)),
                      ' (',sprintf("%0.3f",round(confint(plotDs[[1]])$CI_AUC[index,1]/100,3)),', ',sprintf("%0.3f",round(confint(plotDs[[1]])$CI_AUC[index,2]/100,3)),')'),
               paste0('APOLLO: ',sprintf("%0.3f",round(plotDs[[2]]$AUC[index],3)),
                      ' (',sprintf("%0.3f",round(confint(plotDs[[2]])$CI_AUC[index,1]/100,3)),', ',sprintf("%0.3f",round(confint(plotDs[[2]])$CI_AUC[index,2]/100,3)),')'))
myPlot(plotDs=plotDs,
       colorList=c('darkblue','darkred'),
       smooth=T,
       main="ROC for 36-month OS of Weller",
       legendLab=legendLab,
       legendCol=c('darkblue','darkred'),
       legendLty=c(1,1))

# GSE16011 3
plotDs=list(ROC_GSE16011_C,ROC_GSE16011_A);index=1
legendLab <- c(paste0('Covariates: ',sprintf("%0.3f",round(plotDs[[1]]$AUC[index],3)),
                      ' (',sprintf("%0.3f",round(confint(plotDs[[1]])$CI_AUC[index,1]/100,3)),', ',sprintf("%0.3f",round(confint(plotDs[[1]])$CI_AUC[index,2]/100,3)),')'),
               paste0('APOLLO: ',sprintf("%0.3f",round(plotDs[[2]]$AUC[index],3)),
                      ' (',sprintf("%0.3f",round(confint(plotDs[[2]])$CI_AUC[index,1]/100,3)),', ',sprintf("%0.3f",round(confint(plotDs[[2]])$CI_AUC[index,2]/100,3)),')'))
myPlot(plotDs=plotDs,
       colorList=c('darkblue','darkred'),
       smooth=T,
       main="ROC for 36-month OS of Gravendeel",
       legendLab=legendLab,
       legendCol=c('darkblue','darkred'),
       legendLty=c(1,1))


# TCGA 5
plotDs=list(ROC_TCGA_C,ROC_TCGA_A);index=2
legendLab <- c(paste0('Covariates: ',sprintf("%0.3f",round(plotDs[[1]]$AUC[index],3)),
                      ' (',sprintf("%0.3f",round(confint(plotDs[[1]])$CI_AUC[index,1]/100,3)),', ',sprintf("%0.3f",round(confint(plotDs[[1]])$CI_AUC[index,2]/100,3)),')'),
               paste0('APOLLO: ',sprintf("%0.3f",round(plotDs[[2]]$AUC[index],3)),
                      ' (',sprintf("%0.3f",round(confint(plotDs[[2]])$CI_AUC[index,1]/100,3)),', ',sprintf("%0.3f",round(confint(plotDs[[2]])$CI_AUC[index,2]/100,3)),')'))
myPlot(plotDs=plotDs,
       colorList=c('darkblue','darkred'),
       index=index,
       smooth=T,
       main="ROC for 60-month OS of TCGA",
       legendLab=legendLab,
       legendCol=c('darkblue','darkred'),
       legendLty=c(1,1))

# CGGA-Seq 5
plotDs=list(ROC_CGGASeq_C,ROC_CGGASeq_A);index=2
legendLab <- c(paste0('Covariates: ',sprintf("%0.3f",round(plotDs[[1]]$AUC[index],3)),
                      ' (',sprintf("%0.3f",round(confint(plotDs[[1]])$CI_AUC[index,1]/100,3)),', ',sprintf("%0.3f",round(confint(plotDs[[1]])$CI_AUC[index,2]/100,3)),')'),
               paste0('APOLLO: ',sprintf("%0.3f",round(plotDs[[2]]$AUC[index],3)),
                      ' (',sprintf("%0.3f",round(confint(plotDs[[2]])$CI_AUC[index,1]/100,3)),', ',sprintf("%0.3f",round(confint(plotDs[[2]])$CI_AUC[index,2]/100,3)),')'))
myPlot(plotDs=plotDs,
       colorList=c('darkblue','darkred'),
       smooth=T,
       main="ROC for 60-month OS of CGGA1",
       legendLab=legendLab,
       legendCol=c('darkblue','darkred'),
       legendLty=c(1,1))

# CGGA-Array 5
plotDs=list(ROC_CGGAArray_C,ROC_CGGAArray_A);index=2
legendLab <- c(paste0('Covariates: ',sprintf("%0.3f",round(plotDs[[1]]$AUC[index],3)),
                      ' (',sprintf("%0.3f",round(confint(plotDs[[1]])$CI_AUC[index,1]/100,3)),', ',sprintf("%0.3f",round(confint(plotDs[[1]])$CI_AUC[index,2]/100,3)),')'),
               paste0('APOLLO: ',sprintf("%0.3f",round(plotDs[[2]]$AUC[index],3)),
                      ' (',sprintf("%0.3f",round(confint(plotDs[[2]])$CI_AUC[index,1]/100,3)),', ',sprintf("%0.3f",round(confint(plotDs[[2]])$CI_AUC[index,2]/100,3)),')'))
myPlot(plotDs=plotDs,
       colorList=c('darkblue','darkred'),
       smooth=T,
       main="ROC for 60-month OS of CGGA2",
       legendLab=legendLab,
       legendCol=c('darkblue','darkred'),
       legendLty=c(1,1))

# Rembrandt 5
plotDs=list(ROC_Rembrandt_C,ROC_Rembrandt_A);index=2
legendLab <- c(paste0('Covariates: ',sprintf("%0.3f",round(plotDs[[1]]$AUC[index],3)),
                      ' (',sprintf("%0.3f",round(confint(plotDs[[1]])$CI_AUC[index,1]/100,3)),', ',sprintf("%0.3f",round(confint(plotDs[[1]])$CI_AUC[index,2]/100,3)),')'),
               paste0('APOLLO: ',sprintf("%0.3f",round(plotDs[[2]]$AUC[index],3)),
                      ' (',sprintf("%0.3f",round(confint(plotDs[[2]])$CI_AUC[index,1]/100,3)),', ',sprintf("%0.3f",round(confint(plotDs[[2]])$CI_AUC[index,2]/100,3)),')'))
myPlot(plotDs=plotDs,
       colorList=c('darkblue','darkred'),
       smooth=T,
       main="ROC for 60-month OS of Rembrandt",
       legendLab=legendLab,
       legendCol=c('darkblue','darkred'),
       legendLty=c(1,1))

# GSE61374 5
plotDs=list(ROC_GSE61374_C,ROC_GSE61374_A);index=2
legendLab <- c(paste0('Covariates: ',sprintf("%0.3f",round(plotDs[[1]]$AUC[index],3)),
                      ' (',sprintf("%0.3f",round(confint(plotDs[[1]])$CI_AUC[index,1]/100,3)),', ',sprintf("%0.3f",round(confint(plotDs[[1]])$CI_AUC[index,2]/100,3)),')'),
               paste0('APOLLO: ',sprintf("%0.3f",round(plotDs[[2]]$AUC[index],3)),
                      ' (',sprintf("%0.3f",round(confint(plotDs[[2]])$CI_AUC[index,1]/100,3)),', ',sprintf("%0.3f",round(confint(plotDs[[2]])$CI_AUC[index,2]/100,3)),')'))
myPlot(plotDs=plotDs,
       colorList=c('darkblue','darkred'),
       smooth=T,
       main="ROC for 60-month OS of Weller",
       legendLab=legendLab,
       legendCol=c('darkblue','darkred'),
       legendLty=c(1,1))

# GSE16011 5
plotDs=list(ROC_GSE16011_C,ROC_GSE16011_A);index=2
legendLab <- c(paste0('Covariates: ',sprintf("%0.3f",round(plotDs[[1]]$AUC[index],3)),
                      ' (',sprintf("%0.3f",round(confint(plotDs[[1]])$CI_AUC[index,1]/100,3)),', ',sprintf("%0.3f",round(confint(plotDs[[1]])$CI_AUC[index,2]/100,3)),')'),
               paste0('APOLLO: ',sprintf("%0.3f",round(plotDs[[2]]$AUC[index],3)),
                      ' (',sprintf("%0.3f",round(confint(plotDs[[2]])$CI_AUC[index,1]/100,3)),', ',sprintf("%0.3f",round(confint(plotDs[[2]])$CI_AUC[index,2]/100,3)),')'))
myPlot(plotDs=plotDs,
       colorList=c('darkblue','darkred'),
       smooth=T,
       main="ROC for 36-month OS of Gravendeel",
       legendLab=legendLab,
       legendCol=c('darkblue','darkred'),
       legendLty=c(1,1))


dev.off()



