library(survival)
library(timeROC)
model_A <- coxph(Surv(OS,Censor) ~ Clin_4covar+All_4covar,data = TCGA)

colnames(TCGA)[8] <- 'Histology'
colnames(CGGA1)[1] <- c("UID")

TCGA$ds <- 'TCGA';CGGA1$ds <- 'CGGA-seq';CGGA2$ds <- 'CGGA-array'
Rembrandt$ds <- 'Rembrandt';GSE16011$ds <- 'GSE16011';GSE61374$ds <- 'GSE61374'
colnames(GSE16011)[1] <- 'UID';colnames(GSE61374)[1] <- 'UID'

All_data <- merge(TCGA,CGGA1,all = T)
All_data <- merge(All_data,CGGA2,all = T)
All_data <- merge(All_data,Rembrandt,all = T)
All_data <- merge(All_data,GSE16011,all = T)
All_data <- merge(All_data,GSE61374,all = T)


All_data$Gender[All_data$Gender %in% c('female','F')] <- 'Female'
All_data$Gender[All_data$Gender %in% c('male','M')] <- 'Male'

All_data$Radio_status[All_data$Radio_status == '0'] <- 'NO'
All_data$Radio_status[All_data$Radio_status == '1'] <- 'YES'
All_data$Chemo_status[All_data$Chemo_status == '0'] <- 'NO'
All_data$Chemo_status[All_data$Chemo_status == '1'] <- 'YES'


# Accuracy AUC
# Age
ds <- subset(All_data,Age <= 40)
lp_ds <- predict(model_A,type="lp",newdata = ds)
ROC_ds <- timeROC(T=ds$OS,
                      delta=ds$Censor,marker=lp_ds,
                      cause=1,weighting="marginal",iid=T,
                      times=c(36,60))
ROC_ds$AUC;confint(ROC_ds)$CI_AUC


ds <- subset(All_data,Age > 40)
lp_ds <- predict(model_A,type="lp",newdata = ds)
ROC_ds <- timeROC(T=ds$OS,
                      delta=ds$Censor,marker=lp_ds,
                      cause=1,weighting="marginal",iid=T,
                      times=c(36,60))
ROC_ds$AUC;confint(ROC_ds)$CI_AUC

# Gender
ds <- subset(All_data,Gender == 'Male')
lp_ds <- predict(model_A,type="lp",newdata = ds)
ROC_ds <- timeROC(T=ds$OS,
                      delta=ds$Censor,marker=lp_ds,
                      cause=1,weighting="marginal",iid=T,
                      times=c(36,60))
ROC_ds$AUC;confint(ROC_ds)$CI_AUC

ds <- subset(All_data,Gender == 'Female')
lp_ds <- predict(model_A,type="lp",newdata = ds)
ROC_ds <- timeROC(T=ds$OS,
                      delta=ds$Censor,marker=lp_ds,
                      cause=1,weighting="marginal",iid=T,
                      times=c(36,60))
ROC_ds$AUC;confint(ROC_ds)$CI_AUC

# Grade
ds <- subset(All_data,Grade == 0)
lp_ds <- predict(model_A,type="lp",newdata = ds)
ROC_ds <- timeROC(T=ds$OS,
                      delta=ds$Censor,marker=lp_ds,
                      cause=1,weighting="marginal",iid=T,
                      times=c(36,60))
ROC_ds$AUC;confint(ROC_ds)$CI_AUC

ds <- subset(All_data,Grade == 1)
lp_ds <- predict(model_A,type="lp",newdata = ds)
ROC_ds <- timeROC(T=ds$OS,
                      delta=ds$Censor,marker=lp_ds,
                      cause=1,weighting="marginal",iid=T,
                      times=c(36,60))
ROC_ds$AUC;confint(ROC_ds)$CI_AUC

# IDH_mutation_status
ds <- subset(All_data,IDH_mutation_status == 0)
lp_ds <- predict(model_A,type="lp",newdata = ds)
ROC_ds <- timeROC(T=ds$OS,
                      delta=ds$Censor,marker=lp_ds,
                      cause=1,weighting="marginal",iid=T,
                      times=c(36,60))
ROC_ds$AUC;confint(ROC_ds)$CI_AUC

ds <- subset(All_data,IDH_mutation_status == 1)
lp_ds <- predict(model_A,type="lp",newdata = ds)
ROC_ds <- timeROC(T=ds$OS,
                      delta=ds$Censor,marker=lp_ds,
                      cause=1,weighting="marginal",iid=T,
                      times=c(36,60))
ROC_ds$AUC;confint(ROC_ds)$CI_AUC

# codeletion_1p19q_status
ds <- subset(All_data,codeletion_1p19q_status == 0)
lp_ds <- predict(model_A,type="lp",newdata = ds)
ROC_ds <- timeROC(T=ds$OS,
                      delta=ds$Censor,marker=lp_ds,
                      cause=1,weighting="marginal",iid=T,
                      times=c(36,60))
ROC_ds$AUC;confint(ROC_ds)$CI_AUC

ds <- subset(All_data,codeletion_1p19q_status == 1)
lp_ds <- predict(model_A,type="lp",newdata = ds)
ROC_ds <- timeROC(T=ds$OS,
                      delta=ds$Censor,marker=lp_ds,
                      cause=1,weighting="marginal",iid=T,
                      times=c(36,60))
ROC_ds$AUC;confint(ROC_ds)$CI_AUC

# MGMTp_methylation_status
ds <- subset(All_data,MGMTp_methylation_status == 'methylated')
lp_ds <- predict(model_A,type="lp",newdata = ds)
ROC_ds <- timeROC(T=ds$OS,
                      delta=ds$Censor,marker=lp_ds,
                      cause=1,weighting="marginal",iid=T,
                      times=c(36,60))
ROC_ds$AUC;confint(ROC_ds)$CI_AUC

ds <- subset(All_data,MGMTp_methylation_status == 'un-methylated')
lp_ds <- predict(model_A,type="lp",newdata = ds)
ROC_ds <- timeROC(T=ds$OS,
                      delta=ds$Censor,marker=lp_ds,
                      cause=1,weighting="marginal",iid=T,
                      times=c(36,60))
ROC_ds$AUC;confint(ROC_ds)$CI_AUC

# Radio_status
ds <- subset(All_data,Radio_status == 'YES')
lp_ds <- predict(model_A,type="lp",newdata = ds)
ROC_ds <- timeROC(T=ds$OS,
                      delta=ds$Censor,marker=lp_ds,
                      cause=1,weighting="marginal",iid=T,
                      times=c(36,60))
ROC_ds$AUC;confint(ROC_ds)$CI_AUC

ds <- subset(All_data,Radio_status == 'NO')
lp_ds <- predict(model_A,type="lp",newdata = ds)
ROC_ds <- timeROC(T=ds$OS,
                      delta=ds$Censor,marker=lp_ds,
                      cause=1,weighting="marginal",iid=T,
                      times=c(36,60))
ROC_ds$AUC;confint(ROC_ds)$CI_AUC

# Chemo_status
ds <- subset(All_data,Chemo_status == 'YES')
lp_ds <- predict(model_A,type="lp",newdata = ds)
ROC_ds <- timeROC(T=ds$OS,
                      delta=ds$Censor,marker=lp_ds,
                      cause=1,weighting="marginal",iid=T,
                      times=c(36,60))
ROC_ds$AUC;confint(ROC_ds)$CI_AUC

ds <- subset(All_data,Chemo_status == 'NO')
lp_ds <- predict(model_A,type="lp",newdata = ds)
ROC_ds <- timeROC(T=ds$OS,
                      delta=ds$Censor,marker=lp_ds,
                      cause=1,weighting="marginal",iid=T,
                      times=c(36,60))
ROC_ds$AUC;confint(ROC_ds)$CI_AUC


# Risk HR
All_data$APOLLO <- predict(model_A,type="lp",newdata = All_data)

All_data$APOLLO_Group <- ifelse(All_data$APOLLO <= -0.3595508,'Low-risk','High-risk')
All_data$APOLLO_Group <- relevel(as.factor(All_data$APOLLO_Group),ref = 'Low-risk')
summary(coxph(Surv(OS,Censor)~APOLLO_Group+ds,data=All_data[All_data$Age <= 40,]))
summary(coxph(Surv(OS,Censor)~APOLLO_Group+ds,data=All_data[All_data$Age > 40,]))
summary(coxph(Surv(OS,Censor)~APOLLO_Group+ds,data=All_data[All_data$Gender=='Male',]))
summary(coxph(Surv(OS,Censor)~APOLLO_Group+ds,data=All_data[All_data$Gender=='Female',]))
summary(coxph(Surv(OS,Censor)~APOLLO_Group+ds,data=All_data[All_data$Grade==0,]))
summary(coxph(Surv(OS,Censor)~APOLLO_Group+ds,data=All_data[All_data$Grade==1,]))
summary(coxph(Surv(OS,Censor)~APOLLO_Group+ds,data=All_data[All_data$IDH_mutation_status==1,]))
summary(coxph(Surv(OS,Censor)~APOLLO_Group+ds,data=All_data[All_data$IDH_mutation_status==0,]))
summary(coxph(Surv(OS,Censor)~APOLLO_Group+ds,data=All_data[All_data$codeletion_1p19q_status==1,]))
summary(coxph(Surv(OS,Censor)~APOLLO_Group+ds,data=All_data[All_data$codeletion_1p19q_status==0,]))
summary(coxph(Surv(OS,Censor)~APOLLO_Group+ds,data=All_data[All_data$MGMTp_methylation_status=='methylated',]))
summary(coxph(Surv(OS,Censor)~APOLLO_Group+ds,data=All_data[All_data$MGMTp_methylation_status=='un-methylated',]))
summary(coxph(Surv(OS,Censor)~APOLLO_Group+ds,data=All_data[All_data$Radio_status=='YES',]))
summary(coxph(Surv(OS,Censor)~APOLLO_Group+ds,data=All_data[All_data$Radio_status=='NO',]))
summary(coxph(Surv(OS,Censor)~APOLLO_Group+ds,data=All_data[All_data$Chemo_status=='YES',]))
summary(coxph(Surv(OS,Censor)~APOLLO_Group+ds,data=All_data[All_data$Chemo_status=='NO',]))
