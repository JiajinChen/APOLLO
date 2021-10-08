rm(list=ls())
library(survival)
source("D:\\00 ToolsCode\\DCA\\R code\\stdca.R")

TCGA$ds <- 'TCGA'; CGGA1$ds <- 'CGGA-seq'; CGGA2$ds <- 'CGGA-array'
Rembrandt$ds <- 'Rembrandt'; GSE61374$ds <- 'GSE61374'; GSE16011$ds <- 'GSE16011'
colnames(GSE16011)[1] <- 'UID';colnames(GSE61374)[1] <- 'UID';colnames(CGGA1)[1] <- c("UID")

All_data <- merge(TCGA,CGGA1,all = T)
All_data <- merge(All_data,CGGA2,all = T)
All_data <- merge(All_data,Rembrandt,all = T)
All_data <- merge(All_data,GSE16011,all = T)
All_data <- merge(All_data,GSE61374,all = T)



model_C <- coxph(Surv(OS,Censor) ~ Clin_4covar+ds,data = All_data)
model_A <- coxph(Surv(OS,Censor) ~ Clin_4covar+All_4covar+ds,data = All_data)

All_data <- All_data[complete.cases(All_data[,c('Clin_4covar','All_4covar')]),]


All_data$Covariate3 <- c(1- (summary(survfit(model_C,
                                            newdata=All_data), times=36)$surv))
All_data$APOLLO3 <- c(1- (summary(survfit(model_A,
                                         newdata=All_data), times=36)$surv))

All_data$Covariate5 <- c(1- (summary(survfit(model_C,
                                             newdata=All_data), times=60)$surv))
All_data$APOLLO5 <- c(1- (summary(survfit(model_A,
                                          newdata=All_data), times=60)$surv))

par(mfrow=c(2,2))
fit1 <- stdca(data=All_data, outcome="Censor", ttoutcome="OS", timepoint=36,
              predictors=c('Covariate3','APOLLO3'), xstop=0.88,intervention = F)
fit2 <- stdca(data=All_data, outcome="Censor", ttoutcome="OS", timepoint=36,
              predictors=c('Covariate3','APOLLO3'), xstop=0.88,intervention = T)
fit3 <- stdca(data=All_data, outcome="Censor", ttoutcome="OS", timepoint=60,
              predictors=c('Covariate5','APOLLO5'), xstop=0.97,intervention = F)
fit4 <- stdca(data=All_data, outcome="Censor", ttoutcome="OS", timepoint=60,
              predictors=c('Covariate5','APOLLO5'), xstop=0.97,intervention = T)

# 36-month
a <- fit1$net.benefit
a <- a[a$threshold >=0 & a$threshold <=0.5,]
a$delta <- a$APOLLO- a$Covariate
mean(a$APOLLO);mean(a$Covariate)
t.test(a$APOLLO,a$Covariate,paired = T)$p.value
mean(a$delta)

b <- fit1$interventions.avoided
b <- b[b$threshold >=0 & b$threshold <=0.5,]
b$delta <- b$APOLLO- b$Covariate
mean(b$APOLLO);mean(b$Covariate)
t.test(b$APOLLO,b$Covariate,paired = T)$p.value
mean(b$delta)

# 60-month
a <- fit3$net.benefit
a <- a[a$threshold >=0 & a$threshold <=0.5,]
a$delta <- a$APOLLO- a$Covariate
mean(a$APOLLO);mean(a$Covariate)
t.test(a$APOLLO,a$Covariate,paired = T)$p.value
mean(a$delta)

b <- fit3$interventions.avoided
b <- b[b$threshold >=0 & b$threshold <=0.5,]
b$delta <- b$APOLLO- b$Covariate
mean(b$APOLLO);mean(b$Covariate)
t.test(b$APOLLO,b$Covariate,paired = T)$p.value
mean(b$delta)
