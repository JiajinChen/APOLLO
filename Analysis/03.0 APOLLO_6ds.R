rm(list = ls())
library(survival)
library(data.table)

gene <- c('CHIC2','IGF2BP2','ITGAV','MSN','PLCG1',
          'BCORL1','PRF1','HMGA1','TFG','CTNND2','GOLGA5','FAS','SMAD4')

# TCGA & CGGA1
load("Exp_use_Primary.rdata")
TCGA <- TCGA[TCGA$Grade!="",]
TCGA$Grade[TCGA$Grade==""] <- NA
TCGA$Gender[TCGA$Gender==''] <- NA
TCGA$Grade[TCGA$Grade =="G2"] <- 'WHO II'
TCGA$Grade[TCGA$Grade =="G3"] <- 'WHO III'
TCGA$OS <- TCGA$OS/30
CGGA1$OS <- CGGA1$OS/30
TCGA <- TCGA[,c(colnames(TCGA)[1:11],gene)]
CGGA1 <- CGGA1[,c(colnames(CGGA1)[1:11],gene)]
TCGA$Radio_status[TCGA$Radio_status==''] <- NA
TCGA <- TCGA[complete.cases(TCGA[,c('OS','Censor')]),]
CGGA1 <- CGGA1[complete.cases(CGGA1[,c('OS','Censor')]),]

# CGGA2
load("CGGA_D_clean.rdata")
CGGA2 <- as.data.frame(CGGA_D)
CGGA2 <- CGGA2[CGGA2$Grade %in% c('WHO II','WHO III') &
                   CGGA2$PRS_type == 'Primary',c(colnames(CGGA2)[1:14],gene)]
colnames(CGGA2)[9:13] <- c('Censor','Radio_status','Chemo_status','IDH_mutation_status','codeletion_1p19q_status')
CGGA2$OS <- CGGA2$OS/30
CGGA2 <- CGGA2[complete.cases(CGGA2[,c('OS','Censor')]),]

# Rembrandt
load("Rembrandt_clean.rdata")
Rembrandt <- as.data.frame(Rembrandt)
Age <- strsplit(Rembrandt$Age,"-")
Age <- lapply(Age, as.numeric)
Age <- as.numeric(lapply(Age, mean))
Rembrandt$Age <- Age
Rembrandt$OS <- Rembrandt$OS/30
Rembrandt <- Rembrandt[Rembrandt$Grade %in% c('WHO II','WHO III'),]
Rembrandt <- Rembrandt[,c(colnames(Rembrandt)[1:8],gene)]
Rembrandt <- Rembrandt[complete.cases(Rembrandt[,c('OS','Censor')]),]
rm(Age)

load("GSE16011.rdata")
GSE16011 <- GSE16011[GSE16011$Grade %in% c('G2','G3'),]
GSE16011 <- GSE16011[,c(colnames(GSE16011)[1:28],gene)]
GSE16011 <- GSE16011[complete.cases(GSE16011[,c('OS','Censor')]),]

load("GSE61374.rdata")
GSE61374 <- as.data.frame(GSE61374)
GSE61374 <- GSE61374[GSE61374$Grade %in% c('WHO II','WHO III'),]
GSE61374 <- GSE61374[,c(colnames(GSE61374)[1:18],gene)]
GSE61374$OS <- GSE61374$OS*12


# Scaled
for(k in 1:13){
  TCGA[,k+11] <- as.numeric(scale(TCGA[,k+11], center=T,scale=T))
  CGGA1[,k+11] <- as.numeric(scale(CGGA1[,k+11], center=T,scale=T))
  CGGA2[,k+14] <- as.numeric(scale(CGGA2[,k+14], center=T,scale=T))
  Rembrandt[,k+8] <- as.numeric(scale(Rembrandt[,k+8], center=T,scale=T))
  GSE16011[,k+28] <- as.numeric(scale(GSE16011[,k+28], center=T,scale=T))
  GSE61374[,k+18] <- as.numeric(scale(GSE61374[,k+18], center=T,scale=T))
}

# TCGA
# Main Score
attach(TCGA)
TCGA$Main_4covar <- 0.2616*CHIC2 + 0.3077*IGF2BP2 + 0.2098*ITGAV + 0.4863*MSN + 0.3546*PLCG1
detach(TCGA)
# Transcriptional Score
attach(TCGA)
TCGA$All_4covar <- 1.1376*Main_4covar + 
  0.2361*BCORL1-0.1082*PRF1-0.1674*HMGA1-0.1058*TFG-0.1922*CTNND2-0.1814*GOLGA5-0.0888*FAS-0.2073*SMAD4-
  0.2498*BCORL1*PRF1+0.1930*HMGA1*TFG-0.2340*CTNND2*GOLGA5-0.1724*FAS*SMAD4
detach(TCGA)

# CGGA1
# Main Score
attach(CGGA1)
CGGA1$Main_4covar <- 0.2616*CHIC2 + 0.3077*IGF2BP2 + 0.2098*ITGAV + 0.4863*MSN + 0.3546*PLCG1
detach(CGGA1)
# Transcriptional Score
attach(CGGA1)
CGGA1$All_4covar <- 1.1376*Main_4covar + 
  0.2361*BCORL1-0.1082*PRF1-0.1674*HMGA1-0.1058*TFG-0.1922*CTNND2-0.1814*GOLGA5-0.0888*FAS-0.2073*SMAD4-
  0.2498*BCORL1*PRF1+0.1930*HMGA1*TFG-0.2340*CTNND2*GOLGA5-0.1724*FAS*SMAD4
detach(CGGA1)

# CGGA2
# Main Score
attach(CGGA2)
CGGA2$Main_4covar <- 0.2616*CHIC2 + 0.3077*IGF2BP2 + 0.2098*ITGAV + 0.4863*MSN + 0.3546*PLCG1
detach(CGGA2)
# Transcriptional Score
attach(CGGA2)
CGGA2$All_4covar <- 1.1376*Main_4covar + 
  0.2361*BCORL1-0.1082*PRF1-0.1674*HMGA1-0.1058*TFG-0.1922*CTNND2-0.1814*GOLGA5-0.0888*FAS-0.2073*SMAD4-
  0.2498*BCORL1*PRF1+0.1930*HMGA1*TFG-0.2340*CTNND2*GOLGA5-0.1724*FAS*SMAD4
detach(CGGA2)

# Rembrandt
# Main Score
attach(Rembrandt)
Rembrandt$Main_4covar <- 0.2616*CHIC2 + 0.3077*IGF2BP2 + 0.2098*ITGAV + 0.4863*MSN + 0.3546*PLCG1
detach(Rembrandt)
# Transcriptional Score
attach(Rembrandt)
Rembrandt$All_4covar <- 1.1376*Main_4covar + 
  0.2361*BCORL1-0.1082*PRF1-0.1674*HMGA1-0.1058*TFG-0.1922*CTNND2-0.1814*GOLGA5-0.0888*FAS-0.2073*SMAD4-
  0.2498*BCORL1*PRF1+0.1930*HMGA1*TFG-0.2340*CTNND2*GOLGA5-0.1724*FAS*SMAD4
detach(Rembrandt)

# GSE16011
# Main Score
attach(GSE16011)
GSE16011$Main_4covar <- 0.2616*CHIC2 + 0.3077*IGF2BP2 + 0.2098*ITGAV + 0.4863*MSN + 0.3546*PLCG1
detach(GSE16011)
# Transcriptional Score
attach(GSE16011)
GSE16011$All_4covar <- 1.1376*Main_4covar + 
  0.2361*BCORL1-0.1082*PRF1-0.1674*HMGA1-0.1058*TFG-0.1922*CTNND2-0.1814*GOLGA5-0.0888*FAS-0.2073*SMAD4-
  0.2498*BCORL1*PRF1+0.1930*HMGA1*TFG-0.2340*CTNND2*GOLGA5-0.1724*FAS*SMAD4
detach(GSE16011)


# GSE61374
# Main Score
attach(GSE61374)
GSE61374$Main_4covar <- 0.2616*CHIC2 + 0.3077*IGF2BP2 + 0.2098*ITGAV + 0.4863*MSN + 0.3546*PLCG1
detach(GSE61374)
# Transcriptional Score
attach(GSE61374)
GSE61374$All_4covar <- 1.1376*Main_4covar + 
  0.2361*BCORL1-0.1082*PRF1-0.1674*HMGA1-0.1058*TFG-0.1922*CTNND2-0.1814*GOLGA5-0.0888*FAS-0.2073*SMAD4-
  0.2498*BCORL1*PRF1+0.1930*HMGA1*TFG-0.2340*CTNND2*GOLGA5-0.1724*FAS*SMAD4
detach(GSE61374)


# Clinical
TCGA$Grade<-ifelse(TCGA$Grade=="WHO II",0,1)
TCGA$IDH_mutation_status<-ifelse(TCGA$IDH_mutation_status=="Mutant",1,0)
TCGA$codeletion_1p19q_status<-ifelse(TCGA$codeletion_1p19q_status=="Codel",1,0)


CGGA1$Grade<-ifelse(CGGA1$Grade=="WHO II",0,1)
CGGA1$IDH_mutation_status<-ifelse(CGGA1$IDH_mutation_status=="Mutant",1,0)
CGGA1$codeletion_1p19q_status<-ifelse(CGGA1$codeletion_1p19q_status=="Codel",1,0)

CGGA2$Grade<-ifelse(CGGA2$Grade=="WHO II",0,1)
CGGA2$IDH_mutation_status<-ifelse(CGGA2$IDH_mutation_status=="Mutant",1,0)
CGGA2$codeletion_1p19q_status <-ifelse(CGGA2$codeletion_1p19q_status=="Codel",1,0)

Rembrandt$Grade<-ifelse(Rembrandt$Grade=="WHO II",0,1)
Rembrandt$codeletion_1p19q_status <- ifelse(Rembrandt$codeletion_1p19q_status=="Codel",1,0)

GSE16011$Grade<-ifelse(GSE16011$Grade=="G2",0,1)
GSE16011$IDH_mutation_status[GSE16011$IDH_mutation_status==''] <- NA
GSE16011$IDH_mutation_status<-ifelse(GSE16011$IDH_mutation_status=="mutation",1,0)
GSE16011$codeletion_1p19q_status <- ifelse(GSE16011$codeletion_1p19q_status=="Codel",1,0)

GSE61374$Grade<-ifelse(GSE61374$Grade=='WHO II',0,1)
GSE61374$IDH_mutation_status<-ifelse(GSE61374$IDH_mutation_status=='mut',1,0)
GSE61374$codeletion_1p19q_status <- ifelse(GSE61374$codeletion_1p19q_status=='del',1,0)


TCGA$Clin_4covar <- 0.051046*TCGA$Age+0.861780*TCGA$Grade-0.900118*TCGA$IDH_mutation_status-0.843376*TCGA$codeletion_1p19q_status
CGGA1$Clin_4covar <- 0.051046*CGGA1$Age+0.861780*CGGA1$Grade-0.900118*CGGA1$IDH_mutation_status-0.843376*CGGA1$codeletion_1p19q_status
CGGA2$Clin_4covar <- 0.051046*CGGA2$Age+0.861780*CGGA2$Grade-0.900118*CGGA2$IDH_mutation_status-0.843376*CGGA2$codeletion_1p19q_status
Rembrandt$Clin_4covar <- 0.051046*Rembrandt$Age+0.861780*Rembrandt$Grade-0.843376*Rembrandt$codeletion_1p19q_status
GSE61374$Clin_4covar <- 0.051046*GSE61374$Age+0.861780*GSE61374$Grade-0.900118*GSE61374$IDH_mutation_status-0.843376*GSE61374$codeletion_1p19q_status
GSE16011$Clin_4covar <- 0.051046*GSE16011$Age+0.861780*GSE16011$Grade-0.900118*GSE16011$IDH_mutation_status-0.843376*GSE16011$codeletion_1p19q_status


# APOLLO
model_A <- coxph(Surv(OS,Censor) ~ Clin_4covar+All_4covar,data = TCGA)
TCGA$APOLLO <- predict(model_A,type="lp",newdata = TCGA)
CGGA1$APOLLO <- predict(model_A,type="lp",newdata = CGGA1)
CGGA2$APOLLO <- predict(model_A,type="lp",newdata = CGGA2)
Rembrandt$APOLLO <- predict(model_A,type="lp",newdata = Rembrandt)
GSE61374$APOLLO <- predict(model_A,type="lp",newdata = GSE61374)
GSE16011$APOLLO <- predict(model_A,type="lp",newdata = GSE16011)

save(TCGA,CGGA1,CGGA2,Rembrandt,GSE61374,GSE16011,file='APOLLO.rdata')
