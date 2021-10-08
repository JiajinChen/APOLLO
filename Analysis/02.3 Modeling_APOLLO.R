rm(list=ls())
source('CoxInteraStepwise.R')
load('TCGA_CGGA1.rdata')
library(data.table)
library(survival)

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


# Clinical Score
# > summary(coxph(Surv(OS,Censor)~Age+Grade+IDH_mutation_status+codeletion_1p19q_status,data=TCGA))$coefficient
# coef exp(coef)  se(coef)       z   Pr(>|z|)
# Age                      0.051046   1.05237 0.0080635  6.3305 2.4433e-10
# Grade                    0.861780   2.36737 0.2055866  4.1918 2.7674e-05
# IDH_mutation_status     -0.900118   0.40652 0.2190439 -4.1093 3.9685e-05
# codeletion_1p19q_status -0.843376   0.43026 0.2591314 -3.2546 1.1354e-03
TCGA$Clin_4covar <- 0.051046*TCGA$Age+0.861780*TCGA$Grade-0.900118*TCGA$IDH_mutation_status-0.843376*TCGA$codeletion_1p19q_status


# APOLLO
# > summary(coxph(Surv(OS,Censor) ~ Clin_4covar+All_4covar,data = TCGA))$coefficient
# coef exp(coef) se(coef)      z   Pr(>|z|)
# Clin_4covar 0.61217    1.8444 0.104947 5.8332 5.4374e-09
# All_4covar  0.75284    2.1230 0.084718 8.8864 6.3136e-19
TCGA$APOLLO <- 0.6122*TCGA$Clin_4covar + 0.7528*TCGA$All_4covar
