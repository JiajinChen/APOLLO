source('CoxInteraStepwise.R')
library(data.table)
library(survival)


# Main Effect
# Get stepwise results from SAS: CHIC2 IGF2BP2 ITGAV MSN PLCG1
# > summary(coxph(Surv(OS,Censor)~CHIC2+IGF2BP2+ITGAV+MSN+PLCG1,data=TCGA))$coefficient
# coef exp(coef) se(coef)     z  Pr(>|z|)
# CHIC2   0.2616     1.299  0.08678 3.015 2.571e-03
# IGF2BP2 0.3077     1.360  0.07171 4.291 1.776e-05
# ITGAV   0.2098     1.233  0.08892 2.359 1.831e-02
# MSN     0.4863     1.626  0.11161 4.357 1.319e-05
# PLCG1   0.3546     1.426  0.08960 3.957 7.577e-05
attach(TCGA)
TCGA$Main_4covar <- 0.2616*CHIC2 + 0.3077*IGF2BP2 + 0.2098*ITGAV + 0.4863*MSN + 0.3546*PLCG1
detach(TCGA)


# Interaction Effect
load('LGG_PanGene_Main_Inter.rdata')
LGG_PanGene_Inter1 <- LGG_PanGene_Inter[LGG_PanGene_Inter$FDR_TCGA_4covar <= 0.05 & LGG_PanGene_Inter$P_CGGA_4covar <= 0.05
                                        & LGG_PanGene_Inter$Z_TCGA_4covar*LGG_PanGene_Inter$Z_CGGA_4covar >0,]

varlist <- paste(LGG_PanGene_Inter1$Gene1,LGG_PanGene_Inter1$Gene2,sep = "*");
CoxInteraStepwise(
  variable.list = varlist,
  in.variable = 'Main_4covar',
  Time = "OS",
  Status = "Censor",
  data = TCGA,
  sle = 0.05,
  sls = 0.0501,
  vif.threshold = 10
)

# Result
# summary(coxph(formula = Surv(OS, Censor) ~ Main_4covar + BCORL1 + PRF1 +
#                 HMGA1 + TFG + CTNND2 + GOLGA5 + FAS + SMAD4 + BCORL1:PRF1 +
#                 HMGA1:TFG + CTNND2:GOLGA5 + FAS:SMAD4, data = TCGA))$coefficient
# coef exp(coef) se(coef)       z  Pr(>|z|)
# Main_4covar    1.13756    3.1191  0.11569  9.8328 8.134e-23
# BCORL1         0.23609    1.2663  0.09063  2.6050 9.187e-03
# PRF1          -0.10815    0.8975  0.09992 -1.0823 2.791e-01
# HMGA1         -0.16740    0.8459  0.09919 -1.6877 9.147e-02
# TFG           -0.10578    0.8996  0.11516 -0.9185 3.583e-01
# CTNND2        -0.19216    0.8252  0.10036 -1.9147 5.553e-02
# GOLGA5        -0.18138    0.8341  0.09879 -1.8360 6.636e-02
# FAS           -0.08875    0.9151  0.10790 -0.8225 4.108e-01
# SMAD4         -0.20733    0.8127  0.10646 -1.9475 5.147e-02
# BCORL1:PRF1   -0.24983    0.7789  0.09159 -2.7279 6.375e-03
# HMGA1:TFG      0.19302    1.2129  0.07601  2.5396 1.110e-02
# CTNND2:GOLGA5 -0.23400    0.7914  0.08460 -2.7661 5.673e-03
# FAS:SMAD4     -0.17242    0.8416  0.08001 -2.1550 3.116e-02


# Transcriptional Score
attach(TCGA)
TCGA$All_4covar <- 1.1376*Main_4covar + 
  0.2361*BCORL1-0.1082*PRF1-0.1674*HMGA1-0.1058*TFG-0.1922*CTNND2-0.1814*GOLGA5-0.0888*FAS-0.2073*SMAD4-
  0.2498*BCORL1*PRF1+0.1930*HMGA1*TFG-0.2340*CTNND2*GOLGA5-0.1724*FAS*SMAD4
detach(TCGA)
