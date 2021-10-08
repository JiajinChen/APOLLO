rm(list=ls())
setwd('D:/APOLLO')
library(data.table)
library(survival)

load("Exp_use_Primary.rdata")
TCGA <- TCGA[TCGA$Grade!="",]
TCGA$Grade[TCGA$Grade==""] <- NA
TCGA$Gender[TCGA$Gender==''] <- NA
TCGA$Grade[TCGA$Grade =="G2"] <- 'WHO II'
TCGA$Grade[TCGA$Grade =="G3"] <- 'WHO III'

TCGA$OS <- TCGA$OS/30
CGGA1$OS <- CGGA1$OS/30

Gene_Cancer <-fread("E:/BioTools/Gene_Info/Pan-cancer/Census_allWed Aug 5 09 34 12 2020.csv")
Gene_Cancer$`Gene Symbol` <- gsub("-","_",Gene_Cancer$`Gene Symbol`,fixed=T)
Gene_Cancer$`Gene Symbol` <- gsub(" ","",Gene_Cancer$`Gene Symbol`,fixed=T)

Gene_Cancer <- Gene_Cancer$`Gene Symbol`[Gene_Cancer$`Gene Symbol` %in% colnames(TCGA)]
Gene_Cancer <- Gene_Cancer[Gene_Cancer %in% colnames(CGGA1)]
TCGA <- subset(TCGA,select = c(1:11,which(colnames(TCGA) %in% Gene_Cancer)))
CGGA1 <- subset(CGGA1,select = c(1:11,which(colnames(CGGA1) %in% Gene_Cancer)))
TCGA <- as.data.frame(TCGA)
CGGA1 <- as.data.frame(CGGA1)


# Scaled
for(k in 12:691){
  TCGA[,k] <- as.numeric(scale(TCGA[,k], center=T,scale=T))
  CGGA1[,k] <- as.numeric(scale(CGGA1[,k], center=T,scale=T))
}

TCGA$Grade<-ifelse(TCGA$Grade=="WHO II",0,1)
TCGA$IDH_mutation_status<-ifelse(TCGA$IDH_mutation_status=="Mutant",1,0)
TCGA$codeletion_1p19q_status<-ifelse(TCGA$codeletion_1p19q_status=="Codel",1,0)


CGGA1$Grade<-ifelse(CGGA1$Grade=="WHO II",0,1)
CGGA1$IDH_mutation_status<-ifelse(CGGA1$IDH_mutation_status=="Mutant",1,0)
CGGA1$codeletion_1p19q_status<-ifelse(CGGA1$codeletion_1p19q_status=="Codel",1,0)


# > summary(coxph(Surv(OS,Censor)~Age+Grade+IDH_mutation_status+codeletion_1p19q_status,data=TCGA))$coefficient
# coef exp(coef)  se(coef)       z   Pr(>|z|)
# Age                      0.051046   1.05237 0.0080635  6.3305 2.4433e-10
# Grade                    0.861780   2.36737 0.2055866  4.1918 2.7674e-05
# IDH_mutation_status     -0.900118   0.40652 0.2190439 -4.1093 3.9685e-05
# codeletion_1p19q_status -0.843376   0.43026 0.2591314 -3.2546 1.1354e-03


# Main Effect
LGG_PanGene_Main <- data.frame()
for(k in 12:691){
  Gene <- colnames(TCGA)[k]
  fit_TCGA <- try(coxph(as.formula(paste0("Surv(OS, Censor) ~ Age+as.factor(Grade)+as.factor(IDH_mutation_status)+as.factor(codeletion_1p19q_status)+",Gene)), data=TCGA))
  Z_TCGA <- summary(fit_TCGA)$coefficient[5,4]
  P_TCGA <- summary(fit_TCGA)$coefficient[5,5]
  fit_CGGA <- try(coxph(as.formula(paste0("Surv(OS, Censor) ~ Age+as.factor(Grade)+as.factor(IDH_mutation_status)+as.factor(codeletion_1p19q_status)+",Gene)), data=CGGA1))
  Z_CGGA <- summary(fit_CGGA)$coefficient[5,4]
  P_CGGA <- summary(fit_CGGA)$coefficient[5,5]

  LGG_PanGene_Main <- rbind(LGG_PanGene_Main,data.frame(Gene=Gene,Z_TCGA=Z_TCGA,P_TCGA=P_TCGA,
                                                        Z_CGGA=Z_CGGA,P_CGGA=P_CGGA,
                                                        stringsAsFactors=F))
}


LGG_PanGene_Main$FDR_TCGA <- p.adjust(LGG_PanGene_Main$P_TCGA,method="fdr")


# Interaction Effect
library(parallel)
cl <- makeCluster(32)
params <- expand.grid(
  i = 12:690,
  j = 13:691)
params <- subset(params,j>i)
clusterEvalQ(cl, library(survival))
clusterExport(cl, "params")
clusterExport(cl, "TCGA")
clusterExport(cl, "CGGA1")

system.time(
  paraRes <- parSapply(cl,1L:nrow(params), function(row_id) {
    i <- params$i[row_id]
    j <- params$j[row_id]
    Gene1 <- colnames(TCGA)[i]
    Gene2 <- colnames(TCGA)[j]
    
    fit_TCGA <- try(coxph(as.formula(paste0("Surv(OS, Censor) ~ Age+as.factor(Grade)+as.factor(IDH_mutation_status)+as.factor(codeletion_1p19q_status)+",Gene1,"*",Gene2)), data=TCGA))
    Z_TCGA <- summary(fit_TCGA)$coefficient[7,4]
    P_TCGA <- summary(fit_TCGA)$coefficient[7,5]
    fit_CGGA <- try(coxph(as.formula(paste0("Surv(OS, Censor) ~ Age+as.factor(Grade)+as.factor(IDH_mutation_status)+as.factor(codeletion_1p19q_status)+",Gene1,"*",Gene2)), data=CGGA1))
    Z_CGGA <- summary(fit_CGGA)$coefficient[7,4]
    P_CGGA <- summary(fit_CGGA)$coefficient[7,5]

    result <- c(Gene1=Gene1,Gene2=Gene2,Z_TCGA=Z_TCGA,P_TCGA=P_TCGA,
                Z_CGGA=Z_CGGA,P_CGGA=P_CGGA)
  })
)
stopCluster(cl)
LGG_PanGene_Inter <- as.data.frame(t(paraRes))
LGG_PanGene_Inter$Gene1 <- as.character(LGG_PanGene_Inter$Gene1)
LGG_PanGene_Inter$Gene2 <- as.character(LGG_PanGene_Inter$Gene2)
LGG_PanGene_Inter[,3] <- as.numeric(as.character(LGG_PanGene_Inter[,3]))
LGG_PanGene_Inter[,4] <- as.numeric(as.character(LGG_PanGene_Inter[,4]))
LGG_PanGene_Inter[,5] <- as.numeric(as.character(LGG_PanGene_Inter[,5]))
LGG_PanGene_Inter[,6] <- as.numeric(as.character(LGG_PanGene_Inter[,6]))

LGG_PanGene_Inter$FDR_TCGA <- p.adjust(LGG_PanGene_Inter$P_TCGA,method="fdr")

save(LGG_PanGene_Main,LGG_PanGene_Inter,file='LGG_PanGene_Main_Inter.rdata')


LGG_PanGene_Main1 <- LGG_PanGene_Main[LGG_PanGene_Main$FDR_TCGA <= 0.05 & LGG_PanGene_Main$P_CGGA <= 0.05
                                      & LGG_PanGene_Main$Z_TCGA*LGG_PanGene_Main$Z_CGGA >0,]
LGG_PanGene_Inter1 <- LGG_PanGene_Inter[LGG_PanGene_Inter$FDR_TCGA <= 0.05 & LGG_PanGene_Inter$P_CGGA <= 0.05
                                        & LGG_PanGene_Inter$Z_TCGA*LGG_PanGene_Inter$Z_CGGA >0,]

TCGA_main_FPKM <- TCGA[,c('OS','Censor',LGG_PanGene_Main1$Gene)]

write.table(TCGA_main_FPKM,col.names=T,row.names=F,quote=F,sep='\t',file="TCGA_main_FPKM.txt")
