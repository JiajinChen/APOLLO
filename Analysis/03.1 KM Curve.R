rm(list=ls())
TCGA$ds <- 'TCGA'; CGGA1$ds <- 'CGGA-seq'; CGGA2$ds <- 'CGGA-array'
Rembrandt$ds <- 'Rembrandt'; GSE61374$ds <- 'GSE61374'; GSE16011$ds <- 'GSE16011'
colnames(GSE16011)[1] <- 'UID';colnames(GSE61374)[1] <- 'UID';colnames(CGGA1)[1] <- c("UID")

All_data <- merge(TCGA,CGGA1,all = T)
All_data <- merge(All_data,CGGA2,all = T)
All_data <- merge(All_data,Rembrandt,all = T)
All_data <- merge(All_data,GSE16011,all = T)
All_data <- merge(All_data,GSE61374,all = T)

##6 KM
dev.off()
library(survminer)
library(survival)
library(cowplot)
data<-TCGA
data$Group_SCore<-ifelse(is.na(data$APOLLO),NA,
                         ifelse(data$APOLLO<(0.6944875),"Low-risk","High-risk"))

data$Group_SCore<-as.factor(data$Group_SCore)
data$Group_SCore<-relevel(data$Group_SCore,ref = "Low-risk")
summary(coxph(Surv(OS,Censor)~as.factor(Group_SCore),data=data))
fit<-survfit(Surv(OS,Censor)~Group_SCore,data=data)
a<-ggsurvplot(fit,data = data,title="(A) TCGA",
              
              xlim=c(0,120),
              ylab="Survival probability",
              xlab="Overall Survival Months",#censor=FALSE,
              font.x=12,
              font.y=12,
              legend.title="",
              legend.labs=c("Low-risk","High-risk"),
              legend=c(0.2,0.2),
              #risk.table = TRUE,
              break.x.by = 12,
              break.y.by = 0.2,
              palette = c("#31789F","#9C2351"),
              # risk.table.title=NULL,
              #risk.table.pos="out",
              #tables.height=0.15,
              ggtheme = theme(#text = element_text(family = "Arial"),
                plot.margin = unit(rep(0,4),"cm"),
                panel.spacing = unit(rep(0,4),"cm"),
                panel.border=element_blank(),
                axis.line=element_line(size=0.5),
                legend.background = element_blank(),
                legend.text =element_text(size=11),
                plot.background = element_blank(),
                panel.background = element_blank(),
                legend.box.background=element_blank(),
                legend.key = element_blank(),
                axis.text.y = element_text(color="black"),
                axis.text.x = element_text(color="black"),
                #axis.title.y = element_text(margin=unit(c(0.1,0.1,0,0.1),"cm"),size = 10),
                axis.title.x = element_text(margin=unit(c(0.1,0.1,0,0.1),"cm"))),
              tables.theme  = theme(plot.margin = unit(rep(0,4),"cm"),
                                    panel.spacing = unit(rep(0,4),"cm"),
                                    title = element_blank(),
                                    panel.border=element_blank(),
                                    legend.key = element_blank(),
                                    axis.line=element_blank(),
                                    axis.ticks = element_blank(),
                                    axis.text.y = element_blank(),
                                    axis.text.x = element_blank(),
                                    legend.background = element_blank(),
                                    plot.background = element_blank(),
                                    panel.background = element_blank(),
                                    # axis.text.y = element_text(color="black"),
                                    legend.box.background=element_blank(),
                                    axis.title.x = element_blank()
              ),
              fontsize=4,
              tables.y.text=T,
              tables.y.text.col=TRUE
)
a$plot<-a$plot+
  ggplot2::annotate("text", x =82, y = 1,
                    label = expression(paste('HR = 8.51 ',', ',  italic(P),' = 2.14×10-16')),colour="black", size = 3.5,parse = TRUE)



data<-CGGA1
data$Group_SCore<-ifelse(is.na(data$APOLLO),NA,
                         ifelse(data$APOLLO<(0.6944875),"Low-risk","High-risk"))
data$Group_SCore<-as.factor(data$Group_SCore)
data$Group_SCore<-relevel(data$Group_SCore,ref = "Low-risk")
summary(coxph(Surv(OS,Censor)~as.factor(Group_SCore),data=data))
fit<-survfit(Surv(OS,Censor)~Group_SCore,data=data)
b<-ggsurvplot(fit,data = data,title="(B) CGGA1",
              xlim=c(0,120),
              ylab="Survival probability",
              cex=5,
              xlab="Overall Survival Months",#censor=FALSE,
              font.x=12,
              font.y=12,
              legend.title="",
              legend.labs=c("Low-risk","High-risk"),
              legend=c(0.2,0.2),
              #risk.table = TRUE,
              break.x.by = 12,
              break.y.by = 0.2,
              palette = c("#31789F","#9C2351"),
              # risk.table.title=NULL,
              #risk.table.pos="out",
              #tables.height=0.15,
              ggtheme = theme(plot.margin = unit(rep(0,4),"cm"),
                              panel.spacing = unit(rep(0,4),"cm"),
                              panel.border=element_blank(),
                              axis.line=element_line(size=0.5),
                              legend.background = element_blank(),
                              plot.background = element_blank(),
                              panel.background = element_blank(),
                              legend.box.background=element_blank(),
                              legend.text =element_text(size=11),
                              legend.key = element_blank(),
                              axis.text.y = element_text(color="black"),
                              axis.text.x = element_text(color="black"),
                              #axis.title.y = element_text(size = 10),
                              axis.title.x = element_text(margin=unit(c(0.1,0.1,0,0.1),"cm"))),
              tables.theme  = theme(plot.margin = unit(rep(0,4),"cm"),
                                    panel.spacing = unit(rep(0,4),"cm"),
                                    title = element_blank(),
                                    panel.border=element_blank(),
                                    legend.key = element_blank(),
                                    axis.line=element_blank(),
                                    axis.ticks = element_blank(),
                                    axis.text.y = element_blank(),
                                    axis.text.x = element_blank(),
                                    legend.background = element_blank(),
                                    plot.background = element_blank(),
                                    panel.background = element_blank(),
                                    # axis.text.y = element_text(color="black"),
                                    legend.box.background=element_blank(),
                                    axis.title.x = element_blank()
              ),
              fontsize=4,
              tables.y.text=T,
              tables.y.text.col=TRUE
)
b$plot<-b$plot+
  ggplot2::annotate("text", x =80, y = 1,
                    label = expression(paste('HR = 4.86 ',', ',  italic(P),' = 1.75×10-14')),colour="black", size = 3.5,parse = TRUE)


data<-CGGA2

data$Group_SCore<-ifelse(is.na(data$APOLLO),NA,
                         ifelse(data$APOLLO<(0.6944875),"Low-risk","High-risk"))
data$Group_SCore<-as.factor(data$Group_SCore)
data$Group_SCore<-relevel(data$Group_SCore,ref = "Low-risk")
summary(coxph(Surv(OS,Censor)~as.factor(Group_SCore),data=data))
fit<-survfit(Surv(OS,Censor)~Group_SCore,data=data)
c<-ggsurvplot(fit,data = data,title="(C) CGGA2",
              xlim=c(0,120),
              ylab="Survival probability",
              xlab="Overall Survival Months",#censor=FALSE,
              font.x=12,
              font.y=12,
              legend.title="",
              legend.labs=c("Low-risk","High-risk"),
              legend=c(0.2,0.2),
              break.x.by = 12,
              break.y.by=0.2,
              #risk.table = TRUE,
              palette = c("#31789F","#9C2351"),
              # risk.table.title=NULL,
              #risk.table.pos="out",
              #tables.height=0.15,
              ggtheme = theme(plot.margin = unit(rep(0,4),"cm"),
                              panel.spacing = unit(rep(0,4),"cm"),
                              panel.border=element_blank(),
                              axis.line=element_line(size=0.5),
                              legend.background = element_blank(),
                              legend.text =element_text(size=11),
                              plot.background = element_blank(),
                              panel.background = element_blank(),
                              legend.box.background=element_blank(),
                              
                              legend.key = element_blank(),
                              axis.text.y = element_text(color="black"),
                              axis.text.x = element_text(color="black"),
                              #axis.title.y = element_text(size = 10),
                              axis.title.x = element_text(margin=unit(c(0.1,0.1,0,0.1),"cm"))),
              tables.theme  = theme(plot.margin = unit(rep(0,4),"cm"),
                                    panel.spacing = unit(rep(0,4),"cm"),
                                    title = element_blank(),
                                    panel.border=element_blank(),
                                    legend.key = element_blank(),
                                    axis.line=element_blank(),
                                    axis.ticks = element_blank(),
                                    axis.text.y = element_blank(),
                                    axis.text.x = element_blank(),
                                    legend.background = element_blank(),
                                    plot.background = element_blank(),
                                    panel.background = element_blank(),
                                    # axis.text.y = element_text(color="black"),
                                    legend.box.background=element_blank(),
                                    axis.title.x = element_blank()
              ),
              fontsize=4,
              tables.y.text=T,
              tables.y.text.col=TRUE
)
c$plot<-c$plot+
  ggplot2::annotate("text", x =80, y = 1,
                    label = expression(paste('HR = 6.26 ',', ',  italic(P),' = 4.41×10-06')),colour="black", size = 3.5,parse = TRUE)


data<-Rembrandt
data$Group_SCore<-ifelse(is.na(data$APOLLO),NA,
                         ifelse(data$APOLLO<(0.6944875),"Low-risk","High-risk"))
data$Group_SCore<-as.factor(data$Group_SCore)
data$Group_SCore<-relevel(data$Group_SCore,ref = "Low-risk")
summary(coxph(Surv(OS,Censor)~as.factor(Group_SCore),data=data))
fit<-survfit(Surv(OS,Censor)~Group_SCore,data=data)
d<-ggsurvplot(fit,data = data,title="(D) Rembrandt",
              xlim=c(0,100),
              ylab="Survival probability",
              xlab="Overall Survival Months",#censor=FALSE,
              font.x=12,
              font.y=12,
              legend.title="",
              legend.labs=c("Low-risk","High-risk"),
              legend=c(0.2,0.2),
              #risk.table = TRUE,
              break.x.by = 12,
              break.y.by=0.2,
              palette = c("#31789F","#9C2351"),
              # risk.table.title=NULL,
              #risk.table.pos="out",
              #tables.height=0.15,
              ggtheme = theme(plot.margin = unit(rep(0,4),"cm"),
                              panel.spacing = unit(rep(0,4),"cm"),
                              panel.border=element_blank(),
                              axis.line=element_line(size=0.5),
                              legend.background = element_blank(),
                              legend.text =element_text(size=11),
                              plot.background = element_blank(),
                              panel.background = element_blank(),
                              legend.box.background=element_blank(),
                              
                              legend.key = element_blank(),
                              axis.text.y = element_text(color="black"),
                              axis.text.x = element_text(color="black"),
                              #axis.title.y = element_text(size = 10),
                              axis.title.x = element_text(margin=unit(c(0.1,0.1,0,0.1),"cm"))),
              tables.theme  = theme(plot.margin = unit(rep(0,4),"cm"),
                                    panel.spacing = unit(rep(0,4),"cm"),
                                    title = element_blank(),
                                    panel.border=element_blank(),
                                    legend.key = element_blank(),
                                    axis.line=element_blank(),
                                    axis.ticks = element_blank(),
                                    axis.text.y = element_blank(),
                                    axis.text.x = element_blank(),
                                    legend.background = element_blank(),
                                    plot.background = element_blank(),
                                    panel.background = element_blank(),
                                    # axis.text.y = element_text(color="black"),
                                    legend.box.background=element_blank(),
                                    axis.title.x = element_blank()
              ),
              fontsize=4,
              tables.y.text=T,
              tables.y.text.col=TRUE
)
d$plot<-d$plot+
  ggplot2::annotate("text", x =70, y = 1,
                    label = expression(paste('HR = 3.49 ',', ',  italic(P),' = 3.32×10-06')),colour="black", size = 3.5,parse = TRUE)



data<-GSE61374
data$Group_SCore<-ifelse(is.na(data$APOLLO),NA,
                         ifelse(data$APOLLO<(0.6944875),"Low-risk","High-risk"))
data$Group_SCore<-as.factor(data$Group_SCore)
data$Group_SCore<-relevel(data$Group_SCore,ref = "Low-risk")
summary(coxph(Surv(OS,Censor)~as.factor(Group_SCore),data=data))
fit<-survfit(Surv(OS,Censor)~Group_SCore,data=data)
e<-ggsurvplot(fit,data = data,title="(E) Weller",
              xlim=c(0,120),
              #pval = "HR = 4.05 , P = 2.55×10-12",
              #pval.coord=c(100,1),
              #ggtheme=theme_survminer(),  
              ylab="Survival probability",
              xlab="Overall Survival Months",#censor=FALSE,
              font.x=12,
              font.y=12,
              legend.title="",
              legend.labs=c("Low-risk","High-risk"),
              legend=c(0.2,0.2),
              break.x.by = 12,
              break.y.by =0.2,
              #risk.table = TRUE,
              palette = c("#31789F","#9C2351"),
              # risk.table.title=NULL,
              #risk.table.pos="out",
              #tables.height=0.15,
              ggtheme = theme(plot.margin = unit(rep(0,4),"cm"),
                              panel.spacing = unit(rep(0,4),"cm"),
                              panel.border=element_blank(),
                              axis.line=element_line(size=0.5),
                              legend.background = element_blank(),
                              legend.text =element_text(size=11),
                              plot.background = element_blank(),
                              panel.background = element_blank(),
                              legend.box.background=element_blank(),
                              
                              legend.key = element_blank(),
                              axis.text.y = element_text(color="black"),
                              axis.text.x = element_text(color="black"),
                              #axis.title.y = element_text(size = 10),
                              axis.title.x = element_text(margin=unit(c(0.1,0.1,0,0.1),"cm"))),
              tables.theme  = theme(plot.margin = unit(rep(0,4),"cm"),
                                    panel.spacing = unit(rep(0,4),"cm"),
                                    title = element_blank(),
                                    panel.border=element_blank(),
                                    legend.key = element_blank(),
                                    axis.line=element_blank(),
                                    axis.ticks = element_blank(),
                                    axis.text.y = element_blank(),
                                    axis.text.x = element_blank(),
                                    legend.background = element_blank(),
                                    plot.background = element_blank(),
                                    panel.background = element_blank(),
                                    # axis.text.y = element_text(color="black"),
                                    legend.box.background=element_blank(),
                                    axis.title.x = element_blank()
              ),
              fontsize=4,
              tables.y.text=T,
              tables.y.text.col=TRUE
)
e$plot<-e$plot+
  ggplot2::annotate("text", x =80, y = 1,
                    label = expression(paste('HR = 3.41 ',', ',  italic(P),' = 3.99×10-04')),colour="black", size = 3.5,parse = TRUE)



data<-GSE16011
data$Group_SCore<-ifelse(is.na(data$APOLLO),NA,
                         ifelse(data$APOLLO<(0.6944875),"Low-risk","High-risk"))
data$Group_SCore<-as.factor(data$Group_SCore)
data$Group_SCore<-relevel(data$Group_SCore,ref = "Low-risk")
summary(coxph(Surv(OS,Censor)~as.factor(Group_SCore),data=data))
fit<-survfit(Surv(OS,Censor)~Group_SCore,data=data)
f<-ggsurvplot(fit,data = data,title="(F) Gravendeel",
              xlim=c(0,120),
              #pval = "HR = 4.05 , P = 2.55×10-12",
              #pval.coord=c(100,1),
              #ggtheme=theme_survminer(),  
              ylab="Survival probability",
              xlab="Overall Survival Months",#censor=FALSE,
              font.x=12,
              font.y=12,
              legend.title="",
              legend.labs=c("Low-risk","High-risk"),
              legend=c(0.2,0.2),
              break.x.by = 12,
              break.y.by = 0.2,
              #risk.table = TRUE,
              palette = c("#31789F","#9C2351"),
              # risk.table.title=NULL,
              #risk.table.pos="out",
              #tables.height=0.15,
              ggtheme = theme(plot.margin = unit(rep(0,4),"cm"),
                              panel.spacing = unit(rep(0,4),"cm"),
                              panel.border=element_blank(),
                              axis.line=element_line(size=0.5),
                              legend.background = element_blank(),
                              legend.text =element_text(size=11),
                              plot.background = element_blank(),
                              panel.background = element_blank(),
                              legend.box.background=element_blank(),
                              
                              legend.key = element_blank(),
                              axis.text.y = element_text(color="black"),
                              axis.text.x = element_text(color="black"),
                              #axis.title.y = element_text(size = 10),
                              axis.title.x = element_text(margin=unit(c(0.1,0.1,0,0.1),"cm"))),
              tables.theme  = theme(plot.margin = unit(rep(0,4),"cm"),
                                    panel.spacing = unit(rep(0,4),"cm"),
                                    title = element_blank(),
                                    panel.border=element_blank(),
                                    legend.key = element_blank(),
                                    axis.line=element_blank(),
                                    axis.ticks = element_blank(),
                                    axis.text.y = element_blank(),
                                    axis.text.x = element_blank(),
                                    legend.background = element_blank(),
                                    plot.background = element_blank(),
                                    panel.background = element_blank(),
                                    # axis.text.y = element_text(color="black"),
                                    legend.box.background=element_blank(),
                                    axis.title.x = element_blank()
              ),
              fontsize=4,
              tables.y.text=T,
              tables.y.text.col=TRUE
)
f$plot<-f$plot+
  ggplot2::annotate("text", x =80, y = 1,
                    label = expression(paste('HR = 2.19 ',', ',  italic(P),' = 2.88×10-03')),colour="black", size = 3.5,parse = TRUE)



ggdraw()+draw_plot(a$plot,0.01,0.5,0.32,0.45)+draw_plot(b$plot,0.333,0.5,0.32,0.45)+draw_plot(c$plot,0.666,0.5,0.32,0.45)+
  draw_plot(d$plot,0.01,0.05,0.32,0.45)+draw_plot(e$plot,0.333,0.05,0.32,0.45)+draw_plot(f$plot,0.666,0.05,0.32,0.45)

##6条KM
sigDs<-All_data
sigDs$ds_adj1<-ifelse(sigDs$ds %in% 'CGGA-array',1,ifelse(is.na(sigDs$ds),NA,0))
sigDs$ds_adj2<-ifelse(sigDs$ds %in% 'CGGA-seq',1,ifelse(is.na(sigDs$ds),NA,0))
sigDs$ds_adj3<-ifelse(sigDs$ds %in% 'GSE16011',1,ifelse(is.na(sigDs$ds),NA,0))
sigDs$ds_adj4<-ifelse(sigDs$ds %in% 'GSE61374',1,ifelse(is.na(sigDs$ds),NA,0))
sigDs$ds_adj5<-ifelse(sigDs$ds %in% 'Rembrandt',1,ifelse(is.na(sigDs$ds),NA,0))
sigDs$ds_adj6<-ifelse(sigDs$ds %in% 'TCGA',1,ifelse(is.na(sigDs$ds),NA,0))
sigDs$dose6 <- cut(sigDs$APOLLO,breaks = quantile(sigDs$APOLLO,probs=c(0,0.2,0.4,0.6,0.8,0.9,1),na.rm = T),labels = 1:6)
plotFit <- survfit(coxph(Surv(OS,Censor)~strata(dose6)+sigDs$ds_adj1+sigDs$ds_adj2+sigDs$ds_adj3+sigDs$ds_adj4+sigDs$ds_adj5
                         ,data = sigDs))
sur<-data.frame(plotFit[["time"]],plotFit[["surv"]])

library(survminer)
library(RColorBrewer)
help(RColorBrewer)
color <- c("#31789F","#9C2351","#DBA541","#8FBC8F","#00CED1","#7B68EE")
p <- ggsurvplot(plotFit,data = sigDs,title="Discriminative ability of APOLLO",
                #xlim=c(0,120),
                #pval = "HR = 4.05 , P = 2.55×10-12",
                #pval.coord=c(100,1),
                #ggtheme=theme_survminer(),  
                ylab="3-survival probability",
                xlab="Overall Survival Months",censor=FALSE,
                font.x=12,
                font.y=12,
                legend.title="",
                legend.labs=c("Level 1","Level 2","Level 3","Level 4","Level 5","Level 6"),
                legend="bottom",
                break.x.by = 12,
                #break.y.by = c(0.33,0.50,0.70,0.86,0.93,0.94),
                surv.median.line="v",
                surv.plot.height=0.8,
                palette = color,
                ggtheme = theme(plot.title = element_text(hjust = 0.5),
                                plot.margin = unit(rep(0,4),"cm"),
                                panel.spacing = unit(rep(0,4),"cm"),
                                panel.border=element_rect(color="black",fill=NA),
                                axis.line=element_line(size=0.5),
                                legend.background = element_blank(),
                                legend.text =element_text(size=9),
                                plot.background = element_blank(),
                                panel.background = element_blank(),
                                legend.box.background=element_blank(),
                                #legend.box = TRUE,
                                #legend.box.spacing = unit(rep(0,4),"cm"),
                                legend.key = element_blank(),
                                axis.text.y = element_text(color=c('black',rev(color)[1:2],'black',rev(color)[3:6]),size = 9),
                                axis.text.x = element_text(color="black",size = 9),
                                #axis.title.y = element_text(size = 10),
                                axis.title.x = element_text(margin=unit(c(0.1,0.1,0,0.1),"cm"))))
p$plot<-p$plot+
  ggplot2::scale_y_continuous(expand = c(0.01,0.01),sec.axis = sec_axis(~.,name="5-survival probability",breaks = c(0,0.016,0.204,0.50,0.521,0.723,0.801,0.945),labels = c(0,1.6,20.4,50,52.1,72.3,80.1,94.5)),breaks = c(0,0.079,0.331,0.50,0.715,0.885,0.952,0.986),labels = c(0,7.9,33.1,50,71.5,88.5,95.2,98.6))+
  #ggplot2::scale_y_continuous(expand = c(0.01,0.01),sec.axis = sec_axis(~.,name="5-survival probability",breaks = c(0.23,0.50,0.54,0.68,0.77,0.86),labels = c(0.23,0.50,0.54,0.68,0.77,0.86)))+    
  
  ggplot2::scale_x_continuous(expand = c(0.01,0.01))+
  #ggplot2::geom_abline(intercept=c(36,0.863),linetype=2,size=0.5)+
  #ggplot2::geom_vline(xintercept = c(36,60),linetype=2,size=0.5)+
  ggplot2::geom_hline(yintercept=c(0.5),linetype=2,size=0.5)
p

plotFit <- coxph(Surv(OS,Censor)~as.factor(dose6)+as.factor(ds),data = sigDs)
res_dose <- data.frame(summary(plotFit)$conf.int[1:5,c(1,3:4)])
colnames(res_dose) <- c('HR','LCI','UCI')
res_dose <- rbind(data.frame(HR=1,LCI=1,UCI=1),res_dose)
res_dose[6,3] <- 65
Median_sur <- c(15.7,25.3,61.7,92.0,136.1,192.6)
res_dose <- cbind(Median_sur=rev(Median_sur),res_dose)
PRS_plot <- ggplot(aes(x=Median_sur,y=HR),data=res_dose)+
  #xlab("PRS (%)")+
  #ylab("Hazard Ratio")+
  scale_x_continuous(position="top",breaks=c(15.7,35,61.1,92.0,136.2,192.6),
                     labels=Median_sur,limits=c(0,200))+
  scale_y_continuous(breaks=c(0,10,20,30,40,50,60),trans="reverse",limits=c(70,0))+#expand=c(0.08,0.08))+
  labs(x='Median survival months',y="Hazard Ratio")+
  geom_point(shape='circle',size=1.5,colour=color)+
  geom_errorbar(aes(ymin=LCI,ymax=UCI),width=0,size=0.75,colour=color)+
  geom_text(aes(x=Median_sur+12,y=HR-0.3,label=round(HR,2)),size=3)+
  theme(plot.margin = unit(rep(0,4),"cm"),
        panel.spacing = unit(rep(0,4),"cm"),
        panel.border=element_blank(),
        axis.line=element_line(size=0.5),
        legend.background = element_blank(),
        legend.text =element_text(size=11),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.box.background=element_blank(),
        legend.key = element_blank(),
        axis.text.y = element_text(color="black",size = 9),
        axis.text.x = element_text(color="black",size = 9),
        axis.title.y = element_text(margin=unit(c(0.1,0,0,0.1),"cm"),size = 14,face = "plain"),
        axis.title.x = element_text(margin=unit(c(0.1,0.1,0,0.1),"cm"),size = 14,face = "plain"))

##整合所有图形
ggdraw()+draw_plot(a$plot,0.01,0.5,0.23,0.45)+
  draw_plot(b$plot,0.242,0.5,0.23,0.45)+
  draw_plot(c$plot,0.474,0.5,0.23,0.45)+
  draw_plot(d$plot,0.01,0.02,0.23,0.45)+
  draw_plot(e$plot,0.242,0.02,0.23,0.45)+
  draw_plot(f$plot,0.474,0.02,0.23,0.45)+
  draw_plot(p$plot,0.706,0.5,0.28,0.45)+
  draw_plot(PRS_plot,0.706,0.02,0.28,0.45)
