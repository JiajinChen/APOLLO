# Meta-analysis of AUC
# 3-YEAR
AUC_3 <- data.frame(Study='TCGA',N=505,AUC=0.879,lci=0.831,uci=0.928,stringsAsFactors = F)
AUC_3 <- rbind(AUC_3,c(Study='CGGA1',N=408,AUC=0.832,lci=0.779,uci=0.885))
AUC_3 <- rbind(AUC_3,c(Study='CGGA2',N=143,AUC=0.818,lci=0.679,uci=0.956))
AUC_3 <- rbind(AUC_3,c(Study='Rembrandt',N=121,AUC=0.875,lci=0.813,uci=0.938))
AUC_3 <- rbind(AUC_3,c(Study='Weller',N=137,AUC=0.841,lci=0.735,uci=0.948))
AUC_3 <- rbind(AUC_3,c(Study='Gravendeel',N=106,AUC=0.817,lci=0.718,uci=0.915))
AUC_3$N <- as.numeric(AUC_3$N)
AUC_3$AUC <- as.numeric(AUC_3$AUC)
AUC_3$lci <- as.numeric(AUC_3$lci)
AUC_3$uci <- as.numeric(AUC_3$uci)
AUC_3$CI_char <- paste0("(",AUC_3$lci,', ',AUC_3$uci,')')

# 5-YEAR
AUC_5 <- data.frame(Study='TCGA',N=505,AUC=0.817,lci=0.718,uci=0.915,stringsAsFactors = F)
AUC_5 <- rbind(AUC_5,c(Study='CGGA1',N=408,AUC=0.797,lci=0.743,uci=0.851))
AUC_5 <- rbind(AUC_5,c(Study='CGGA2',N=143,AUC=0.789,lci=0.662,uci=0.916))
AUC_5 <- rbind(AUC_5,c(Study='Rembrandt',N=121,AUC=0.834,lci=0.762,uci=0.907))
AUC_5 <- rbind(AUC_5,c(Study='Weller',N=137,AUC=0.780,lci=0.683,uci=0.876))
AUC_5 <- rbind(AUC_5,c(Study='Gravendeel',N=106,AUC=0.710,lci=0.596,uci=0.824))
AUC_5$N <- as.numeric(AUC_5$N)
AUC_5$AUC <- as.numeric(AUC_5$AUC)
AUC_5$lci <- as.numeric(AUC_5$lci)
AUC_5$uci <- as.numeric(AUC_5$uci)
AUC_5$CI_char <- paste0("(",AUC_5$lci,', ',AUC_5$uci,')')



library(meta)
m3 = metagen(TE = AUC_3$AUC,lower = AUC_3$lci,upper = AUC_3$uci,n.e=AUC_3$N)
m5 = metagen(TE = AUC_5$AUC,lower = AUC_5$lci,upper = AUC_5$uci,n.e=AUC_5$N)

pdf("Meta_AUC.Clin3_CIRpackage.pdf")
forest(m3,studlab = AUC_3$Study,comb.fixed = T,comb.random = F,ref=0.5,
       xlim=c(0.5,1),digits = 3,
       leftcols=c("studlab",'n.e',"effect","ci"),leftlabs = c("Study","N","AUC","95% CI"),
       rightcols = F,
       print.tau2 = F,col.diamond.fixed = "#B70031",col.diamond.lines.fixed="#B70031",col.square="black" )
dev.off()

pdf("Meta_AUC.Clin5_CIRpackage.pdf")
forest(m5,studlab = AUC_5$Study,comb.fixed = T,comb.random = F,ref=0.5,
       xlim=c(0.5,1),digits = 3,
       leftcols=c("studlab",'n.e',"effect","ci"),leftlabs = c("Study","N","AUC","95% CI"),
       rightcols = F,
       print.tau2 = F,col.diamond.fixed = "#B70031",col.diamond.lines.fixed="#B70031",col.square="black" )
dev.off()
