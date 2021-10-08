# Meta-analysis of AUC
# 3-YEAR
AUC_3 <- data.frame(Study='TCGA',N=505,AUC=0.933,lci=0.897,uci=0.970,stringsAsFactors = F)
AUC_3 <- rbind(AUC_3,c(Study='CGGA1',N=408,AUC=0.888,lci=0.847,uci=0.928))
AUC_3 <- rbind(AUC_3,c(Study='CGGA2',N=143,AUC=0.898,lci=0.806,uci=0.991))
AUC_3 <- rbind(AUC_3,c(Study='Rembrandt',N=121,AUC=0.893,lci=0.837,uci=0.950))
AUC_3 <- rbind(AUC_3,c(Study='Weller',N=137,AUC=0.844,lci=0.752,uci=0.936))
AUC_3 <- rbind(AUC_3,c(Study='Gravendeel',N=106,AUC=0.861,lci=0.779,uci=0.943))
AUC_3$N <- as.numeric(AUC_3$N)
AUC_3$AUC <- as.numeric(AUC_3$AUC)
AUC_3$lci <- as.numeric(AUC_3$lci)
AUC_3$uci <- as.numeric(AUC_3$uci)
AUC_3$CI_char <- paste0("(",AUC_3$lci,', ',AUC_3$uci,')')

# 5-YEAR
AUC_5 <- data.frame(Study='TCGA',N=505,AUC=0.854,lci=0.788,uci=0.921,stringsAsFactors = F)
AUC_5 <- rbind(AUC_5,c(Study='CGGA1',N=408,AUC=0.851,lci=0.806,uci=0.896))
AUC_5 <- rbind(AUC_5,c(Study='CGGA2',N=143,AUC=0.896,lci=0.814,uci=0.978))
AUC_5 <- rbind(AUC_5,c(Study='Rembrandt',N=121,AUC=0.817,lci=0.741,uci=0.892))
AUC_5 <- rbind(AUC_5,c(Study='Weller',N=137,AUC=0.806,lci=0.717,uci=0.895))
AUC_5 <- rbind(AUC_5,c(Study='Gravendeel',N=106,AUC=0.790,lci=0.687,uci=0.893))
AUC_5$N <- as.numeric(AUC_5$N)
AUC_5$AUC <- as.numeric(AUC_5$AUC)
AUC_5$lci <- as.numeric(AUC_5$lci)
AUC_5$uci <- as.numeric(AUC_5$uci)
AUC_5$CI_char <- paste0("(",AUC_5$lci,', ',AUC_5$uci,')')



library(meta)
m3 = metagen(TE = AUC_3$AUC,lower = AUC_3$lci,upper = AUC_3$uci,n.e=AUC_3$N)
m5 = metagen(TE = AUC_5$AUC,lower = AUC_5$lci,upper = AUC_5$uci,n.e=AUC_5$N)

pdf("meta_AUC3.pdf")
forest(m3,studlab = AUC_3$Study,comb.fixed = T,comb.random = F,ref=0.5,
       xlim=c(0.5,1),digits = 3,
       leftcols=c("studlab",'n.e',"effect","ci"),leftlabs = c("Study","N","AUC","95% CI"),
       rightcols = F,
       print.tau2 = F,col.diamond.fixed = "#B70031",col.diamond.lines.fixed="#B70031",col.square="black" )
dev.off()

pdf("meta_AUC5.pdf")
forest(m5,studlab = AUC_5$Study,comb.fixed = T,comb.random = F,ref=0.5,
       xlim=c(0.5,1),digits = 3,
       leftcols=c("studlab",'n.e',"effect","ci"),leftlabs = c("Study","N","AUC","95% CI"),
       rightcols = F,
       print.tau2 = F,col.diamond.fixed = "#B70031",col.diamond.lines.fixed="#B70031",col.square="black" )
dev.off()
