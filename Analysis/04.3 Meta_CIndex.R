library(meta)
CIndex <- data.frame(Study='TCGA',N=505,CIndex=0.874,lci=0.844,uci=0.904,stringsAsFactors = F)
CIndex <- rbind(CIndex,c(Study='CGGA1',N=408,CIndex=0.804,lci=0.769,uci=0.839))
CIndex <- rbind(CIndex,c(Study='CGGA2',N=143,CIndex=0.807,lci=0.731,uci=0.884))
CIndex <- rbind(CIndex,c(Study='Rembrandt',N=121,CIndex=0.772,lci=0.729,uci=0.815))
CIndex <- rbind(CIndex,c(Study='Weller',N=137,CIndex=0.787,lci=0.718,uci=0.856))
CIndex <- rbind(CIndex,c(Study='Gravendeel',N=106,CIndex=0.759,lci=0.704,uci=0.815))
CIndex$N <- as.numeric(CIndex$N)
CIndex$CIndex <- as.numeric(CIndex$CIndex)
CIndex$lci <- as.numeric(CIndex$lci)
CIndex$uci <- as.numeric(CIndex$uci)
CIndex$CI_char <- paste0("(",CIndex$lci,', ',CIndex$uci,')')

library(meta)
m3 = metagen(TE = CIndex$CIndex,lower = CIndex$lci,upper = CIndex$uci,n.e=CIndex$N)

pdf("meta_CIndex.pdf")
forest(m3,studlab = CIndex$Study,comb.fixed = F,comb.random = T,ref=0.5,
       xlim=c(0.5,1),digits = 3,
       leftcols=c("studlab",'n.e',"effect","ci"),leftlabs = c("Study","N","AUC","95% CI"),
       rightcols = F,
       print.tau2 = F,col.diamond.random = "#B70031",col.diamond.lines.random = "#B70031",col.square="black" )
dev.off()
