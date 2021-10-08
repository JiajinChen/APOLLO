rm(list=ls())
library(forestplot)
data <- read.csv("Subgroup.csv",stringsAsFactors=FALSE,header = FALSE)


forestplot(data[,c(1:2)],data$V9,data$V10,data$V11,
           graphwidth=unit(35,"mm"),
           col =fpColors(box = "brown4",lines = "grey0",zero = "grey0"),
           colgap=unit(7,"mm"),
           lineheight = "auto",
           txt_gp = fpTxtGp(ticks = gpar(height=0.2,cex=0.8),
                            summary = gpar(cex=1.0),cex = 0.9),
           zero=1,
           clip=c(0.8,1.1),
           xticks=seq(1,9.3,by=1)
)

forestplot(labeltext=rep("",25),data$V3,data$V4,data$V5,
           #graph.pos="left",
           graphwidth=unit(35,"mm"),
           col =fpColors(box = c("brown4","red"),lines = "grey0",zero = "grey0"),
           colgap=unit(3,"mm"),
           lineheight = "auto",
           txt_gp = fpTxtGp(ticks = gpar(height=0.2,cex=0.8),
                            summary = gpar(cex=1.0),cex = 0.9),
           zero=0.7,
           clip=c(0.7,1.1),
           xticks=seq(0.7,1,by=0.1)
)

forestplot(labeltext=rep("",25),data$V6,data$V7,data$V8,
           graphwidth=unit(35,"mm"),
           col =fpColors(box = "brown4",lines = "grey0",zero = "grey0"),
           colgap=unit(3,"mm"),#xlog = TRUE,
           lineheight = "auto",
           txt_gp = fpTxtGp(ticks = gpar(height=0.2,cex=0.8),
                            summary = gpar(cex=1.0),cex = 0.9),
           zero=0.7,
           clip=c(0.7,1.1),
           xticks=seq(0.7,1,by=0.1)
)
