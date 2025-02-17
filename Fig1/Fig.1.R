wd="F:\\文章\\果皮微生物风土\\投稿3\\revision\\data and code\\Fig1\\"

###### Fig 1

setwd(wd)
library(sciplot)
library(doBy)
library(plotrix)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(tidyr)
library(reshape2)  
library(Rmisc)


don=read.csv("climate_soil_data.csv",check.names = F,row.names = 1)
colnames(don)

Var=don[,3:33]
Var = na.omit(Var)

df_pca <- prcomp(Var,scale. = TRUE) 
df_pcs <-data.frame(df_pca$x, Site=don$vineyard,Variety=don$Variety)

df_pcs$Site = factor(df_pcs$Site,levels=c("HD","YSGD","CHFC","SSB","XG","HSP"))
df_pcs$Variety = factor(df_pcs$Variety,levels = c('CS','ME','CH'))
##1-2
percentage<-round(df_pca$sdev^2 / sum(df_pca$sdev^2) * 100,2)
percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))

p1 <-ggplot(df_pcs,aes(x=PC1,y=PC2,color =Site,shape=Variety)) + 
  geom_point(size=7) + 
  xlab(percentage[1]) + ylab(percentage[2]) +
  theme_bw(base_size = 15) +
  scale_color_manual(values=c ( HD='#20a391', 
                                YSGD= 'dodgerblue4', 
                                CHFC= '#a00528', 
                                SSB= '#fc7839', 
                                XG= '#580346', 
                                HSP= 'palegreen3' )) +
  theme(plot.title = element_text(hjust = 0.5))
  # geom_text_repel(aes(x = PC1,y=PC2,label = rownames(df_pcs)),size=4,adj = 1.2)
p1


features0 = as.data.frame(df_pca$rotation)
features0$features = row.names(features0)

features0$features = factor(features0$features,levels = row.names(features0))
rownames(features0)

features0$group <- c("climate","berry",rep("climate", 7), rep("soil", 22))

p2 <- ggplot(features0,aes(x=PC1,y=PC2,label=features,color=group)) + geom_point()+
  geom_segment(data = features0, aes(x = 0, xend = PC1, y=0, yend = PC2, color = features), 
               size = 1, arrow = arrow(length = unit(0.02, "npc")), alpha = 1)+
  scale_color_manual(values = c(climate="slateblue1", soil="gold4", berry='red1'))+
  xlab(percentage[1]) + ylab(percentage[2]) +
  geom_text_repel(size = 5) +
  geom_hline(yintercept = 0, linetype = "dotted",size =1.0) +
  geom_vline(xintercept = 0, linetype = "dotted",size =1.0) +
  theme_bw(base_size = 15) +theme(plot.title = element_text(hjust = 0.8))+
  theme(legend.position = "none")
p2



pdf("Fig.1.pdf",height =6,width = 6)
p1
p2
dev.off()

