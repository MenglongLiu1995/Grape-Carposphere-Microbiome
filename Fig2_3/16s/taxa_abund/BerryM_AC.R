wd="E:\\Data\\BerryMicromeGene2021\\WGCNA"
setwd(wd)
library(progress)
library(ggplot2)
library("rlang")
library("dplyr")
library(ggpmisc)
library(ggpubr)
library(xlsx)


# datExpr=read.csv(datExpr,"2021_NX_berry_Stage2_Transcripts.csv")
gene_lst=read.xlsx("gene_module_ann_v3.xlsx",sheetIndex=2,header = T)
gene_lst1=gene_lst[,1]
dim(datExpr)
rownames(datExpr)

don=read.csv("BerryTraitMicrobe.csv",header = T,row.names = 1) 

rownames(don)
colnames(don)
formula <- y ~ x
don$site <- factor(don$site)
don$variety=factor(don$variety)
##
# i="VIT_200s0199g00100"
lst <- list() ##用于批量画图的排版
n <- 0 ##  用于批量画图的排版
nr_ann
rownames(nr_ann)
colnames(nr_ann)
filter_gene <- c()
r <- c()
p_val <- c()
for (i in gene_lst1) {
n=n+1

data0=cbind(don[,1:15],datExpr[,i])
colnames(data0)[16]=i
data0$site <- factor(data0$site)
data0$variety=factor(data0$variety)
colnames(data0)
p <- ggplot(data0,aes_string('AC_final',i)) + 
  geom_point(size=4, aes(color=site,shape=variety))+
  stat_smooth(method = "lm") + 
  scale_shape_manual(values = c("CS"=17,"ME"=6,"CH"=8))+
  stat_cor(label.y = max(data0[,i])*0.95,size=5) + 
  stat_poly_eq(aes(label = ..eq.label..),size=5,
               formula = formula, parse = TRUE) +
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5,size=10),
        axis.text.y = element_text(size = 12),
        axis.title =  element_text(size = 12),
        legend.title = element_text(size=10),
        legend.text  = element_text(size=10))+
  theme(plot.title = element_text(hjust = 0.5,size=15))+ 
  labs(x="Total anthocyanin (mg/g)", y= paste0("Relative expression of ",i)) +
  # labs(x="Climatic distance (Euclidean)", y= "Community similarity (%)") +  # 
  ggtitle(nr_ann[i,9])+
   scale_color_manual(values=c (HD='cyan3', YSGD='dodgerblue4', CHFC='red4', SSB='sandybrown', XG='orangered2', HSP='palegreen3',XGS = '#fff143'))
  # scale_y_continuous(limits = c(0,6))
 print(p)
lst[[n]] <- p ### 把画的


filter_gene=append(filter_gene,i)
r <- append(r,ggplot_build(p)$data[[3]]$r)
p_val <- append(p_val, ggplot_build(p)$data[[3]]$p.value)
}
pdf("Module_BerryM_cor_2.pdf",width=20,height=15)
do.call("ggarrange", c(lst,ncol=3, nrow=3))
dev.off()

relation_lst=data.frame(filter_gene=filter_gene, r=r, p_val=p_val)
relation_lst1=nr_ann[relation_lst[,1],]
relation_lst=cbind(relation_lst,relation_lst1)
write.csv(relation_lst,"salmon_AC_gene.csv")



### 计算Burkholdria的种水平与AC的相关性
don=
metadat1=read.xlsx("sample_lst.xlsx",sheetIndex = 1,header = T)
rownames(metadat1)=metadat1[,1]
row_dat=read.csv('E:\\Data\\BerryMicromeGene2021\\allsampleOTU_16s.csv',header = T,check.names = F,row.names = 1)
OTU_annotation=read.csv("E:\\Data\\BerryMicromeGene2021\\OTU_annotation_16s.csv",header = T,check.names = F,row.names = 1)
BerryOTU0=row_dat[,rownames(metadat1)]
BerryOTU0=cbind(BerryOTU0,OTU_annotation)
berryOTU1=read.csv("E:\\Data\\BerryMicromeGene2021\\berry_OTU_16s_56sample.csv",header = T,row.names = 1,check.names = F)
colnames(BerryOTU0)

berryOTU2=as.data.frame(t(BerryOTU0))
berryOTU2=berryOTU2[1:54,]
berryOTU2=as.data.frame(apply(berryOTU2,2,as.numeric))
rownames(berryOTU2)=colnames(BerryOTU0)[1:54]

Burk_lst1=OTU_annotation[which(OTU_annotation$f=="f__Xanthomonadaceae"),]

formula <- y ~ x

lst <- list() ##用于批量画图的排版
n <- 0 ##  用于批量画图的排版
filter_gene <- c()
r <- c()
p_val <- c()
for (j in rownames(Burk_lst1)) {
  n=n+1
  # j="OTU33840"
  data0=cbind(don[,1:15],berryOTU2[,j])
  colnames(data0)[16]=j
  data0$site <- factor(data0$site)
  data0$variety=factor(data0$variety)
  colnames(data0)
  p <- ggplot(data0,aes_string('AC_final',j)) + 
    geom_point(size=4, aes(color=site,shape=variety))+
    stat_smooth(method = "lm") + 
    scale_shape_manual(values = c("CS"=17,"ME"=6,"CH"=8))+
    stat_cor(label.y = max(data0[,j])*0.95,size=5) + 
    stat_poly_eq(aes(label = ..eq.label..),size=5,
                 formula = formula, parse = TRUE) +
    theme_bw() + 
    theme(axis.text.x = element_text(hjust = 0.5,size=10),
          axis.text.y = element_text(size = 12),
          axis.title =  element_text(size = 12),
          legend.title = element_text(size=10),
          legend.text  = element_text(size=10))+
    theme(plot.title = element_text(hjust = 0.5,size=15))+ 
    labs(x="Total anthocyanin (mg/g)", y= paste0("Abundance of ",j)) +
    # labs(x="Climatic distance (Euclidean)", y= "Community similarity (%)") +  # 
    ggtitle(paste0(OTU_annotation[j,6]," ",OTU_annotation[j,7]))+
    scale_color_manual(values=c (HD='cyan3', YSGD='dodgerblue4', CHFC='red4', SSB='sandybrown', XG='orangered2', HSP='palegreen3',XGS = '#fff143'))
  # scale_y_continuous(limits = c(0,6))
  print(p)
  lst[[n]] <- p ### 把画的
  
  
  filter_gene=append(filter_gene,j)
  r <- append(r,ggplot_build(p)$data[[3]]$r)
  p_val <- append(p_val, ggplot_build(p)$data[[3]]$p.value)
}
pdf("Burk_s_cor_AC1.pdf",width=20,height=15)
do.call("ggarrange", c(lst,ncol=3, nrow=3))
dev.off()

relation_lst1=data.frame(filter_gene=filter_gene, r=r, p_val=p_val)
write.csv(relation_lst1,"Burk_s_cor_AC.csv")


########### 筛选OTU

Genus=read.csv("E:\\Data\\BerryMicromeGene2021\\microeco\\taxa_abund\\Genus_abund.csv", row.names = 1, header = T,check.names = F)
head(Genus)

tax_lst=read.csv("E:\\Data\\BerryMicromeGene2021\\microeco\\taxa_abund\\tax_table1.csv", row.names = 1, header = T,check.names = F)

gene_filter1 <- data.frame(Genus[which(apply(Genus, 1, function(x){mean(x)})
                                       >0.01),], check.names=F)

write.csv(gene_filter1,"E:\\Data\\BerryMicromeGene2021\\microeco\\taxa_abund\\GenusAbun_0.1.csv")

lst=c()
for (i in rownames(gene_filter1)) {
  
  extract_otu=rownames(tax_lst[which(tax_lst$Taxa_g==i),])
  lst=append(lst,extract_otu)
  
}

lst #### otu list with their Genus taxa >1%
otu=read.csv("E:\\Data\\BerryMicromeGene2021\\microeco\\basic_files\\feature_table.csv",
             row.names = 1,header = T,check.names = F )
otu_reabunde=otu/rowSums(otu)
otu_Genus0.01=otu_reabunde[,lst]
write.csv(otu_Genus0.01,'otu_Genus0.01.csv')


######## ### TOP10 OTU~ berry quality
library(progress)
library(ggplot2)
library("rlang")
library("dplyr")
library(ggpmisc)
library(ggpubr)
library(xlsx)

wd="E:\\Data\\BerryMicromeGene2021\\ggclusterNet\\"
setwd(wd)

otu_lst=read.csv("X_loading.csv",
                row.names = 1,header = T,check.names = F )
otu_tab=read.csv("otu_tab_filter1.csv",
                 row.names = 1,header = T,check.names = F )
otu_tax=read.csv("otu_tax_filter.csv",
                 row.names = 1,header = T,check.names = F )
berry_var=read.csv("berry_var.csv",
                 row.names = 1,header = T,check.names = F )
berry_var1=as.data.frame(t(berry_var))

otu=read.csv("E:\\Data\\BerryMicromeGene2021\\microeco\\basic_files\\feature_table.csv",
                 row.names = 1,header = T,check.names = F )
otu_reabunde=otu/rowSums(otu)

sample=read.csv("E:\\Data\\BerryMicromeGene2021\\microeco\\basic_files\\sample_table.csv",row.names = 1,header = T,check.names = F )

quality_lst=colnames(berry_var1)
otu_lst1=rownames(otu_lst)[1:10]

formula <- y ~ x

lst <- list() ##用于批量画图的排版
n <- 0 ##  用于批量画图的排版
filter_gene <- c()
r <- c()
p_val <- c()
for (i in quality_lst) {
  for (j in otu_lst1) {
    n=n+1
    # i="beta_Ionone"
    # j='b_Erysipelotrichaceae'
    data0=cbind(otu_reabunde[,j],berry_var1[,i])
    data0=as.data.frame(data0)
    colnames(data0)[1]=j
    colnames(data0)[2]=i
    data0[,j]=data0[,j]*100
    data0$site <- factor(sample$site)
    data0$cultivar=factor(sample$cultivar,levels = c("CS","ME","CH"))
    colnames(data0)
    p <- ggplot(data0,aes_string(j,i)) + 
      geom_point(size=4, aes(color=site,shape=cultivar))+
      stat_smooth(method = "lm") + 
      scale_shape_manual(values = c("CS"=16,"ME"=17,"CH"=15))+
      stat_cor(label.y = max(data0[,j])*0.95,size=5) + 
      stat_poly_eq(aes(label = ..eq.label..),size=5,
                   formula = formula, parse = TRUE) +
      theme_bw() + 
      theme(axis.text.x = element_text(hjust = 0.5,size=10),
            axis.text.y = element_text(size = 12),
            axis.title =  element_text(size = 12),
            legend.title = element_text(size=10),
            legend.text  = element_text(size=10))+
      theme(plot.title = element_text(hjust = 0.5,size=15))+ 
      labs(x=paste0("Relative Abundance of ",j," (%)"), y=i ) +
      # labs(x="Climatic distance (Euclidean)", y= "Community similarity (%)") +  # 
      ggtitle(j)+
      scale_color_manual(values=c (HD='cyan3', YSGD='dodgerblue4', CHFC='red4', SSB='sandybrown', XG='orangered2', HSP='palegreen3',XGS = '#fff143'))
    # scale_y_continuous(limits = c(0,6))
    print(p)
    
    lst[[n]] <- p ### 
    filter_gene=append(filter_gene,j)
    r <- append(r,ggplot_build(p)$data[[3]]$r)
    p_val <- append(p_val, ggplot_build(p)$data[[3]]$p.value)
  }
}
pdf("g_Q_cor_v2.pdf",width=20,height=15)
do.call("ggarrange", c(lst,ncol=3, nrow=3))
dev.off()

relation_lst1=data.frame(filter_gene=filter_gene, r=r, p_val=p_val)
write.csv(relation_lst1,"g_Q_cor_1.csv")



############### Genus0.1~berry quality
library(ggplot2)
gene_filter1=read.csv("E:\\Data\\BerryMicromeGene2021\\microeco\\taxa_abund\\GenusAbun_0.1.csv",row.names = 1,check.names = F,header = T)
rownames(gene_filter1)
g_lst = c("f_Moraxellaceae_g_Cavicella","f_Atopobiaceae_g_","f_Gemmatimonadaceae_g_Gemmatimonas",
        "f_Sphingomonadaceae_g_Sphingomonas","c_Subgroup_6_g_","f_Muribaculaceae_g_","f_Lachnospiraceae_g_")
gene_filter2=as.data.frame(t(gene_filter1))
gene_filter2 = gene_filter2[,g_lst]


berry_var=read.csv("E:\\Data\\BerryMicromeGene2021\\ggclusterNet\\berry_var.csv",
                   row.names = 1,header = T,check.names = F )
colnames(berry_var)
aroma_lst = c("cis_3_Hexenal","trans_2_Hexenal","trans_2_Pentenal")
berry_var1 = berry_var[,aroma_lst]

sample=read.csv("E:\\Data\\BerryMicromeGene2021\\microeco\\basic_files\\sample_table.csv",row.names = 1,header = T,check.names = F )

formula <- y ~ x

lst <- list() ##用于批量画图的排版
n <- 0 ##  用于批量画图的排版
filter_gene <- c()
r <- c()
p_val <- c()
for (i in colnames(berry_var1)) {
  for (j in colnames(gene_filter2)) {
    n=n+1
    # i="beta_Ionone"
    # j='b_Erysipelotrichaceae'
    data0=cbind(gene_filter2[,j],berry_var[,i])
    data0=as.data.frame(data0)
    colnames(data0)[1]=j
    colnames(data0)[2]=i
    data0[,j]=data0[,j]*100
    data0$site <- factor(sample$site)
    data0$cultivar=factor(sample$cultivar,levels = c("CS","ME","CH"))
    colnames(data0)
    p <- ggplot(data0,aes_string(j,i)) + 
      geom_point(size=4, aes(color=site,shape=cultivar))+
      stat_smooth(method = "lm") + 
      scale_shape_manual(values = c("CS"=16,"ME"=17,"CH"=15))+
      stat_cor(label.y = max(data0[,j])*0.95,size=5) + 
      stat_poly_eq(aes(label = ..eq.label..),size=5,
                   formula = formula, parse = TRUE) +
      theme_bw() + 
      theme(axis.text.x = element_text(hjust = 0.5,size=10),
            axis.text.y = element_text(size = 12),
            axis.title =  element_text(size = 12),
            legend.title = element_text(size=10),
            legend.text  = element_text(size=10))+
      theme(plot.title = element_text(hjust = 0.5,size=15))+ 
      labs(x=paste0("Relative Abundance of ",j," (%)"), y=i ) +
      # labs(x="Climatic distance (Euclidean)", y= "Community similarity (%)") +  # 
      ggtitle(j)+
      scale_color_manual(values=c (HD='cyan3', YSGD='dodgerblue4', CHFC='red4', SSB='sandybrown', XG='orangered2', HSP='palegreen3',XGS = '#fff143'))
    # scale_y_continuous(limits = c(0,6))
    print(p)
    
    lst[[n]] <- p ### 
    filter_gene=append(filter_gene,j)
    r <- append(r,ggplot_build(p)$data[[3]]$r)
    p_val <- append(p_val, ggplot_build(p)$data[[3]]$p.value)
  }
}
pdf("genus_aldehydes.pdf",width=20,height=15)
do.call("ggarrange", c(lst,ncol=3, nrow=3))
dev.off()

relation_lst1=data.frame(filter_gene=filter_gene, r=r, p_val=p_val)
write.csv(relation_lst1,"g_Q_cor_1.csv")


str(otu_tab)
sample=read.csv("E:\\Data\\BerryMicromeGene2021\\microeco\\basic_files\\sample_table.csv",
             row.names = 1,header = T,check.names = F )

tax=read.csv("E:\\Data\\BerryMicromeGene2021\\microeco\\basic_files\\tax_table.csv",
             row.names = 1,header = T,check.names = F )
g_lst1=unique(tax[which(tax$Family=="f__Sphingomonadaceae"),"Genus"])

Family=read.csv("E:\\Data\\BerryMicromeGene2021\\microeco\\taxa_abund\\Family_abund1.csv",
               row.names = 1,header = T,check.names = F )
Family=as.data.frame(t(Family))
colnames(Family)

Genus=read.csv("E:\\Data\\BerryMicromeGene2021\\microeco\\taxa_abund\\Genus_abund1.csv",
             row.names = 1,header = T,check.names = F )
Genus=as.data.frame(t(Genus))

Berry=read.csv("E:\\Data\\BerryMicromeGene2021\\microeco\\taxa_abund\\BerryTrait_comp_short.csv",
             row.names = 1,header = T,check.names = F )
Berry=as.data.frame(t(Berry))



quality_lst=colnames(Berry)
f_lst1=colnames(Family)

formula <- y ~ x

lst <- list() ##用于批量画图的排版
n <- 0 ##  用于批量画图的排版
filter_gene <- c()
r <- c()
p_val <- c()
for (i in quality_lst) {
for (j in g_lst1) {
  n=n+1
  i="beta_Ionone"
     # j='b_Erysipelotrichaceae'
  data0=cbind(Genus[,j],Berry[,i])
  data0=as.data.frame(data0)
   colnames(data0)[1]=j
   colnames(data0)[2]=i
  data0[,j]=data0[,j]*100
  data0$site <- factor(sample$site)
  data0$cultivar=factor(sample$cultivar,levels = c("CS","ME","CH"))
  colnames(data0)
  p <- ggplot(data0,aes_string(j,i)) + 
    geom_point(size=4, aes(color=site,shape=cultivar))+
    stat_smooth(method = "lm") + 
    scale_shape_manual(values = c("CS"=16,"ME"=17,"CH"=15))+
    stat_cor(label.y = max(data0[,j])*0.95,size=5) + 
    stat_poly_eq(aes(label = ..eq.label..),size=5,
                 formula = formula, parse = TRUE) +
    theme_bw() + 
    theme(axis.text.x = element_text(hjust = 0.5,size=10),
          axis.text.y = element_text(size = 12),
          axis.title =  element_text(size = 12),
          legend.title = element_text(size=10),
          legend.text  = element_text(size=10))+
    theme(plot.title = element_text(hjust = 0.5,size=15))+ 
    labs(x=paste0("Relative Abundance of ",j," (%)"), y=i ) +
    # labs(x="Climatic distance (Euclidean)", y= "Community similarity (%)") +  # 
    ggtitle(j)+
    scale_color_manual(values=c (HD='cyan3', YSGD='dodgerblue4', CHFC='red4', SSB='sandybrown', XG='orangered2', HSP='palegreen3',XGS = '#fff143'))
  # scale_y_continuous(limits = c(0,6))
  print(p)
  
  lst[[n]] <- p ### 
  filter_gene=append(filter_gene,j)
  r <- append(r,ggplot_build(p)$data[[3]]$r)
  p_val <- append(p_val, ggplot_build(p)$data[[3]]$p.value)
}
}
pdf("g_Q_cor_v1.pdf",width=20,height=15)
do.call("ggarrange", c(lst,ncol=3, nrow=3))
dev.off()

relation_lst1=data.frame(filter_gene=filter_gene, r=r, p_val=p_val)
write.csv(relation_lst1,"g_Q_cor_1.csv")


#################### 
install.packages("remotes")
remotes::install_github("chentianlu/gramm4R")

library(devtools)   
install_github("PhDMeiwp/basicTrendline@master", force = TRUE)
library(basicTrendline)

library("gramm4R")

gene_filter1=read.csv("E:\\Data\\BerryMicromeGene2021\\microeco\\taxa_abund\\GenusAbun_0.1.csv",row.names = 1,check.names = F,header = T)
gene_filter2=as.data.frame(t(gene_filter1))

gene_filter1=read.csv("E:\\Data\\BerryMicromeGene2021\\microeco\\basic_files\\feature_table.csv",row.names = 1,check.names = F,header = T)
# 假设df是你的数据框  
df <- data.frame(A = c(0, 1, 0, 2), B = c(0, 0, 3, 0), C = c(4, 0, 0, 5))  

# 使用lapply替换所有0为0.01，然后直接转换为数据框  
df_new <- data.frame(lapply(gene_filter1, function(x) replace(x, x == 0, 0.0000000001)))  

# 查看结果  
print(df_new)

write.csv(df_new, 'feature_clr.csv')



gene_filter2=as.data.frame(t(gene_filter1))

berry_var=read.csv("E:\\Data\\BerryMicromeGene2021\\ggclusterNet\\berry_var.csv",
                   row.names = 1,header = T,check.names = F )
berry_var1=as.data.frame(t(berry_var))

micro_pro=microbesPrepro()




preGramm(gene_filter2,berry_var1)

Gramm(gene_filter2,berry_var1,metaNor=TRUE,rarefaction = FALSE,r = 0.5,alpha = 0.05)

load("E:\\文献\\AAAAA微生物\\microbes.rda")
load("E:\\文献\\AAAAA微生物\\covariates.rda")
load("E:\\文献\\AAAAA微生物\\metabolites.rda")
Gramm(microbes,metabolites,covariates,metaNor=TRUE,rarefaction = FALSE,r = 0.5,alpha = 0.05)

preGramm(metabolites,microbes)
naiveGramm(metabolites,microbes,covariates)
data("metabolites")
data("microbes")
nlfitGramm(metabolites,microbes)
Gramm(metabolites,microbes,covariates)

metabolites[[1]]

setClass("Mymetabolites",  
         slots = list(  
           data = "data.frame"  # 假设我们想要存储一个数据框  
         )  
)
Mymetabolites <- new("Mymetabolites", data = berry_var1[1:54,1:10])  

setClass("Mymicrobe",  
         slots = list(  
           data = "data.frame"  # 假设我们想要存储一个数据框  
         )  
)
Mymicrobe <- new("Mymicrobe", data = gene_filter2[1:54,1:17])  

preGramm(Mymetabolites,Mymicrobe)

#####
install.packages("imputeLCMD")
install.packages("pcaMethods")
BiocManager::install("pcaMethods")
library("pcaMethods")
library("imputeLCMD")
library(BiocManager)
library("phyloseq")


soil_data=read.csv("E:\\Data\\BerryMicromeGene2021\\ggclusterNet\\soil_data_20_depth.csv",row.names = 1,header = T)
mebolite_Prepro = metabolitesPrepro(berry_var1, missPro = TRUE, missMethod = "QRILC",
                                    scalingPro = TRUE, transPro = TRUE, kValue = 3,
                                    phenoData = NA)

mebolite_Prepro 

gene_filter2[] <- lapply(gene_filter2, function(x) ifelse(x == 0, 1e-10, x))
microbe_Prepro = microbesPrepro(gene_filter2,missPro = TRUE,missMethod = "mean",
                                rarePro = F, scalingPro = F, transPro = T, 
                                kValue = 3)


microbe_Prepro
res = naivegramm(mebolite_Prepro,microbe_Prepro,covdata = NA, r = 0.3,alpha = 0.99)

res = as.data.frame(res)



write.csv(res,"E:\\Data\\BerryMicromeGene2021\\ggclusterNet\\res.csv")

res1 = naivegramm(berry_var1,microbe_Prepro,covdata = NA, r = 0.3,alpha = 0.99)
res3 = as.data.frame(res1)
write.csv(res3,"E:\\Data\\BerryMicromeGene2021\\ggclusterNet\\res3.csv")

#######
install.packages("chemometrics")
install.packages("compositions")
library("compositions")
library(psych)
library(dplyr)
genus_abun = read.csv("E:\\Data\\BerryMicromeGene2021\\ggclusterNet\\GenusAbun_0.1_CTR.csv",
                      row.names = 1,header = T,check.names = F)
genus_abun = clr(df_new)
head(genus_abun)

soil_data=read.csv("E:\\Data\\BerryMicromeGene2021\\ggclusterNet\\soil_data_20_depth.csv",row.names = 1,header = T)
mebolite_Prepro = metabolitesPrepro(berry_var, missPro = TRUE, missMethod = "QRILC",
                                    scalingPro = TRUE, transPro = TRUE, kValue = 3,
                                    phenoData = NA)

mebolite_Prepro 

gene_filter2[] <- lapply(gene_filter2, function(x) ifelse(x == 0, 1e-10, x))
microbe_Prepro = microbesPrepro(df_new,missPro = TRUE,missMethod = "mean",
                                rarePro = F, scalingPro = F, transPro = T, 
                                kValue = 3)


res=naivegramm(mebolite_Prepro,microbe_Prepro, covdata  = NA)

cor.result <- pcor.test(genus_abun,berry_var1,soil_data)
r <- cor.result$r     #相关性
r

p <- cor.result$p     #相关性的pvalue值
p
write.csv(p,"p.csv")



#############
# function - naivegramm #
########################################################################
# File: naivegramm.R
# Aim : Generalized Correlation Analysis for Metabolome and Microbiome (GRaMM), for inter-correlation pairs discovery 
#       among metabolome and microbiome.
#---------------------------------------------------------------------------------------------------------------------
# Author : Tianlu Chen, Tao Sun, Dandan Liang, Mengci Li
# Email  : chentianlu@sjtu.edu.cn
# Date   : 2020-08
# Version: 1.0
#---------------------------------------------------------------------------------------------------------------------
#
#
######################################################################################
## Input:                                                                           ##
##    metaData --- A dataframe of metabolites                                       ##
##                       (rows: samples,columns: metabolites)                       ##                    
##    micData --- A dataframe of microbes                                           ##
##                       (rows: samples,columns: microbes)                          ##  
##    covdata --- A dataframe of confounder                                         ##
##                       (rows:samples,columns:confounder)                          ##
##    r --- Correlation coefficient threshold                                       ##
##                Default: TRUE                                                     ##
##    alpha --- FDR threshold                                                       ##
##                Default: 0.99                                                     ##
##    pheno --- Phenotype name                                                      ##
##                                                                                  ##
######################################################################################
## Output:                                                                          ## 
##    res --- A list included coefficient, p value and type of correlation          ##
##                                                                                  ##
######################################################################################
#---------------------------------------------------------------------------------------------------------------------

#####
# function - microbesPrepro #
########################################################################
# File: microbesPrepro.R
# Aim : Preprocessing for microbiome data
#---------------------------------------------------------------------------------------------------------------------
# Author : Tianlu Chen, Tao Sun, Dandan Liang, Mengci Li
# Email  : chentianlu@sjtu.edu.cn
# Date   : 2020-08
# Version: 1.0
#---------------------------------------------------------------------------------------------------------------------
#
#
######################################################################################
## Input:                                                                           ##
##    micData --- A dataframe of microbes                                           ##
##                       (rows: samples,columns: microbes)                          ##                    
##    missPro --- Flag of whether to complete the missing value                     ##
##                Default: TRUE                                                     ##
##    missMethod --- Method for completing missing values                           ##
##                Default: mean                                                     ##
##    rarePro --- Flag of whether to rarefy                                         ##
##                Default: TRUE                                                     ##
##    scalingPro --- Flag of whether to scale                                       ##
##                Default: TRUE                                                     ##
##    transPro --- Flag of whether to transform                                     ##
##                Default: TRUE                                                     ##
##    kValue --- The number of nearest neighbours to use in knn imputation          ##
##                Default: 3                                                        ##
##    phenoData --- A dataframe of phenotype data                                   ##
##                Default: NA  (rows: samples,columns: phenotype)                   ##
##    phenoDataType --- Phenotype data type                                         ##
##                Default: continuous                                               ##
######################################################################################
## Output:                                                                          ## 
##     microbes_preproRes --- The dataframe of microbes after preprocessing         ##
##                       (rows: samples,columns: microbes)                          ##
##                                                                                  ##
######################################################################################
#---------------------------------------------------------------------------------------------------------------------

microbesPrepro <- function(micData, missPro = TRUE,missMethod = "mean",rarePro = T, scalingPro = T, transPro = T, kValue = 3,phenoData = NA,phenoDataType = "continuous")
{
  if(!file.exists("./results/prepro/")){
    dir.create("./results/prepro/",recursive = T)
  }
  
  microbes_preproRes <- micData
  
  microbes <- micData
  
  if(missPro){
    ######################################
    ############ Step 1 - knn ############
    ######################################
    if(missMethod == "knn"){
      
      microbes_preproRes <- DMwR::knnImputation(microbes_preproRes,k = kValue,scale = TRUE)
    }else if(missMethod == "QRILC"){
      obj.QRILC = imputeLCMD::impute.QRILC(microbes_preproRes)
      microbes_preproRes = obj.QRILC[[1]]
    }else if(missMethod == "mean"){
      for (i in colnames(microbes_preproRes)) {
        microbes_preproRes[,i] <- Hmisc::impute(microbes_preproRes[,i], fun = mean)
      }
    }else if(missMethod == "min"){
      for (i in colnames(microbes_preproRes)) {
        microbes_preproRes[,i] <- Hmisc::impute(microbes_preproRes[,i], fun = min)
      }
    }else if(missMethod == "median"){
      for (i in colnames(microbes_preproRes)) {
        microbes_preproRes[,i] <- Hmisc::impute(microbes_preproRes[,i], fun = median)
      }
    } 
  }
  
  if(phenoDataType != "categorical"){
    # If the value of 0 is more than 70%, the microbe is deleted 
    index <- c()
    for (i in 1:ncol(microbes_preproRes)) {
      sum_count <- sum(microbes_preproRes[,i] == 0)
      if(sum_count > (nrow(microbes_preproRes)*0.7) ){
        index <- append(index,i)
      }
    }
    if(!is.null(index)){
      microbes_preproRes <- microbes_preproRes[,-index]
    }
  }else if(phenoDataType == "categorical"){
    if(is.na(phenoData)){
      stop("Please input your phenoData!")
    }else{
      index <- c()
      for (i in unique(phenoData[,1])) {
        tempmicData <- micData[rownames(phenoData)[which(phenoData[,1] == i)],]
        # If the value of 0 is more than 70% in each group, the microbe is deleted
        for (j in 1:ncol(tempmicData)) {
          sum_count <- sum(tempmicData[,j] == 0)
          if(sum_count > (nrow(tempmicData)*0.7) ){
            index <- append(index,j)
          }
        }
      }
      if(!is.null(index)){
        index <- unique(index)
        microbes_preproRes <- microbes_preproRes[,-index]
      }
    }
  }
  
  if(rarePro){
    
    ######################################
    ######## Step 1 - rarefaction ########
    ######################################
    
    microbes_preproRes <- apply(microbes,2,as.integer)
    microbes_preproRes <- t(microbes_preproRes)
    colnames(microbes_preproRes) <- rownames(microbes)
    microbes_number <- nrow(microbes_preproRes) * ncol(microbes_preproRes)
    
    # Create a pretend taxonomy table
    taxmat <- matrix(sample(letters, microbes_number, replace = TRUE), nrow = nrow(microbes_preproRes),
                     ncol = 7)
    
    rownames(taxmat) <- rownames(microbes_preproRes)
    colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family",
                          "Genus", "Species")
    OTU = phyloseq::otu_table(microbes_preproRes, taxa_are_rows = TRUE)
    TAX = phyloseq::tax_table(taxmat)
    physeq = phyloseq(OTU, TAX)
    ## set.seed(100)
    raretest <- phyloseq::rarefy_even_depth(physeq,rngseed=100)
    microbes_preproRes <- raretest@otu_table@.Data
    microbes_preproRes <- t(microbes_preproRes)
    colnames(microbes_preproRes) <- colnames(microbes)
    rownames(microbes_preproRes) <- rownames(microbes)
    
    ##  deal with 0 value
    microbes_preproRes[microbes_preproRes == 0] <- 1
  }
  
  if(scalingPro){
    ######################################
    ####### Step 2 - normalization #######
    ######################################
    
    microbes_preproRes<- apply(microbes_preproRes,2,function(each_col){
      col_sum  <-  sum(each_col,na.rm = TRUE)
      each_col <-  each_col / col_sum * 30000
      return(each_col)
    })
  }
  
  if(transPro){
    ######################################
    ############ Step 3 - clr ############
    ######################################
    
    microbes_preproRes <- apply(microbes_preproRes,2,function(each_col){
      log_col <- psych::geometric.mean(each_col)
      each_col <- abs(log(each_col) / log(log_col))
      return(each_col)
    })
  } 
  
  microbes_preproRes <- as.data.frame(microbes_preproRes)
  write.csv(microbes_preproRes,"./results/prepro/microbes_prepro_res.csv")
  return(microbes_preproRes)
  
}
############
# function - metabolitesPrepro #
########################################################################
# File: metabolitesPrepro.R
# Aim : Preprocessing for metabolome data
#---------------------------------------------------------------------------------------------------------------------
# Author : Tianlu Chen, Tao Sun, Dandan Liang, Mengci Li
# Email  : chentianlu@sjtu.edu.cn
# Date   : 2020-08
# Version: 1.0
#---------------------------------------------------------------------------------------------------------------------
#
#
######################################################################################
## Input:                                                                           ##
##    metaData --- A dataframe of metabolites                                       ##
##                       (rows: samples,columns: metabolites)                       ##
##    missPro --- Flag of whether to complete the missing value                     ##
##                Default: TRUE                                                     ##
##    missMethod --- Method for completing missing values                           ##
##                Default: QRILC                                                    ##
##    scalingPro --- Flag of whether to scale                                       ##
##                Default: TRUE                                                     ##
##    transPro --- Flag of whether to transform                                     ##
##                Default: TRUE                                                     ##
##    kValue --- The number of nearest neighbours to use in knn imputation          ##
##                Default: 3                                                        ##
##    phenoData --- A dataframe of phenotype data                                   ##
##                Default: NA  (rows: samples,columns: phenotype)                   ##
##    phenoDataType --- Phenotype data type                                         ##
##                Default: continuous                                               ##
######################################################################################
## Output:                                                                          ## 
##    metabolites_preproRes --- The dataframe of metabolites after preprocessing    ##
##                       (rows: samples,columns: metabolites)                       ##
##                                                                                  ##
######################################################################################
#---------------------------------------------------------------------------------------------------------------------
metabolitesPrepro <- function(metaData, missPro = TRUE, missMethod = "QRILC", scalingPro = TRUE, transPro = TRUE, kValue = 3,phenoData = NA,phenoDataType = "continuous")  
{
  if(!file.exists("./results/prepro")){
    dir.create("./results/prepro",recursive = T)
  }
  
  metabolites_preproRes <- metaData
  
  ######################################
  ############ missing value ###########
  ######################################
  if(missPro){
    if(missMethod == "knn"){
      # metabolites_preproRes[metabolites_preproRes == 0] <- NA
      metabolites_preproRes <- DMwR::knnImputation(metabolites_preproRes,k = kValue,scale = F)
    }else if(missMethod == "QRILC"){
      obj.QRILC = imputeLCMD::impute.QRILC(metabolites_preproRes)
      metabolites_preproRes = obj.QRILC[[1]]
    }else if(missMethod == "mean"){
      for (i in colnames(metabolites_preproRes)) {
        metabolites_preproRes[,i] <- Hmisc::impute(metabolites_preproRes[,i], fun = mean)
      }
    }else if(missMethod == "min"){
      for (i in colnames(metabolites_preproRes)) {
        metabolites_preproRes[,i] <- Hmisc::impute(metabolites_preproRes[,i], fun = min)
      }
    }else if(missMethod == "median"){
      for (i in colnames(metabolites_preproRes)) {
        metabolites_preproRes[,i] <- Hmisc::impute(metabolites_preproRes[,i], fun = median)
      }
    } 
  }
  
  if(phenoDataType != "categorical"){
    # If the value of 0 is more than 70%, the metabolite is deleted
    index <- c()
    for (i in 1:ncol(metabolites_preproRes)) {
      sum_count <- sum(metabolites_preproRes[,i] == 0)
      if(sum_count > (nrow(metabolites_preproRes)*0.7) ){
        index <- append(index,i)
      }
    }
    if(!is.null(index)){
      metabolites_preproRes <- metabolites_preproRes[,-index]
    }
  }else if(phenoDataType == "categorical"){
    if(is.na(phenoData)){
      stop("Please input your phenoData!")
    }else{
      index <- c()
      for (i in unique(phenoData[,1])) {
        tempmetaData <- metaData[rownames(phenoData)[which(phenoData[,1] == i)],]
        # If the value of 0 is more than 70% in each group, the metabolite is deleted
        for (j in 1:ncol(tempmetaData)) {
          sum_count <- sum(tempmetaData[,j] == 0)
          if(sum_count > (nrow(tempmetaData)*0.7) ){
            index <- append(index,j)
          }
        }
      }
      if(!is.null(index)){
        index <- unique(index)
        metabolites_preproRes <- metabolites_preproRes[,-index]
      }
    }
  }
  
  ######################################
  ############ normalization ###########
  ######################################
  if(scalingPro){
    metabolites_preproRes<- apply(metabolites_preproRes,2,function(each_col){
      col_sum  <-  sum(each_col,na.rm = TRUE)
      each_col <-  each_col / col_sum * 30000
      return(each_col)
    })
  }
  
  ######################################
  ########## transformation ############
  ######################################
  if(transPro){
    metabolites_preproRes <- log(metabolites_preproRes + 1)
  }
  
  metabolites_preproRes <- as.data.frame(metabolites_preproRes) 
  write.csv(metabolites_preproRes,"./results/prepro/metabolites_prepro_res.csv")
  return(metabolites_preproRes)
}

######
naivegramm<-function(metaData, micData, covdata, r = 0.3,alpha = 0.99,pheno="All")
{
  metaData <- as.data.frame(t(metaData))
  micData <- as.data.frame(t(micData))
  covdata <- as.data.frame(t(covdata))
  ## Exist confounders
  if(is.data.frame(covdata) && is.na(covdata[1,1]) != T){
    options(warn = -1)
    
    ## Generate result table
    r_result <- p_result <- r2_result <- type_result <- matrix(0, nrow = nrow(metaData),ncol = nrow(micData))
    rownames(p_result) <- rownames(r_result) <- rownames(type_result) <- rownames(metaData)
    colnames(p_result) <- colnames(r_result) <- colnames(type_result) <- rownames(micData)
    
    ## lm function was used to obtain the adjusted data and fit the linear model include linear regression and multiple linear regression
    for(i in seq_len(nrow(metaData))){
      for(j in seq_len(nrow(micData))){
        options(warn = -1)
        x1 <- t(metaData[i,])
        y1 <- t(rbind(micData[j,], covdata))
        lmx <- lm(x1 ~ y1)
        x11 <- y1[,1] - lmx$coefficients[3] * y1[,2]
        x2 <- t(rbind(metaData[i,], covdata))
        y2 <- t(micData[j,])
        lmy <- lm(y2 ~ x2)
        y22 <- x2[,1] - lmy$coefficients[3] * x2[,2]
        summ <- summary(lmx)
        pvalue <- summ$coefficients[2,4]
        rvalue <- lmx$coefficients[2] * sd(y1[,1]) / sd(x1)
        r2value <- rvalue^2
        
        ## linear: multiple linear regression
        if(pvalue < alpha & rvalue > r & !is.na(pvalue) & !is.na(rvalue)){
          p_result[i,j] <- pvalue
          r_result[i,j] <- rvalue
          r2_result[i,j] <- r2value
          type_result[i,j] <- "linear"
        }
        
        ## nonlinear: MIC
        else{
          micxy <- minerva::mine(y22,x11)
          micr <- micxy$MIC
          micp <- 0
          for(k in seq_len(101))
          {
            #bootstrap
            options(warn = -1)
            bootx <- matrix(sample(x11, replace = TRUE))
            booty <- matrix(sample(y22, replace = TRUE))
            if(sd(bootx) == 0){
              bootx <- matrix(sample(x11, replace = FALSE))
            }
            if(sd(booty) == 0){
              booty <- matrix(sample(y22, replace = FALSE))
            }
            tmp <- minerva::mine(booty, bootx)
            MICtp <- tmp$MIC
            tempM <- ifelse(micr <= MICtp,1,0)
            micp <- tempM + micp
          }
          micp <- micp / 101
          prsmic <- cor.test(y22, x11, method="pearson",use = "pairwise.complete.obs")
          prsr <- prsmic$estimate
          ltmic <- 1-micr + prsr^2
          p_result[i,j] <- micp
          r_result[i,j] <- micr
          #r2_result[i,j]<-ltmic
          type_result[i,j] <- "nonlinear"
        }
        
      }
      
    }
  }
  
  ## No confounders
  if(is.na(covdata[1,1]) == T){
    options(warn = -1)
    r_result <- p_result <- matrix(0, nrow = nrow(metaData),ncol = nrow(micData))
    r2_result <- type_result <- matrix(0, nrow = nrow(metaData),ncol = nrow(micData))
    rownames(p_result) <- rownames(r_result) <- rownames(type_result) <- rownames(metaData)
    colnames(p_result) <- colnames(r_result) <- colnames(type_result) <- rownames(micData)
    
    ## lm function was used to obtain the adjusted data and fit the linear model include linear regression and multiple linear regression
    for(i in seq_len(nrow(metaData))){
      for(j in seq_len(nrow(micData))){
        x1 <- t(metaData[i,])
        y1 <- t(micData[j,])
        lmx <- lm(x1 ~ y1)
        summ <- summary(lmx)
        pvalue <- summ$coefficients[2,4]
        rvalue <- lmx$coefficients[2] * sd(y1) / sd(x1)
        
        ## linear: linear regression
        if(pvalue < alpha & rvalue > r){
          p_result[i,j] <- pvalue
          r_result[i,j] <- rvalue
          r2value <- rvalue^2
          r2_result[i,j] <- r2value
          type_result[i,j] <- "linear"
        }
        
        ## nonlinear: MIC
        else{
          warnings('off')
          micxy <- minerva::mine(y1,x1)
          micr <- micxy$MIC
          micp <- 0
          for(k in seq_len(101))
          {
            bootx <- matrix(sample(x1, replace = FALSE))
            booty <- matrix(sample(y1, replace = FALSE))
            tmp <- minerva::mine(bootx, booty)
            MICtp <- tmp$MIC
            tempM <- ifelse(micr <= MICtp,1,0)
            micp <- tempM + micp
          }
          micp <- micp/101
          prsmic <- cor.test(y1,x1,method="pearson",use = "pairwise.complete.obs")
          prsr <- prsmic$estimate
          ltmic <- 1 - micr + prsr^2
          p_result[i,j] <- micp
          r_result[i,j] <- micr
          type_result[i,j] <- "nonlinear"
        }
      }
    }
  }
  res <- list()
  res[[paste0("r_",pheno)]] <- as.data.frame(t(r_result))
  res[[paste0("p_",pheno)]] <- as.data.frame(t(p_result))
  res[[paste0("type_",pheno)]] <- as.data.frame(t(type_result))
  return(res)
}
