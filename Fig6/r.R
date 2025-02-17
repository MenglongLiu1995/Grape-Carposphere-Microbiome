######### Fig 6b
library(microeco) # Microbial Community Ecology Data Analysis
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics
library(dplyr)
library(magrittr)

otu <- read.csv("featureTable.csv",header=T,check.names=FALSE ,row.names=1)

sample <- read.csv("group1.csv",header=T,check.names=FALSE ,row.names=1)

sample_cs <- sample[which(sample$cultivar=="CS"),]
sample_cH <- sample[which(sample$cultivar=="CH"),]

otu_CS = otu[,rownames(sample_cs)]
otu_CH = otu[,rownames(sample_cH)]

tax <- read.csv("taxa.csv", header=T,check.names=FALSE ,row.names=1)

tax %<>% tidy_taxonomy

df <- microtable$new(sample_table = sample,
                     otu_table = otu,
                     tax_table = tax,
                     auto_tidy = F)
df

dfCS <- microtable$new(sample_table = sample_cs,
                       otu_table = otu_CS,
                       tax_table = tax,
                       auto_tidy = F)


dfCH <- microtable$new(sample_table = sample_cH,
                       otu_table = otu_CH,
                       tax_table = tax,
                       auto_tidy = F)

###

dfCS$filter_pollution(taxa = c("mitochondria", "chloroplast"))
dfCS

dfCS$tidy_dataset()
dfCS
# 

dfCS$sample_sums() %>% range
# 

dfCS$rarefy_samples(sample.size = 21651)
dfCS$sample_sums() %>% range


dfCS$cal_abund()

dfCS$cal_alphadiv(PD = FALSE)

dfCS$cal_betadiv(unifrac = FALSE)

dfCH$filter_pollution(taxa = c("mitochondria", "chloroplast"))
dfCH

dfCH$tidy_dataset()
dfCH

dfCH$sample_sums() %>% range
# 

dfCH$rarefy_samples(sample.size = 22259)
dfCH$sample_sums() %>% range

dfCH$cal_abund()

dfCH$cal_alphadiv(PD = FALSE)


dfCH$cal_betadiv(unifrac = FALSE)

df$sample_table$group %<>% factor(., levels = c("Control", "inoculation"))
df$sample_table$cultivar %<>% factor(., levels = c("CS", "CH"))


t1 <- trans_abund$new(dataset = df, taxrank = "Phylum", ntaxa = 10)
t1 <- trans_abund$new(dataset = df, taxrank = "Family", ntaxa = 10)
t1 <- trans_abund$new(dataset = df, taxrank = "Genus", ntaxa = 10)
t1$plot_bar(others_color = "grey70", facet = c("cultivar","group"), xtext_keep = FALSE, legend_text_italic = FALSE)


### Fig 6c
t1 <- trans_diff$new(dataset = dfCS, method = "anova", group = "group", taxa_level = "Genus", filter_thres = 0.001)
P1 = t1$plot_diff_abund(use_number = 1:5, add_sig = T, coord_flip = F,color_values = c(Control="#7F7F7F", inoculation="chartreuse3"))+
  theme(text = element_text(size = 15))+
  ggtitle("CS")+
  lims(y = c(0,1))+
  theme_bw() + theme(panel.grid=element_blank(),axis.text.x = element_text(angle = 30,vjust = 0.5,hjust = 0.5))
#####
t1 <- trans_diff$new(dataset = dfCH, method = "anova", group = "group", taxa_level = "Genus", filter_thres = 0.001)
P2 = t1$plot_diff_abund(use_number = 1:5, add_sig = T, coord_flip = F,color_values = c(Control="#7F7F7F", inoculation="chartreuse3"))+
  theme(text = element_text(size = 15))+
  ggtitle("CH")+
  lims(y = c(0,1))+
  theme_bw() + theme(panel.grid=element_blank(),axis.text.x = element_text(angle = 30,vjust = 0.5,hjust = 0.5))

pdf("Genus_abund.pdf", width = 4, height = 4)
P1
P2
dev.off()

### Fig 6d

#####################################################################################
library(tidyr)
library(Rmisc)
library(ggplot2)
library(ggpubr)
library(multcompView)
library(tidyverse)
library(xlsx)

raw_fpkm <- read.csv('aroma_res.csv',check.names = F,header = T,row.names = 1)
raw_fpkm=as.data.frame(t(raw_fpkm))
rownames(raw_fpkm)

metadata <- read.csv('aroma_res_metadata.csv',check.names = F,header = T,row.names = 1)

rownames(metadata)

CH_lst = metadata[which(metadata$Cultivar=="CH"),]
extract_fpkm <- raw_fpkm[,rownames(CH_lst)]

CS_lst = metadata[which(metadata$Cultivar=="CS"),]
extract_fpkm <- raw_fpkm[,rownames(CS_lst)]
######################################################################################
keyword <- 'CH_res1'
oppath <- paste0(keyword,'.txt')
pcsvpath <- paste0(keyword,'.pdf')
annsvpath <- paste0(keyword,'.csv')
######################################################################################
gene_lstdata <- read.table('aroma_lst.txt',header = F,sep='\t')

#####################################################################################
p_lst <- list()

gene_lst <- gene_lstdata$V1 
for (gene_id in gene_lst) {
  extract_data <- extract_fpkm[gene_id,]
  new_data <- gather(extract_data,id,fpkm)
  new_data$DAT=factor(CH_lst$DAT, levels = c('DAT7','DAT14','DAT28'))
  new_data$Group=factor(CH_lst$Group)
  new_data$sample=factor(CH_lst$sample)
  new_data$fpkm=as.numeric(new_data$fpkm)

  str(new_data$fpkm)
  new_data0 <- summarySE(new_data,measurevar = 'fpkm',
                         groupvars=c('DAT','Group',"sample"))
 
  new_data0$DAT=factor(new_data0$DAT, levels = c('DAT7','DAT14','DAT28'))
  new_data0$Group=factor(new_data0$Group, levels = c('Control','Treat'))

  p <- ggplot(new_data0,aes(DAT,fpkm,fill=Group)) +
    geom_bar(stat="identity",position=position_dodge(0.9),color='black') +
    geom_errorbar(aes(ymin=fpkm, ymax=fpkm+se,group=Group), position = position_dodge(0.9),width=.1) +
    theme_bw() + theme(panel.grid=element_blank()) + theme_classic() +
    ggtitle(gene_id) + 
    ylab("Aroma content (ng/g FW)")+
    scale_fill_manual(name='Group',
                      values=c ( Control = "#7F7F7F", Treat = 'chartreuse3')) +
    theme(axis.title.x=element_text(size=15,vjust=-1,face = "bold",colour='black'),
          axis.title.y=element_text(size=15,vjust=3,face = "bold",colour='black'),
          axis.text.x=element_text(size=15,face = "bold",colour='black'),
          axis.text.y=element_text(size=15,face = "bold",colour='black'),
          axis.ticks.length = unit(6.5,'pt')) 
  p

  p_lst[[gene_id]] <- p
}
pdf(pcsvpath,
    height = 9,width=12)
do.call("ggarrange", c(p_lst, ncol=3, nrow=3))
dev.off()
###################################
### Fig 6e
library(ggplot2); library(ggpubr)

data=read.csv("20241021_lml_36_FPKM.csv",header = T,row.names = 1)
group=read.csv("group.csv",header = T,row.names = 1)
data1 = data[c('VIT_205s0020g03170','VIT_201s0011g04170','VIT_205s0049g00010'),]
data2 = as.data.frame(t(data1))
rownames(data2)
rownames(group)
data2 = cbind(data2,group)
data2$stage = factor(data2$stage, levels = c("DAT7","DAT14","DAT28"))
CSdata = data2[which(data2$cultivar=="CS"),]
CHdata = data2[which(data2$cultivar=="CH"),]
### for CS
p<- ggplot(CSdata, aes(x = stage, y = VIT_205s0020g03170, fill = group)) +  ### change the  y manually 
  geom_bar(stat = "summary", fun = mean, color = "black", position = position_dodge()) +
  scale_fill_manual(values=c(Control = "#7F7F7F", inoculation = "#66CD00"))+
  stat_summary(fun.data = 'mean_sd', geom = "errorbar", colour = "black",
               width = 0.25,position = position_dodge( .9))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p
### for CH
p<- ggplot(CHdata, aes(x = stage, y = VIT_205s0020g03170, fill = group)) +  ### change the  y manually 
  geom_bar(stat = "summary", fun = mean, color = "black", position = position_dodge()) +
  scale_fill_manual(values=c(Control = "#7F7F7F", inoculation = "#66CD00"))+
  stat_summary(fun.data = 'mean_sd', geom = "errorbar", colour = "black",
               width = 0.25,position = position_dodge( .9))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p

