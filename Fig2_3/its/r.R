wd="microeco\\fungi"
setwd(wd)

library(igraph)
library(microeco) # Microbial Community Ecology Data Analysis
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics
library(dplyr)
library(magrittr)
library("ggh4x")

####

otu <- read.csv("otuits_berry.csv",header=T,check.names=FALSE ,row.names=1)

sample <- read.csv("metadata_its.csv",header=T,check.names=FALSE ,row.names=1)

tax <- read.csv("its_ann.csv", header=T,check.names=FALSE ,row.names=1)

tax %<>% tidy_taxonomy

df <- microtable$new(sample_table = sample,
                     otu_table = otu,
                     tax_table = tax,
                     auto_tidy = F)
df
###

df$filter_pollution(taxa = c("mitochondria", "chloroplast"))
df
# 
df$tidy_dataset()
df
# 
df$sample_sums() %>% range
# 
df$rarefy_samples(sample.size = 56353)
df$sample_sums() %>% range
 
df$save_table(dirpath = "basic_files", sep = ",")

df$cal_abund()
# return dataset$taxa_abund

df$cal_alphadiv(PD = FALSE)

# return dataset$alpha_diversity
class(df$alpha_diversity)
# 
df$cal_betadiv(unifrac = FALSE)
# return dataset$beta_diversity
# class(df$beta_diversity)
# 
# # save dataset$beta_diversity to a directory
dir.create("beta_diversity")
df$save_betadiv(dirpath = "beta_diversity")


# create trans_abund object
# use 12 Phyla with the highest abundance in the dataset.
df$sample_table$site %<>% factor(., levels = c("HD", "YSGD", "CHFC","SSB", "XG", "HSP"))
df$sample_table$cultivar %<>% factor(., levels = c("CH", "CS", "ME"))

####3 
t1 <- trans_abund$new(dataset = df, taxrank = "Phylum", ntaxa = 10)

### Fig 2e
cultivars_cycle <- c("CH", "CS", "ME")  

cultivar <- rep(cultivars_cycle, times = 66)  

t1$data_abund$cultivar <- cultivar
t1$data_abund$cultivar = factor(t1$data_abund$cultivar, levels = c("CS","ME","CH"))
site_cycle <- c("CHFC", "CHFC", "CHFC","HD","HD","HD","HSP","HSP","HSP", "SSB","SSB","SSB",'XG','XG','XG','YSGD','YSGD','YSGD')  
site = rep(site_cycle, times = 11)  
t1$data_abund$site <- site
t1$data_abund$site = factor(t1$data_abund$site, levels = c("HD",'YSGD','CHFC','SSB','XG','HSP'))

p1=t1$plot_bar(others_color = "grey70",facet = c("site","cultivar"),legend_text_italic = FALSE,
               xtext_keep = FALSE,barwidth = 1.5)
pdf("bar_p_f.pdf",height = 4,width = 7.5)
p1
dev.off()



### Fig 2f
# Then we show the heatmap with the high abundant genera.

t1 <- trans_abund$new(dataset = df, taxrank = "Family", ntaxa = 10,groupmean = "sample")
t1$data_abund$Sample=factor(t1$data_abund$Sample,levels = c("HD_CS","HD_ME","HD_CH",
                                                            "YSGD_CS","YSGD_ME","YSGD_CH",
                                                            "CHFC_CS","CHFC_ME","CHFC_CH",
                                                            "SSB_CS","SSB_ME","SSB_CH",
                                                            "XG_CS","XG_ME","XG_CH",
                                                            "HSP_CS","HSP_ME","HSP_CH"))

###
cultivars_cycle <- c("CH", "CS", "ME")  

cultivar <- rep(cultivars_cycle, times = 966)  

t1$data_abund$cultivar <- cultivar
t1$data_abund$cultivar = factor(t1$data_abund$cultivar, levels = c("CS","ME","CH"))
site_cycle <- c("CHFC", "CHFC", "CHFC","HD","HD","HD","HSP","HSP","HSP", "SSB","SSB","SSB",'XG','XG','XG','YSGD','YSGD','YSGD')  
site = rep(site_cycle, times = 161)  
t1$data_abund$site <- site
t1$data_abund$site = factor(t1$data_abund$site, levels = c("HD",'YSGD','CHFC','SSB','XG','HSP'))
###
p2=t1$plot_heatmap(facet = c("site","cultivar"), withmargin = FALSE,
                   ytext_size = 10,strip_text = 10,  xtext_keep = FALSE)
p2
pdf("Heat_f_f.pdf",height = 3,width = 7.9)
p2
dev.off()

#####  shannon index   Fig 2g
t1 <- trans_alpha$new(dataset = df, group = "cultivar")
t1$cal_diff(method = "KW")
# return t1$res_diff
head(t1$res_diff)
# p.adi=0.77524834  ns

t1 <- trans_alpha$new(dataset = df, group = "site")
t1$cal_diff(method = "KW")
# return t1$res_diff
head(t1$res_diff)
# 1.307321e-06 

library(ggpubr)

###Fig 3g 

data=read.csv("alpha_diversity\\alpha_diversity.csv",header = T,
              row.names = 1,encoding="UTF-8",check.names = F)
colnames(data)
rownames(data)
data$site = factor(df$sample_table$site, levels = c("HD",'YSGD','CHFC','SSB','XG','HSP'))
data$cultivar = factor(df$sample_table$cultivar, levels = c("CS","ME","CH"))

  p<- ggplot(data, aes_string(x = "site", y = 'Shannon',fill="cultivar")) + 
    geom_bar(stat = "summary", fun = mean, color = "black", position = position_dodge()) +
    scale_fill_manual(values=c(CH = "gold2", CS = "darkorchid3", ME = "dodgerblue1"))+
    stat_summary(fun.data = 'mean_sd', geom = "errorbar", colour = "black",
                 width = 0.25,position = position_dodge( .9))+
    labs(y=paste("Shannon"))+
    ggtitle(("site: p.adj=1.31 e-6 *** ; cultivar: p.adj=0.78  ns"))
  p

 
  
pdf("Shannon_its.pdf",width=6,height=4)
p
dev.off()
#############################
### Fig 3h
tmp <- df$merge_samples("site")
tmp
t1$data_summary

t1 <- trans_venn$new(dataset = tmp)

# only show some sets with large intersection numbers
t1$data_summary %<>% .[.[, 1] > 20, ]
g1 <- t1$plot_bar(left_plot = T, bottom_height = 0.5, left_width = 0.15, up_bar_fill = "grey50", left_bar_fill = "grey50", bottom_point_color = "black")
g1
# g1 is aplot class and can be saved with ggplot2::ggsave, aplot::ggsave or cowplot::save_plot function
# as g1 is comprised of several sub-plots, please adjust the details for each sub-plot

pdf('Upset_f.pdf',height = 4,width = 6)
g1
dev.off()

#############################
# g1 is aplot class and can be saved with ggplot2::ggsave, aplot::ggsave or cowplot::save_plot function
# as g1 is comprised of several sub-plots, please adjust the details for each sub-plot
dataset1 <- df$merge_samples("site")
t1 <- trans_venn$new(dataset1)

# transform venn results to the sample-species table, here do not consider abundance, only use presence/absence.
t2 <- t1$trans_comm(use_frequency = TRUE)
# t2 is a new microtable class, each part is considered a sample

# calculate taxa abundance, that is, the frequency
t2$cal_abund()
# transform and plot
t3 <- trans_abund$new(dataset = t2, taxrank = "Family", ntaxa = 10) 
unique(t3$data_abund$Group)

dat_venn=as.data.frame(t3$data_abund)
dat_venn1=dat_venn[which(dat_venn$Group=="CHFC&HD&HSP&SSB&XG&YSGD"),]
write.csv(dat_venn1,'dat_venn2_f.csv')  ### data for dount plot


## beta 
### Fig 3f
# create an trans_beta object
# measure parameter must be one of names(dataset$beta_diversity)
t1 <- trans_beta$new(dataset = df, group = "site", measure = "bray")

# PCoA, PCA and NMDS are available
t1$cal_ordination(ordination = "NMDS")
# t1$res_ordination is the ordination result list
class(t1$res_ordination)
# plot the PCoA result with confidence ellipse
p1= t1$plot_ordination(plot_color = "site", plot_shape = "cultivar",
                   shape_values = c(CS=16,ME=17,CH=15),
                   plot_type = c("point"),
                   # group_order = c("HD","YSGD","CHFC","SSB", "XG","HSP"),
                   color_values=c(HD='#20a391', 
                                  YSGD= 'dodgerblue4', 
                                  CHFC= '#a00528', 
                                  SSB= '#fc7839', 
                                  XG= '#580346', 
                                  HSP= 'palegreen3' ))+
  ggtitle("Stress=0.13")

pdf("NMDS_f.pdf",height = 4,width = 5)
p1
dev.off()
####  manova result
t1 <- trans_beta$new(dataset = df, group = "site", measure = "bray")
# manova for all groups when manova_all = TRUE
t1$cal_manova(manova_all = TRUE)
t1$res_manova
# adonis2(formula = use_formula, data = metadata)
# Df SumOfSqs      R2      F Pr(>F)    
# site      5   5.1384 0.40931 6.6523  0.001 ***
##
t1 <- trans_beta$new(dataset = df, group = "cultivar", measure = "bray")
# manova for all groups when manova_all = TRUE
t1$cal_manova(manova_all = TRUE)
t1$res_manova
# adonis2(formula = use_formula, data = metadata)
# Df SumOfSqs      R2      F Pr(>F)
# cultivar  2   0.5454 0.04345 1.1582  0.252

#########################
# Fig.S3 different abundance of OTUs between site and cultivar

### Fig 3g
otu = read.csv("basic_files_its\\feature_table.csv",row.names = 1, header = T,check.names = F)
tax = read.csv("basic_files_its\\tax_table.csv",row.names = 1, header = T,check.names = F)
sample = read.csv("basic_files_its\\sample_table.csv",row.names = 1, header = T,check.names = F)

otu <- otu[, rownames(sample)]

dds <- DESeqDataSetFromMatrix(
  countData = otu,
  colData = sample,
  design = ~ cultivar 
)

keep <- rowSums(counts(dds) >= 10) >= 5
dds <- dds[keep, ]

dds <- DESeq(dds)

dds_lrt <- DESeq(dds, test="LRT", reduced=~1)
res_lrt <- results(dds_lrt, alpha=0.05)

sig_otus <- rownames(res_lrt)[which(res_lrt$padj < 0.05)]

#### the propotion of RA.otu for diff_cultivar
otu_RA <- otu/colSums(otu)
mean_otu = as.data.frame(rowMeans(otu_RA))
sum(mean_otu)

sum(mean_otu[sig_otus,])
# [1] 0.205817
#################
### for  sites

dds <- DESeqDataSetFromMatrix(
  countData = otu,
  colData = sample,
  design = ~ site  
)

keep <- rowSums(counts(dds) >= 10) >= 5
dds <- dds[keep, ]

dds <- DESeq(dds)

dds_lrt <- DESeq(dds, test="LRT", reduced=~1)
res_lrt <- results(dds_lrt, alpha=0.05)

sig_otus <- rownames(res_lrt)[which(res_lrt$padj < 0.05)]

otu_RA <- otu/colSums(otu)
mean_otu = as.data.frame(rowMeans(otu_RA))
sum(mean_otu)

sum(mean_otu[sig_otus,])
# [1] 0.6222026

otu_RA_ave = otu_RA_ave[significant_otu,]

#  
# Fig.3g Tern plot
otu_RA = as.data.frame(t(otu_RA))
value_cols <- colnames(otu_RA)  ### calculate each OTU


result <- otu_RA %>%  
  group_by(sample$cultivar) %>%  
  summarise(across(all_of(value_cols), mean, .names = "{col}"))
result = as.data.frame(result)
rownames(result) = result[,1]
result = result[,c(-1)]

tern_res = as.data.frame(t(result))
tern_res = tern_res[sig_otus,]

otu_RA_ave = as.data.frame(otu_RA_ave)
tern_res$RA = otu_RA_ave$otu_RA_ave

tax = tax[sig_otus,]

tern_res$Phylum = tax$Phylum
tern_res$Family = tax$Family
tern_res$Class = tax$Class
tern_res$Mean = otu_RA_ave[sig_otus,]
str(tern_res)

###
library(ggtern)
ggtern(data = tern_res, aes(x = CS, y = ME, z = CH)) +  
  geom_mask() +  
  geom_point(aes(size = Mean, color = Family)) +  # 
  theme_showarrows() +  
  scale_color_brewer(palette = "Paired")+
  scale_size_area(max_size = 10)+
  theme(
    text = element_text(size = 14),       
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12),  
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )

# Fig S2b
# use Genus level for parameter taxa_level, if you want to use all taxa, change to "all"
# install.packages("randomForest")
library("randomForest")
##
# nresam = 1 and boots = 1 represent no bootstrapping and use all samples directly
t1 <- trans_diff$new(dataset = df, method = "rf", group = "cultivar", taxa_level = "Family",p_adjust_method = "none")
# plot the MeanDecreaseGini bar
g1 <- t1$plot_diff_bar(use_number = 1:10, group_order = c("CS", "ME", "CH"),
                       color_values=c(CH = "gold2", CS = "darkorchid3", ME = "dodgerblue1"))

# plot the abundance using same taxa in g1
g2 <- t1$plot_diff_abund(group_order = c("CS", "ME", "CH"), 
                         select_taxa = t1$plot_diff_bar_taxa,
                         color_values=c(CH = "gold2", CS = "darkorchid3", ME = "dodgerblue1"))

g1 <- g1 + theme(legend.position = "none")
g2 <- g2 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
gridExtra::grid.arrange(g1, g2, ncol = 2, nrow = 1, widths = c(2, 1.7))

pdf("rm_f_cultivar.pdf",height = 5,width = 7)
gridExtra::grid.arrange(g1, g2, ncol = 2, nrow = 1, widths = c(2, 1.7))
dev.off()

# Fig S2d
# nresam = 1 and boots = 1 represent no bootstrapping and use all samples directly
t1 <- trans_diff$new(dataset = df, method = "rf", group = "site", taxa_level = "Family",p_adjust_method = "fdr")
# plot the MeanDecreaseGini bar
# group_order is designed to sort the groups
g1 <- t1$plot_diff_bar(use_number = 1:10, group_order = c("HD","YSGD","CHFC","SSB", "XG","HSP"),
                       color_values=c(HD='#20a391', 
                                      YSGD= 'dodgerblue4', 
                                      CHFC= '#a00528', 
                                      SSB= '#fc7839', 
                                      XG= '#580346', 
                                      HSP= 'palegreen3' ))
# plot the abundance using same taxa in g1
g2 <- t1$plot_diff_abund(group_order = c("HD","YSGD","CHFC","SSB", "XG","HSP"), 
                         select_taxa = t1$plot_diff_bar_taxa,
                         color_values=c(HD='#20a391', 
                                        YSGD= 'dodgerblue4', 
                                        CHFC= '#a00528', 
                                        SSB= '#fc7839', 
                                        XG= '#580346', 
                                        HSP= 'palegreen3' ))

# now the y axis in g1 and g2 is same, so we can merge them
# remove g1 legend; remove g2 y axis text and ticks
g1 <- g1 + theme(legend.position = "none")
g2 <- g2 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
gridExtra::grid.arrange(g1, g2, ncol = 2, nrow = 1, widths = c(2, 1.7))

pdf("rm_b_site2.pdf",height = 5,width = 7)
gridExtra::grid.arrange(g1, g2, ncol = 2, nrow = 1, widths = c(2, 1.7))

dev.off()

####
#### RDA  ANALYSIS
 ### Fig 3h
## Explainable class env
env_data_16S=read.csv('soil_data_20_depth_Climate2.csv', row.names = 1,header = T)
str(env_data_16S)
rownames(env_data_16S)
new_test <- clone(df)
new_test$sample_table <- data.frame(new_test$sample_table, env_data_16S[rownames(new_test$sample_table), ])
# now new_test$sample_table has the whole data
colnames(new_test$sample_table)
# let’s use env_cols to select the required columns from sample_table in the microtable object.
t1 <- trans_env$new(dataset = new_test, env_cols = 4:34,standardize = FALSE)

### RDA analysis
t1$cal_ordination(method = "RDA", taxa_level = "Family")
# select 10 features and adjust the arrow length
t1$trans_ordination(show_taxa = 10, adjust_arrow_length = TRUE, max_perc_env = 1.5, max_perc_tax = 1.5, min_perc_env = 0.2, min_perc_tax = 0.2)
# t1$res_rda_trans is the transformed result for plot
p1=t1$plot_ordination(plot_color = "site",plot_shape = "cultivar",
                      shape_values = c(CS=16,ME=17,CH=15),
                      color_values=c ( HD='#20a391', 
                                       YSGD= 'dodgerblue4', 
                                       CHFC= '#a00528', 
                                       SSB= '#fc7839', 
                                       XG= '#580346', 
                                       HSP= 'palegreen3' ))
pdf("RDA_f.pdf",height = 4,width = 6)
p1
dev.off()

### Fig 3d 
# Decay_curve
# first cal the bray similarity 
library(ggplot2)
library(ggvegan)
dat_f=read.csv("\\basic_files_its\\feature_table.csv",row.names = 1)
dat_f=as.data.frame(t(dat_f))
distance<-vegdist(dat_f,method='bray')
bray_similarity=1-distance
bray_similarity <- as.matrix(bray_similarity)  
write.table(bray_similarity,"bray_similarity_fungi.txt",sep = "\t")

### read the bray-dist data of bacteria
bray_dist0=read.csv("bray_similarity_fungi.csv",row.names = 1,header = T)
var_names=rownames(bray_dist0)

## adjust the data frame, for cal the cor between bray_bactertia and other dataset
# 
result_df <- data.frame()  

#  
for (i in 1:(nrow(bray_dist0) - 1)) {  
  for (j in (i + 1):nrow(bray_dist0)) {  
    # 
    result_df <- rbind(result_df, data.frame(  
      Var1 = var_names[i],  
      Var2 = var_names[j],  
      Correlation = bray_dist0[i, j]  
    ))  
  }  
}  

# 
print(result_df)
### read the bray-dist data of env
env=read.csv("RRC curve\\soil_data_20_depth_Climate2.csv",row.names = 1,header = T)
colnames(env)
climate=env[,1:10]
soil=env[,11:30]

dist_Climate<-vegdist(climate,method='euclidean')
dist_Climate=as.matrix(dist_Climate)

## adjust the data frame, for cal the cor between bray_fungi and other dataset
# 
result_climate <- data.frame()  

# 
for (i in 1:(nrow(dist_Climate) - 1)) {  
  for (j in (i + 1):nrow(dist_Climate)) {  
    # 
    result_climate <- rbind(result_climate, data.frame(  
      Var1 = rownames(dist_Climate)[i],  
      Var2 = rownames(dist_Climate)[j],  
      Correlation = dist_Climate[i, j]  
    ))  
  }  
}  

#   
print(result_climate)

###
## adjust the data frame, for cal the cor between bray_fungi and other dataset
# 
result_climate <- data.frame()  

#  
for (i in 1:(nrow(dist_Climate) - 1)) {  
  for (j in (i + 1):nrow(dist_Climate)) {  
    # 
    result_climate <- rbind(result_climate, data.frame(  
      Var1 = rownames(dist_Climate)[i],  
      Var2 = rownames(dist_Climate)[j],  
      Correlation = dist_Climate[i, j]  
    ))  
  }  
}  

# 查看结果  
print(result_climate)

#### geo_dist
library("geosphere")
lonlat = read.csv("NX_site_lon_lat1.csv",row.names = 1,header = T);
# matrix
muer.lonlat = cbind(lonlat$Lon, lonlat$Lat);
# 精确计算，椭圆
muer.dists = distm(muer.lonlat, fun=distVincentyEllipsoid);

rownames(muer.dists) = rownames(lonlat);
colnames(muer.dists) = rownames(lonlat);

### 
result_geo <- data.frame()  

# 
for (i in 1:(nrow(muer.dists) - 1)) {  
  for (j in (i + 1):nrow(muer.dists)) {  
    # 
    result_geo <- rbind(result_geo, data.frame(  
      Var1 = rownames(muer.dists)[i],  
      Var2 = rownames(muer.dists)[j],  
      Correlation = muer.dists[i, j]  
    ))  
  }  
}  

# 
print(result_geo)


#### 
## adjust the data frame, for cal the cor between bray_fungi and other dataset

dist_Soil<-vegdist(soil,method='euclidean')
dist_Soil=as.matrix(dist_Soil)

# 
result_soil <- data.frame()  

# 
for (i in 1:(nrow(dist_Soil) - 1)) {  
  for (j in (i + 1):nrow(dist_Soil)) {  
    # 
    result_soil <- rbind(result_soil, data.frame(  
      Var1 = rownames(dist_Soil)[i],  
      Var2 = rownames(dist_Soil)[j],  
      Correlation = dist_Soil[i, j]  
    ))  
  }  
}  

# 
print(result_soil)

### merge data 
RC_data=cbind(result_df,result_geo[,3],result_climate[,3],result_soil[,3])
colnames(RC_data)
colnames(RC_data)[3]="dist_f"
colnames(RC_data)[4]="dist_geo"
colnames(RC_data)[5]="dist_climate"
colnames(RC_data)[6]="dist_soil"


RC_data$dist_f=RC_data$dist_f*100  # change to the %
RC_data$dist_geo=RC_data$dist_geo*0.001  # from unit (m) to (km)

write.csv(RC_data,'RC_data_fungi.csv')

RC_data = read.csv("RC_data_fungi.csv", row.names = 1, header = T, check.names = F)
#####
install.packages("ggExtra")
install.packages("Matrix",type="binary")
library(ggpmisc)
library(ggplot2)
library(ggpubr)
library(ggExtra)

formula <- y ~ x 
formula=y ~ poly(x, 2) ### for regression between fungi and soil 

colnames(RC_data)
p = ggplot(RC_data,aes(x = dist_climate,y = dist_f)) + 
  geom_point(size=4) + ## aes(color=site,shape=variety),
  # stat_smooth( method = "lm") +
  geom_smooth(method="lm", formula=y ~ poly(x, 2), se=TRUE, color="blue")+  ### for regression between fungi and soil
  stat_cor(label.y = max(RC_data$dist_f)*0.95,size=10) + 
  stat_poly_eq(aes(label = ..eq.label..),size=10,
               formula = formula, parse = TRUE) +
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5,size=15),
        axis.text.y = element_text(size = 15),
        axis.title =  element_text(size = 20),
        legend.title = element_text(size=10),
        legend.text  = element_text(size=10))+
  theme(plot.title = element_text(hjust = 0.5,size=15))+ 
  labs(x="Geographic distance (km)", y= "Community similarity (%)") +
  # labs(x="Climatic distance (Euclidean)", y= "Community similarity (%)") +  #
  # labs(x="Edaphic distance (Euclidean)", y= "Community similarity (%)") +  #
  scale_y_continuous(limits = c(0,100))
p
p1 <- ggMarginal(p, type="histogram")
p1
# p2 <- ggMarginal(p, type="density")
# p2
pdf('RC_dist_geo_f1.pdf',width = 6,height =5.5)
print(p1)
dev.off()
############
### residual plot
colnames(RC_data)
# 
model <- lm(dist_f ~ dist_geo, data = RC_data)     ####  geo, climate
model <- lm(dist_f ~ poly(dist_soil, 2), data = RC_data)  ####  soil
 
RC_data$fitted_values <- fitted(model) 
RC_data$residuals <- residuals(model)  


ggplot(RC_data, aes(x = fitted_values, y = residuals)) +  
  geom_point(color = "blue") +  
  geom_hline(yintercept = 0, color = "red") +  
  labs(x = "Fitted Values",  
       y = "Residuals",  
       title = "fungi~dist_soil") +  
  theme_minimal()

