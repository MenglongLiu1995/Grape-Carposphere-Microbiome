wd="E:\\Data\\BerryMicromeGene2021\\rare_curve"

# github安装包需要devtools，检测是否存在，不存在则安装
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
library(devtools)
# 检测amplicon包是否安装，没有从源码安装
if (!requireNamespace("amplicon", quietly = TRUE))
  install_github("microbiota/amplicon")
# library加载包，suppress不显示消息和警告信息
suppressWarnings(suppressMessages(library(amplicon)))

# Biconductor包安装，需要BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library("BiocManager")
# 检测amplicon包是否安装，没有从源码安装
p_list = c("phyloseq", "microbiome")
for(p in p_list){
  if (!requireNamespace(p, quietly = TRUE))
    BiocManager::install(p)}

#######
install(microbiome)
library(phyloseq)
library(microbiome)
library(amplicon)

otu <- read.csv("E:\\Data\\BerryMicromeGene2021\\microeco\\fungi\\basic_files_its\\feature_table.csv",header=T,check.names=FALSE ,row.names=1)
#样本信息
sample <- read.csv("E:\\Data\\BerryMicromeGene2021\\microeco\\fungi\\basic_files_its\\sample_table.csv",header=T,check.names=FALSE ,row.names=1)

# 构造phyloseq对象
ps = phyloseq(otu_table(otu, taxa_are_rows=TRUE), sample_data(sample))  
# 输入为Phyloseq的绘图
result = alpha_rare_all(ps = ps, group = "sample", method = "chao1",
                        start = 1000, step = 1000)
result[[3]]

(p = result[[4]])

# 按组均值绘图
(p = result[[3]])
ggsave(paste0("p4.rare_curve_group.pdf"), p, width=89*1.5, height=56, units="mm")
ggsave(paste0("p4.rare_curve_group.png"), p, width=89*1.5, height=56, units="mm")

