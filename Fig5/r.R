####
## Fig 5a
library(sciplot)
library(doBy)
library(plotrix)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(tidyr)
library(reshape2)  
library(Rmisc)
library(xlsx)

fpkm1 <- read.csv('Berry_transcriptome_54sample.csv',row.names = 1)
don=fpkm1[, which(apply(fpkm1, 2, var) != 0)]

df_pca <- prcomp(don,scale. = TRUE) 
df_pcs <-data.frame(df_pca$x, Site=metadata$site, Cultivar=metadata$cultivar)
head(df_pcs)
df_pcs$Site = factor(df_pcs$Site,levels=c("HD","YSGD","CHFC","SSB","XG","HSP"))
df_pcs$Variety = factor(df_pcs$Cultivar,levels = c('CS','ME','CH'))

percentage<-round(df_pca$sdev^2 / sum(df_pca$sdev^2) * 100,2)
percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))


p <-ggplot(df_pcs,aes(x=PC1,y=PC2,color =Site,shape=Cultivar)) + 
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
p

######### WGCNA analysis
### Fig 5b

############ 分步法WGCNA分析
setwd('E:\\Data\\BerryMicromeGene2021\\WGCNA')
# install.packages("Rmisc")
library(xlsx)
library(ggplot2)
library(xlsx)
library(ggplot2)
library(ggpubr)
library(Rmisc) 
library(tidyr)
library(reshape2)  
library('WGCNA')

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

fpkm <- read.csv('Berry_transcriptome_54sample.csv',row.names = 1)
# metadata <- read.xlsx('E:\\Data\\BerryMicromeGene2021\\sample_lst.xlsx',row.names = 1,sheetIndex = 1)

dim(fpkm)

gsg_all = goodSamplesGenes(fpkm, verbose = 3)
gsg_all$allOK

if (!gsg_all$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg_all$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(fpkm)[!gsg_all$goodGenes], collapse = ", ")));
  if (sum(!gsg_all$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(fpkm)[!gsg_all$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  fpkm = fpkm[gsg_all$goodSamples, gsg_all$goodGenes]
}
gsg_all = goodSamplesGenes(fpkm, verbose = 3)
gsg_all$allOK
dim(fpkm)

sampleTree = hclust(dist(fpkm), method = "average")
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
### # 
abline(h = 45000, col = "red")

clust = cutreeStatic(sampleTree, cutHeight = 45000, minSize = 10)
table(clust)

keepSamples = (clust==1|clust==2)
datExpr = datExpr0[keepSamples, ]
dim(datExpr)

powers = c(c(1:10), seq(from = 12, to=20, by=2)) 
powers

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
##
pdf(" Analysis of network topology for various soft-thresholding powers.pdf",width = 9, height=5)

par(mfrow = c(1,2))

cex1 = 0.7
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red") 

sft$powerEstimate
# [1] 7 
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

softPower = 8
adjacency = adjacency(datExpr, power = softPower)

TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM


geneTree = hclust(as.dist(dissTOM), method = "average")

pdf(" Gene clustering on TOM-based dissimilarity.pdf",width = 12, height=9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

minModuleSize = 500
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 3, pamRespectsDendro = FALSE, minClusterSize = minModuleSize) 
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
pdf("Gene dendrogram and module colors .pdf",width = 8, height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()

MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes

MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")
plot(METree, main = "Clustering of module eigengenes0921", xlab = "", sub = "")

MEDissThres = 0
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
pdf("PlotsgeneDendro.pdf",width = 12, height=9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
save(MEs, moduleLabels, moduleColors, geneTree, file = "networkConstruction-stepByStep0913.RData")


datTraits0 <- read.csv("Berry_Microbe_Trait.csv",header=TRUE,row.names = 1)
dim(datTraits0)

datTraits=datTraits0[rownames(datExpr),]

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

table(moduleColors)

### Fig 4b
pdf("Module-trait associations_module.pdf",width = 7, height=5) 

textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(datTraits), yLabels = names(MEs), ySymbols = names(MEs), colorLabels = FALSE, colors = greenWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1), main = paste("Module-trait relationships"))

dev.off()
# 
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

genes_module_ann=as.data.frame(cbind(rownames(geneModuleMembership),moduleLabels,moduleColors))
write.csv(genes_module_ann,"gene_module_ann.csv")

##### Fig 4d
library(tidyverse)
library(igraph)
library(psych)
library(stringr)

correlate = function(other, metabo, route)
{

  result=data.frame(print(corr.test(other, metabo, use="pairwise", method="spearman", adjust="fdr", alpha=.05, ci=TRUE, minlength=100), short=FALSE, digits=5))

  result_raw=data.frame(print(corr.test(other, metabo, use="pairwise", method="spearman", adjust="none", alpha=.05, ci=TRUE, minlength=100), short=FALSE, digits=5))

  pair=rownames(result) 
  result2=data.frame(pair, result[, c(2, 4)]) 

  result4=data.frame(str_split_fixed(result2$pair, "-", 2), result2[, c(2, 3)], p_value=result_raw[, 4])
  colnames(result4)=c("feature_1", "feature_2", "r_value", "fdr_p_value", "raw_p_value")
  
  write.table(result4, file=paste(route, "Correlation_result.txt", sep="/"), sep="\t", row.names=F, quote=F)
}

dir.create("Result_v1") 

dat=read.csv("network_data.csv", header = T, check.names = F, row.names = 1) 
group=read.csv("group.csv", header = T, check.names = F) 


correlate(dat, dat, "Result_v1")

data = read.table("Result_v1/Correlation_result.txt", sep="\t", header=T)

data = data[data$r_value != 1,]

rownames(data) = 1:nrow(data)
box = paste(data$feature_1, data$feature_2, sep="-")
delete = c()
for(i in 1:nrow(data))
{
  tmp = paste(data[i, 2], data[i, 1], sep="-")
  box = box[-1]
  if(tmp %in% box)
  {
    delete = c(delete, i)
  }
}

data = data[-delete,]

data = data[data$raw_p_value <= 0.05,]
r_label = c()
for(i in 1:nrow(data))
{
  if(data[i, 3] < 0)
  {
    r_label = c(r_label, "neg")
  }
  else
  {
    r_label = c(r_label, "pos")
  }
}

data$r_label = r_label
write.table(data, file="Result_v1/input_pre.txt", sep="\t", row.names=F, quote=F)


data$r_value = abs(data$r_value)
input_network = data[, c(1,2,3,6)]
write.table(input_network, file="Result_v1/input_network.txt", sep="\t", row.names=F, quote=F)
degree=as.data.frame(table(input_network$feature_1))
degree1=as.data.frame(table(input_network$feature_2))
write.csv(degree,'node_degree.csv')

