################### berry metabolites pca 
## Fig 4a
don1=read.csv("BerryTrait.csv",row.names = 1,header = T)
colnames(don2)
colnames(don1)
don1$site=don2$site
don1$cultivar=don2$cultivar


Var=don1[,1:38]
Var = na.omit(Var)

df_pca <- prcomp(Var,scale. = TRUE) 
df_pcs <-data.frame(df_pca$x, Site=don1$site,Variety=don1$cultivar)

df_pcs$Site = factor(df_pcs$Site,levels=c("HD","YSGD","CHFC","SSB","XG","HSP"))
df_pcs$Variety = factor(df_pcs$Variety,levels = c('CS','ME','CH'))
##1-2
percentage<-round(df_pca$sdev^2 / sum(df_pca$sdev^2) * 100,2)
percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))


p <-ggplot(df_pcs,aes(x=PC1,y=PC2,color =Site,shape=Variety)) + 
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

features0 = as.data.frame(df_pca$rotation)
features0$features = row.names(features0)
features0$features = factor(features0$features,levels = row.names(features0))
write.csv(features0,"features0.csv") 
features0=read.csv("features0.csv",row.names = 1, header = T) ### added the group label for each variables
unique(features0$group)

## Fig 4b
p <- ggplot(features0,aes(x=PC1,y=PC2,label=features,color=group)) + geom_point()+
  geom_segment(data = features0, aes(x = 0, xend = PC1, y=0, yend = PC2, color = features), 
               size = 1, arrow = arrow(length = unit(0.02, "npc")), alpha = 1)+
  scale_color_manual(values = c(Maturity="#96C37D", Totalanthocyanins="#FA7F6F", Elementcomposition='gold4',Aroma='#82B0D2'))+
  xlab(percentage[1]) + ylab(percentage[2]) +
  geom_text_repel(size = 4) +
  geom_hline(yintercept = 0, linetype = "dotted",size =1.0) +
  geom_vline(xintercept = 0, linetype = "dotted",size =1.0) +
  theme_bw(base_size = 15) +theme(plot.title = element_text(hjust = 0.8))+
  theme(legend.position = "none")
p

#########  VPA analysis
## Fig 4d
library(vegan)

#
b_otu <- read.csv("feature_table_b.csv",header = T, check.names = F,row.names = 1)
f_otu <- read.csv("feature_table_f.csv",header = T, check.names = F,row.names = 1)
berry=read.csv('BerryTrait.csv',header = T,check.names = F,row.names = 1)
env <- read.csv("soil_data_20_depth_Climate2.CSV",row.names = 1,header = T)

rownames(berry)
rownames(b_otu)
#
climate=env[,1:10]
soil=env[,11:30]

climate=prcomp(climate, scale=T,center = T)
summary(climate) #PC3 0.9078
climate_pca=climate$x[,3]

soil=prcomp(soil, scale=T,center = T)
summary(soil) #PC8 0.92259
soil_pca=soil$x[,1:8]

b_otu_pca=prcomp(b_otu,scale=T,center = T)
summary(b_otu_pca)## PC37 0.90294
b_otu_pca1=b_otu_pca$x[,1:37] 

f_otu_pca=prcomp(f_otu,scale=T,center = T)
summary(f_otu_pca)## PC32 0.90455
f_otu_pca1=f_otu_pca$x[,1:32]

colnames(berry[,c(6,29,2,17,31,32,24,33,15,18)])
berry_var=berry[,c(6,29,2,17,31,32,24,33,15,18)]

#### climate soil and bacteria
rda.vpa <- varpart(berry[,c(6,29,2,17,31,32,24,33,15,18)], 
                   climate_pca,soil_pca,b_otu_pca1,transfo = "hellinger",chisquare = FALSE) 
rda.vpa
plot(rda.vpa, cex = 1, bg = c("skyblue","goldenrod1","hotpink"),
     Xnames=c("Climate","Soil","Bacteria"))  


#### climate soil and fungi 
rda.vpa <- varpart(berry[,c(6,29,2,17,31,32,24,33,15,18)], 
                   climate_pca,soil_pca,f_otu_pca1,transfo = "hellinger",chisquare = FALSE) 

rda.vpa
plot(rda.vpa, cex = 1, bg = c("skyblue","goldenrod1","hotpink"),
     Xnames=c("Climate","Soil","Fungi"))  

################   Heatmap
### Fig 4e
berry=read.csv('berry_var.csv', row.names = 1,header = T)
Microbe=read.csv('GenusAbun.csv', row.names = 1,header = T)
library(psych)
library(stringr)

correlate = function(other, metabo, route)
{
  
  result=data.frame(print(corr.test(other, metabo, use="pairwise", method="spearman", adjust="fdr", alpha=.05, ci=TRUE, minlength=50), short=FALSE, digits=5))
  

  pair=rownames(result)  
  result2=data.frame(pair, result[, c(2, 4)]) 
  

  result3=data.frame(result2[order(result2[,"raw.p"], decreasing=F),])
  
  result4=data.frame(str_split_fixed(result3$pair, "-", 2), result3[, c(2, 3)])
  colnames(result4)=c("feature_1", "feature_2", "r_value", "p_value")
  
  write.table(result4, file=paste(route, "Correlation_result.txt", sep="/"), sep="\t", row.names=F, quote=F)
}

dir.create("Result")

correlate(berry,Microbe ,"Result")
###
library(reshape2)
library(pheatmap)

correlate_pheatmap = function(infile,route)
{
  data=read.table(paste(route, infile, sep="/"), sep="\t", header=T)
  
  data_r=dcast(data, feature_1 ~ feature_2, value.var="r_value")
  data_p=dcast(data, feature_1 ~ feature_2, value.var="p_value")
  rownames(data_r)=data_r[,1]
  data_r=data_r[,-1]
  rownames(data_p)=data_p[,1]
  data_p=data_p[,-1]
  
  data_mark=data_p
  for(i in 1:length(data_p[,1])){
    for(j in 1:length(data_p[1,])){
      data_mark[i,j]=ifelse(data_p[i,j] <= 0.05, "*", "")
    }
  }
  
  pheatmap(data_r,display_numbers=data_mark, cluster_rows = TRUE, cluster_cols = TRUE, cellwidth=20, cellheight=20, fontsize_number=18, filename=paste(route, "Correlation_result.pdf", sep="/"))
  pheatmap(data_r, display_numbers=data_mark,cluster_rows = TRUE, cluster_cols = TRUE, cellwidth=20, cellheight=20, fontsize_number=18, filename=paste(route, "Correlation_result.png", sep="/"))
}

#
correlate_pheatmap("Correlation_result.txt", "Result") 

####### GRaMM
##Fig 4f


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
###########
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
##########
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
          for(k in seq_len(200))
          {
            bootx <- matrix(sample(x1, replace = FALSE))
            booty <- matrix(sample(y1, replace = FALSE))
            tmp <- minerva::mine(bootx, booty)
            MICtp <- tmp$MIC
            tempM <- ifelse(micr <= MICtp,1,0)
            micp <- tempM + micp
          }
          micp <- micp/200
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
##########
library(minerva)

otu_raw=read.csv("feature_table.csv",row.names = 1,check.names = F)
head(otu_raw)

tax_lst=read.csv("tax_table.csv", row.names = 1, header = T)
otu_lst=tax_lst[which(tax_lst$Genus=="g__Sphingomonas"),]

berry_var=read.csv("berry_var.csv",
                   row.names = 1,header = T,check.names = F )

sample_sums <- rowSums(otu_raw) 

relative_abundance <- otu_raw / sample_sums *100
rela_otu_ave <- apply(relative_abundance, MARGIN = 2, FUN = mean)
otu_filter = rela_otu_ave[rela_otu_ave > 0.02]

sphingmonas_lst = which(rela_otu_ave[rownames(otu_lst)]>0.02) ### OTU in Sphingmoas 

relative_abundance_filter = relative_abundance[,names(otu_filter)] ## relative OTU abund for the next CLR transformation

otu_clr <- lapply(relative_abundance_filter, function(x) ifelse(x == 0, 1e-10, x))

otu_clr = as.data.frame(otu_clr)

rownames(otu_clr) = rownames(relative_abundance_filter)

write.csv(otu_clr,"otu_filtered_clr.csv")
microbe_Prepro = microbesPrepro(otu_clr,missPro = TRUE,missMethod = "mean",
                                rarePro = F, scalingPro = F, transPro = T, 
                                kValue = 3,phenoData = NA,phenoDataType = "continuous")

otu_RA_CLR_sphingmoas = microbe_Prepro[,names(sphingmonas_lst)]

write.csv(otu_RA_CLR_sphingmoas,"otu_RA_CLR_sphingmoas.csv")

mebolite_Prepro = metabolitesPrepro(berry_var, missPro = TRUE, missMethod = "QRILC",
                                    scalingPro = F, transPro = T, kValue = 3,
                                    phenoData = NA)

res = naivegramm(mebolite_Prepro, otu_RA_CLR_sphingmoas, covdata = NA, r = 0.3,alpha = 0.99,pheno="All")
res_r = res$r_All
res_p = res$p_All

res_type = res$type_All

otu_abun = colMeans(otu_RA_CLR_sphingmoas, na.rm = T)

otu_abun = as.data.frame(otu_abun)
unique(otu_abun$otu_abun)
otu_abun = otu_abun[rownames(res_p),]

res_p_abun = cbind(otu_abun,res_p)
# write.csv(res_p_abun,'res_p_abun.csv')

library(ggpubr)
library(ggrepel)
library(ggbreak)
colnames(res_p_abun)

#### otu~ cis_3_Hexenal
p1 = ggplot(res_p_abun,aes(x=otu_abun,y=cis_3_Hexenal))+
  geom_point(aes(size=otu_abun,color=res_p_abun[,'cis_3_Hexenal']),alpha=0.6)+
  scale_size(range=c(1,6))+
  scale_color_gradient(low = "red",high = "green") +
  theme_bw()+
  scale_y_continuous(limits = c(0,0.075))+
  geom_text_repel(data = subset(res_p_abun, res_p_abun[,'cis_3_Hexenal'] <= 0.05),label = rownames(subset(res_p_abun, res_p_abun[,'cis_3_Hexenal'] <= 0.05)) ,size = 3,segment.color = "black",max.overlaps=20, 
                  show.legend = FALSE )+
  labs(title ="cis_3_Hexenal", colour=expression("GRaMM"),size="Abundance",x="OTU abundance (CLR transformed)",
       y="P-value")+
  geom_hline(yintercept=c(0.05),lty = 2, lwd = 1)
p1


pdf('Gramm_1_cis_3_Hexenal.pdf', height = 4, width = 4)
p1
dev.off()
#### otu~ "trans_2_Hexenal"
p1= ggplot(res_p_abun,aes(x=otu_abun,y=trans_2_Hexenal))+
  geom_point(aes(size=otu_abun,color=res_p_abun[,'trans_2_Hexenal']),alpha=0.6)+
  scale_size(range=c(1,6))+
  scale_color_gradient(low = "red",high = "green") +
  theme_bw()+
  scale_y_continuous(limits = c(0,0.075))+
  geom_text_repel(data = subset(res_p_abun, res_p_abun[,'trans_2_Hexenal'] <= 0.05),label = rownames(subset(res_p_abun, res_p_abun[,'trans_2_Hexenal'] <= 0.05)) ,size = 3,segment.color = "black",max.overlaps=10, 
                  show.legend = FALSE )+
  labs(title ="trans_2_Hexenal", colour=expression("GRaMM"),size="Abundance",x="OTU abundance (CLR transformed)",
       y="P-value")+
  geom_hline(yintercept=c(0.05),lty = 2, lwd = 1)
p1
# 
pdf('Gramm_2_trans_2_Hexenal.pdf', height = 4, width = 4)
p1
dev.off()

#### otu~ "trans_2_Pentenal"
p1= ggplot(res_p_abun,aes(x=otu_abun,y=trans_2_Pentenal))+
  geom_point(aes(size=otu_abun,color=res_p_abun[,'trans_2_Pentenal']),alpha=0.6)+
  scale_size(range=c(1,6))+
  scale_color_gradient(low = "red",high = "green") +
  theme_bw()+
  scale_y_continuous(limits = c(0,0.075))+
  geom_text_repel(data = subset(res_p_abun, res_p_abun[,'trans_2_Pentenal'] <= 0.05),label = rownames(subset(res_p_abun, res_p_abun[,'trans_2_Pentenal'] <= 0.05)) ,size = 3,segment.color = "black",max.overlaps=10, 
                  show.legend = FALSE )+
  labs(title ="trans_2_Pentenal", colour=expression("GRaMM"),size="Abundance",x="OTU abundance (CLR transformed)",
       y="P-value")+
  geom_hline(yintercept=c(0.05),lty = 2, lwd = 1)
p1
pdf('Gramm_3_trans_2_Pentenal.pdf', height = 4, width = 4)
p1
dev.off()

####### Supplemental S4
######### Supplemental Fig S4
############### Genus0.1~berry quality
library(ggplot2)
gene_filter1=read.csv("GenusAbun.csv",row.names = 1,check.names = F,header = T)
rownames(gene_filter1)
g_lst = c("f_Moraxellaceae_g_Cavicella","f_Atopobiaceae_g_","f_Gemmatimonadaceae_g_Gemmatimonas",
          "f_Sphingomonadaceae_g_Sphingomonas","c_Subgroup_6_g_","f_Muribaculaceae_g_","f_Lachnospiraceae_g_")
gene_filter2=as.data.frame(t(gene_filter1))
gene_filter2 = gene_filter2[,g_lst]


berry_var=read.csv("berry_var.csv",row.names = 1,header = T,check.names = F )
colnames(berry_var)
aroma_lst = c("cis_3_Hexenal","trans_2_Hexenal","trans_2_Pentenal")
berry_var1 = berry_var[,aroma_lst]

sample=read.csv("sample_table.csv",row.names = 1,header = T,check.names = F )

formula <- y ~ x

lst <- list() ##用于批量画图的排版
n <- 0 ##  用于批量画图的排版
filter_gene <- c()
r <- c()
p_val <- c()
for (i in colnames(berry_var1)) {
  for (j in colnames(gene_filter2)) {
    n=n+1
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


