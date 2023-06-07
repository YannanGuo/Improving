#####################################################
################  UTF-8 @ R 4.2.0  ##################
###############  Author: YannanGuo  #################
#############  From DEGSeq2 Analysis  ###############
############  Created in May 10, 2023  ##############
############  Pre-install Packages:  ################
## Bioconductor version 3.16 (BiocManager 1.30.20)###
##Pheatmap version 1.0.12############################
##RColorBrewer version 1.1-3#########################
##glmpca version 0.2.0###############################
##ggplot2 version 3.4.1##############################
##gggorce version 0.4.1##############################
##concaveman version 1.1.0###########################
##reshape2 version 1.4.4#############################
#####################################################

#Clean the environment
rm(list=ls())

#Import the data
setwd("/home/data/t150327/HT22-RNAseq/DGE/DGEAnalysis")
data_matrix <- read.table("/home/data/t150327/HT22-RNAseq/DGE/Exmatrix/counts.matrix", sep = "\t", header = TRUE)

#Checking the matrix
dim(data_matrix)
head(data_matrix)
summary(data_matrix)

#Package Loading
if (length(find.package("DESeq2", quiet = TRUE)) == 0) {
  library(BiocManager)
  BiocManager::install("DESeq2")
  library(DESeq2)
} else {
  library(DESeq2)
  packageVersion("DESeq2")
}
library(pheatmap)
library(RColorBrewer)
library(glmpca)
library(ggplot2)
library(ggforce)
library(concaveman)
library(reshape2)

#Generating coldata
# Sample names
sample_names <- colnames(data_matrix)[-1]

# Group labels for each sample
#Regarding the design of experiment as two factors interaction
AD <- c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
SIRT6_upregulation <- c(0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1)

# Creating the coldata dataframe
coldata <- data.frame(Sample = sample_names,
                      Disease = factor(AD),
                      Upregulation = factor(SIRT6_upregulation))

# Make sure the row names match the sample names
rownames(coldata) <- sample_names
data_matrix <- data_matrix[,-1]

#Round the non-integer in count matrix
count_matrix_rounded <- round(data_matrix)

#Generate dds dataset
dds <- DESeqDataSetFromMatrix(countData = count_matrix_rounded,
                              colData = coldata,
                              design = ~ Disease+Upregulation+Upregulation:Disease)  

#Normalize the dds matrix
dds_nor <- DESeq(dds)
resultsNames(dds_nor)

#Counts perfilter
keep <- rowSums(counts(dds_nor)>=10) >=3
dds_nor <- dds_nor[keep,]
write.csv(assay(dds_nor),"dds_nor-1.csv")

#Reviewing sample
#Refer to https://mp.weixin.qq.com/s?src=11&timestamp=1684222487&ver=4531&signature=gusl-COjxMpKhoDVniOWpqxo0vskjuk*NugiXKVkiiy*x4SbeVDyR2jioOdN5-PnNc*Go4bWM56zSswfP666XUvZ1UFK1WZLwZPyeJh2PspYZmHKj6eJZmCFZb8*l6HV&new=1
#sample similarity heatmap
sampleDists = dist(t(assay(dds_nor)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(dds_nor$Sample)
colnames(sampleDistMatrix) <- paste(dds_nor$Sample)
color=colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
pdf("sample_pheatmap.pdf",width=5,height=4)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col=color)
dev.off()
#sample PCA
gpca = glmpca(counts(dds_nor),L=2)
gpca.dat=gpca$factors
gpca.dat$Group <- c("AD","AD","AD","ADTG","ADTG","ADTG","Ctrl","Ctrl","Ctrl","TG","TG","TG")
pdf("sample_PCA.pdf",width=5, height =4)
ggplot(gpca.dat,aes(x=dim1,y=dim2,color=Group))+geom_point(size=1.5)+
  geom_mark_ellipse(aes(fill=Group),alpha=0.1)+theme_bw()+
  coord_cartesian(clip="on")+labs(x="PC1",y="PC2")+
  theme(text=element_text(family="serif"))+
  xlim(-60,60)+
  ylim(-45,45)
dev.off()
#sample density plot
data <- log10(counts(dds_nor))
gene_num <- ncol(data)
data <- melt(data)
pdf("sample_density.pdf",width=5,height=4)
ggplot(data,aes(data$value,color=data$Var2))+
  xlab("log10(counts)")+
  geom_density(alpha=0.6)+
  geom_rug() + theme_bw()
dev.off()

