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
##dplyr version 1.1.0################################
#####################################################

#Clean the environment
rm(list=ls())
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
library(dplyr)


#After reviewing the sample, S6_3 might need to exclude
#The code below just go through the process before
#Excluded S6_3
data_matrix_ad <- read.table("/home/data/t150327/HT22-RNAseq/DGE/Exmatrix/counts.matrix", sep = "\t", header = TRUE)
row.names(data_matrix_ad) <- data_matrix_ad[,1]
data_matrix_ad <- data_matrix_ad[,c(-1,-13)]
sample_names_ad <- colnames(data_matrix_ad)
AD <- c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0)
SIRT6_upregulation <- c(0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1)
coldata_ad <- data.frame(Sample = sample_names_ad,
                         Disease = factor(AD),
                         Upregulation = factor(SIRT6_upregulation))
count_matrix_ad <- round(data_matrix_ad)
dds_ad <- DESeqDataSetFromMatrix(countData = count_matrix_ad,
                                 colData = coldata_ad,
                                 design = ~ Disease * Upregulation) 
dds_ad <- DESeq(dds_ad)
res <- results(dds_ad)

#Get the DESeq result of disease vs non-disease
resultsNames(dds_ad)
res_disease <- results(dds_ad, name="Disease_1_vs_0")
summary(results(dds_ad,name = "Disease_1_vs_0"))
res_disease_LFC <- lfcShrink(dds_ad,coef="Disease_1_vs_0",type="apeglm")
res_disease_LFC <- res_disease_LFC[order(res_disease_LFC$pvalue),]

#Generate DE table
DE_disease_sub <- subset(res_disease_LFC, padj <0.05 & abs(log2FoldChange)>1)
write.csv(DE_disease_sub,file="DE_disease.csv")

#Volcano plot
res_disease_LFC$Change <- ifelse(res_disease_LFC$pvalue<0.05 & abs(res_disease_LFC$log2FoldChange)>=1,
                            ifelse(res_disease_LFC$log2FoldChange>1,"Up","Down"), "Stable")
res_disease_LFC$neg_log10_padj <- -log10(res_disease_LFC$padj)
DE_disease <- as.data.frame(res_disease_LFC)
DE_disease <- DE_disease[complete.cases(DE_disease),]
p <- ggplot(DE_disease, aes(x=log2FoldChange, y=neg_log10_padj, color=Change)) +  
    geom_point() + 
    theme_minimal() + 
    labs(x="Log2 Fold Change", y="-Log10 Adjusted P-value", title="Differential AD vs non-AD") + 
    scale_color_manual(values=c("Up"="red", "Down"="blue", "Stable"="grey")) +
    xlim(c(-5,5)) +
    geom_vline(xintercept=c(-1,1), lty=2, col="black", lwd=0.8) +
    geom_hline(yintercept = -log10(0.01), lty=2, col="black", lwd=0.8) +
    theme(axis.line=element_line(color="black"),plot.title = element_text(hjust = 0.5),legend.title = element_blank()) 
for_label <- DE_disease %>% filter(abs(log2FoldChange) >2 & -log10(pvalue)> -log10(0.001))
symbol <- row.names(for_label)
p + geom_point(size = 3.5, shape = 1, data = for_label) +
  ggrepel::geom_text_repel(
    aes(label = symbol),
    data = for_label,max.overlaps =30,force = 2,
    color="black")


#Get the DESeq result of Intercept
resultsNames(dds_nor)
res_inter <- results(dds_nor, name="Intercept")
summary(results(dds_nor,name = "Intercept"))
res_inter_LFC <- lfcShrink(dds_nor,type="apeglm")
res_inter_LFC <- res_inter_LFC[order(res_inter_LFC$pvalue),]

#Generate DE table
DE_inter_sub <- subset(res_inter, padj <0.05 & abs(log2FoldChange)>1)
write.csv(DE_inter_sub,file="DE_inter.csv")

#Volcano plot
res_inter$Change <- ifelse(res_inter$pvalue<0.05 & abs(res_inter$log2FoldChange)>=1,
                                 ifelse(res_inter$log2FoldChange>1,"Up","Down"), "Stable")
res_inter$neg_log10_padj <- -log10(res_inter$padj)
DE_inter <- as.data.frame(res_inter)
DE_inter <- DE_inter[complete.cases(DE_inter),]
p <- ggplot(DE_inter, aes(x=log2FoldChange, y=neg_log10_padj, color=Change)) +  
  geom_point() + 
  theme_minimal() + 
  labs(x="Log2 Fold Change", y="-Log10 Adjusted P-value", title="Differential Analysis") + 
  scale_color_manual(values=c("Up"="red", "Down"="blue", "Stable"="grey")) +
  xlim(c(-5,5)) +
  geom_vline(xintercept=c(-1,1), lty=2, col="black", lwd=0.8) +
  geom_hline(yintercept = -log10(0.01), lty=2, col="black", lwd=0.8) +
  theme(axis.line=element_line(color="black"),plot.title = element_text(hjust = 0.5),legend.title = element_blank()) 
for_label <- DE_inter %>% filter(abs(log2FoldChange) >2 & -log10(pvalue)> -log10(0.001))
symbol <- row.names(for_label)
p + geom_point(size = 3.5, shape = 1, data = for_label) +
  ggrepel::geom_text_repel(
    aes(label = symbol),
    data = for_label,max.overlaps =30,force = 2,
    color="black")


res_disease.upregulation <- summary(results(dds_nor,name = "Disease1.Upregulation1"))

