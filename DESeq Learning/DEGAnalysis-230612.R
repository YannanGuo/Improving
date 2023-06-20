###################################
######utf-8########################
#####Author Yannan Guo#############
###################################
##Created @ 202306012##############
###################################
##Combination of DEGAnalysis and##
##Practical uses###################
###################################

#Preparation
rm(list=ls())
library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(tidyverse)

#Data input and pre-cleaning
data_matrix <- read.table("/home/data/t150327/HT22-RNAseq/DGE/Exmatrix/counts.matrix",sep="\t",header = TRUE)
summary(data_matrix)
dim(data_matrix)
rownames(data_matrix) <- data_matrix[,1]
data_matrix <- round(data_matrix[,c(-1,-13)])

#Edit DESeq input data
sample_name <- colnames(data_matrix)
Disease <- c(1,1,1,1,1,1,0,0,0,0,0)
Gene <- c(0,0,0,1,1,1,0,0,0,1,1)
coldata <- data.frame(Sample = sample_name,
                      Disease = factor(Disease),
                      Gene = factor(Gene))

#Generating DESeq DataSet
dds <- DESeqDataSetFromMatrix(countData = data_matrix,
                              colData = coldata,
                              design = ~ 1 + Disease + Gene + Disease:Gene)
dds <- DESeq(dds)

#get the model matrix
mod_mat <- model.matrix(design(dds),colData(dds))
#Define coefficient vectors for each condition
AD_Sirt6 <- colMeans(mod_mat[dds$Disease == 1 & dds$Gene == 1, ])
AD_WT <- colMeans(mod_mat[dds$Disease == 1 & dds$Gene == 0, ])
Ctrl_WT <- colMeans(mod_mat[dds$Disease == 0 & dds$Gene == 0, ])
Ctrl_Sirt6 <- colMeans(mod_mat[dds$Disease == 0 & dds$Gene == 1, ])

#Get the results of different groups
res.AD_Sirt6vsAD <- results(dds, contrast = AD_Sirt6 - AD_WT)
res.ADvsWT <- results(dds, contrast = AD_WT - Ctrl_WT)
res.Sirt6vsWT <- results(dds, contrast = Ctrl_Sirt6 - Ctrl_WT)
res.AD_Sirt6vsWT <- results(dds, contrast = AD_Sirt6 - Ctrl_WT)
res.cross <- results(dds, contrast = (AD_Sirt6 - AD_WT) - (Ctrl_Sirt6 - Ctrl_WT))

#Generate DE table
DE.AD_Sirt6vsAD <- subset(res.AD_Sirt6vsAD, padj<0.05 & abs(log2FoldChange)>0.3)
write.csv(DE.AD_Sirt6vsAD, file="./DE_table/DE.AD_Sirt6vsAD.csv")
DE.ADvsWT <- subset(res.ADvsWT, padj<0.05 & abs(log2FoldChange)>0.3)
write.csv(DE.ADvsWT, file="./DE_table/DE.ADvsWT.csv")
DE.Sirt6vsWT <- subset(res.Sirt6vsWT, padj<0.05 & abs(log2FoldChange)>0.3)
write.csv(DE.Sirt6vsWT, file="./DE_table/DE.Sirt6vsWT.csv")
DE.AD_Sirt6vsWT <- subset(res.AD_Sirt6vsWT, padj<0.05 & abs(log2FoldChange)>0.3)
write.csv(DE.AD_Sirt6vsWT, file="./DE_table/DE.AD_Sirt6vsWT.csv")
DE.cross <- subset(res.cross, padj<0.05 & abs(log2FoldChange)>0.3)
write.csv(DE.cross, file="./DE_table/DE.cross.csv")

#Volcano plot
res.AD_Sirt6vsAD$Change <- ifelse(res.AD_Sirt6vsAD$pvalue < 0.05 & abs(res.AD_Sirt6vsAD$log2FoldChange) >=1,
                                  ifelse(res.AD_Sirt6vsAD$log2FoldChange >1, "Up", "Down"), "Stable")
res.AD_Sirt6vsAD$symbol <- row.names(res)

