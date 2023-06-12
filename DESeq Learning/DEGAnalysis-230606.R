###################################
######utf-8########################
#####Author Yannan Guo#############
###################################
##Created @ 20230606###############
###################################
##Learning the design of DESeq2####
###################################


rm(list = ls())
###Learning material from ####
##https://github.com/tavareshugo/tutorial_DESeq2_contrasts/blob/main/DESeq2_contrasts.md
##The github is refered to https://www.seqanswers.com/
##Two Factors with interaction DESeq2 design######

#Loading library
library(DESeq2)

##Stimulate DESeqDataSet of 1000genes(rows) and 12samples(columns) the standard deviation of the 
##log2 fold changes(effect size) across the genes is 4
dds_stimu <- makeExampleDESeqDataSet(n=1000, m=12, betaSD=4)
#create new variable"disease alternates between "disease" and "Ctrl" and assigned to "disease" and as factor
dds_stimu$disease <- factor(rep(c("AD","Ctrl"),each=6))
#This line changes the base/reference level of the factor "colour" to "white".
#This is useful when performing statistical analyses, because the first level of a factor is usually treated as the reference category against which other levels are compared.
dds_stimu$disease <- relevel(dds_stimu$disease, "Ctrl")
dds_stimu$condition <-  factor(rep(c("Sirt 6","WT"),6))
#This line rearranges the columns of the DESeqDataSet according to the levels of the factors "colour" and "condition". 
#order gives the indices that would sort the input; when used like this, 
#it rearranges the columns of the data set according to the sorted order of "colour" and "condition".
dds_stimu <- dds_stimu[,order(dds_stimu$disease, dds_stimu$condition)]
# renames the column names of the DESeqDataSet to be "sample1", "sample2", ..., "sampleN" where N is the number of columns in dds. 
#paste0 concatenates the strings without a separator, 
#and 1:ncol(dds) generates a sequence from 1 to the number of columns in dds.
colnames(dds_stimu) <- paste0("sample",1:ncol(dds_stimu))
#See the sample information
colData(dds_stimu)

#Set the DESeq analysis design
design(dds_stimu) <- ~ 1 + disease + condition + disease:condition
dds_stimu <-  DESeq(dds_stimu)
resultsNames(dds_stimu)

#get the model matrix
mod_mat <- model.matrix(design(dds_stimu),colData(dds_stimu))
#Define coefficient vectors for each condition
AD_Sirt6 <- colMeans(mod_mat[dds_stimu$disease == "AD" & dds_stimu$condition == "Sirt 6", ])
AD_WT <- colMeans(mod_mat[dds_stimu$disease == "AD" & dds_stimu$condition == "WT", ])
Ctrl_WT <- colMeans(mod_mat[dds_stimu$condition == "WT" & dds_stimu$disease == "Ctrl", ])
Ctrl_Sirt6 <- colMeans(mod_mat[dds_stimu$condition == "Sirt 6" & dds_stimu$disease == "Ctrl", ])


#Comparing the results of different groups
res1 <- results(dds_stimu, contrast = AD_WT - Ctrl_WT)
#res1 <- results(dds_stimu, contrast = list("disease_AD_vs_Ctrl")) # equivalently to the previous one

res2 <- results(dds_stimu, contrast = AD_Sirt6 - Ctrl_Sirt6)
#res2 <- results(dds_stimu, contrast = list(c("disease_AD_vs_Ctrl","diseaseAD.conditionWT")))

res3 <- results(dds_stimu, contrast= Ctrl_WT - Ctrl_Sirt6)
#res3 <- results(dds_stimu, contrast = list(c("condition_WT_vs_Sirt.6")))

res4 <- results(dds_stimu,contrast = AD_WT - AD_Sirt6)
#res4 <- results(dds_stimu, contrast = list(c("disease_AD_vs_Ctrl","diseaseAD.conditionWT")))

res5 <- results(dds_stimu, contrast = (AD_Sirt6 - AD_WT) - (Ctrl_Sirt6 - Ctrl_WT))
#res5 <- results(dds_stimu, contrast = list("diseaseAD.conditionWT"))


