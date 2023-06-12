library(clusterProfiler)
library(DESeq2)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(enrichplot)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(stats)
setwd('C:/Users/YuanJS/Desktop/gyn/')

#DESeq2
RNAseq.0809.counts <- read.delim("C:/Users/YuanJS/Desktop/gyn/counts.matrix")
sample_group <- read.csv("C:/Users/YuanJS/Desktop/gyn/sample_group.csv")
rownames(RNAseq.0809.counts)<-RNAseq.0809.counts[,1]
count <- RNAseq.0809.counts[,-1]
count <- round(count)
condition <- factor(sample_group$condition)
colData <- data.frame(row.names=colnames(count), condition)
dds <- DESeqDataSetFromMatrix(count, colData, design= ~ condition)
dds <- DESeq(dds)

#DE gene save 
DE1 <- results(dds, contrast=c("condition", "AD", "Ctrl"))
DE1 = DE1[order(DE1$log2FoldChange),]
DE_sub1 <-subset(DE1, pvalue < 0.05 & abs(log2FoldChange) > 1)
dim(DE_sub1)
write.csv(DE_sub1,file = "AD_vs_Ctrl.csv")

#volcano plot
res <- results(dds, contrast=c("condition", "AD", "Ctrl"))
res = res[order(res$pvalue),]#Set the pvalue ascending order
res$change = ifelse(res$pvalue < 0.01 & abs(res$log2FoldChange) >= 1, 
                     ifelse(res$log2FoldChange> 1 ,'Up','Down'),
                     'Stable')
res$symbol <- row.names(res)
res1 <- as.data.frame(res)
res2 <- res1[complete.cases(res1),]#remove the row with empty value
# method1
p1 <- ggplot(data = res2, 
            aes(x = log2FoldChange, 
                y = -log10(pvalue), 
                colour=change,
                label = symbol)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("blue", "grey","red"))+
  xlim(c(-5, 5)) +
  geom_vline(xintercept=c(-1,1),lty=2,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=2,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (p-value)",
       title="Differential AD_vs_Ctrl")  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())
for_label <- res2 %>% 
  filter(abs(log2FoldChange) >2 & -log10(pvalue)> -log10(0.001))
p1 + geom_point(size = 3.5, shape = 1, data = for_label) +
  ggrepel::geom_text_repel(
    aes(label = symbol),
    data = for_label,max.overlaps =30,force = 2,
    color="black")

# add certain gene on volcano plot
qPCR <- read.csv("C:/Users/YuanJS/Desktop/lwq220113/qPCR.csv", sep="")
csv <- as.list(qPCR[,1])
labellist <- res2[which(rownames(res2)%in%csv),]
p1 + geom_point(size = 3.5, shape = 1, data = labellist) +
  ggrepel::geom_text_repel(
    aes(label = symbol),
    data = labellist,max.overlaps =50,force = 8,
    color="black")

#over representation analysis,ORA   set threshhold pvalue&log2FoldChange
res <- results(dds, contrast=c("condition", "sg6", "N2"))
res = res[order(res$pvalue),]
diff_gene_deseq2 <-subset(res, pvalue < 0.05 & abs(log2FoldChange) > 1)
dim(diff_gene_deseq2)
write.csv(diff_gene_deseq2,file= "sg6_vs_N2.csv")

df <- read.csv("sg6_vs_N2.csv")
colnames(df) <- c('SYMBOL', 'baseMean', 'log2FoldChange', 'ifcSE', 'stat', 'pvalue', 'padj')
 #org.Mm.eg.db
dim(df.id)
easy.df<-merge(df,df.id,by= "SYMBOL",all=F)
head(easy.df)
sortdf<-easy.df[order(easy.df $log2FoldChange, decreasing = T),]
gene.expr = sortdf $log2FoldChange
names(gene.expr) <- sortdf $ENTREZID
ego <- enrichGO(
  gene = names(gene.expr),
  OrgDb = "org.Hs.eg.db",
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05)
#One of "BP", "MF", and "CC" subontologies, or "ALL" for all three
dotplot(ego,showCategory = 30,font.size = 7,title = "sg6 vs N2, GO_BP p<0.05 & abs(log2FC)>1") 

ekegg <- enrichKEGG(
  gene = names(gene.expr),
  organism     = 'hsa',
  pvalueCutoff = 0.05)

dotplot(ekegg,showCategory = 30,font.size = 7,title = "sg6_vs_N2,KEGG p<0.05 & abs(log2FC)>1")

#GSEA 
res <- results(dds, contrast=c("condition", "AD", "Ctrl"))
res = res[order(res$pvalue),]
diff_gene_deseq2 <-subset(res, pvalue < 0.05 )
dim(diff_gene_deseq2)
write.csv(diff_gene_deseq2,file= "AD_vs_Ctrl_0.05.csv")

df <- read.csv("AD_vs_Ctrl_0.05.csv")
colnames(df) <- c('SYMBOL', 'baseMean', 'log2FoldChange', 'ifcSE', 'stat', 'pvalue', 'padj')
df.id <- bitr(df$SYMBOL, fromType = "SYMBOL", toType =c("ENTREZID","ENSEMBL"), OrgDb = "org.Mm.eg.db")  #org.Mm.eg.db
head(df.id)
dim(df.id)
easy.df<-merge(df,df.id,by= "SYMBOL",all=F)
head(easy.df)
sortdf<-easy.df[order(easy.df $log2FoldChange, decreasing = T),]

gseadf<- sortdf[!duplicated(sortdf$ENTREZID),]%>% filter(log2FoldChange!="NA")
geneList = gseadf$log2FoldChange
names(geneList) <- gseadf$ENTREZID
length(geneList)

GSEAgo  <- gseGO(geneList  = geneList,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "BP", 
                 nPerm        = 1000,
                 minGSSize    = 10,
                 maxGSSize    = 500,
                 pvalueCutoff = 1)
#One of "BP", "MF", and "CC" subontologies, or "ALL" for all three
dotplot(GSEAgo,color="pvalue", showCategory = 30,font.size = 8,
        title = "sg1_vs_N2_GSEA_GO_BP p<0.05 & abs(log2FC)>1,pvalueCutoff = 1")
GSEAgo2 = setReadable(GSEAgo, OrgDb = org.Hs.eg.db) #transfer gene name
GSEAgo2 <-as.data.frame(GSEAgo2@result)
write.csv(GSEAgo2,file= "GSEAgo2_ALL_p_0.05_log2FC_0.5.csv")

GSEAkk<- gseKEGG(geneList   = geneList,
                 organism = 'hsa',
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 1,
                 verbose      = FALSE)
dotplot(GSEAkk,color="pvalue",showCategory = 30,font.size = 7,
        title = "sg1_vs_N2_GSEA_KEGG")

dotplot(GSEAkk,color="pvalue",showCategory = 30,font.size = 7,
        title = "sg1_vs_N2_GSEA_KEGG p<0.05 & abs(log2FC)>1,pvalueCutoff = 1")

GSEAkk2 = setReadable(GSEAkk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID") #transfer gene name
GSEAkk2 <-as.data.frame(GSEAkk2@result)
write.csv(GSEAkk2,file= "GSEA_KEGG_p_0.05log2FC_0.5.csv")

#GSEA down and up
GSEAgo2 = setReadable(GSEAgo, OrgDb = org.Hs.eg.db) #transfer gene name
GSEAgo2 <-as.data.frame(GSEAgo2@result)
GSEAgo3<-GSEAgo2[order(GSEAgo2 $NES, decreasing = T),]
GSEAgo3sub <- GSEAgo3[c(1:20),]
GSEAgo3sub1 <- tail(GSEAgo3,n=20)
GSEAgo3submer <- rbind(GSEAgo3sub,GSEAgo3sub1)
write.csv(GSEAgo2,file= "GSEAgo2_BP_p_0.05_log2FC_1.csv")

GSEAkk1 = setReadable(GSEAkk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID") #transfer gene name
GSEAkk1 <-as.data.frame(GSEAkk1@result)
GSEAkk2 <- subset(GSEAkk1,qvalues<0.25)
GSEAkk3<-GSEAkk2[order(GSEAkk2 $NES, decreasing = T),]

GSEAkk3sub <- GSEAkk3[c(1:20),]
GSEAkk3sub1 <- tail(GSEAkk3,n=20)
GSEAkk3submer <- rbind(GSEAkk3sub,GSEAkk3sub1)

library('forcats')
ggplot(GSEAgo3submer, aes(NES, fct_reorder(Description, NES), fill=pvalue, showCategory=30)) + 
  geom_bar(stat='identity') + 
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  theme_minimal() + ylab(NULL)

ggplot(GSEAkk3submer, aes(NES, fct_reorder(Description, NES), fill=pvalue, showCategory=30)) + 
  geom_bar(stat='identity') + 
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  theme_minimal() + ylab(NULL)

ggplot(GSEAkk3, aes(NES, fct_reorder(Description, NES), fill=NES, showCategory=30)) + 
  geom_bar(stat='identity') + 
  scale_fill_continuous(high='red', low='blue', guide=guide_colorbar(reverse=TRUE)) + 
  theme_minimal() + ylab(NULL)+
  theme(axis.text.x = element_text(face = "bold", size=14),
        axis.text.y = element_text(face = "bold", size=14))

ggplot(GSEAkk3, aes(NES, fct_reorder(Description, NES), fill=qvalues, showCategory=30)) + 
  geom_bar(stat='identity') + 
  scale_fill_continuous(high='red', low='blue', guide=guide_colorbar(reverse=TRUE)) + 
  theme_minimal() + ylab(NULL)+
  theme(axis.text.x = element_text(face = "bold", size=14),
        axis.text.y = element_text(face = "bold", size=14))








