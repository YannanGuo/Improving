library(ggplot2)

data <- read.csv("DEG.csv",row.names = 1)

#################
# ggplot2绘制火山图
data$label <- c(rownames(data)[1:10],rep(NA,nrow(data) - 10))

ggplot(data,aes(log2FoldChange,-log10(pvalue),color = regulate)) + 
  xlab("log2FC") + 
  geom_point(size = 0.6) + 
  scale_color_manual(values=c("#00AFBB","#999999","#FC4E07")) + 
  geom_vline(xintercept = c(-1,1), linetype ="dashed") +
  geom_hline(yintercept = -log10(0.05), linetype ="dashed") + 
  theme(title = element_text(size = 15), text = element_text(size = 15)) + 
  theme_classic() + 
  geom_text(aes(label = label),size = 3, vjust = 1,hjust = -0.1)



#############
# 使用ggpurb绘制火山图：
library(ggpubr)
data$pvalue <- -log10(data$pvalue)
ggscatter(data, 
          x = "log2FoldChange", 
          y = "pvalue", 
          ylab="-log10(P.value)", 
          size=0.6, 
          color = "regulate", 
          label = rownames(data),
          label.select = rownames(data)[1:10],
          repel = T,
          palette = c("#00AFBB", "#999999", "#FC4E07")) +
          geom_hline(yintercept = 1.30,linetype ="dashed") +
          geom_vline(xintercept = c(-1,1),linetype ="dashed")
