##############
data <- read.csv("DEG.csv",row.names = 1)

color <- rep("#999999",nrow(data))

color[data$pvalue <0.05 & data$log2FoldChange > 1] <- "#FC4E07"
color[data$pvalue <0.05 & data$log2FoldChange < -1] <- "#00AFBB"

par(oma = c(0,2,0,0))

plot(data$log2FoldChange,-log10(data$pvalue),pch = 16,cex = 0.5,
     xlim = c(-4,4), ylim = c(0,32), col = color, frame.plot = F,
     xlab = "log2FC", ylab = "-log10(Pvalue)", cex.axis = 1, cex.lab = 1.3)

# 添加参考线：
abline(h = -log10(0.05),lwd = 2, lty = 3)  # lwd设置线的宽度，lty设置线的类型；
abline(v = c(-1,1),lwd = 2, lty = 3)  # lwd设置线的宽度，lty设置线的类型；

# 添加图例
legend(x = 3, y = 32, legend = c("Up","Normal","Down"), 
       bty = "n", # 去除边框
       pch = 19,cex = 1, # 设置点的样式和大小
       x.intersp = 0.3, # 设置字与点之间的距离；
       y.intersp = 0.3, # 设置点与点的高度差，相当于行距；
       col = c("#FC4E07","#999999","#00AFBB"))

# 添加标签：
color = c()
color[which(data[1:10,]$regulate == "Up")] = "#FC4E07"
color[which(data[1:10,]$regulate != "Up")] = "#00AFBB"
text(data$log2FoldChange[1:10],-log10(data$pvalue)[1:10],
     labels = data$row[1:10],
     adj = c(0,1.5),
     cex = 0.6,
     col = color)


# 包装函数：
# 调整1: xlim和ylim得去掉
# 调整2: 修改图例的位置
plotVoc <- function(data){
  color <- rep("#999999",nrow(data))
  
  color[data$pvalue <0.05 & data$log2FoldChange > 1] <- "#FC4E07"
  color[data$pvalue <0.05 & data$log2FoldChange < -1] <- "#00AFBB"
  
  par(oma = c(0,2,0,0))
  
  plot(data$log2FoldChange,-log10(data$pvalue),pch = 16,cex = 0.5,
       col = color, frame.plot = F,
       xlab = "log2FC", ylab = "-log10(Pvalue)", cex.axis = 1, cex.lab = 1.3)
  
  # 添加参考线：
  abline(h = -log10(0.05),lwd = 2, lty = 3)  # lwd设置线的宽度，lty设置线的类型；
  abline(v = c(-1,1),lwd = 2, lty = 3)  # lwd设置线的宽度，lty设置线的类型；
  
  # 添加图例
  legend(x = 3, y = max(-log10(data$pvalue)), legend = c("Up","Normal","Down"), 
         bty = "n", # 去除边框
         pch = 19,cex = 1, # 设置点的样式和大小
         x.intersp = 0.3, # 设置字与点之间的距离；
         y.intersp = 0.3, # 设置点与点的高度差，相当于行距；
         col = c("#999999", "#FC4E07","#00AFBB"))
  
  # 添加标签：
  color = c()
  color[which(data[1:10,]$regulate == "Up")] = "#FC4E07"
  color[which(data[1:10,]$regulate != "Up")] = "#00AFBB"
  text(data$log2FoldChange[1:10],-log10(data$pvalue)[1:10],
       labels = data$row[1:10],
       adj = c(0,1.5),
       cex = 0.6,
       col = color)
}


data <- read.csv("DEG2.csv",row.names = 1)

plotVoc(data)

