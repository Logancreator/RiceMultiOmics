# ====RStudio中设置当前源文件目录作为工作目录================
# RStudio中设置当前源文件目录作为工作目录
# 路径中中文字符也可以识别。
rm(list=ls())
library(rstudioapi)
library(pheatmap)
library(xlsx)
library(rJava)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
pcr_1 = na.omit(read.xlsx("qpcr.xlsx",sheetIndex = 1,head=T));pcr_1
pcr_3 = na.omit(read.xlsx("qpcr.xlsx",sheetIndex = 2,head=T));pcr_3
rownames(pcr_1)=pcr_1[,1];pcr_1=pcr_1[,-1]
rownames(pcr_3)=pcr_3[,1];pcr_3=pcr_3[,-1]
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

p1=pheatmap(pcr_1,scale = 'column',
         cluster_rows = F,
         cluster_cols = F,
         main = "qPCR 1cm",
         border_color = "white",
         color = colorRampPalette(colors = c("blue","white","red"))(100))

save_pheatmap_pdf(p1, "qpcr_1cm.pdf",12,3.5)

p2=pheatmap(pcr_3,scale = 'column',
            cluster_rows = F,
            cluster_cols = F,
            main = "qPCR 3cm",
            border_color = "white",
            color = colorRampPalette(colors = c("blue","white","red"))(100))
save_pheatmap_pdf(p2, "qpcr_3cm.pdf",12,3.5)

