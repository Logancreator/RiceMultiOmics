library(ggplot2)
library(reshape2)
library(gg.gap)
library(readxl)
library(openxlsx)
library(ggbreak)
setwd("D:\\hzh\\RiceMultiOmics")
##########File1#############
#file1=read.csv("file1.csv",header=T,row.names = 1,check.names = F)
file1 = as.data.frame(read_excel("2.xlsx",range = "A1:F6"))

rownames(file1) = file1$`Sequence Length`;file1 = as.data.frame(file1[,c(2,3,4,5,6)])

file1$Algorithm <- rownames(file1) 

file1=melt(file1,id.vars = c("Algorithm"))

file1$Algorithm = factor(file1$Algorithm,levels = c("align_bench_seq","align_bench_par(30t)","align_bench_wave(30t)","bsalign(no band)","TSTA"))

p1=ggplot(data = file1, mapping = aes(x = variable, y = value,fill= Algorithm))+
  geom_bar(stat="identity",position=position_dodge(1))+
  scale_fill_manual(values = c("darkorchid3","darkolivegreen4","goldenrod1","dodgerblue4","red"))+
  theme_classic()+
  xlab("Sequence length (bp)")+
  ylab("Time (ms)")+scale_y_break(c(35000,80000),#截断位置及范围
                                  space = 1,#间距大小
                                  scales = 0.5);p1
ggsave(filename = "file1.pdf",plot = p1, width = 8, height =4)

##########File2#############
file2 = as.data.frame(read_excel("1.xlsx",range = "A2:E6"))

rownames(file2) = file2$`Number of threads/Block`;file2 = as.data.frame(file2[,c(2,3,4,5)])


#file2=read.csv("file2.csv",header=T,row.names = 1,check.names = F)

file2$Threads <- rownames(file2) 

file2=melt(file2,id.vars = c("Threads"))
colnames(file2)[2]="Block"
file2$Threads=factor(file2$Threads,levels=c(5,10,20,30))
p2=ggplot(data = file2, mapping = aes(x = factor(Threads,levels = c(5,10,20,30)), y = value, group=Block,linetype = Block, colour = Block, shape =Block, fill = Block))+
  geom_line(size=0.8)+ 
  scale_linetype_manual(values = c(1,1,1,1))+ 
  #scale_color_manual(values = c('steelblue','darkred',"gray"))+ 
  scale_shape_manual(values = c(20,20,20,20))+ 
  scale_color_manual(values = c("darkorchid3","darkolivegreen4","goldenrod1","dodgerblue4"))+
  theme_classic()+
  xlab("Number of threads")+
  ylab("Time (ms)")+ggtitle("Sequence Length: 200k (bp) - psa");p2
ggsave(filename = "file2.pdf",plot = p2, width = 6, height =3)

##########File3#############
#file3=read.csv("file3.csv",header=T,row.names = 1,check.names = F)

file3 = as.data.frame(read_excel("3.xlsx",range = "A2:E6"))

rownames(file3) = file3$`Number of threads/Block`;file3 = as.data.frame(file3[,c(2,3,4,5)])


file3$Block <- rownames(file3) 

file3=as.data.frame(melt(file3,id.vars = c("Block")))

colnames(file3)[2]="Threads"

file3$Block=factor(file3$Block,levels=c(5,10,20,30))

p3=ggplot(data=file3,mapping=aes(x=Threads,y=value,color=Block,group = Block))+geom_line()+
  scale_color_manual(values = c("darkorchid3","darkolivegreen4","goldenrod1","dodgerblue4"))+
  theme_classic()+
  xlab("Number of threads")+
  ylab("Time (ms)")+ggtitle("Sequence Length: 50k *10 (bp)");p3
ggsave(filename = "file3.pdf",plot = p3, width = 4, height =3)

##########File4#############
#file4=read.csv("file4.csv",header=T,row.names = 1,check.names = F)
file4 = as.data.frame(read_excel("4.xlsx",range = "A1:I5"))

rownames(file4) = file4$`Sequence Length`;file4 = as.data.frame(file4[,c(2:9)])

file4$Algorithm <- rownames(file4) 

file4=melt(file4,id.vars = c("Algorithm"))

file4$Algorithm = factor(file4$Algorithm,levels = c("SPOA","abPOA","bsalign(no band)","TSTA"))


p4=ggplot(data = file4, mapping = aes(x = variable, y = value, group=Algorithm,linetype =Algorithm,  shape =Algorithm, fill = Algorithm))+
  geom_bar(stat="identity",, position = 'dodge')+ 
  scale_linetype_manual(values = c(1,1,1,1,1))+ 
  scale_shape_manual(values = c(20,20,20,20,20))+ 
  scale_fill_manual(values = c("darkorchid3","darkolivegreen4","goldenrod1","dodgerblue4","orangered2"))+
  theme_classic()+
  xlab("Sequence Length * Quantity")+
  ylab("Time (ms)")+scale_y_break(c(50000,100000),#截断位置及范围
                                   space = 0.2,#间距大小
                                   scales = 0.4);p4
ggsave(filename = "file4.pdf",plot = p4, width = 7, height =4)

variables <- c("5K*5", "10K*5", "20K*5", "50K*5", "5K*10", "10K*10", "20K*10", "50K*10")
for (variable in variables) {
  subset <- file4[which(file4$variable == variable),]
  subset$Algorithm = factor(subset$Algorithm,levels = c("SPOA","abPOA","bsalign(no band)","TSTA"))
  
  p <- ggplot(data = subset, mapping = aes(x = variable, y = value, group = Algorithm, linetype = Algorithm, shape = Algorithm, fill = Algorithm)) +
    geom_bar(stat = "identity", position = 'dodge') +
    scale_linetype_manual(values = c(1, 1, 1, 1, 1)) +
    scale_shape_manual(values = c(20, 20, 20, 20, 20)) +
    scale_fill_manual(values = c("darkorchid3", "darkolivegreen4", "goldenrod1", "dodgerblue4", "orangered2")) +
    theme_classic() + xlab(NULL) + ylab(NULL) + theme(legend.position = 'none')
  filename <- paste0("file4_", gsub("\\*", "", variable), ".pdf")
  ggsave(filename = filename, plot = p, width = 2, height = 2)
}
