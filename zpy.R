library(ggplot2)
library(reshape2)
library(gg.gap)
setwd("C:\\Users\\jiany\\Desktop")
##########File1#############
file1=read.csv("file1.csv",header=T,row.names = 1,check.names = F)

file1$Algorithm <- rownames(file1) 

file1=melt(file1,id.vars = c("Algorithm"))

p1=ggplot(data = file1, mapping = aes(x = variable, y = value, group=Algorithm,linetype = Algorithm, shape =Algorithm, fill = Algorithm))+
  geom_bar(stat="identity",, position = 'dodge')+ 
  scale_linetype_manual(values = c(1,1,1))+ 
  #scale_color_manual(values = c('steelblue','darkred',"gray"))+ 
  scale_shape_manual(values = c(20,20,20))+ 
  scale_fill_manual(values = c("darkorchid3","darkolivegreen4","goldenrod1"))+
  theme_classic()+
  xlab("Sequence length (bp)")+
  ylab("Time (ms)");p1
ggsave(filename = "file1.pdf",plot = p1, width = 6, height =3)

##########File2#############
file2=read.csv("file2.csv",header=T,row.names = 1,check.names = F)

file2$Block <- rownames(file2) 

file2=as.data.frame(melt(file2,id.vars = c("Block")))
colnames(file2)[2]="Threads"
file2$Block=factor(file2$Block,levels=c(5,10,20,30))
p2=ggplot(data = file2, mapping = aes(x = factor(Threads,levels = c(5,10,20,30)), y = value, group=Block,linetype = Block, colour = Block, shape =Block, fill = Block))+
  geom_line(size=0.8)+ 
  scale_linetype_manual(values = c(1,1,1,1))+ 
  #scale_color_manual(values = c('steelblue','darkred',"gray"))+ 
  scale_shape_manual(values = c(20,20,20,20))+ 
  scale_color_manual(values = c("darkorchid3","darkolivegreen4","goldenrod1","dodgerblue4"))+
  theme_classic()+
  xlab("Number of threads")+
  ylab("Time (ms)")+ggtitle("Sequence Length: 200k (bp)");p2
ggsave(filename = "file2.pdf",plot = p2, width = 6, height =3)

##########File3#############
file3=read.csv("file3.csv",header=T,row.names = 1,check.names = F)

file3$Block <- rownames(file3) 

file3=as.data.frame(melt(file3,id.vars = c("Block")))
colnames(file3)[2]="Threads"
file3$Block=factor(file2$Block,levels=c(5,10,20,30))

p3=ggplot(data=file3,mapping=aes(x=Threads,y=value,color=Block,group=Block))+geom_line(size=.8)+ 
  scale_linetype_manual(values = c(1,1,1,1))+ 
  #scale_color_manual(values = c('steelblue','darkred',"gray"))+ 
  scale_shape_manual(values = c(20,20,20,20))+ 
  scale_color_manual(values = c("darkorchid3","darkolivegreen4","goldenrod1","dodgerblue4"))+
  theme_classic()+
  xlab("Number of threads")+
  ylab("Time (ms)")+ggtitle("Sequence Length: 50k *10 (bp)");p3
ggsave(filename = "file3.pdf",plot = p3, width = 4, height =3)

##########File4#############
file4=read.csv("file4.csv",header=T,row.names = 1,check.names = F)

file4$Algorithm <- rownames(file4) 

file4=melt(file4,id.vars = c("Algorithm"))

p4=ggplot(data = file4, mapping = aes(x = variable, y = value, group=Algorithm,linetype =Algorithm,  shape =Algorithm, fill = Algorithm))+
  geom_bar(stat="identity",, position = 'dodge')+ 
  scale_linetype_manual(values = c(1,1,1,1,1))+ 
  #scale_color_manual(values = c('steelblue','darkred',"gray"))+ 
  scale_shape_manual(values = c(20,20,20,20,20))+ 
  scale_fill_manual(values = c("darkorchid3","darkolivegreen4","goldenrod1","dodgerblue4","orangered2"))+
  theme_classic()+
  xlab("Sequence Length * Quantity")+
  ylab("Time (ms)")+scale_y_break(c(35000,80000),#截断位置及范围
                                   space = 0.8,#间距大小
                                   scales = 0.3);p4
# gg.gap(plot = p4,
#        ylim = c(0,80000),
#          segments = c(10000,40000))
# add.legend(plot = p4,
#            margin = c(top=1,right=1,bottom=1,left=460))
file4_1=file4[which(file4$variable=="5K*5"),]
file4_2=file4[which(file4$variable=="10K*5"),]
file4_3=file4[which(file4$variable=="5K*10"),]
file4_4=file4[which(file4$variable=="10K*10"),]
p5=ggplot(data = file4_1, mapping = aes(x = variable, y = value, group=Algorithm,linetype =Algorithm,  shape =Algorithm, fill = Algorithm))+
  geom_bar(stat="identity",, position = 'dodge')+ 
  scale_linetype_manual(values = c(1,1,1,1,1))+ 
  #scale_color_manual(values = c('steelblue','darkred',"gray"))+ 
  scale_shape_manual(values = c(20,20,20,20,20))+ 
  scale_fill_manual(values = c("darkorchid3","darkolivegreen4","goldenrod1","dodgerblue4","orangered2"))+
  theme_classic()+xlab(NULL)+ylab(NULL)+theme(legend.position = 'none');p5
ggsave(filename = "file4_1.pdf",plot = p5,width = 2,height = 2)
p6=ggplot(data = file4_2, mapping = aes(x = variable, y = value, group=Algorithm,linetype =Algorithm,  shape =Algorithm, fill = Algorithm))+
  geom_bar(stat="identity",, position = 'dodge')+ 
  scale_linetype_manual(values = c(1,1,1,1,1))+ 
  #scale_color_manual(values = c('steelblue','darkred',"gray"))+ 
  scale_shape_manual(values = c(20,20,20,20,20))+ 
  scale_fill_manual(values = c("darkorchid3","darkolivegreen4","goldenrod1","dodgerblue4","orangered2"))+
  theme_classic()+xlab(NULL)+ylab(NULL)+theme(legend.position = 'none');p6
ggsave(filename = "file4_2.pdf",plot = p6,width = 2,height = 2)
p7=ggplot(data = file4_3, mapping = aes(x = variable, y = value, group=Algorithm,linetype =Algorithm,  shape =Algorithm, fill = Algorithm))+
  geom_bar(stat="identity",, position = 'dodge')+ 
  scale_linetype_manual(values = c(1,1,1,1,1))+ 
  #scale_color_manual(values = c('steelblue','darkred',"gray"))+ 
  scale_shape_manual(values = c(20,20,20,20,20))+ 
  scale_fill_manual(values = c("darkorchid3","darkolivegreen4","goldenrod1","dodgerblue4","orangered2"))+
  theme_classic()+xlab(NULL)+ylab(NULL)+theme(legend.position = 'none');p7
ggsave(filename = "file4_3.pdf",plot = p7,width = 2,height = 2)
p8=ggplot(data = file4_4, mapping = aes(x = variable, y = value, group=Algorithm,linetype =Algorithm,  shape =Algorithm, fill = Algorithm))+
  geom_bar(stat="identity",, position = 'dodge')+ 
  scale_linetype_manual(values = c(1,1,1,1,1))+ 
  #scale_color_manual(values = c('steelblue','darkred',"gray"))+ 
  scale_shape_manual(values = c(20,20,20,20,20))+ 
  scale_fill_manual(values = c("darkorchid3","darkolivegreen4","goldenrod1","dodgerblue4","orangered2"))+
  theme_classic()+xlab(NULL)+ylab(NULL)+theme(legend.position = 'none');p8
ggsave(filename = "file4_4.pdf",plot = p8,width = 2,height = 2)
