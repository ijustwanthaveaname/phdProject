rm(list = ls())
setwd("~/phdProject/Omics/metagenome")
library(vegan) 
######第一步，利用vegan包实现抽平分析的方法##################
otu <- read.csv('otu_table.csv',row.names = 1)

colSums(otu)#求和查看每个样品的和

otu_Flattening <- as.data.frame(t(rrarefy(t(otu),min(colSums(otu)))))#使用代码进行抽平

colSums(otu_Flattening)#查看抽平后每个样品的和

write.csv(otu_Flattening,"otu_Flattening.csv")#将抽平后的otu表格保存到目录下，准备后续的工作

############第二布，计算各种alpha指数#####################
#install.packages("picante")#没有这个包就先安装，有的话就省略
library("picante")

otu <- read.csv("otu_Flattening.csv", row.names = 1)#导入第一步中得出的抽平后的otu丰度表

otu <- t(otu)#将otu表格转置

alpha_index <- function(x, base = exp(1)) {    #为了方便快捷，这里制作一个函数，用来一次性计算所有alpha指数
  observed_species <- estimateR(x)[1,]#observed_species
  richness <- rowSums(x > 0)	#丰富度指数
  chao1 <- estimateR(x)[2, ]	#Chao1 指数
  ace <- estimateR(x)[4, ]	#ACE 指数
  shannon <- diversity(x, index = 'shannon', base = base)	#Shannon 指数
  simpson <- diversity(x, index = 'simpson')	#Gini-Simpson 指数
  pielou <- diversity(x, index = 'shannon', base = base) / log(estimateR(x)[1, ], base)	#Pielou 均匀度
  goods_covers <- 1 - rowSums(x == 1) / rowSums(x)	#goods_coverage
  result <- data.frame(observed_species, richness, chao1, ace, shannon, simpson, pielou, goods_covers)
}
alpha <- alpha_index(otu)
alpha$name <- rownames(alpha)#这一步和下一步给已经算出来的数据进行分组，以便第三步作图
alpha$group <- c(rep('X', 12), rep('Z', 12), rep('C', 12))
write.csv(alpha, "YBC_alpha.csv",row.names = T )#将所有alpha指数写入一个文档中

##########第三步，作箱线图图##########
library(ggplot2)

p <- ggplot(alpha, aes(x = group, y = chao1, fill = group))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.title=element_blank(),axis.title.x=element_blank(),
        strip.background=element_rect(fill="black",size=0.1),
        strip.text=element_text(colour="white"))
p
ggsave("YBC_chao1.PDF", p, width = 5.5, height = 5.5)
  
pp<- ggplot(alpha, aes(x = group, y = shannon, fill = group))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.title=element_blank(),axis.title.x=element_blank(),
        strip.background=element_rect(fill="black",size=0.1),
        strip.text=element_text(colour="white"))
pp
ggsave("YBC_shannon.PDF", pp, width = 5.5, height = 5.5)
