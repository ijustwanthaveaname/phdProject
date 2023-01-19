rvcheck::update_all(check_R=F)

setwd("/home/tuyx/Apc_metab/")
org_mix <- read.delim("score50/metab_abund.txt",row.names = 1)
dim(org_mix)
# 保留样本分组表
Label<-sapply(strsplit(colnames(org_mix),split = "G|0"), function(x){
  return(x[[1]])
})
meta_data<-data.frame(Label=Label,row.names = colnames(org_mix))


# 缺失值处理
#####将缺失值替换为数值0，方便后续我们进行数值判断
## 如果没有缺失值，就把0换成组内最小值或者是组内最小值的一半
# 计算组间最小值
Replace_Min<-function(mt_group,mt_df){
  mt_group=unique(mt_group)  ## 不同组
  mt_list=list()
  for (i in mt_group) {
    ## 不同组的矩阵
    group_df=mt_df[,str_detect(colnames(mt_df),i)]  
    kk = apply(group_df,1,function(x){
      min_replace = min(x,na.rm = T) ## 行:代谢物的最小值
      # half_replace=min_replace/2  ## 也可以用最小值的一半替代
      x[is.na(x)] = min_replace   ## 该行中NA值替换成最小值,如果都是NA，则最小值是Inf
      return(x)
    })%>%t()
    
    mt_list[[i]]<-kk
  }
  return(mt_list)
}
a <- org_mix[,!str_detect(colnames(org_mix),pattern = "^QC")]
a[a==0]<-NA  ## 假设这里的0还没有处理好，还都是NA值

c<-Replace_Min(unique(meta_data$Label)[1:4],a)
c<-do.call(cbind,c)
c[c==Inf]=0





# score50
setwd("/home/tuyx/Apc_metab/")
library(reshape2)
library(dplyr)
library(mixOmics)
library(limma)

org_mix <- read.delim("score50/metab_abund.txt",row.names = 1)
dim(org_mix)
colnames(org_mix)
#QC <- org_mix[,grep("QC",colnames(org_mix))]
org_mix <- org_mix[,colnames(org_mix)[order(colnames(org_mix))]%>%.[c(1:17,23:43,18:22)]]
group <- paste0(strsplit2(colnames(org_mix),"_")[,1],"_",strsplit2(colnames(org_mix),"_")[,2])
group[39:43] <- "QC" 
Group <- group[1:38]
dim(org_mix)

# 缺失值处理
#####将缺失值替换为数值0，方便后续我们进行数值判断
str(org_mix)
a <- org_mix[,1:38]
str(a)
a[a==0] <- NA

dim(a)
#View(a)
result <- data.frame(rowmname=rownames(a),missing=rowSums(is.na(a)))
head(result)
result$missing_percent<-result$missing/length(colnames(a)) ## 34个样本
res_filter<-filter(result,missing_percent>0.5 | missing_percent==0.5)  ## 0.5筛选

dim(res_filter)
res_keep<-filter(result,missing_percent<0.5)  ## 0.5筛选
dim(res_keep)
c=a[rownames(a)%in%res_keep$rowmname,]
dim(c)
# min 填充
min <- as.data.frame(t(apply(c,1,function(x){
  tapply(x,Group,min,na.rm=T)})))
min
dim(min)
#t <- c(NA,NA,NA)
#min(t,na.rm=T)
colnames(c)
colnames(min)
table(Group)
for (i in 1:nrow(a)){
  c[i,1:9][is.na(c[i,1:9])]=min[i,1]}
for (i in 1:nrow(a)){
  c[i,10:17][is.na(c[i,10:17])]=min[i,2]}
for (i in 1:nrow(a)){
  c[i,18:28][is.na(c[i,18:28])]=min[i,3]}
for (i in 1:nrow(a)){
  c[i,29:38][is.na(c[i,29:38])]=min[i,4]}
c[c==Inf]=0

# RSD
QC<- org_mix[,39:43]
head(QC)
RSD <- data.frame(rsd=apply(QC,1,sd)/apply(QC,1,mean))
head(RSD,3)
dim(RSD)
RSD0.2 <- subset(RSD,RSD$rsd< 0.2)
dim(RSD0.2)
c <- c[rownames(c)%in%rownames(RSD0.2),]
dim(c)
# log 转换
logc <- log10(c+1)
head(logc)
#View(logc)
name <- read.delim("score50/metab_desc.txt",encoding = "utf-8")

rownames(names) <- name$X
name=name[grep("metab_",name$Metabolite,invert = T),]
View(name)

#View(name)
dim(name)
colnames(name)
name <- name[,c(1:3,12:13)]
#write.table(name,"score50/name.txt", sep = '\t', col.names = NA, quote = FALSE)
library(ropls)
dim(logc)
intersect(name$metab_id,rownames(logc))
logc <- logc[rownames(logc)%in%name$metab_id,]
logc
dim(logc)
dim(name)
colnames(logc)

logc_scale <- t(scale(t(c[,1:38])))
dim(logc_scale)
# VIP value
pos <- logc_scale[rownames(logc_scale)%in%name[grep("pos",name$ID),]$metab_id,]
neg <- logc_scale[rownames(logc_scale)%in%name[grep("neg",name$ID),]$metab_id,]
dim(pos)
#dim(neg)



#########################  Apc_M vs Apc_F ###################################
dim(logc)
d <- logc[,c(1:17)]
colnames(d)
e <- c[,c(1:17)][rownames(c[,c(1:17)])%in%rownames(logc),]
dim(e)
rownames(logc)
Pvalue<-c(rep(0,nrow(d))) 
FC<-c(rep(0,nrow(d))) 
colnames(d)
for(i in 1:nrow(d)){
  if(sd(d[i,10:17])==0&&sd(d[i,1:9])==0){
    Pvalue[i] <-"NA"
    FC[i]<-"NA"
  }else{
    y=t.test(as.numeric(d[i,10:17]),as.numeric(d[i,1:9]))
    Pvalue[i]<-y$p.value
    FC[i]<-mean(as.numeric(e[i,10:17]))/mean(as.numeric(e[i,1:9])) 
  }
}
dim(d)
FDR=p.adjust(Pvalue, "BH")
mix <- cbind(d,FC)%>%cbind(.,FDR)%>%as.data.frame(.)
dim(mix)
l <- list(pos[,c(1:17)],neg[,c(1:17)])
dim(pos)
dim(neg)
dim(l[[2]])
dat <- list()
group <- list()
data <- list()
out <- list()
oplsda <- purrr::map(1:2,function(x){
  dat[[x]] <- as.data.frame(t(l[[x]]),stringsAsFactors = F)
  group[[x]] <- substr(colnames(d),5,5) %>% factor()
  data[[x]]= as.data.frame(apply(dat[[x]],2,as.numeric))
  rownames(data[[x]]) <- rownames(dat[[x]])
  out[[x]] <- opls(data[[x]], group[[x]], orthoI = 1) 
})
vipVn <- list()
vipVn_select <- list()
vip<- do.call(rbind,purrr::map(1:2,function(x){
  vipVn[[x]] <- getVipVn(oplsda[[x]])%>%as.data.frame(.)
}))
dim(vip)
#View(vip)
dim(mix)
mix$log2FC <- log2(as.numeric(mix$FC))
dim(vip)

mix <- mix[rownames(mix)%in%rownames(vip),]
mix <- cbind(mix,vip)
dim(mix)

# fdr vipmix
mix_select1=subset(mix,FDR<0.05 & .>1 )
# fc
#mix_select2=subset(mix_select1,FC< 0.7| FC >1.5 )
subset(mix_select1,FC > 1)%>%dim()
subset(mix_select1,FC< 1)%>%dim()

dim(mix_select1)
mix_select1$metab_id <- rownames(mix_select1)
mix_select <- merge(mix_select1,name,all.x=T,by="metab_id")
mix_select <- mix_select[order(mix_select$FC,decreasing = T),]
#View(mix_select)
mix_select$group <- c(rep("Apc_M",147),rep("Apc_F",139))
Apc_MvsApc_F_select <- mix_select

#View(mix_select)
write.table(mix_select,"score50/Apc_MvsApc_F_select20220725.txt", sep = '\t', col.names = NA, quote = FALSE)
res_MF<- as.data.frame(mix)
dim(res_MF)
res_MF$type <- ifelse(res_MF$FDR < 0.05,
                      ifelse(abs(res_MF$log2FC) > 0 ,
                             ifelse(res_MF$log2FC < 0 ,'down','up'),'noSig'),'noSig')

class(res_MF)
res_MF[which(res_MF$.<1),]$type <- "noSig"
res_MF$VIP <- res_MF$. 
res_MF$. <- NULL
# check
head(res_MF,5)
res_MF$metab_id <- rownames(res_MF)
res_MF<- merge(res_MF,name,all.x=T,by="metab_id")

# filter gene
#metab_14325
#metab_17992
Metabolites_Origin <- read.csv("score50/Apc_MetOrigin_20220728200453/02_Origin_Analysis/Metabolites_Origin.csv")
meta <- res_MF[res_MF$Metabolite%in%c(Metabolites_Origin[Metabolites_Origin$Origin%in%c("Microbiota","Co-Metabolism"),]$Name ,
                                      "N-[(3a,5b,7b)-7-hydroxy-24-oxo-3-(sulfooxy)cholan-24-yl]-Glycine",
                                      "Aldosterone 18-glucuronide","Dehydroisoandrosterone 3-glucuronide",
                                      "Palmitoyl-L-carnitine","Pantothenic Acid",
                                      "Alpha-Linolenic acid","4-Guanidinobutanoic acid",
                                      "N-[(3a,5b,7b)-7-hydroxy-24-oxo-3-(sulfooxy)cholan-24-yl]-Glycine",
                                      "N-Methylglutamic acid"),]
meta
meta=meta[meta$Metabolite%in%c("Niacinamide","LysoPC(16:1(9Z))","LysoPC(16:0)","PC(16:0/0:0)","Glycerophosphocholine",
                               "4-Guanidinobutanoic acid","Alpha-Linolenic acid"),]

p=ggplot(res_MF,aes(x = log2FC,y = -log10(FDR))) +
  geom_point(aes(color = type,size=VIP),size = 3) +
  scale_color_manual(name = '',
                     # color or three types
                     values = c('up'='#DA1212','noSig'='grey','down'='#3E7C17'),
                     # legend labels
                     label = c('up'='up (num=147)','noSig'='noSig (num=579)','down'='down (num=139)')) +
  # 阈值分割线
  geom_hline(yintercept = -log10(0.05),lty = 'dashed',size = 0.8) +
  geom_vline(xintercept = c(-1,1),lty = 'dashed',size = 0.8) +
  ggtitle('Apc-M vs Apc-F')+theme_classic()
p + geom_text_repel(data = meta,aes(x = log2FC,y = -log10(FDR),label = meta$Metabolite),
                    force=20,color="black",size=6,point.padding = 0.5,hjust = 0.5,
                    arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
                    segment.color="black",segment.size=0.2,nudge_y=1)+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "top",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+
 theme(plot.title = element_text(size = 24,face = "bold"))



############################# WT_M vs WT_F###############################################
colnames(logc)
d <- logc[,c(18:38)]

e <- c[,c(18:38)][rownames(c[,c(18:38)])%in%rownames(logc),]
Pvalue<-c(rep(0,nrow(d))) 
FC<-c(rep(0,nrow(d))) 
colnames(d)
for(i in 1:nrow(d)){
  if(sd(d[i,12:21])==0&&sd(d[i,1:11])==0){
    Pvalue[i] <-"NA"
    FC[i]<-"NA"
  }else{
    y=t.test(as.numeric(d[i,12:21]),as.numeric(d[i,1:11]))
    Pvalue[i]<-y$p.value
    FC[i]<-mean(as.numeric(e[i,12:21]))/mean(as.numeric(e[i,1:11])) 
  }
}
FDR=p.adjust(Pvalue, "BH")
mix <- cbind(d,FC)%>%cbind(.,FDR)
l <- list(pos[,c(18:38)],neg[,c(18:38)])
dat <- list()
group <- list()
data <- list()
out <- list()
oplsda <- purrr::map(1:2,function(x){
  dat[[x]] <- as.data.frame(t(l[[x]]),stringsAsFactors = F)
  group[[x]] <- substr(colnames(d),4,4) %>% factor()
  data[[x]]= as.data.frame(apply(dat[[x]],2,as.numeric))
  rownames(data[[x]]) <- rownames(dat[[x]])
  out[[x]] <- opls(data[[x]], group[[x]], orthoI = 1) 
})
vipVn <- list()
vipVn_select <- list()
vip<- do.call(rbind,purrr::map(1:2,function(x){
  vipVn[[x]] <- getVipVn(oplsda[[x]])%>%as.data.frame(.)
}))
head(vip)
mix$log2FC <- log2(as.numeric(mix$FC))
dim(vip)

mix <- mix[rownames(mix)%in%rownames(vip),]
mix <- cbind(mix,vip)
dim(mix)

# fdr vipmix
mix_select1=subset(mix,FDR<0.05 & .>1 )
# fc
#mix_select2=subset(mix_select1,FC< 0.7| FC >1.5 )
subset(mix_select1,FC > 1)%>%dim()
subset(mix_select1,FC< 1)%>%dim()

dim(mix_select1)
mix_select1$metab_id <- rownames(mix_select1)
mix_select <- merge(mix_select1,name,all.x=T,by="metab_id")
mix_select <- mix_select[order(mix_select$FC,decreasing = T),]
#View(mix_select)
mix_select$group <- c(rep("WT_M",137),rep("WT_F",130))

#View(mix_select)
#write.table(mix_select,"score50/WT_MvsWT_F_select20220725.txt", sep = '\t', col.names = NA, quote = FALSE)
res_MF<- as.data.frame(mix)
dim(res_MF)

res_MF$type <- ifelse(res_MF$FDR < 0.05,
                      ifelse(abs(res_MF$log2FC) > 0 ,
                             ifelse(res_MF$log2FC < 0 ,'down','up'),'noSig'),'noSig')

class(res_MF)
res_MF[which(res_MF$.<1),]$type <- "noSig"
res_MF$VIP <- res_MF$. 
res_MF$. <- NULL
# check
head(res_MF,5)
res_MF$metab_id <- rownames(res_MF)
res_MF<- merge(res_MF,name,all.x=T,by="metab_id")

# filter gene
#metab_14325
#metab_17992
Metabolites_Origin <- read.csv("score50/WT_MetOrigin_20220728200453/02_Origin_Analysis/Metabolites_Origin.csv")
meta <- res_MF[res_MF$Metabolite%in%Metabolites_Origin[Metabolites_Origin$Origin%in%c("Microbiota","Co-Metabolism"),]$Name ,]
meta <- meta[meta$Metabolite%in%c("Indole-3-carboxaldehyde","Benzaldehyde",
                                  "4-Hydroxybenzaldehyde","Isomaltose",
                                  "S-Adenosylmethionine","Deoxyinosine"),]
p=ggplot(res_MF,aes(x = log2FC,y = -log10(FDR))) +
  geom_point(aes(color = type,size=VIP),size = 3) +
  theme_classic(base_size = 16)  +
  scale_color_manual(name = '',
                     # color or three types
                     values = c('up'='#DA1212','noSig'='grey','down'='#3E7C17'),
                     # legend labels
                     label = c('up'='up (num=137)','noSig'='noSig (num=596)','down'='down (num=130)')) +
  # 阈值分割线
  geom_hline(yintercept = -log10(0.05),lty = 'dashed',size = 0.8) +
  geom_vline(xintercept = c(-1,1),lty = 'dashed',size = 0.8) +
  ggtitle('WT_M vs WT_F')
p +geom_text_repel(data = meta,aes(x = log2FC,y = -log10(FDR),label = meta$Metabolite),
                   force=20,color="black",size=6,point.padding = 0.5,hjust = 0.5,
                   arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
                   segment.color="black",segment.size=0.2,nudge_y=1)+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "top",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+
  theme(plot.title = element_text(size = 24,face = "bold"))
# fdr vipmix
mix_select1=subset(mix,FDR<0.05 & .>1 )
# fc
#mix_select2=subset(mix_select1,FC< 0.7| FC >1.5 )
subset(mix_select1,FC< 1)%>%dim()
subset(mix_select1,FC> 1)%>%dim()
dim(mix_select)
mix_select1$metab_id <- rownames(mix_select1)
mix_select <- merge(mix_select1,name,all.x=T,by="metab_id")
mix_select <- mix_select[order(mix_select$FC,decreasing = T),]
View(mix_select)
mix_select$group <- c(rep("WT_M",137),rep("WT_F",130))
WT_MvsWT_F_select <- mix_select
write.table(mix_select,"score50/WT_MvsWT_F_select20220727.txt", sep = '\t', col.names = NA, quote = FALSE)
Apc_MvsApc_F_select <- read.delim("score50/Apc_MvsApc_F_select20220725.txt",row.names = 1)
WT_MvsWT_F_select <- read.delim("score50/WT_MvsWT_F_select20220727.txt",row.names = 1)
diff <- setdiff(Apc_MvsApc_F_select$Metabolite,WT_MvsWT_F_select$Metabolite)
diff
diff_metab <- Apc_MvsApc_F_select[Apc_MvsApc_F_select$Metabolite%in%diff,]
View(diff_metab)
write.table(diff_metab,"score50/diff_metab.txt20220727", sep = '\t', col.names = NA, quote = FALSE)

#diff_metab <- read.delim("score50/diff_metab.txt",row.names = 1)
#diff_Apc_name <- data.frame(ID=strsplit2(diff_metab$Library.ID,";")[grep("HMDB",strsplit2(diff_metab$Library.ID,";"))])
#write.table(diff_Apc_name,"score50/diff_Apc_name.txt", sep = '\t', col.names = NA, quote = FALSE)


#PCA（主成分分析）
#BiocManager::install("PCAtools")
library(PCAtools)
meta_data <- data.frame(group=Group)
rownames(meta_data) <- colnames(logc_scale)
pca <- pca(logc, metadata =meta_data )
biplot(pca, x = 'PC1', y = 'PC2') #可以看到这两个成分对样品的解释度
#screen(p) #可以看到所有成分对样品的解释
##将pca与关联样本信息表 
pca$variance
pca_rotated_plus <- cbind(pca$rotated,meta_data)
library(ggsci)
pca_rotated_plus$group <- factor(pca_rotated_plus$group,levels=c("WT_M","WT_F","Apc_M","Apc_F"))
ggplot(pca_rotated_plus,aes(x = PC1 , y = PC2))+
  geom_point(size = 8,aes(shape = group, fill = group)) +
  stat_ellipse(aes(color = group,
                   fill=pca_rotated_plus$group),
               linetype = 'dashed',
               size = 1, show.legend = FALSE) + #添加分组椭圆
  labs(x = 'PC1 (22.12%)',y = 'PC2 (14.73%)') + 
  scale_shape_manual(values = c(21,22,23,24))+ 
  scale_fill_manual(values = c("#4575B4","#FCB271","#689F93","#CEB855"))+
  scale_color_manual(values = c("#4575B4","#FCB271","#689F93","#CEB855"))+
  theme_classic()+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),	
         legend.text = element_text(size=12,face = "bold"),legend.position = "top",legend.title =NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+
  theme(plot.title = element_text(size = 24,face = "bold"))


install.packages("scatterplot3d")
library(scatterplot3d)
my_color = c("#FCB271","#CEB855","#4575B4","#689F93")
library(FactoMineR)
tdata <- as.data.frame(t(as.data.frame(logc_scale)))
tdata$type <- factor(Group,levels = c("Apc_M", "Apc_F","WT_M","WT_F"))
dim(tdata)
dat.pca = PCA(tdata[,-906],graph = F)
dat.pca$eig
a = dat.pca[["ind"]][["coord"]] 
my_pch = 21:24
pchs = my_pch[as.numeric(tdata$type)]
colors = my_color[as.numeric(tdata$type)]
?scatterplot3d
p = scatterplot3d(a[,1:3],color = "black",pch = pchs,bg = colors,cex.lab=2,cex.axis=1.5,cex.symbols = 2,
                  xlab = "Dim.1(21.57%)",ylab = "Dim.2(13.99%)",zlab = "Dim.3(12.44%)")
legend("bottom",col = "black", legend = levels(tdata$type),pt.bg =  my_color, pch = my_pch,
       inset = -0.2, xpd = TRUE, horiz = TRUE)

tdata <- as.data.frame(t(as.data.frame(logc_scale)))
tdata$type <- factor(Group,levels = c("Apc_M", "Apc_F","WT_M","WT_F"))
dim(tdata)
dat.pca = PCA(tdata[,-906],graph = F)
dat.pca$eig
a = dat.pca[["ind"]][["coord"]] 
my_pch = 21:24
pchs = my_pch[as.numeric(tdata$type)]
colors = my_color[as.numeric(tdata$type)]
?scatterplot3d


p = scatterplot3d(a[,1:3],color = "black",pch = pchs,bg = colors,cex.lab=2,cex.axis=1.5,cex.symbols = 2,
                  xlab = "Dim.1(21.57%)",ylab = "Dim.2(13.99%)",zlab = "Dim.3(12.44%)"#,main = "PCA of metablites"
                  )
legend("bottom",col = "black", legend = levels(tdata$type),pt.bg =  my_color, pch = my_pch,
       inset = -0.2, xpd = T, horiz = T)

tdata <- as.data.frame(t(as.data.frame(logc)))
tdata$type <- factor(Group,levels = c("Apc_M", "Apc_F","WT_M","WT_F"))
dim(tdata)
dat.pca = PCA(tdata[,-906],graph = F)
dat.pca$eig
a = dat.pca[["ind"]][["coord"]] 
my_pch = 21:24
pchs = my_pch[as.numeric(tdata$type)]
colors = my_color[as.numeric(tdata$type)]
?scatterplot3d
p = scatterplot3d(a[,1:3],color = "black",pch = pchs,bg = colors,cex.lab=2,cex.axis=1.5,cex.symbols = 2,
                  xlab = "Dim.1(21.57%)",ylab = "Dim.2(13.99%)",zlab = "Dim.3(12.44%)")
legend("bottom",col = "black", legend = levels(tdata$type),pt.bg =  my_color, pch = my_pch,
       inset = -0.2, xpd = TRUE, horiz = TRUE)


#metab_origin
#metab_origin
Metabolites_Origin <- read.csv("score50/Diff_MetOrigin_20220728200453//02_Origin_Analysis/Metabolites_Origin.csv")
#Metabolites_Origin <- read.csv("score50/MetOrigin_20220522192255/02_Origin_Analysis/Metabolites_Origin.csv")

View(Metabolites_Origin)
#Metabolites_Origin$diff <- gsub(1,"SigDiff",Metabolites_Origin$Diff)%>%gsub("0","UnSigDiff",.)
dim(Metabolites_Origin)
Metabolites_Origin$value <- 1
Metabolites_Origin$Origin <- gsub("Drug related", "Others",Metabolites_Origin$Origin)%>%gsub("Unknown", "Others",.)
Metabolites_Origin$Origin <- factor(Metabolites_Origin$Origin,levels = c("Microbiota","Co-Metabolism","Host","Food related","Others"  ))
unique(Metabolites_Origin$Origin)
Metabolites_Origin$sample <- "ApcMF"
unique(Metabolites_Origin$Origin)

ggplot(Metabolites_Origin) +
  geom_bar(aes(x =sample,y=value,fill = Origin),stat = "identity",position = "fill") +
  scale_fill_manual(values  = c("#4D7587","#85AEC9","#976276","#E1AB8D"))+
  theme_classic()+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=14,face = "bold"),legend.position = "right",legend.title = element_text(size=14,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+
  theme(plot.title = element_text(size = 24,face = "bold"))+xlab("")+ylab("% of Origin")+
  scale_fill_manual(name = '',
                     values = c('Microbiota'='#4D7587','Co-Metabolism'='#85AEC9','Food related'='#976276','Others'='#E1AB8D'),
                     label = c('Microbiota'='Microbiota (num=6)','Co-Metabolism'='Co-Metabolism (num=17)',
                               'Food related'='Food related (num=54)','Others'='Others (num=94)'))


ggdata = list(c(Metabolites_Origin$Name[Metabolites_Origin$Origin%in%c("Microbiota","Co-Metabolism")]),c(Metabolites_Origin$Name[Metabolites_Origin$Origin%in%c("Host","Co-Metabolism")]))
names(ggdata) <- c("Microbiota","Host")
ggdata
BiocManager::install("ggvenn")
library(ggvenn)
p=ggvenn(ggdata,       
         show_percentage = F,show_elements=F,text_size = 8,
         stroke_color = "black",
         fill_color = c("#85AEC9","#A0B3A1"),
         set_name_color = c("#85AEC9","#A0B3A1")) 
p

MPEA_Results_Co.Metabolism <- read.csv("score50/Diff_MetOrigin_20220728200453/03_Function_Analysis/MPEA_Results_Co-Metabolism.csv")
MPEA_Results_Mico.Metabolism <- read.csv("score50/Diff_MetOrigin_20220728200453//03_Function_Analysis/MPEA_Results_Microbiota.csv")
MPEA_Results_H.Metabolism <- read.csv("score50/Diff_MetOrigin_20220728200453//03_Function_Analysis/MPEA_Results_Host.csv")
MPEA_Results_Co.Metabolism$group <- "Co_Metabolism"
MPEA_Results_Mico.Metabolism$group <- "Micro_Metabolism"
MPEA_Results_H.Metabolism$group <- "H_Metabolism"
MPEA_Results.Metabolism <- rbind(MPEA_Results_Co.Metabolism,MPEA_Results_Mico.Metabolism)%>%rbind(.,MPEA_Results_H.Metabolism)
MPEA_Results.Metabolism$logPvalue <- -log(MPEA_Results.Metabolism$Pvalue)
data <- MPEA_Results.Metabolism[order(MPEA_Results.Metabolism$Pvalue,decreasing = T),]
data$group <- factor(data$group,levels = c("Micro_Metabolism","H_Metabolism","Co_Metabolism"))

ggplot(data, aes(x=logPvalue,y=Name,fill=group)) +
  geom_bar(position="dodge",stat="identity") +scale_fill_manual(values=c("#03989E","#FF9800"))+
  ylab("")+xlab("log Pvalue")+
  scale_y_discrete(limits=unique(data$Name))+theme_classic()+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=14,face = "bold"),legend.position = "top",legend.title =element_text(size=14,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+
  theme(plot.title = element_text(size = 24,face = "bold"))


View(diff_ApcMF_metab)
Metabolites_Origin[Metabolites_Origin$Origin%in%c("Microbiota","Co-Metabolism"),]$Name

diff_ApcMF_micrometab <- diff_ApcMF_metab[diff_ApcMF_metab$Metabolite%in%c(Metabolites_Origin[Metabolites_Origin$Origin%in%c("Microbiota","Co-Metabolism"),]$Name,
                                                  "N-[(3a,5b,7b)-7-hydroxy-24-oxo-3-(sulfooxy)cholan-24-yl]-Glycine",
                                                  "Aldosterone 18-glucuronide","Dehydroisoandrosterone 3-glucuronide",
                                                  "Palmitoyl-L-carnitine","Pantothenic Acid",
                                                  "Alpha-Linolenic acid","4-Guanidinobutanoic acid",
                                                  "N-[(3a,5b,7b)-7-hydroxy-24-oxo-3-(sulfooxy)cholan-24-yl]-Glycine",
                                                  "N-Methylglutamic acid"),]
dim(diff_ApcMF_micrometab)
write.table(diff_ApcMF_micrometab,"score50/diff_ApcMF_micrometab20220729.txt", sep = '\t', col. <-  = NA, quote = FALSE)




# pvalue FC Apc/WT
Pvalue<-c(rep(0,nrow(logc))) 
FC<-c(rep(0,nrow(logc))) 
for(i in 1:nrow(logc)){
  if(sd(logc[i,1:17])==0&&sd(logc[i,18:38])==0){
    Pvalue[i] <-"NA"
    FC[i]<-"NA"
  }else{
    y=t.test(as.numeric(logc[i,1:20]),as.numeric(logc[i,21:38]))
    Pvalue[i]<-y$p.value
    FC[i]<-mean(as.numeric(c[i,1:20]))/mean(as.numeric(c[i,21:38])) 
  }
}
FC
group
FDR=p.adjust(Pvalue, "BH")
mix <- cbind(logc,FC)%>%cbind(.,FDR)
# z-score
dim(mix)
dim(logc)
dim(c)

l <- list(pos,neg)
dat <- list()
group <- list()
data <- list()
out <- list()
oplsda <- purrr::map(1:2,function(x){
  dat[[x]] <- as.data.frame(t(l[[x]]),stringsAsFactors = F)
  group[[x]] <- Group %>% factor()
  data[[x]]= as.data.frame(apply(dat[[x]],2,as.numeric))
  rownames(data[[x]]) <- rownames(dat[[x]])
  out[[x]] <- opls(data[[x]], group[[x]], orthoI = 0) 
})
vipVn <- list()
vipVn_select <- list()
vip<- do.call(rbind,purrr::map(1:2,function(x){
  vipVn[[x]] <- getVipVn(oplsda[[x]])%>%as.data.frame(.)
}))
dim(mix)
mix <- cbind(mix,vip)
#head(mix)
mix$log2FC <- log2(mix$FC)
res_MF<- mix
res_MF$type <- ifelse(res_MF$FDR < 0.05,
                      ifelse(abs(res_MF$log2FC) > 0 ,
                             ifelse(res_MF$log2FC < 0 ,'down','up'),'noSig'),'noSig')
class(res_MF)
res_MF$VIP <- res_MF$. 
res_MF$. <- NULL
# check
View(res_MF)

res_MF$metab_id <- rownames(res_MF)
res_MF<- merge(res_MF,name,all.x=T,by="metab_id")
rownames(res_MF) <- res_MF$Metabolite
# filter gene
#meta <- res_MF[res_MF$type%in%c("up","down"),] %>% subset(.,log2FC >0.5 |log2FC< -0.5)
#View(meta)
#View(res_MF)
table(res_MF$type)
p=ggplot(res_MF,aes(x = log2FC,y = -log10(FDR))) +
  geom_point(aes(color = type),alpha = 0.5,size = 3) +
  theme_bw(base_size = 16) +
  theme(aspect.ratio = 1,
        # 标题居中
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(color = 'black',size=20),
        axis.title = element_text(color = 'black')) +
  scale_color_manual(name = '',
                     # color or three types
                     values = c('up'='#DA1212','noSig'='grey','down'='#3E7C17'),
                     # legend labels
                     label = c('up'='up (29)','noSig'='noSig (819)','down'='down (17)')) +
  # 阈值分割线
  geom_hline(yintercept = -log10(0.05),lty = 'dashed',size = 0.8) +
  geom_vline(xintercept = c(-1,1),lty = 'dashed',size = 0.8) +
  ggtitle('Apc vs WT')
p
#View(meta)
#View(res_MF)
meta=res_MF[res_MF$Metabolite%in%c("LysoPE(0:0/20:3(8Z,11Z,14Z))","LysoPA(0:0/18:0)",
                                   "4,7-Megastigmadien-9-ol","4,5-Dihydrovomifoliol","(-)-Stercobilin"),]
head(meta)
p + geom_text_repel(data = meta,aes(x = log2FC,y = -log10(FDR),label = meta$Metabolite))


# fdr vipmix
mix_select1=subset(mix,FDR<0.05 & .>1 )
# fc
#mix_select2=subset(mix_select1,FC< 0.7| FC >1.5 )
subset(mix_select1,FC< 1)%>%dim()
subset(mix_select1,FC> 1)%>%dim()
#dim(mix_select)
mix_select1$metab_id <- rownames(mix_select1)
mix_select <- merge(mix_select1,name,all.x=T,by="metab_id")
mix_select <- mix_select[order(mix_select$FC,decreasing = T),]
#View(mix_select)
mix_select$group <- c(rep("Apc",29),rep("WT",16))
View(mix_select)
write.table(mix_select,"score50/ApcvsWT_select.txt", sep = '\t', col.names = NA, quote = FALSE)
write.table(res_MF,"score50/ApcvsWT_res.txt", sep = '\t', col.names = NA, quote = FALSE)


 