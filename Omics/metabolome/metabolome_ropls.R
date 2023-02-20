setwd("/home/weishen/phdProject/Omics/metabolome")
## 方法一 使用mixOmics进行PLS-DA
library(mixOmics)
data(breast.tumors)
X <- breast.tumors$gene.exp
Y <- breast.tumors$sample$treatment
plsda.breast <- plsda(X, Y, ncomp = 2)
pdf("plsda_breast.pdf")
plotIndiv(plsda.breast, ind.names = TRUE, ellipse = TRUE, legend = TRUE)
dev.off()
# 得到VIP值
library(RVAideMemoire)
# pdf("plsda_breast_vip.pdf")
# PLSDA.VIP(plsda.breast, graph = TRUE)
# dev.off()
# 对于两组比较，VIP >= 1认为有显著差异
vip <- PLSDA.VIP(plsda.breast, graph = FALSE)

## 方法二 使用ropls进行PCA, PLS, PLS-DA, OPLS, or OPLS-DA
library(ropls)
library(ggplot2)
data(sacurine)
# 提各指标读数
dataMatrix <- sacurine$dataMatrix
# 提取性别数据列
genderFc <- sacurine$sampleMetadata[, "gender"]
### PLS-DA
# 图解释
# 结果中，R2X和R2Y分别表示所建模型对X和Y矩阵的解释率，Q2表示模型的预测能力，它们的值越接近于1表明模型的拟合度越好，训练集的样本越能够被准确划分到其原始归属中。
# 显著性诊断(左上)：实际和模拟模型的R2Y和Q2Y值经随机排列后的散点图，模型R2Y和Q2Y(散点)大于真实值时(横线)，表明产生过拟合2。
# Inertia(惯量)柱形图(右上)：通过展示累计解释率评估正交组分是否足够
# 离群点展示(左下)：通过scoreMN和loadingMN计算出各样本在投影平面及正交平面的坐标，并标明相互差异较大的样本。
# x-score plot(右下)：各样本在PLS-DA轴中的坐标；R2X、R2Y等值展示在下方，
pdf("plsda_ropls.pdf")
# 不指定或orthoI = 0，执行PLS-DA
sacurine.plsda <- opls(dataMatrix, genderFc, orthoI = 0) 
dev.off()
# 根据VIP筛选, 一般选择VIP>=1的代谢物
# 根据fold change筛选，一般选择fold change >= 2，fold change <= .5的
# 若生物学重复 < 3，根据fold change进行筛选，若 > 3则根据两者同时筛选
vipVn <- sacurine.plsda@vipVn  # getVipVn()
vipVn_select <- vipVn[vipVn > 1]

# cbind Metadata and vipVn
vipVn_select <- cbind(sacurine$variableMetadata[names(vipVn_select), ], vipVn_select)

# rename
names(vipVn_select)[4] <- 'VIP'

# reorder
vipVn_select <- vipVn_select[order(vipVn_select$VIP, decreasing = TRUE), ]

pdf("plsda_ropls_xloading.pdf")
plot(sacurine.plsda, typeVc = "x-loading") #展示前10个
dev.off()
head(vipVn_select)

### OPLS-DA
# Orthogonal partial least squares(OPLS) 
# 将观测值矩阵X的差异分为两个部分：第一部分代表与Y相关的差异，第二部分代表与Y不相关（正交垂直）的差异，
# 结果展示时需要结合起来讨论；由于OPLS区分了无关变量数据，从而使模型更加容易解读。
# 通过orthoI指定正交组分数目
# orthoI = NA时，执行OPLS，并通过交叉验证自动计算适合的正交组分数
# predI 响应值组数：默认设定为1(for OPLS)
pdf("oplsda_ropls.pdf")
sacurine.oplsda <- opls(dataMatrix, genderFc, predI = 1, orthoI = NA)
dev.off()
# plot(sacurine.oplsda, typeVc = 'x-loading')
# VIP提取分析同上
# 执行OPLS后的数据提取，与PLS和PCA略有不同，需要同时考虑得分矩阵和正交矩阵。
# Score矩阵
scoreMN <- cbind(sacurine.oplsda@scoreMN[, 1], sacurine.oplsda@orthoScoreMN[, 1])

# rename 
colnames(scoreMN) <- c("h1", paste0("o", 1))
x_lab <- paste0("t1(", sacurine.oplsda@modelDF[1, "R2X"] * 100, "%)")
y_lab <- "to1"
head(scoreMN)

# Loading矩阵
loadMN <- cbind(sacurine.oplsda@loadingMN[, 1], sacurine.oplsda@orthoLoadingMN[, 1])

# rename
colnames(loadMN) <- c("h1", paste0("o", 1))

x_lab <- sacurine.plsda@modelDF[1, "R2X"] * 100
y_lab <- paste0("pOrtho1(", sacurine.oplsda@modelDF[1, "R2X"] * 100, "%)")
head(loadMN)
