setwd("~/phdProject/population_genetics")
library(qtl)
# write data to local
data(listeria)
write.cross(listeria, "csv", "listeria")
data(hyper)
summary(hyper)
hyper
geno.table(hyper)
geno.table(fill.geno(hyper))
# read from local
hyper <- read.cross("csv", file="./hyperexample.csv", genotypes = c("BA", "BB", "AA"), alleles = c("A", "B"))
calc.genoprob(hyper)
# read 4-way-cross from local
sim.4way <- read.cross("csv", file="./4-way-cross.csv", genotypes = c("AC", "AD", "BC", "BD"), crosstype = "4way", alleles = c("A", "B", "C", "D"))
calc.genoprob(sim.4way)
fill.geno(sim.4way)
geno.table(sim.4way)
est.map(fake.4way, verbose=TRUE, map.function="haldane")
calc.errorlod(fake.4way, err=0.01)
# 4-way-cross example
data(fake.4way)
fake.4way
# estimate recombination fractions
fake.4way <- est.rf(fake.4way)
plotRF(fake.4way)

# estimate genetic maps
ssmap <- est.map(fake.4way, verbose=TRUE)
samap <- est.map(fake.4way, sex.sp=FALSE, verbose=TRUE)
plot(ssmap, samap)

# error lod scores
fake.4way <- calc.genoprob(fake.4way, err=0.01)
fake.4way <- calc.errorlod(fake.4way, err=0.01)
top.errorlod(fake.4way, cutoff=2.5)

# genome scan
fake.4way <- calc.genoprob(fake.4way, step=1)
out.em <- scanone(fake.4way, method="em")
## 使用Haley-Knott回归方法进行全基因组扫描
out.hk <- scanone(fake.4way, method="hk")
## 使用Multiple imputation法进行全基因组扫描
fake.4way <- sim.geno(fake.4way, step=1, n.draws=64)
out.imp <- scanone(fake.4way, method="imp")
## 比较三种方法结果的差异
plot(out.em, out.hk, out.imp, col=c("blue", "red", "green"))  
## 进行1000次Permutation test
operm <- scanone(fake.4way, method="hk", n.perm=1000)
## 获得显著性阈值
summary(operm, alpha=c(0.05, 0.2))
## 从扫描结果中挑选显著的位点
summary(out.hk, perms=operm, alpha=0.2, pvalues=TRUE)
## 获得7号染色体1.5倍LOD区间和95%贝叶斯区间
# 第一行和第三行是区间的范围，第二行是预测QTL的位置。
lodint(out.hk, chr=7, drop=1.5)
bayesint(out.hk, chr=7, prob=0.95)
## 获得区间两侧最近的标记
lodint(out.hk, chr=7, expandtomarkers=TRUE)
bayesint(out.hk, chr=7, expandtomarkers=TRUE)
## 获得离QTL最近的标记
max(out.hk)
mar <- find.marker(fake.4way, chr=7, pos=47.7)
## 统计不同基因型个体的表型
plotPXG(fake.4way, marker=mar)
## 统计不同基因型个体的效应
effectplot(fake.4way, mname1=mar)
