setwd("~/phdProject/meta-analysis")
rm(list = ls())
library(metafor)
library(meta)
library(metap)

twodata <- read.csv("test_bi.csv")
contidata <- read.csv("test_con.csv")
Rate <- read.csv("test_rate1.csv")
hrdata <- read.csv("test_hr.csv")
# combine pvalue, stuffer method
# data(dat.metap)
# teachexpect <- dat.metap$teachexpect
p.test <- runif(10)
sample.size <- sample(1:100, 10)
weight.root <- sample.size^(1/2)
sumz(p.test, weight.root)
# binary variable
m1 <- metabin(
    Tevent, Ttotal, Cevent, Ctotal,
    data = twodata, sm = "OR",
    studlab = paste(Study, Year),
    label.e = "Intervention",
    label.c = "Control")
# continuous variable
m2 <- metacont(
    Tsample, Tmean, Tsd, Csample, Cmean, Csd,
    data = contidata, sm = "SMD",
    studlab = paste(Author, Year),
    label.e = "Intervention",
    label.c = "Control")
# 
m3 <- metaprop(
    Case, Number,
    data = Rate, sm = "PFT",
    studlab = paste(Study, Year)
)
# survival data
m4 <- metagen(
    LogHR, SE,
    data = hrdata, sm = "HR",
    studlab = paste(Study, Year)
)
summary(m1)
summary(m2)
summary(m3)
# draw forest plot
forest(m1, col.square = "blue", col.diamond = "red") 
forest(m2, col.square = "blue", col.diamond = "red")
forest(m3)

# draw funnel plot
funnel(m1)
funnel(m2)

# trimfill
trimfill(m1)
funnel(trimfill(m1))
funnel(trimfill(m2))

#begg and egger test
metabias(m1, method.bias = "rank")
metabias(m1, method.bias = "linreg")

# labbe plot
labbe(m1)
labbe(m2)

# Galbraith plot
radial(m1)
radial(m2)

# cumulative meta analysis
metacum(m1, pooled = "random") # or common
forest(metacum(m1, pooled = "random"))

# sensitive analysis
metainf(m1, pooled = "random")
forest(metainf(m1, pooled = "random"))

# subgroup analysis
m1 <- metabin(
    Tevent, Ttotal, Cevent, Ctotal,
    data = twodata, sm = "OR",
    studlab = paste(Study, Year),
    label.e = "Intervention",
    label.c = "Control",
    byvar = Dose,
    comb.fixed = FALSE)
forest(m1)

# meta regression
m <- metareg(m1, latitude+allocation)
bubble(m)


#########################################
# Using metafor to do meta-analysis
#########################################
# Combine regression coefficients
png("combine_beta2.png")
metabirth3 <- data.frame(
    author = c("Fenster", "Callan", "Lopez", "Lauritzen", "Lauritzen", "Our study"),
    year = c(2006, 2016, 2011, 2017, 2017, 2018),
    country = c("Australia", "Australia", "Spain"),
    beta = rnorm(6),
    se = rnorm(6)
)
# DL means random effect model，fiexd effect model need to set `method="FE"`
metamod <- rma(yi = beta, data = metabirth3, sei = se, method = "DL")
summary(metamod)
forestplot <- forest(metamod, refline = 1, mlab = "Random-effect Model for All Studies", slab = paste(metabirth3$author, metabirth3$year, sep = ","), xlab = "β", showweights = T)
# text(-500, 6:1, pos = 2, metabirth3$country)
# text(c(-1600, -500, 300, 800), 8, pos = c(4, 2, 4, 4), c("Author(s) and Year", "Location", "Weight", "β[95%CI]"), cex = 1, font = 2)
dev.off()
#######################
# Publication bias
#######################
# Calculate Fail-Safe N
data(dat.bcg)
### calculate log relative risks and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
fsn(yi, vi, data=dat)
fsn(yi, vi, data=dat, type="Orwin")
fsn(yi, vi, data=dat, type="Rosenberg")
