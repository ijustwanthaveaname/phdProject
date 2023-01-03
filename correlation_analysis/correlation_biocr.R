rm(list=ls())
setwd("./correlation_analysis")
library(correlation)
library(tidyverse)
library(see)  # It is a very useful r package to visualize model results
bi.cor <- correlation(iris, method = "biweight")
correlation(iris, method = "spearman")
correlation(iris, method = "pearson")
library(WGCNA)
bicor(iris)

png("cor_see.png")
bi.cor %>%
  summary(redundant = TRUE) %>%
  plot()
dev.off()
summary(bi.cor)
