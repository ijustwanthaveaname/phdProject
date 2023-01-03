rm(list = ls())
library(gee)
table(warpbreaks$wool, warpbreaks$tension)
xtabs(~ wool + tension, data = warpbreaks)
## 等相关
summary(gee(breaks ~ tension, id = wool, data = warpbreaks, corstr = "exchangeable"))
# 自相关
summary(gee(breaks ~ tension, id = wool, data = warpbreaks, corstr = "AR-M", Mv = 1))
str(warpbreaks)

