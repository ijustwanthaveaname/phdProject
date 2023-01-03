rm(list = ls())
setwd("./Project")
library(nlme)
library(lmerTest)
library(tidyverse)

####################
######## Example 1
####################
fm1 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
fm2 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
fm3 <- lmer(Reaction ~ Days + (Days || Subject), sleepstudy)
coef(fm1)
coef(fm2)
anova(fm1, fm2, fm3)
summary(fm1)
summary(fm2)
summary(fm3)

#########################################
# Example 2: Cross and nested effect
#########################################
dt <- read.table("http://bayes.acs.unt.edu:8083/BayesContent/class/Jon/R_SC/Module9/lmm.data.txt",
    header = TRUE, sep = ",", na.strings = "NA", dec = ".", strip.white = TRUE
)
xtabs(~ school + class, dt)
# nested model
m0 <- lmer(extro ~ open + agree + social + (1 | school / class), data = dt)
summary(m0)
# cross model
m1 <- lmer(extro ~ open + agree + social + (1 | school) + (1 | class), data = dt)
summary(m1)
# other way to show nested model
dt$classID <- paste(dt$school, dt$class, sep = ".")
xtabs(~ school + classID, dt)
# Both nested and cross modell will generate same result with m0
m2 <- lmer(extro ~ open + agree + social + (1 | school / classID), data = dt)
summary(m2)
m3 <- lmer(extro ~ open + agree + social + (1 | school) + (1 | classID), data = dt)
summary(m3)
