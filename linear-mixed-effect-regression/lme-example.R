rm(list = ls())
setwd("./Project")
library(nlme)
library(lmerTest)
library(tidyverse)
library(lme4)

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
# nested model,1 | school / class equal to 1|school + 1|school:class
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

### Example 3 Binomial Generalized Linear Mixed Models
# Generate a dataset, id (indentifying number of subject), 
# day (the day of observation), trt (treatment status: control or treat), sex (male or female)
n <- 250
id <- seq(n)
day <- 1:14
d <- expand.grid(id = id, day = day)
set.seed(1)
trt <- sample(c("control", "treat"), size = n, replace = TRUE)
sex <- sample(c("female", "male"), size = n, replace = TRUE)
d$trt <- trt[d$id]
d$sex <- sex[d$id]
d <- d[order(d$id, d$day),]
rownames(d) <- NULL
head(d, n = 15)
d$trtsex <- interaction(d$trt, d$sex)
probs <- c(0.40, 0.85, 0.30, 0.50)
names(probs) <- levels(d$trtsex)
d$p <- probs[d$trtsex]
set.seed(3)
r_probs <- rnorm(n = n, mean = 0, sd = 0.03)
d$random_p <- r_probs[d$id]
d$p <- d$p + d$random_p
d$y <- rbinom(n = nrow(d), size = 1, prob = d$p)
head(d[c("id", "day", "trt", "sex", "p", "y")])
head(d[d$id == 5, c("id", "day", "trt", "sex", "p", "y")])
# fit model
m <- glmer(y ~ trt * sex + (1|id), data = d, family = binomial)
# The fixed effect coefficients are not on the probability scale but on the log-odds, or logit, scale.
summary(m)
# The predicted log-odds a male in the control group eats vegetables is the intercept plus the coefficient for sexmale: -0.30840 + -0.63440 = -0.9428. The inverse logit is
plogis(-0.9428)
