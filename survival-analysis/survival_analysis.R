setwd("~/Project/survival-analysis")
library(survival)
library(survminer)
library(eha)
library(coxme) 
### run PH model using eha package
head(oldmort)
# enter: start time, exit: end time, event: true means death, false means others
fit.g <- phreg(Surv(enter - 60, exit - 60, event) ~ sex + civ + region, 
             dist = "gompertz", data = oldmort)
summary(fit.g)
### run cox model using eha package
fit.c <- coxreg(Surv(enter - 60, exit - 60, event) ~ sex + civ + region, data = oldmort)
summary(fit.c)
### Fit a mixed effects Cox model
# A non-significant institution effect
fit1 <- coxph(Surv(time, status) ~ ph.ecog + age, data=lung,
              subset=(!is.na(inst)))
fit2 <- coxme(Surv(time, status) ~ ph.ecog + age + (1|inst), lung)
anova(fit1, fit2)
# Shrinkage effects (equivalent to ridge regression)
temp <- with(lung, scale(cbind(age, wt.loss, meal.cal)))
rfit <- coxme(Surv(time, status) ~ ph.ecog + (temp | 1), data=lung)

### fit the km model for the data
# set ~1 because no x variables
# default type is km
head(lung)
km.model <- survfit(Surv(time, status) ~ 1, type = "kaplan-meier", data = lung)
summary(km.model)
# include 1 categorical variable with in KM model
km.model2 <- survfit(Surv(time, status) ~ sex, data = lung)
summary(km.model2)
# plot
png("base-kmplot.png")
plot(
    km.model, conf.int = T, xlab = "Times", 
    ylab = "%ALIVE = S(t)", main = "KM-Model",
    las = 1, mark.time = TRUE
    )
abline(h = 0.5, col = "red")
dev.off()
# plot with sex
png("base-2-kmplot.png")
plot(
    km.model2, conf.int = F, xlab = "Times", 
    ylab = "%ALIVE = S(t)", main = "KM-Model",
    col = c("red", "blue"), 
    las = 1, mark.time = TRUE,
    lwd = 2
    )
legend(
    800, 0.8, legen = c("Male","Female"), lty = 1, lwd = 2, 
    col = c("red", "blue", bty = "", cex = 0.6))
dev.off()
### log rank test which used for testing if two survival curve are significantly different
# other methods: Cochran–Mantel–Haenszel, wilcoxon test
survdiff(Surv(time, status) ~ sex, data = lung)
# other ways to plot
ggsurvplot(km.model, data = lung)
ggsave("km.png")
dev.off()
ggsurvplot(km.model2,
           conf.int=TRUE, # add confidence intervals
           pval=TRUE, # show the p-value for the log-rank test
        #    risk.table=TRUE, # show a risk table below the plot
           legend.labs=c("Male", "Female"), # change group labels
           legend.title="Sex",  # add legend title
           palette=c("dodgerblue4", "orchid2"), # change colors of the groups
           title="Kaplan-Meier Curve for Lung Cancer Survival", # add title to plot
           risk.table.height=.2)
ggsave("km2.png")
dev.off()
# Cox regression
res.cox <- coxph(Surv(time, status) ~ age + sex + ph.ecog, data = lung)
res.cox2 <- coxph(Surv(time, status) ~ age + sex, data = lung)
res.cox3 <- coxph(Surv(time, status) ~ age, data = lung)
hr <- predict(res.cox2, lung, type = "risk")
hr[1]/hr[2]
# Calculate baseline hazard
bhest <- basehaz(res.cox)
# do the LRT
anova(res.cox2, res.cox3, test = "LRT")
fit.km <- survfit(Surv(time, status) ~ sex, data = lung)
fit.cox <- survfit(res.cox, data = lung)
summary(res.cox)
# plot
ggsurvplot(fit.cox, data = lung)
ggsave("coxplot.png")
dev.off()
# elastic cox regression
library(glmnet)
data(CoxExample)
x <- CoxExample$x
y <- CoxExample$y
x[1:5, 1:5]
y[1:5, ]
fit <- glmnet(x, y, family = "cox")
# use "response" to get relative risk
fit.pred <- predict(fit, x[1:5, ], type = "response", s = 0.05)
fit$lambda
cvfit <- cv.glmnet(x, y, family = "cox")
# lambda.min is the value of λ that gives minimum mean cross-validated error
# lambda.1se is the value of λ that gives the most regularized model 
## such that the cross-validated error is within one standard error of the minimum.
predict(cvfit, x, type = "response", s = "lambda.min")
##########################################
####Checking COX PH MODEL ASSUMPTIONS
##########################################
# Using MARTINGALE residuals
png("cox-martingale2.png")
plot(
    predict(res.cox3), residuals(res.cox3, type = "martingale"), 
    xlab = "fitted values", ylab = "Martingale residuals",
    main = "Residual Plot", las = 1
    )
abline(h = 0)
lines(smooth.spline(predict(res.cox3),
      residuals(res.cox3, type = "martingale")), col = "red")
dev.off()
# Using deviance plot
png("cox-deviance.png")
plot(
    predict(res.cox3), residuals(res.cox3, type = "deviance"), 
    xlab = "fitted values", ylab = "Deviance residuals",
    main = "Residual Plot", las = 1
    )
abline(h = 0)
lines(smooth.spline(predict(res.cox3),
      residuals(res.cox3, type = "deviance")), col = "red")
dev.off()
#############################################
###Checking PROPORTIONAL HAZARDS ASSUMPTION
#############################################
test.ph <- cox.zph(res.cox2)
ggcoxzph(test.ph)
ggsave("Schoenfeld-test.png")
dev.off()
# See if 95percent confidence interval contains zeros most of the time
png("beta-change-over-time.png")
par(mfrow = c(2, 1))
plot(test.ph)
dev.off()
