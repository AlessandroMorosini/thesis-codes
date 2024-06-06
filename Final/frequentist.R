library(survival)
library(survminer)
library(gridExtra)
library(smoothSurv)


############################### DATA ##############################

armA <- read.csv("ArmAB_data/ArmA.csv")
armB <- read.csv("ArmAB_data/ArmB.csv")

head(armA)
head(armB)

# Combine to one dataframe
armA$group <- "A"
armB$group <- "B"
armA$treatment <- 1
armB$treatment <- 0
armA$status <- 1-armA$cens
armB$status <- 1-armB$cens

combined_data <- rbind(armA, armB)
head(combined_data)

# Shuffle the final data
set.seed(123)  
df <- combined_data[sample(nrow(combined_data)), ]
head(df)


# Color to use throughout the thesis
bocconi_blue <- "#0062B1"


########################### KAPLAN-MEIER ############################

# Fit KM both on all data and stratified by treatment
km_fit_all <- survfit(Surv(time=time, event=status) ~ 1, data = df)
km_fit_treatment <- survfit(Surv(time=time, event=status) ~ treatment, data = df)


# Plot the two KM estimators
par(mfrow=c(1,2))

plot1 <- ggsurvplot(km_fit_all, 
                    conf.int = TRUE, 
                    ggtheme = theme_minimal() + theme(panel.grid = element_blank()),
                    title = "Kaplan-Meier Estimator",
                    xlab = "Time",
                    ylab = "Survival Probability",
                    palette = "black")

plot2 <- ggsurvplot(km_fit_treatment, 
                    conf.int = TRUE, 
                    ggtheme = theme_minimal() + theme(panel.grid = element_blank()),
                    title = "Kaplan-Meier Estimator by Treatment",
                    legend.title = "Treatment",
                    legend.labs = c("No", "Yes"),
                    xlab = "Time",
                    ylab = "Survival Probability",
                    palette = c("black", bocconi_blue))

grid.arrange(plot1$plot, plot2$plot, ncol = 2)


# LOG RANK TEST to check significance of treatment

log_rank_test <- survdiff(Surv(time=time, event=status) ~ treatment, data = df)
print(log_rank_test)

# p-value = 0.008 -> significant effect


######################### COX REGRESSION ##########################

# Fit Cox proportional hazards model
cox_model <- coxph(Surv(time=time, event=status) ~ treatment + age, data = df)

cox_summary <- summary(cox_model)
print(cox_summary)

# Extract the coefficients and confidence intervals
cox_coefficients <- cox_summary$coefficients
cox_confint <- cox_summary$conf.int

print(cox_coefficients)
print(cox_confint)


# TEST FOR PROPORTIONALITY

# Schoenfeld residuals
cox_zph <- cox.zph(cox_model)
print(cox_zph)


# Plot residuals for treatment and age
par(mfrow = c(1, 2))

plot(cox_zph, var = "treatment", main = "Schoenfeld Residuals for Treatment")

plot(cox_zph, var = "age", main = "Schoenfeld Residuals for Age")


############################# WEIBULL AFT ###############################

# Fit a Weibull AFT model
weibull_model <- survreg(Surv(time, status) ~ treatment + age, data = df, dist = "weibull")
summary(weibull_model)

# Extract coefficients and scale
coefficients <- coef(weibull_model)
scale <- 1 / weibull_model$scale

print(coefficients)
print(shape)

# Plot log(-log(KM)) against log(time)
par(mfrow = c(1, 1))

plot(log(km_fit_all$time), log(-log(km_fit_all$surv)),
     xlab = "log(Time)", ylab = "log(-log(Survival Probability))",
     main = "Log-Log Plot for Weibull Model")

# Fit a survival regression model using the smoothSurvReg to plot
smooth.lung <- smoothSurvReg(Surv(time=time, event=status) ~ treatment, data = df, init.dist = 'weibull')
cov <- matrix(c(0, 1), ncol = 1, byrow = FALSE)

par(mfrow = c(2, 2))

# Plot survival functions
surv_fit <- survfit(smooth.lung, cov = cov)
cdf_fit <- survfit(smooth.lung, cdf = TRUE, cov = cov)
hazard_fit <- hazard(smooth.lung, cov = cov)
density_fit <- fdensity(smooth.lung, cov = cov)

