library(survival)
library(survminer)
library(gridExtra)
library(smoothSurv)

armA <- read.csv("/Users/alessandromorosini/Desktop/Thesis/Codes/ArmAB_data/ArmA.csv")
armB <- read.csv("/Users/alessandromorosini/Desktop/Thesis/Codes/ArmAB_data/ArmB.csv")

# Display the first few rows of the data to confirm it's imported correctly
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

########################### KEPLEN MEIER ############################
km_fit_all <- survfit(Surv(time=time, event=status) ~ 1, data = df)
km_fit_treatment <- survfit(Surv(time=time, event=status) ~ treatment, data = df)

blue_shades <- c(
  "#003A70", # Bocconi color
  "#004080", # Slightly lighter and more vibrant
  "#004488", # Adds a bit more brightness
  "#004C99", # Even lighter and more vivid
  "#0054A3", # A noticeable step lighter with a hint of turquoise
  "#005CAA", # Lighter still, with a balanced tone
  "#0062B1", # Adds more brightness and is clearly lighter
  "#0068B8", # Approaching a sky blue shade
  "#0070BF"  # Light and vibrant, nearly a cyan shade
)
bocconi_blue <- "#0062B1"

plot1 <- ggsurvplot(km_fit_all, 
                    conf.int = TRUE, 
                    ggtheme = theme_minimal(),
                    title = "Kaplan-Meier Estimator",
                    xlab = "Time",
                    ylab = "Survival Probability",
                    palette = "black")

plot2 <- ggsurvplot(km_fit_treatment, 
                    conf.int = TRUE, 
                    ggtheme = theme_minimal(),
                    title = "Kaplan-Meier Estimator by Treatment",
                    legend.title = "Treatment",
                    legend.labs = c("No", "Yes"),
                    xlab = "Time",
                    ylab = "Survival Probability",
                    palette = c("black", bocconi_blue))

# Arrange plots side by side
grid.arrange(plot1$plot, plot2$plot, ncol = 2)

############ LOG RANK TEST

log_rank_test <- survdiff(Surv(time=time, event=status) ~ treatment, data = df)
print(log_rank_test)


######################### COX REGRESSION ##########################

# Fit Cox proportional hazards model
cox_model <- coxph(Surv(time=time, event=status) ~ treatment + age, data = df)

# Summarize the model to get the coefficients and their significance
cox_summary <- summary(cox_model)

# Print the summary of the Cox model
print(cox_summary)

# Extract the coefficients, p-values, and confidence intervals
cox_coefficients <- cox_summary$coefficients
cox_confint <- cox_summary$conf.int

# Print the coefficients, p-values, and confidence intervals
print(cox_coefficients)
print(cox_confint)

############ TEST FOR PROPORTIONALITY

# Fit Cox proportional hazards model with treatment and age as predictors
cox_model <- coxph(Surv(time=time, event=status) ~ treatment + age, data = df)

# Test proportional hazards assumption using Schoenfeld residuals
cox_zph <- cox.zph(cox_model)

# Print the test results
print(cox_zph)

par(mfrow = c(1, 2))

# Plot Schoenfeld residuals for treatment
plot(cox_zph, var = "treatment", main = "Schoenfeld Residuals for Treatment")

# Plot Schoenfeld residuals for age
plot(cox_zph, var = "age", main = "Schoenfeld Residuals for Age")

# Reset the plotting area
par(mfrow = c(1, 1))

######################### WEIBULL ###############################

# Fit a Weibull AFT model
weibull_model <- survreg(Surv(time, status) ~ treatment + age, data = df, dist = "weibull")

# Summarize the model
summary(weibull_model)

# Extract coefficients
coefficients <- coef(weibull_model)

# Extract the scale parameter (sigma)
scale <- weibull_model$scale

# Convert scale to Weibull shape parameter
shape <- 1 / scale

# Print coefficients and shape parameter
print(coefficients)
print(shape)

# Plot log(-log(KM)) against log(time)
plot(log(km_fit_all$time), log(-log(km_fit_all$surv)),
     xlab = "log(Time)", ylab = "log(-log(Survival Probability))",
     main = "Log-Log Plot for Weibull Model")

# Fit a survival regression model using the smoothSurvReg function
smooth.lung <- smoothSurvReg(Surv(time=time, event=status) ~ treatment, data = df, init.dist = 'weibull')
cov <- matrix(c(0, 1), ncol = 1, byrow = FALSE)

# Plot the survival function with covariates
par(mfrow = c(2, 2))

# Plot survival function
surv_fit <- survfit(smooth.lung, cov = cov)
cdf_fit <- survfit(smooth.lung, cdf = TRUE, cov = cov)
hazard_fit <- hazard(smooth.lung, cov = cov)
density_fit <- fdensity(smooth.lung, cov = cov)

