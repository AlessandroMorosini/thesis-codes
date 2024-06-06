library(survival)
library(DPWeibull)
library(ggplot2)

set.seed(1)

########################### SYNTHETIC DATA ########################

# Define the parameters for log normals
w1 <- 0.8
mu1 <- 0
sigma1 <- sqrt(0.25)
w2 <- 0.2
mu2 <- 1.2
sigma2 <- sqrt(0.02)

# probability of death -> P[censoring] = 1-prob
prob <- 0.95

# function to generate data
generate_data <- function(n_samples, mu1, mu2, sigma1, sigma2, prob) {
  samples1 <- rlnorm(n = round(w1 * n_samples), meanlog = mu1, sdlog = sigma1)
  samples2 <- rlnorm(n = round(w2 * n_samples), meanlog = mu2, sdlog = sigma2)
  synthetic_data <- c(samples1, samples2)
  data <- data.frame(time = synthetic_data)
  data$status <- rbinom(n_samples, 1, prob)
  return(data)
}

# create dataset
n_samples <- 200
df <- generate_data(n_samples, mu1, mu2, sigma1, sigma2, prob)
df$time <- as.numeric(df$time)
df$status <- as.numeric(df$status)

min_time <- min(df$time)
max_time <- max(df$time)

xs <- seq(min_time, max_time, length.out = 50)

# f(t)
pdf_mix <- w1 * dlnorm(xs, meanlog = mu1, sdlog = sigma1) +
  w2 * dlnorm(xs, meanlog = mu2, sdlog = sigma2)

# S(t)
survival_mix <- 1 - (
  w1 * plnorm(xs, meanlog = mu1, sdlog = sigma1) + 
    w2 * plnorm(xs, meanlog = mu2, sdlog = sigma2)
)

# h(t)
hazard_mix <- pdf_mix / survival_mix


##### KM FOR COMPARISON

km_fit <- survfit(Surv(time=time, event=status) ~ 1, data = df)

# Extract KM data for manual plotting
km_data <- data.frame(
  time = km_fit$time,
  surv = km_fit$surv,
  upper = km_fit$upper,
  lower = km_fit$lower
)

true_surv_data <- data.frame(
  time = xs,
  surv = survival_mix
)

p <- ggplot() +
  geom_step(data = km_data, aes(x = time, y = surv, color = "KM Estimator")) +
  geom_ribbon(data = km_data, aes(x = time, ymin = lower, ymax = upper), alpha = 0.1, color = NA, fill = "#0062B1", show.legend = FALSE) +
  geom_line(data = true_surv_data, aes(x = time, y = surv, color = "True Survival"), size = 0.5) +
  scale_color_manual(values = c("KM Estimator" = "#0062B1", "True Survival" = "black")) +
  labs(
    title = "Kaplan-Meier Estimator",
    x = "Time",
    y = "S(t)"
  ) +
  theme_minimal() +
  theme(legend.title = element_blank())

print(p)

####################################################################
############################### DPMM ###############################
############################# SYTHETHIC ############################
####################################################################

q95 <- quantile(subset(df, status == 1)$time, 0.95)

dpw <- dpweib(Surv(time=time, event=status)~1,
              data=df,
              high.pct=q95, 
              predtime=xs,
              iteration=10000)

par(mfrow=c(1,3))

# Plot density
plot(xs, dpw$dpred, type="l", lwd=2, col="#0062B1", xlab="t", ylab="f(t)", main="Estimated Density", axes=TRUE, ylim=c(0, 1))
lines(xs, dpw$dpredu, lty=2, col="#0062B1")
lines(xs, dpw$dpredl, lty=2, col="#0062B1")
lines(xs, pdf_mix,  lwd=1.5, col="black")
legend("topright", legend=c("Estimated", "95% CI", "True"), col=c("#0062B1", "#0062B1", "black"), lty=c(1, 2, 1), lwd=c(2, 1, 2))

# Survival density
plot(xs, dpw$Spred, type="l", lwd=2, col="#0062B1", xlab="t", ylab="S(t)", main="Estimated Survival", axes=TRUE, ylim=c(0, 1))
lines(xs, dpw$Spredu, lty=2, col="#0062B1")
lines(xs, dpw$Spredl, lty=2, col="#0062B1")
lines(xs, survival_mix,  lwd=1.5, col="black")
legend("topright", legend=c("Estimated", "95% CI", "True"), col=c("#0062B1", "#0062B1", "black"), lty=c(1, 2, 1), lwd=c(2, 1, 2))

# Hazard rate
plot(xs, dpw$hpred, type="l", lwd=2, col="#0062B1", xlab="t", ylab="h(t)", main="Estimated Hazard", axes=TRUE, ylim=c(0, 5))
lines(xs, dpw$hpredu, lty=2, col="#0062B1")
lines(xs, dpw$hpredl, lty=2, col="#0062B1")
lines(xs, hazard_mix,  lwd=1.5, col="black")
legend("topright", legend=c("Estimated", "95% CI", "True"), col=c("#0062B1", "#0062B1", "black"), lty=c(1, 2, 1), lwd=c(2, 1, 2))


####################################################################
############################### DPMM ###############################
############################# REAL DATA ###############################
####################################################################

armA <- read.csv("ArmAB_data/ArmA.csv")
armB <- read.csv("ArmAB_data/ArmB.csv")

armA$status <- 1-armA$cens
armB$status <- 1-armB$cens

head(armA)
head(armB)

q95A <- quantile(armA$time, 0.95)
q95B <- quantile(armB$time, 0.95)

max_A <- max(armA$time)
max_B <- max(armB$time)
max_AB <- max(max_A, max_B)

dpwA <- dpweib(Surv(time=time, event=status)~1,
               data=armA,
               predtime=seq(0, max_AB, 50),
               burnin = 5000,
               a=2,
               b=0.5,
               iteration=10000)

dpwB <- dpweib(Surv(time=time, event=status)~1,
               data=armB,
               predtime=seq(0, max_AB, 50),
               burnin = 5000,
               a=2,
               b=0.5,
               iteration=10000)


######################### PLOTTING ######################

# Set layout for the plots and the shared legend
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), heights = c(4, 1))
par(mar=c(4,4,2,1), oma=c(0,0,2,0))

# Survival density plot
plot(dpwA$predtime, dpwA$Spred, type="l", lwd=2, col="#0062B1", xlab="t", ylab="S(t)", main="Survival Function", ylim=c(0, 1))
lines(dpwA$predtime, dpwA$Spredu, lty=2, col="#0062B1")
lines(dpwA$predtime, dpwA$Spredl, lty=2, col="#0062B1")

lines(dpwB$predtime, dpwB$Spred, lwd=2, col="black")
lines(dpwB$predtime, dpwB$Spredu, lty=2, col="black")
lines(dpwB$predtime, dpwB$Spredl, lty=2, col="black")

# Hazard rate plot
plot(dpwA$predtime, dpwA$hpred, type="l", lwd=2, col="#0062B1", xlab="t", ylab="h(t)", main="Hazard Rate", ylim=c(0, 0.01))
lines(dpwA$predtime, dpwA$hpredu, lty=2, col="#0062B1")
lines(dpwA$predtime, dpwA$hpredl, lty=2, col="#0062B1")

lines(dpwB$predtime, dpwB$hpred, lwd=2, col="black")
lines(dpwB$predtime, dpwB$hpredu, lty=2, col="black")
lines(dpwB$predtime, dpwB$hpredl, lty=2, col="black")

# Add the shared legend at the bottom
par(mar=c(0,0,0,0))
plot.new()
legend("center", legend=c("Group A", "Group B"), 
       col=c("black", "#0062B1"), lty=1, lwd=2, ncol=2, bty="n")

