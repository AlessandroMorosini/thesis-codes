library(dirichletprocess)

# Set a seed for reproducibility
set.seed(6)  

########################### DATA #############################

# Define the parameters for log normals
w1 <- 0.8
mu1 <- 0
sigma1 <- sqrt(0.25)
w2 <- 0.2
mu2 <- 1.2
sigma2 <- sqrt(0.02)

# probability of death -> P[censoring] = 1-prob
prob <- 0.8

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

# f(t)
min_data = min(df$time)
max_data = max(df$time)
t <- seq(min_data, max_data, length.out = 200)
pdf1 <- dlnorm(t, meanlog = mu1, sdlog = sigma1)
pdf2 <- dlnorm(t, meanlog = mu2, sdlog = sigma2)
pdf_mix <- w1 * pdf1 + w2 * pdf2


# S(t)
survival1 <- 1 - plnorm(t, meanlog = mu1, sdlog = sigma1)
survival2 <- 1 - plnorm(t, meanlog = mu2, sdlog = sigma2)
survival_mix <- w1 * survival1 + w2 * survival2

# h(t)
hazard_mix = pdf_mix / survival_mix

# function to infer parameters of fitted weibull on data
estimate_weibull_params <- function(samples){
  # Calculate sample median and IQR
  median <- median(samples)
  q1 <- quantile(samples, 0.25)
  q3 <- quantile(samples, 0.75)
  iqr <- q3 - q1
  
  # Define equations to solve for alpha and lambda
  equations <- function(params) {
    alpha <- params[1]
    lambda <- params[2]
    eq1 <- median - (lambda * -log(0.5))^(1 / alpha)
    eq2 <- iqr - ((lambda * -log(0.25))^(1 / alpha) - (lambda * -log(0.75))^(1 / alpha))
    return(c(eq1, eq2))
  }
  
  # Initial guess for alpha and lambda
  initial_guess <- c(1, 1)
  
  # Solve the equations
  solution <- nleqslv::nleqslv(initial_guess, equations)
  alpha <- solution$x[1]
  lambda <- solution$x[2]
  
  return(c(alpha, lambda))
}

# eestimated parameters
estimated_params <- estimate_weibull_params(df$time)
estimated_alpha <- estimated_params[1]
estimated_lambda <- estimated_params[2]

# check fit of weibull on data, meeeh
ggplot(df, aes(x = time)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "blue", alpha = 0.5) +
  stat_function(fun = function(x) dweibull(x, shape = estimated_alpha, scale = estimated_lambda),
                color = "red", size = 1) +
  labs(title = "Histogram of Times with Weibull PDF",
       x = "Time",
       y = "Density") +
  theme_minimal()

###### plot the true density and survival

# Create data frame for plotting
data <- data.frame(
  t = t,
  pdf1 = pdf1,
  pdf2 = pdf2,
  pdf_mix = pdf_mix,
  survival1 = survival1,
  survival2 = survival2,
  survival_mix = survival_mix
)

# Plot the PDF
ggplot(data, aes(x = t)) +
  # geom_line(aes(y = pdf1, color = "PDF1")) +
  # geom_line(aes(y = pdf2, color = "PDF2")) +
  geom_line(aes(y = pdf_mix, color = "Combined PDF")) +
  labs(title = "Probability Density Function (PDF)",
       y = "Density",
       color = "Legend") +
  theme_minimal()

# Plot the Survival Function
ggplot(data, aes(x = t)) +
  # geom_line(aes(y = survival1, color = "Survival1")) +
  # geom_line(aes(y = survival2, color = "Survival2")) +
  geom_line(aes(y = survival_mix, color = "Combined Survival")) +
  labs(title = "Survival Function",
       y = "Survival Probability",
       color = "Legend") +
  theme_minimal()

km_fit <- survfit(Surv(time=time) ~ 1, data = df)
# Plot the true survival function
plot(t, survival_mix, type = "l", col = "#0062B1", lwd = 2,
     xlab = "Time", ylab = "Survival Probability", main = "True vs. Kaplan-Meier Survival Function")
# Plot the Kaplan-Meier survival function
lines(km_fit, col = "black", lwd = 2, lty = 1)
# Add a legend
legend("topright", legend = c("True Survival", "Kaplan-Meier Survival"),
       col = c("#0062B1", "black"), lwd = 2, lty = c(1, 1))


################# DPMM WITH NO CENSORING #################

set.seed(3)  

g0Priors <- c(12, 3, 10)
alphaPriors = c(2, 0.1)
mhStepSize = c(0.1, 0.1)
hyperPriorParameters = c(20, 2, 1, 0.25)

dpw <- DirichletProcessWeibull(df$time, 
                               g0Priors = g0Priors,
                               alphaPriors = alphaPriors,
                               mhStepSize = mhStepSize,
                               hyperPriorParameters = hyperPriorParameters)
its <- 500
dpw <- Fit(dpw, its)
plot(dpw, data_method="hist")

########## DENSITY COMPARISON ###########

compute_density_matrix <- function(t, num_samples = 20) {
  density_matrix <- sapply(1:num_samples, function(i) {
    sapply(t, PosteriorFunction(dpw))
  })
  return(density_matrix)
}


plot_dpmm_density <- function(t, pdf_mix, density_matrix, df, title = "Density Comparison", show_legend = TRUE) {
  
  # Calculate the mean and the confidence intervals (e.g., 95% CI)
  mean_density <- rowMeans(density_matrix)
  lower_ci <- apply(density_matrix, 1, quantile, probs = 0.025)
  upper_ci <- apply(density_matrix, 1, quantile, probs = 0.975)
  
  # First plot: Density Comparison
  plot(t, pdf_mix, col = "black", lwd = 2, lty = 1, type = "l", ylim = c(0, 1), 
       xlab = "Time", ylab = "Density", main = title)
  lines(t, mean_density, col = "#0062B1", lwd = 2)
  lines(t, lower_ci, col = "#0062B1", lwd = 2, lty = 2)
  lines(t, upper_ci, col = "#0062B1", lwd = 2, lty = 2)
  hist(df$time, probability = TRUE, breaks = 30, col = rgb(0.7, 0.7, 0.7, 0.3), border = "white", add = TRUE)
  if (show_legend) {
    legend("topright", legend = c("True Density", "Estimated Density", "Confidence Interval"), 
           col = c("black", "#0062B1", "#0062B1"), lwd = 2, lty = c(1, 1, 2))
  }
}

density_matrix = compute_density_matrix(t)
plot_dpmm_density(t, pdf_mix, density_matrix, df)

########### COMPARE SURVIVAL ###########

plot_dpmm_survival <- function(t, survival_mix, density_matrix) {
  # Estimate the CDF from the density values
  cdf_matrix <- apply(density_matrix, 2, cumsum) * diff(t)[1]
  
  # Calculate survival function as 1 - CDF
  survival_matrix <- 1 - cdf_matrix
  
  mean_survival <- rowMeans(survival_matrix)
  lower_survival_ci <- apply(survival_matrix, 1, quantile, probs = 0.025)
  upper_survival_ci <- apply(survival_matrix, 1, quantile, probs = 0.975)
  
  # Plot the Estimated Survival Function with Confidence Intervals
  plot(t, mean_survival, type = "l", col = "#0062B1", lwd = 2, ylim = c(0, 1), 
       ylab = "Survival Probability", xlab = "Time", main = "Survival Curve with Confidence Intervals")
  lines(t, lower_survival_ci, col = "#0062B1", lwd = 2, lty = 2)
  lines(t, upper_survival_ci, col = "#0062B1", lwd = 2, lty = 2)
  lines(t, survival_mix, col = "black", lwd = 2, lty = 1)
  legend("topright", legend = c("Estimated Survival", "Confidence Interval", "True Survival"), 
         col = c("#0062B1", "#0062B1", "black"), lwd = 2, lty = c(1, 2, 1))
}

plot_dpmm_survival(t, survival_mix, density_matrix)


########### COMPARE HAZARDS #################

plot_dpmm_hazard <- function(t, hazard_mix, density_matrix) {
  
  cdf_matrix <- apply(density_matrix, 2, cumsum) * diff(t)[1]
  survival_matrix <- 1 - cdf_matrix
  
  # Calculate the hazard rate for each posterior sample
  hazard_matrix <- density_matrix / survival_matrix
  
  # Calculate the mean and the confidence intervals (e.g., 95% CI)
  mean_hazard <- rowMeans(hazard_matrix)
  lower_hazard_ci <- apply(hazard_matrix, 1, quantile, probs = 0.025)
  upper_hazard_ci <- apply(hazard_matrix, 1, quantile, probs = 0.975)
  
  # Plot the Hazard Rate with estimated hazard and confidence intervals
  plot(t, hazard_mix, col = "black", lwd = 2, lty = 1, type = "l", ylim = c(0, max(c(hazard_mix, upper_hazard_ci))), 
       xlab = "Time", ylab = "Hazard Rate", main = "Hazard Rate Comparison")
  lines(t, mean_hazard, col = "#0062B1", lwd = 2)
  lines(t, lower_hazard_ci, col = "#0062B1", lwd = 2, lty = 2)
  lines(t, upper_hazard_ci, col = "#0062B1", lwd = 2, lty = 2)
  legend("topright", legend = c("True Hazard", "Estimated Hazard", "Confidence Interval"), 
         col = c("black", "#0062B1", "#0062B1"), lwd = 2, lty = c(1, 1, 2))
}

plot_dpmm_hazard(t, hazard_mix, density_matrix)

############### COMPARE SURVIVAL UNDER DIFFERNT PRIORS ##############
par(mfrow = c(1, 3))
its <- 1000

alphaPriors = c(1.5, 0.5)
dpw <- DirichletProcessWeibull(df$time, 
                               g0Priors = g0Priors,
                               alphaPriors = alphaPriors,
                               mhStepSize = mhStepSize,
                               hyperPriorParameters = hyperPriorParameters)
dpw <- Fit(dpw, its)
density_matrix1 <- compute_density_matrix(t)
plot_dpmm_density(t, pdf_mix, density_matrix1, df, title=("(1.5, 0.5)"),
                  show_legend=FALSE)


alphaPriors = c(2, 0.1)
dpw <- DirichletProcessWeibull(df$time, 
                               g0Priors = g0Priors,
                               alphaPriors = alphaPriors,
                               mhStepSize = mhStepSize,
                               hyperPriorParameters = hyperPriorParameters)
dpw <- Fit(dpw, its)
density_matrix2 <- compute_density_matrix(t)
plot_dpmm_density(t, pdf_mix, density_matrix2, df, title=("(2, 0.1)"),
                  show_legend=FALSE)


alphaPriors = c(3, 0.05)
dpw <- DirichletProcessWeibull(df$time, 
                               g0Priors = g0Priors,
                               alphaPriors = alphaPriors,
                               mhStepSize = mhStepSize,
                               hyperPriorParameters = hyperPriorParameters)
dpw <- Fit(dpw, its)
density_matrix3 <- compute_density_matrix(t)
plot_dpmm_density(t, pdf_mix, density_matrix3, df, title=("(3, 0.05)"),
                  show_legend=FALSE)







#####################################################################
######################### CENSORING HERE ############################

Likelihood.weibullcens <- function(mdobj, x, theta){
  a = theta[[1]][,,,drop=TRUE]
  b = theta[[2]][,,,drop=TRUE]
  y <- as.numeric(
    b^(-1) * a * x[,1]^(a-1) * exp(-b^(-1) * x[, 1]^a)
  )
  y_cens <- as.numeric(1 - exp(-x[,1]^a / b))
  if(nrow(x) == 1){
    if(x[,2] == 0) return(y)
    if(x[,2] == 1) return(y_cens)
  } else{
    y_ret <- y
    y_ret[x[, 2] == 1] <- y_cens[x[, 2]==1]
    return(y_ret)
  } }

mdobjA <- MixingDistribution("weibullcens", 
                             c(1,2,0.5), 
                             "nonconjugate",
                             mhStepSize=c(0.11,0.11),
                             hyperPriorParameters=c(2.222, 2, 1, 0.05)
)

class(mdobjA) <- c("list", "weibullcens", "weibull", "nonconjugate")

dpA <- DirichletProcessCreate(df, mdobjA, c(1, 1))
dpA <- Initialise(dpA)