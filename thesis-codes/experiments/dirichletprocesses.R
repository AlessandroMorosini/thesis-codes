# Load necessary libraries
library(dirichletprocess)
library(ggplot2)

# Define the remission and censored data for Data A
remission_A <- c(1, 3, 3, 6, 7, 7, 10, 12, 14, 15, 18, 19, 22, 26, 28, 29, 34, 40, 48, 49)
censored_A <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1)
data_A <- data.frame(remission_A, censored_A)
data_A <- as.matrix(data_A)

# Define the remission and censored data for Data B
remission_B <- c(1, 1, 2, 2, 3, 4, 5, 8, 8, 9, 11, 12, 14, 16, 18, 21, 27, 31, 38, 44)
censored_B <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1)
data_B <- data.frame(remission_B, censored_B)
data_B <- as.matrix(data_B)

# Define the likelihood function for Weibull distribution with censoring
Likelihood.weibullcens <- function(mdobj, x, theta) {
  a <- theta[[1]]
  b <- theta[[2]]
  
  y <- as.numeric(x[, 1])
  b_inv <- 1 / b
  
  y_uncens <- b_inv * a * y^(a - 1) * exp(-b_inv * y^a)
  y_cens <- 1 - exp(-y^a / b)
  
  if (nrow(x) == 1) {
    if (x[, 2] == 0) return(y_uncens)
    if (x[, 2] == 1) return(y_cens)
  } else {
    y_ret <- y_uncens
    y_ret[x[, 2] == 1] <- y_cens[x[, 2] == 1]
    return(y_ret)
  }
}

# Create the mixingDistribution object for Data A
mobjA <- MixingDistribution(
  "weibullcens",         
  c(1, 2, 0.5),          
  "nonconjugate",        
  mhStepSize = c(0.11, 0.11),             
  hyperPriorParameters = c(2.222, 2, 1, 0.05)
)

class(mobjA) <- c("list", "weibullcens", "weibull", "nonconjugate")

# Create the mixingDistribution object for Data B
mobjB <- MixingDistribution(
  "weibullcens",         
  c(1, 2, 0.5),          
  "nonconjugate",        
  mhStepSize = c(0.11, 0.11),             
  hyperPriorParameters = c(2.069, 2, 1, 0.08)
)

class(mobjB) <- c("list", "weibullcens", "weibull", "nonconjugate")

# Create DirichletProcess objects
dpA <- DirichletProcessCreate(
  data_A,  
  mobjA,   
  c(2, 0.99)
)

dpB <- DirichletProcessCreate(
  data_B,  
  mobjB,   
  c(2, 0.99)
)

# Initialize and fit the Dirichlet Process objects
dpA <- Initialise(dpA)
dpA <- Fit(dpA, 500, TRUE)

dpB <- Initialise(dpB)
dpB <- Fit(dpB, 500, TRUE)

# Extract posterior samples from the fitted models
shape_params_A <- dpA$clusterParameters[[1]]
scale_params_A <- dpA$clusterParameters[[2]]

shape_params_B <- dpB$clusterParameters[[1]]
scale_params_B <- dpB$clusterParameters[[2]]

# Define functions to compute Weibull density and survival
weibull_density <- function(x, shape, scale) {
  (shape / scale) * (x / scale)^(shape - 1) * exp(-(x / scale)^shape)
}

weibull_survival <- function(x, shape, scale) {
  exp(-(x / scale)^shape)
}

# Create a sequence of x values for the plots
x_values <- seq(0, 60, length.out = 100)

# Compute the density and survival estimates for each x value for Data A
density_estimates_A <- sapply(x_values, function(x) {
  mean(sapply(1:length(shape_params_A), function(i) {
    weibull_density(x, shape_params_A[i], scale_params_A[i])
  }))
})

survival_estimates_A <- sapply(x_values, function(x) {
  mean(sapply(1:length(shape_params_A), function(i) {
    weibull_survival(x, shape_params_A[i], scale_params_A[i])
  }))
})

# Compute the density and survival estimates for each x value for Data B
density_estimates_B <- sapply(x_values, function(x) {
  mean(sapply(1:length(shape_params_B), function(i) {
    weibull_density(x, shape_params_B[i], scale_params_B[i])
  }))
})

survival_estimates_B <- sapply(x_values, function(x) {
  mean(sapply(1:length(shape_params_B), function(i) {
    weibull_survival(x, shape_params_B[i], scale_params_B[i])
  }))
})

# Create data frames for plotting
density_df_A <- data.frame(x = x_values, y = density_estimates_A, dataset = "A")
density_df_B <- data.frame(x = x_values, y = density_estimates_B, dataset = "B")

survival_df_A <- data.frame(x = x_values, y = survival_estimates_A, dataset = "A")
survival_df_B <- data.frame(x = x_values, y = survival_estimates_B, dataset = "B")

# Combine data frames for ggplot
density_df <- rbind(density_df_A, density_df_B)
survival_df <- rbind(survival_df_A, survival_df_B)

# Plot the density function using ggplot2
density_plot <- ggplot(density_df, aes(x = x, y = y, color = dataset)) +
  geom_line() +
  labs(title = "Density Function", x = "x", y = "Density") +
  theme_minimal()

# Plot the survival function using ggplot2
survival_plot <- ggplot(survival_df, aes(x = x, y = y, color = dataset)) +
  geom_line() +
  labs(title = "Survival Function", x = "x", y = "Survival Probability") +
  theme_minimal()

# Print the plots
print(density_plot)
print(survival_plot)
