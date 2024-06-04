# Load necessary library
library(MASS)
library(ggplot2)
library(dirichletprocess)
library(survival)

######################## DATA GENERATION ###########################

set.seed(123)  # Set a seed for reproducibility

# Define the parameters
mu1 <- 0
sigma1 <- sqrt(0.25)
mu2 <- 1.2
sigma2 <- sqrt(0.02)
n_samples <- 60  # Define the number of samples you want to generate

# Generate samples from the two distributions
samples1 <- rlnorm(n = round(0.8 * n_samples), meanlog = mu1, sdlog = sigma1)
samples2 <- rlnorm(n = round(0.2 * n_samples), meanlog = mu2, sdlog = sigma2)

# Combine the samples
synthetic_data <- c(samples1, samples2)

# Convert synthetic data to a data frame for ggplot2
data <- data.frame(values = synthetic_data)
data$status <- 0

size <- 1
prob <- 0.1
dead <- rbinom(n_samples, size, prob)

data$status = dead
data

# Combined Histogram and Density Plot
ggplot(data, aes(x = values)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.2, fill = "blue", color = "black", alpha = 0.7) +
  geom_density(fill = "green", alpha = 0.5) +
  labs(title = "Histogram and Density Plot of Synthetic Data", x = "Value", y = "Density") +
  theme_minimal()


###################### REAL DISTRIBUTIONS #####################

x_values <- seq(min(synthetic_data), max(synthetic_data), length.out = 1000)

# Calculate the density for each distribution and combine them
density1 <- dlnorm(x_values, meanlog = mu1, sdlog = sigma1)
density2 <- dlnorm(x_values, meanlog = mu2, sdlog = sigma2)
combined_density <- 0.8 * density1 + 0.2 * density2

# Calculate the CDF for each distribution and combine them
cdf1 <- plnorm(x_values, meanlog = mu1, sdlog = sigma1)
cdf2 <- plnorm(x_values, meanlog = mu2, sdlog = sigma2)
combined_cdf <- 0.8 * cdf1 + 0.2 * cdf2

# Calculate the survival function
survival_function <- 1 - combined_cdf
# Hazard function
hazard_function <- combined_density / survival_function

# Create data frames for plotting
density_data <- data.frame(x = x_values, density = combined_density)
cdf_data <- data.frame(x = x_values, cdf = combined_cdf)
survival_data <- data.frame(x = x_values, survival = survival_function)
hazard_data <- data.frame(x = x_values, hazard = hazard_function)

# Plot Density
ggplot(data, aes(x = values)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.15, fill = "blue", alpha = 0.3) +
  geom_line(data = density_data, aes(x = x, y = density), color = "blue", size = 0.5) +
  labs(title = "Histogram and Theoretical Density of Synthetic Data", x = "Value", y = "Density") +
  theme_minimal()

# Plot CDF
ggplot(cdf_data, aes(x = x, y = cdf)) +
  geom_line(color = "green") +
  labs(title = "CDF of Synthetic Data", x = "Value", y = "CDF") +
  theme_minimal()

# Plot Survival Function
ggplot(survival_data, aes(x = x, y = survival)) +
  geom_line(color = "red") +
  labs(title = "Survival Function of Synthetic Data", x = "Value", y = "Survival Probability") +
  theme_minimal()

# Plot Hazard Function
ggplot(hazard_data, aes(x = x, y = hazard)) +
  geom_line(color = "pink") +
  labs(title = "Hazard Function of Synthetic Data", x = "Value", y = "Hazard") +
  theme_minimal()


################# KEPLEN-MEIER ################### 

# Create a survival object
surv_object <- Surv(time = synthetic_data)

# Fit the Kaplan-Meier estimator
km_fit <- survfit(surv_object ~ 1)

# Extract the Kaplan-Meier estimate with confidence intervals
km_data <- data.frame(
  time = km_fit$time,
  survival = km_fit$surv,
  lower = km_fit$lower,
  upper = km_fit$upper
)

# Plot the Kaplan-Meier estimate with confidence intervals and the theoretical survival function
ggplot() +
  geom_step(data = km_data, aes(x = time, y = survival, color = "Kaplan-Meier"), size = 1) +
  geom_ribbon(data = km_data, aes(x = time, ymin = lower, ymax = upper), fill = "blue", alpha = 0.2) +
  geom_line(data = survival_data, aes(x = x, y = survival, color = "Theoretical"), size = 1, linetype = "dashed") +
  labs(title = "Kaplan-Meier vs True Survival Function", x = "Time", y = "Survival Probability") +
  theme_minimal() +
  scale_color_manual(name = "Type", values = c("Kaplan-Meier" = "blue", "True" = "red")) +
  theme(legend.title = element_blank(), legend.position = "bottom")


################# COX REGRESSION ####################

# Fit the Cox proportional hazards model (since we have no covariates, the model is fitted to the intercept)
cox_fit <- coxph(surv_object ~ 1)

# Estimate the baseline survival function from the Cox model
cox_baseline <- survfit(cox_fit)

# Extract the Cox model baseline survival estimate
cox_data <- data.frame(
  time = cox_baseline$time,
  survival = cox_baseline$surv
)

# Combine Cox model and theoretical survival data for ggplot
cox_data$Type <- "Cox Model"
survival_data$Type <- "Theoretical"

# Plot the Cox model estimate and the theoretical survival function
ggplot() +
  geom_step(data = cox_data, aes(x = time, y = survival, color = "Cox Model"), size = 1) +
  geom_line(data = survival_data, aes(x = x, y = survival, color = "Theoretical"), size = 1, linetype = "dashed") +
  labs(title = "Cox Model Estimate vs. Theoretical Survival Function", x = "Time", y = "Survival Probability") +
  theme_minimal() +
  scale_color_manual(name = "Type", values = c("Cox Model" = "blue", "True" = "red")) +
  theme(legend.title = element_blank(), legend.position = "bottom")


############### PARAMETRIC - WEIBULL ##################

# Fit a Weibull parametric model
weibull_fit <- survreg(surv_object ~ 1, dist = "weibull")

# Extract the scale and shape parameters from the Weibull fit
scale_weibull <- 1 / weibull_fit$scale
shape_weibull <- exp(weibull_fit$coefficients)

# Define the survival function for the Weibull distribution
weibull_survival <- function(t, shape, scale) {
  exp(- (t / scale) ^ shape)
}

# Generate survival estimates from the Weibull model
weibull_data <- data.frame(
  time = x_values,
  survival = weibull_survival(x_values, shape_weibull, scale_weibull),
  Type = "Weibull"
)

# Calculate the theoretical survival function
survival_function <- 1 - combined_cdf
survival_data <- data.frame(
  time = x_values,
  survival = survival_function,
  Type = "Theoretical"
)

# Fit the Kaplan-Meier estimator for comparison
km_fit <- survfit(surv_object ~ 1)
km_data <- data.frame(
  time = km_fit$time,
  survival = km_fit$surv,
  Type = "Kaplan-Meier"
)

# Combine all data into one dataframe with consistent column names
combined_data <- rbind(km_data, weibull_data, survival_data)

# Plot the Weibull model estimate, Kaplan-Meier estimate, and theoretical survival function
ggplot() +
  geom_step(data = combined_data[combined_data$Type == "Kaplan-Meier", ], aes(x = time, y = survival, color = Type), size = 1) +
  geom_line(data = combined_data[combined_data$Type == "Weibull", ], aes(x = time, y = survival, color = Type), size = 1, linetype = "dotted") +
  geom_line(data = combined_data[combined_data$Type == "Theoretical", ], aes(x = time, y = survival, color = Type), size = 1, linetype = "dashed") +
  labs(title = "Weibull Model Estimate vs. Kaplan-Meier and Theoretical Survival Functions", x = "Time", y = "Survival Probability") +
  theme_minimal() +
  scale_color_manual(name = "Type", values = c("Kaplan-Meier" = "blue", "Weibull" = "green", "Theoretical" = "red")) +
  theme(legend.title = element_blank(), legend.position = "bottom")


################ dirichlet processes  #################

# Fit a Dirichlet Process Mixture Model with Weibull kernel
g0Priors <- c(1, 2, 0.5)
alphaPriors <- c(2, 1)
dp <- DirichletProcessWeibull(synthetic_data,
                              g0Priors=g0Priors,
                              alphaPriors=alphaPriors)
dp <- Fit(dp, 1000)  # Fit the model with 1000 iterations

plot(dp)
plot(dp, data_method="hist")

# Extract posterior samples from the fitted models
shape_params <- dp$clusterParameters[[1]]
scale_params <- dp$clusterParameters[[2]]

# Define functions to compute Weibull density and survival
weibull_density <- function(x, shape, scale) {
  (shape / scale) * (x / scale)^(shape - 1) * exp(-(x / scale)^shape)
}

weibull_survival <- function(x, shape, scale) {
  exp(-(x / scale)^shape)
}

# Compute the density and survival estimates for each x value for Data A
density_estimates <- sapply(x_values, function(x) {
  mean(sapply(1:length(shape_params), function(i) {
    weibull_density(x, shape_params[i], scale_params[i])
  }))
})

survival_estimates <- sapply(x_values, function(x) {
  mean(sapply(1:length(shape_params), function(i) {
    weibull_survival(x, shape_params[i], scale_params[i])
  }))
})

# Create data frames for plotting
density_df <- data.frame(x = x_values, y = density_estimates, dataset = "A")
survival_df <- data.frame(x = x_values, y = survival_estimates, dataset = "A")

# Combine data frames for ggplot
density_df <- rbind(density_df)
survival_df <- rbind(survival_df)

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



######################### DPWEIBULL PACKAGE ###########################

install.packages("truncdist")
install.packages("prodlim")
install.packages("RcppArmadillo")
install.packages("/Users/alessandromorosini/Desktop/Thesis/Codes/packages/binaryLogic_0.3.9.tar.gz", repos = NULL, type = "source")
install.packages("/Users/alessandromorosini/Desktop/Thesis/Codes/packages/DPWeibull_1.8.tar.gz", repos = NULL, type = "source")
library(DPWeibull)
library(Rsolnp)


############ LOG-NORMAL MIXTURE GENERATED DATA ###################

MSE.times = seq(0,4,by=0.1) 
msim=1
nsim=200
set.seed(2024)

##### generation of data ##### 
full.norm.simulate = function(n_samples, mu1, mu2, sigma1, sigma2, prob) {    

  # Generate samples from the two distributions
  samples1 <- rlnorm(n = round(0.8 * n_samples), meanlog = mu1, sdlog = sigma1)
  samples2 <- rlnorm(n = round(0.2 * n_samples), meanlog = mu2, sdlog = sigma2)
  
  # Combine the samples
  synthetic_data <- c(samples1, samples2)
  
  # Convert synthetic data to a data frame for ggplot2
  data <- data.frame(time = synthetic_data)

  dead <- rbinom(n_samples, 1, prob)
  data$indicator = dead
  
  return(data)
}

mu1 <- 0
sigma1 <- sqrt(0.25)
mu2 <- 1.2
sigma2 <- sqrt(0.02)
prob <- 0.1

norm_data <- full.norm.simulate(nsim, mu1, mu2, sigma1, sigma2, prob)

full.norm.sim.data=list()
for(i in 1:msim){
  full.norm.sim.data[[i]]=full.norm.simulate(nsim, mu1, mu2, sigma1, sigma2, prob)
}
sum(full.norm.sim.data[[1]]$indicator==1)/nsim

burnin = 5000
iterations = 10000
DP.lst=list() 

for(i in 1:msim){
  DP.lst[[i]]=dpweib(Surv(time=time, event=indicator)~1,
                     data=full.norm.sim.data[[i]],
                     predtime=MSE.times,
                     burnin=burnin, 
                     iteration=iterations)
}

DPM.est.matrix= matrix(NA,msim,length(MSE.times)) 
DPM.MSE = rep(0, msim)

# underlying truth
truth_norms = 1 - (0.8*plnorm(MSE.times, meanlog = mu1, sdlog = sigma1) 
                   + 0.2*plnorm(MSE.times, meanlog = mu2, sdlog = sigma2))

for (i in 1:msim) {
  DPM.est = DP.lst[[i]]$Spred
  DPM.est.matrix[i,]=DPM.est
  DPM.MSE[[i]] = sum((MSE.times[2] - MSE.times[1])*( DPM.est -truth_norms)^2)
}


# empirical 95% confidence intervals
DPM.est.mean=apply(DPM.est.matrix,2,mean)

DPM.est.low=apply(DPM.est.matrix,2,function(x){
  quantile(x,0.025)
})
DPM.est.high=apply(DPM.est.matrix,2,function(x){
  quantile(x,0.975)
})


plot(MSE.times,truth_norms,type='l',xlab='Time',ylab='Survival Probability',main="10% Censoring")
lines(MSE.times,DPM.est.mean,type='l',col='red');lines(MSE.times,DPM.est.low, lty=2,col='red');lines(MSE.times,DPM.est.high, lty=2,col='red')
legend('topright',inset=0.02, legend=c("Underlying Truth", "Dirichlet Process Weibull Mixture Model","95%  Confidence Interval"),
       col=c("black", "red","red"), lty=c(1,1,2), cex=0.5,box.lty=0)



################# WEIBULL MIXTURE GENERATED DATA ###################

MSE.times = seq(0,4,by=0.1) 
msim=1
nsim=200
set.seed(2024)

##### generation of data ##### 
full.weib.mixed.simulate = function(n, coeff, weibpar1,weibpar2,weibpar3,cpar) {    
  count1=0;count2=0;count3=0
  j = 0
  output = NULL
  while (j < n) {                                   # iterate until a sample of size n is obtained
    mod.sel = which(rmultinom(1,1,coeff)==1)
    if (mod.sel == 1){
      count1=count1+1
      t = rweibull(1, shape=weibpar1[1], scale=weibpar1[2])             
    }
    else if (mod.sel == 2) {
      count2=count2+1
      t = rweibull(1, shape=weibpar2[1], scale=weibpar2[2])     
    }
    else{ #mod.sel == 2
      count3=count3+1
      t = rweibull(1, shape=weibpar3[1], scale=weibpar3[2]) 
    }
    cens.t = rexp(1,cpar)                  # generate a forward censoring time (this form of ltr + cens, ensures the censoring time is included in the cohort)
    surv.time = ifelse(t <= cens.t, t, cens.t)  # take the minimum to
    delta = ifelse(t <= cens.t, 1, 0)           # generate status indicator
    output = rbind(output, c(t, delta)) 
    
    j=j+1
  }
  output = data.frame(output)
  
  dimnames(output)[[2]] = c("time","indicator") # give names to output columns for easy access
  
  return(output)
}

coeff=rep(1/3,3)
weibpar1=c(2,2);weibpar2=c(1.5,3);weibpar3=c(3,3.2)
cpar = 0.15

weib_data <- full.weib.mixed.simulate(nsim,coeff =coeff,weibpar1 = weibpar1,weibpar2 = weibpar2,weibpar3=weibpar3,cpar=cpar)

full.weib.mixed.sim.data=list()
for(i in 1:msim){
  full.weib.mixed.sim.data[[i]] =full.weib.mixed.simulate(nsim,coeff =coeff,weibpar1 = weibpar1,weibpar2 = weibpar2,weibpar3=weibpar3,cpar=cpar)
}
sum(full.weib.mixed.sim.data[[1]]$indicator==1)/nsim

burnin = 5000
iterations = 100000

DP.lst=list() 
for(i in 1:msim){
  DP.lst[[i]]=dpweib(Surv(time=time, event=indicator)~1,
                     data=full.weib.mixed.sim.data[[i]],
                     predtime=MSE.times,
                     burnin=burnin, 
                     iteration=iterations)
}

### MSE and 95% CI 

# a matrix to record the estimates by DPM for each simulated dataset
DPM.est.matrix= matrix(NA,msim,length(MSE.times)) 
DPM.MSE = rep(0, msim)

# truth
truth_weib = 1-(coeff[1]*pweibull(MSE.times, weibpar1[1],weibpar1[2])+coeff[2]*pweibull(MSE.times, weibpar2[1],weibpar2[2])+coeff[3]*pweibull(MSE.times, weibpar3[1],weibpar3[2]))


for (i in 1:msim) {
  DPM.est = DP.lst[[i]]$Spred
  DPM.est.matrix[i,]=DPM.est
  DPM.MSE[[i]] = sum((MSE.times[2] - MSE.times[1])*( DPM.est -truth_weib)^2)
}


# empirical 95% confidence intervals
DPM.est.mean=apply(DPM.est.matrix,2,mean)

DPM.est.low=apply(DPM.est.matrix,2,function(x){
  quantile(x,0.025)
})
DPM.est.high=apply(DPM.est.matrix,2,function(x){
  quantile(x,0.975)
})


plot(MSE.times,truth_weib,type='l',xlab='Time',ylab='Survival Probability',main="10% Censoring")
lines(MSE.times,DPM.est.mean,type='l',col='red');lines(MSE.times,DPM.est.low, lty=2,col='red');lines(MSE.times,DPM.est.high, lty=2,col='red')
legend('topright',inset=0.02, legend=c("Underlying Truth", "Dirichlet Process Weibull Mixture Model","95%  Confidence Interval"),
       col=c("black", "red","red"), lty=c(1,1,2), cex=0.5,box.lty=0)


DP.lst[[1]]
