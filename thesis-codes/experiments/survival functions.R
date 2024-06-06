# Load necessary libraries
install.packages("cowplot")

library(ggplot2)
library(dplyr)
library(gridExtra)
library(cowplot)
# Define the data using Weibull distribution
shape <- 2
scale <- 4
time <- seq(0, 10, by=0.01)
density <- dweibull(time, shape=shape, scale=scale)
cdf <- pweibull(time, shape=shape, scale=scale)
survival <- 1 - cdf

data <- data.frame(time, density, cdf, survival)

# Create the density plot without the legend
density_plot <- ggplot(data, aes(x=time)) +
  geom_area(data = subset(data, time <= 4), aes(y=density, fill="F(t)"), alpha=0.5) +
  geom_area(data = subset(data, time > 4), aes(y=density, fill="S(t)"), alpha=0.5) +
  geom_line(aes(y=density, color="f(t)"), size=1) +
  geom_vline(xintercept=4, linetype="dashed") +
  annotate("text", x=4.2, y=0.2, label="T", size=4) +
  labs(title="",
       x="Time",
       y="Density",
       color=NULL,
       fill=NULL) +
  scale_fill_manual(values=c("F(t)"="#E6E6E6", "S(t)"="darkgrey")) +
  scale_color_manual(values=c("f(t)"="black")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_y_continuous(limits = c(0, 1)) +
  guides(fill="none", color="none")

# Create the survival, CDF, and density plot without the legend
sc_plot <- ggplot(data, aes(x=time)) +
  geom_line(aes(y=density, color="f(t)"), size=1) +
  geom_line(aes(y=survival, color="S(t)"), size=1) +
  geom_line(aes(y=cdf, color="F(t)"), size=1) +
  geom_vline(xintercept=4, linetype="dashed") +
  annotate("text", x=4.2, y=0.2, label="T", size=4) +
  labs(title="",
       x="Time",
       y="Probability",
       color=NULL,
       fill=NULL) +
  scale_color_manual(values=c("f(t)"="black", "S(t)"="#E6E6E6", "F(t)"="darkgrey")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_y_continuous(limits = c(0, 1)) +
  guides(color=guide_legend(override.aes = list(linetype = c(1, 1, 1), 
                                                shape = c(NA, NA, NA),
                                                size = 1),
                            title=NULL))

# Extract the legend from one of the plots
legend <- get_legend(
  ggplot(data, aes(x=time)) +
    geom_line(aes(y=density, color="f(t)"), size=1) +
    geom_line(aes(y=survival, color="S(t)"), size=1) +
    geom_line(aes(y=cdf, color="F(t)"), size=1) +
    scale_color_manual(values=c("f(t)"="black", "S(t)"="#E6E6E6", "F(t)"="darkgrey")) +
    theme_minimal() +
    theme(legend.position = "right") +
    guides(color=guide_legend(override.aes = list(linetype = c(1, 1, 1), 
                                                  shape = c(NA, NA, NA),
                                                  size = 1),
                              title=NULL))
)

# Arrange the plots and the shared legend
plot_grid(density_plot, sc_plot, legend, ncol = 3, rel_widths = c(1, 1, 0.3))
