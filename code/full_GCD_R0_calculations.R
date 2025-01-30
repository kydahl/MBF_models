# Full GCD-R0 comparison

# Pseudo-code ----

# Four model types:
# 1. Exponential
# 2. Empirical
# 3. Phenomenological
# 4. Mechanistic

# A range of GCD values
# These should be on the order of 

# For each model type,
#   i) Calculate parameters from specified GCD values
#      - See specific instructions for mechanistic model
#  ii) Calculate R0
# iii) Calculate a scaled version of R0 so that 
#      min(R0_scaled) = 0, max(R0_scaled) = 1
#  iv) Plots of R0 vs GCD, R0 vs 1/GCD, R0_scaled vs GCD, R0_scaled vs 1/GCD
#   v) Sensitivity analysis of R0 wrt GCD parameter(s)

# Load libraries ----
library(tidyverse)
library(matlib)

# General functions ----

# Calculate R0
R0_calc <- function(A, v_alpha, LambdaB, LambdaH, BetaB, BetaH) {
  # Calculates R0 as a function of the phase-type distribution parameters
  # and specific transmission formulation
  
}

# Calculate theta, biting process duration
theta_calc <- function(A, v_alpha) {
  A_dim = dim(A)[1]
  green_A = inv(-A)
  v_one = matrix(rep(1, A_dim), 1)
  theta = v_one %*% green_A %*% v_alpha
}

# Calculate GCD, gonotrophic cycle duration
GCD_calc <- function(A, v_alpha) {
  # Calculate gonotrophic cycle duration as a function of phase-type distribution parameters
  GCD = theta_calc(A, v_alpha) + (1/gammaV) + (1/gammaR)
}


# Fixed Parameter values ----
# Oviposition and resting
gammaV = 1 / (2 * 1440) #1/(5 * 1440) # exit rate from oviposition to resting, including bloodmeal digestion and site search (5 days)
gammaR = 1/(2 * 1440) # exit rate from resting to return to blood-feeding (2 days)
# Transmission parameters
betaH = betaB = 1
eta = 1/(6 * 1440) # 6 days for infection to develop in vector
mu = 1/(20 * 1440) # 20 day lifespan for vector
gammaH = 1/(7 * 1440) # rate of recovery in hosts (7 days)
muH = 1/(365.25 * 65 * 1440) # host mortality rate (65 years)
KH = 1E8 # host population density

KJ = 3 * 1E8
rhoJ = 1 / (12 * 1440) # 12 day larval development period
muJ = 1 / (20 * 1440) # 62.5% probability of larval survival (20 / (20 + 12))
varPhi = 300 / 1440 # on average 3 eggs per female per day

# GCD range ----
# Set range of GCD values to consider
# Note: GCD = biting duration + oviposition duration + resting duration
#       Only theta = biting duration will be varied
resolution = 1001

theta_min = 8 * 60 # minimum biting duration of 8 hours
theta_max = 20 * 24 * 60 # maximum biting duration of 20 days
theta_vec = seq(theta_min, theta_max, resolution)

GCD_vec = theta_vec + (1/gammaR) + (1/gammaV)


# Exponential model ----
b_vec = 1/theta_vec

## Calculate R0...

# Empirical model ----

## For each theta value...

## Generate 100 samples from a lognormal distribution with mean theta, variance 1...
samples <- tibble(theta = theta_vec) %>% 
  expand_grid(tibble(n = 1:resolution)) %>% 
  rowwise() %>% 
  mutate(sample = rlnorm(1, meanlog = log(theta), sdlog = 1))

## Specify initial guess A and v_alpha...
## Should be dimension 3. Use the same for all theta values (for now)
set.seed(24060)
v_alpha_in = rgamma(3, 1)
v_alpha_in = v_alpha_in/sum(v_alpha_in)

# Initialize data frame
out_df = as_tibble(matrix(
  nrow = 0, 
  ncol = 3, 
  dimnames = list(NULL, c("n", "alpha", "A"))
))


## Run EM algorithm to obtain best-fit A and v_alpha...
library(mapfit)
out_df = data.frame()
for (theta_val in theta_vec) {
  sample_in = filter(samples, theta == theta_val)$sample
  set.seed(24060) # to ensure same initial state
  out = phfit.point(ph = ph(3), 
                    x = sample_in)
  
  out_df = rbind(out_df,
                 tibble(theta = theta_val, alpha = list(out$alpha), A = c(out$Q)))
}
saveRDS(out_df, "data/EM_fits_to_lognormal.rds")

out_df = readRDS("data/EM_fits_to_lognormal.rds")


## [Sanity check] Calculate theta_calc(A, v_alpha) and compare to theta

# For some chosen values of theta, plot pdfs of lognormal and fitted PH
library(PhaseTypeR)
library(gridExtra)
library(cowplot)
library(latex2exp)

test_theta_vals = c(min(samples$theta), median(samples$theta), max(samples$theta))

test_df = tibble(test_theta = test_theta_vals) %>% 
  rowwise() %>% 
  mutate(ph_model = list(phfit.point(ph = ph(3), 
                                     x = filter(samples, theta == test_theta)$sample)$model),
         row_num = row_number()) %>% 
  group_by(row_num)

x_maxs = c(2750, 75000, 200000)
scale_factors = c(1, 10, 10)
legend_list = c("none", "none", "legend")

plot_list = list()

for (i in 1:length(test_df)) {
  test_theta_val = test_df$test_theta[i]
  x_max = x_maxs[i]
  temp_plot = ggplot(data.frame(x_vals = seq(0, x_max, length.out = 1000)), aes(x = x_vals)) +
    stat_function(fun = dphase, args = list(ph = test_df$ph_model[[i]]),
                  aes(color = "PH"),
                  lwd = 1
    ) +
    stat_function(fun = dlnorm, args = list(mean = log(test_theta_val), sd = 1), 
                  aes(color = "lognormal"),
                  lwd = 1
    ) +
    scale_x_continuous(name = "Days", breaks = seq(0, x_max, by = scale_factors[i]*1440), labels = seq(0, x_max / (scale_factors[i]*1440), by = 1)) + 
    # scale_x_continuous(name = "Days", trans = ~ . / 1440) +
    scale_y_continuous(name = "") +
    guides(
      y = "none"
      ) +
    scale_color_discrete(name = "Distribution") +
    ggtitle(TeX(paste0("$\\theta \\approx ", round(test_theta_val/1440, 2)," days"))) +
    theme_minimal(18) +
    theme(
      legend.position = "right",  # Legend positioned to the right
      legend.key.size = unit(1, "cm")  # Adjust legend key size
    )
  
  plotgrob = ggplotGrob(temp_plot)
  legend_index <- which(sapply(plotgrob$grobs, function(x) x$name) == "guide-box")
  legend <- plotgrob$grobs[[legend_index]]
  
  # panel_index <- which(plotgrob$layout$name == "panel")
  # plot_panel <- gtable::gtable_filter(plotgrob, "panel")
  
  plotgrob$grobs[[legend_index]] <- nullGrob()
  
  plot_list[[i]] = temp_plot + theme(legend.position = "none")
}

# plot_only = 

all_plots <- grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], legend, nrow = 1, widths = (3, 1))

## Determine LambdaH, LambdaB, BetaH, BetaB from A...

## Calculate R0...

# Phenomenological model ----

## For each theta value...

## Assemble A and v_alpha following structure in paper...

## Determine LambdaH, LambdaB, BetaH, BetaB from A...

## Calculate R0....

# Mechanistic model ----

## !!! Cannibalize most of this from old code

## Define ranges for each of the parameters...
## Might consider flighty vs persistent

## Set up parameter grid(s)...

## Calculate GCD over parameter grid(s)...

## Determine LambdaH, LambdaB, BetaH, BetaB over parameter grid(s)...

## Calculate R0 over parameter grid(s)...