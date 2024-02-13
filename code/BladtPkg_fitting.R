# Fitting phase type distributions to data
# Fit (structured) phase type distributions to simulated mosquito blood-feeding 
# data using the matrixdist package
# Initialized: February 2024

# Load libraries and define helper functions ----
library(matrixdist)
library(tidyverse)
library(matlib) # for a bunch of matrix operations
library(cowplot)
library(PhaseTypeR) # to create and simulate phase type distributions

# Function: Get the number of parameters for a phase-type distribution of a given class and dimension
param_count_func <- function(class, dimension) {
  p = dimension
  out <- case_when(
    class == "general" ~ p^2 + p,
    class == "coxian" ~ max(1, 2*p - 2),
    class == "hyperexponential" ~ 2*p,
    class == "gcoxian" ~ 3*p - 2,
    class == "gerlang" ~ p
  )
} 

# Get data ----
# Simulated data set
simulated_data = read_csv("data/sample_data.csv")

simulated_outs = sample(simulated_data$Out_time/1440, 100)

# Data set from fixed Generalized Coxian distribution
set_A = matrix(c(
  -1.5, 0, 0,
  1.5, -1, 0,
  0, 1, -0.5), ncol = 3) 
set_alpha = c(0.9, 0.1, 0) # initial probability vector
ph = PH(set_A, set_alpha)

# Get 100 random samples from the distribution
direct_samps = rPH(100, ph)

# Define initial distributions for fitting ----

dimensions = seq(1, 20)
classes = c("general", "coxian", "hyperexponential", "gcoxian")#, "gerlang")

data_sets = c("DTMC", "direct")

# Fit the data ----

# Estimate the AIC for each structure and dimension
# AIC_df = data.frame(class = as.character(), dimension = as.integer(), data_set = as.character(),
#                     logLik = as.double(), param_count = as.integer(), AIC = as.double())

for (dimension in dimensions) {
  for (class in classes) {
    for (data_set in data_sets) {
      if (data_set == "DTMC") {out_times = simulated_outs
      } else {out_times = direct_samps}
    print(paste0("Fitting PH distribution:"))
    print(paste0("Structure: ", class))
    print(paste0("Dimension: ", dimension))
    print(paste0("Data set: ", data_set ))
    temp_ph = ph(structure = class, dimension = dimension)
    
    temp_fit = fit(temp_ph, y = out_times, stepsEM = 1500)
    
    temp_logLik = temp_fit@fit$logLik
    temp_param_count = param_count_func(class, dimension)
    temp_AIC = 2 * temp_param_count - 2 * temp_logLik
    
    AIC_df = add_row(AIC_df, class = class, dimension = dimension, data_set = data_set,
                     logLik = temp_logLik, param_count = temp_param_count, AIC = temp_AIC)
    }
  }
}

AIC_compare_plot <- AIC_df %>% 
  mutate(actual_dim = ifelse(data_set == "direct", 3, 12)) %>% 
  group_by(data_set) %>% 
  mutate(min_AIC = min(AIC)) %>% 
  mutate(delta_AIC = abs(AIC - min_AIC)) %>% 
  # mutate(low_AIC = min_AIC - 2) %>% 
  # mutate(high_AIC = min_AIC + 2) %>% 
  ggplot() +
  geom_point(aes(x = dimension, y = delta_AIC, color = class)) +
  geom_line(aes(x = dimension, y = delta_AIC, color = class)) +
  geom_ribbon(aes(x = dimension, ymin = 0, ymax = 2), alpha = 0.3) +
  # geom_hline(aes(yintercept =min_AIC)) +
  geom_vline(aes(xintercept = actual_dim)) +
  scale_x_continuous(breaks = dimensions, expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle("Comparing AIC") +
  facet_wrap( ~ data_set, nrow = 2, scales = "free") +
  theme_cowplot() +
  coord_cartesian(xlim = range(dimensions), y = c(0, 100))

logLik_compare_plot <- AIC_df %>% 
  mutate(actual_dim = ifelse(data_set == "direct", 3, 12)) %>% 
  ggplot() +
  geom_point(aes(x = dimension, y = logLik, color = class)) +
  geom_line(aes(x = dimension, y = logLik, color = class)) +
  geom_vline(aes(xintercept = actual_dim)) +
  scale_x_continuous(breaks = dimensions) +
  ggtitle("Comparing logLikelihood") +
  facet_wrap( ~ data_set, nrow = 2, scales = "free") +
  theme_cowplot()

# Can we get confidence intervals around each of the parameters?
# Not from this method since the EM algorithm will always give the same results
# Could try boot-strapping or methods like that

# To get initial probability vector alpha
alpha = coef(ph_fit)$alpha

# To get subintensity matrix A
A = coef(ph_fit)$S

# To get logLikelihood
# ph_fit@fit$logLik


# Visualizations ----

# QQ plots to assess fit visually

