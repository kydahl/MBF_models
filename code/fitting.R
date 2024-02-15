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
library(PhaseType) # to fit structured phase type distributions

# Helper function: Get the number of parameters for a phase-type distribution of a given class and dimension
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

# # Data set from fixed Generalized Coxian distribution
# set_A = matrix(c(
#   -1.5, 0, 0,
#   1.5, -1, 0,
#   0, 1, -0.5), ncol = 3) 
# set_alpha = c(0.9, 0.1, 0) # initial probability vector
# ph = PH(set_A, set_alpha)
# 
# # Get 100 random samples from the distribution
# direct_samps = rPH(100, ph)

# Define initial distributions for fitting ----

# Function: fitting function 
init_func <- function(data, dist_class, dist_dimension, n_runs) {
  # Fits a phase-type distribution of class = dist_class and 
  # dimension = dist_dimension to the data
  
  # Runs the fitting process n_runs times with different initial conditions to
  # try to avoid getting stuck in local minima
  fits = list()
  temp_logLik = double()
  
  for (run_num in 1:n_runs) {
    temp_ph = ph(structure = dist_class, dimension = dist_dimension)
    
    # KD: using capture.output to suppress output from "fit"
    capture.output(
      temp_fit <- fit(temp_ph, y = data, stepsEM = 1500), 
      file = nullfile()
    )
    fits[[run_num]] = temp_fit
    temp_logLik[run_num] = temp_fit@fit$logLik
    temp_param_count = param_count_func(dist_class, dist_dimension)
    temp_AIC = 2 * temp_param_count - 2 * temp_logLik
  }
  
  final_fit = fits[temp_logLik == max(temp_logLik)][[1]]
  
}



dimensions = seq(1, 5)
classes = c("general", "coxian", "hyperexponential", "gcoxian", "gerlang")

data_sets = c("CTMC") #, "direct")

# Fit the data ----

# # Estimate the AIC for each structure and dimension
# AIC_df = data.frame(class = as.character(), dimension = as.integer(), data_set = as.character(),
#                     logLik = as.double(), param_count = as.integer(), AIC = as.double())
# 
# fit_df = tibble(class = as.character(), dimension = as.integer(), data_set = as.character(),
#                 logLik = as.double(), AIC = as.double(), A = list(), alpha = list())
# 
# for (dimension in dimensions) {
#   for (class in classes) {
#     for (data_set in data_sets) {
#       if (data_set == "CTMC") {
#         out_times = simulated_outs
#       } else {
#         out_times = direct_samps
#       }
#       print(paste0("Fitting PH distribution:"))
#       print(paste0("Structure: ", class))
#       print(paste0("Dimension: ", dimension))
#       print(paste0("Data set: ", data_set ))
#       
#       # Get fitted alpha and A of given class and dimension
#       # Use 10 different initial conditions to ensure we're not stuck in 
#       # a local minimum
#       temp_fit <-init_func(out_times, class, dimension, 10)
#       
#       temp_logLik = temp_fit@fit$logLik
#       
#       print(paste0("Likelihood = ", temp_logLik ))
#       temp_param_count = param_count_func(class, dimension)
#       temp_AIC = 2 * temp_param_count - 2 * temp_logLik
#       
#       AIC_df = add_row(AIC_df, 
#                        class = class, dimension = dimension, data_set = data_set,
#                        logLik = temp_logLik, param_count = temp_param_count, 
#                        AIC = temp_AIC)
#       
#       fit_df = add_row(fit_df, class = class, dimension = dimension,
#                        data_set = data_set, logLik = temp_logLik,
#                        AIC = temp_AIC, A = list(temp_fit@pars$S), 
#                        alpha = list(temp_fit@pars$alpha))
#     }
#   }
# }

# Set up data to analyze
dimensions = seq(1, 5)
# The three structures we consider for now are:
# 1) Empirical: the general phase distribution, which best fits the data
# 2) Phenomenological: Coxian - equivalent to our "disruption" model or an approximation of the mechanistic model
# 3) Mechanistic: assumes a certain form given below
classes = c("general", "coxian") #, "mechanistic") # "hyperexponential", "gcoxian", "gerlang")

total_sample_size = length(simulated_data$Out_time)
num_samples = 100
sample_size = 100

full_df <- expand_grid(dimension = dimensions, class = classes)

out_df = tibble(class = as.character(), dimension = as.integer(), data_set = as.integer(),
                # logLik = as.double(), AIC = as.double(), 
                A = list(), alpha = list())


for (i in 1:num_samples) {
  
  data_sample = sample(simulated_data$out_days, sample_size)
  
  for (j in 1:dim(full_df)[1]) {
    dist_class = full_df$class[j]
    dist_dimension = full_df$dimension[j]
    print(paste0("Fitting PH distribution:"))
    print(paste0("Structure: ", dist_class))
    print(paste0("Dimension: ", dist_dimension))
    print(paste0("Data set: ", i, " out of ", num_samples))
    
    temp_fit <- init_func(data_sample, dist_class, dist_dimension, 10)
    
    # print(paste0("Likelihood = ", temp_logLik ))
    # temp_param_count = param_count_func(class, dimension)
    # temp_AIC = 2 * temp_param_count - 2 * temp_logLik
    
    out_df = add_row(out_df, class = dist_class, dimension = dist_dimension,
                     data_set = i, 
                     # logLik = temp_logLik, AIC = temp_AIC, 
                     A = list(temp_fit@pars$S), alpha = list(temp_fit@pars$alpha))
  }
  
}


write_csv(out_df, "data/fit_dists.csv")

# Deal with the mechanistic fit separately
# We'll make use of the PhaseType package which makes it easier to generate
# random phase-type subintensity matrices of a given structure

mech_mat_test = matrix(c(
  0,"R","R",0,
  "F",0,0,0,
  "F",0,0,0,
  0,"F","F",0
),4)

# Define the structure of the phase-type generator
mech_mat = matrix(c(
  0,  "L1","P1", "G1", 0,
  "Q",0,   "P2", "G2", 0,
  0,  "L2",0,    0,    0,
  0,  0,   "P3", 0,    0,
  0,  0,   0,    0,    0
),5)
# Gamma priors for shape hyperparameters of model parameters
nu <- list(
  "Q" = 1,
  "L1" = 1, "L2" = 1,
  "P1" = 1, "P2" = 1, "P3" = 1,
  "G1" = 1, "G2" = 1
  )
# Gamma priors for reciprocal scale hyperparamters of model parameters
zeta <- c(
  "Q" = 1,
  "L1" = 1, "L2" = 1,
  "P1" = 1, "P2" = 1, "P3" = 1,
  "G1" = 1, "G2" = 1
)

# alpha is fixed 
alpha = c(1,0,0,0)

# perform N MCMC iterations to fit data structured phase-type distribution to data "x"
N = 1
res <- phtMCMC2(x, mech_mat, alpha, nu, zeta, N)

# Next steps:
# For each sampled data set (indexed by "data_set" in the data frame), we 
# calculate R0 (using "repnum_func" from get_outputs.R). We then will have 
# obtained distributions of R0 across the data sets. 
# Then 


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

