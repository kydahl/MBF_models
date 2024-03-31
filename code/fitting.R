############# Fitting phase type distributions to data #########################
# Fit (structured) phase type distributions to simulated mosquito blood-feeding 
# data using the matrixdist, PhaseTypeR, and PhaseType packages
# Initialized: February 2024

# Outline of script: 
# 0) Load libraries, necessary functions, and data
# 1) Set up data sets
# 2) Set up structured matrices for fitting
# 3) Fit matrices to the data
# 4) Transform matrix entries to model parameters
# 5) Create visualizations of parameter and model fits
# ---------------------------------------------------------------------------- #

# 0) Load libraries, necessary functions, and data ------------------------
library(matrixdist)
library(tidyverse)
library(matlib) # for a bunch of matrix operations
library(cowplot)
library(PhaseTypeR) # to create and simulate phase type distributions
library(PhaseType) # to fit structured phase type distributions
library(mapfit) # alternative way to fit phase type distributions
library(reshape2)

# Function: Get the number of parameters for a phase-type distribution of
#           a given class and dimension
param_count_func <- function(class, dimension) {
  p = dimension
  out <- case_when(
    class == "general" ~ p^2 + p,
    class == "coxian" ~ max(1, 2*p - 2),
    class == "hyperexponential" ~ 2*p,
    class == "gcoxian" ~ 3*p - 2,
    class == "gerlang" ~ p,
    class == "mechanistic" ~ 9 * p / 4 # 9 parameters per 4x4 motif. p must be a multiple of 4
  )
} 

# Function: calculate the likelihood of phase type parameters given  data
likelihood_func <- function(parameters, data) {}

# Function: put together likelihood and AIC for each model and data set
model_selection_func <- function(model_class, dimension, parameters, data) {
  #       temp_fit <-init_func(out_times, class, dimension, 10)
  #       
  #       temp_logLik = temp_fit@fit$logLik
  #       
  #       temp_param_count = param_count_func(class, dimension)
  #       temp_AIC = 2 * temp_param_count - 2 * temp_logLik
}

# OLD Function: initialize parameters for a structured phase-type distribution then fit to data 
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

# Function: Plot phtMCMC outputs nicely (from PhaseType package)
plot.phtMCMC <- function(x, ...) {
  # Get chain info
  StartEndThin <- attr(x$samples, "mcpar")
  
  # Trace plots
  print(
    ggplot(melt(as.data.frame(x$samples[seq(StartEndThin[1], StartEndThin[2], StartEndThin[3]),]), id.vars = NULL)) +
      geom_line(aes_string(x = "1:length(value)", y = "value")) +
      geom_smooth(aes_string(x = "1:length(value)", y = "value"), method = "glm") +
      #geom_hline(aes(yintercept = value), data = truth, colour = "red", linetype = "dashed") +
      facet_wrap(~variable, scales = "free") +
      #theme_grey(base_family = "serif", base_size = 11) +
      ggtitle("Parameter Traces") + #, plot.title = theme_text(size = 14, face = "bold", family = "serif")) +
      xlab("Iteration") + ylab("Parameter Value")
  )
  
  # Marginal posterior densities
  print(
    ggplot() +
      geom_density(aes_string(x = "value"), melt(as.data.frame(x$samples[seq(StartEndThin[1], StartEndThin[2], StartEndThin[3]),]), id.vars = NULL)) +
      #geom_vline(aes(xintercept = value), data = truth, colour = "red") +
      facet_wrap(~variable, scales = "free") +
      #theme_grey(base_family = "serif", base_size = 11) +
      ggtitle("Marginal Posterior Densities") + #, plot.title = theme_text(size = 14, face = "bold", family = "serif")) +
      xlab("Parameter Value") + ylab("Density")
  )
}


# 1) Set up data sets with increasing levels of data ----------------------

# TODO: 
# [] Determine appropriate levels to use
# [] Figure out how to set up submatrices to fit to lower levels of data
# []

#### Get data ----
# Simulated data set
simulated_data = read_csv("data/sample_data_continuous.csv")
total_sample_size = length(simulated_data$Out_time)
sample_size = 100



# 2) Set up structured matrices for fitting -------------------------------

# TODO: 
# [] Define each distribution type and its structure in comments
# [] Determine reasonable uninformed priors to use for each type
# [] Write function that creates a structured matrix of the given type
# [] 

### Empirical distribution matrix ----

# alpha is fixed 
test_data = simulated_data %>% 
  mutate(not_Q = time_G+time_L+time_P+time_G) %>%  select(not_Q)
test_data = test_data$not_Q[1:100]
num_states = 4
alpha = rep(0, num_states); alpha[num_states] = 1
nu <- rep(1, num_states^2)#rep(1/mean(test_data), num_states)
# Gamma prior: reciprocal scale hyperparameters (one per matrix row)
zeta <- rep(1, num_states)
res <- phtMCMC(test_data, states = num_states, alpha, nu, zeta, n = 100, method = "DCS")
plot.phtMCMC(res)

wsample <- rweibull(n=100, shape=2, scale=1)

## PH fitting for general PH
(result_gen <- phfit.point(ph=mapfit::ph(4), x=test_data))
# Canonical form 1: mixture of Erlangs with alpha = 0,0,...,1
(result_cf1 <- phfit.point(ph=mapfit::cf1(4), x=test_data))
(result_erlang <- phfit.point(ph=mapfit::herlang(5), x=test_data))

# Using matrixdist #
# Make a ph distribution
# matrixdist::ph()
# Fit a ph distribution
# matrixdist::fit(ph(), data)

dimensions = seq(1, 5)
classes = c("general", "coxian", "hyperexponential", "gcoxian", "gerlang")

# The three structures we consider for now are:
# 1) Empirical: the general phase distribution, which best fits the data 
# 2) Phenomenological: Coxian - equivalent to our "disruption" model or an approximation of the mechanistic model
# 3) Mechanistic: assumes a certain form given below
classes = c("general", "coxian") #, "mechanistic") # "hyperexponential", "gcoxian", "gerlang")

write_rds(out_df, "data/fit_dists.rds")

# MCMC fitting

# Deal with the mechanistic fit separately
# We'll make use of the PhaseType package which makes it easier to generate
# random phase-type subintensity matrices of a given structure


# Define the structure of the phase-type generator
mech_mat = matrix(c(
  0,  "L1","P1", "G1", 0,
  "Q",0,   "P2", "G2", 0,
  0,  "L2",0,    0,    0,
  0,  0,   "P3", 0,    0,
  0,  0,   0,    "G3", 0
),5)
# Gamma priors for shape hyperparameters of model parameters
nu_ref = 0.01
nu <- list(
  "Q" = nu_ref,
  "L1" = nu_ref, "L2" = nu_ref,
  "P1" = nu_ref, "P2" = nu_ref, "P3" = nu_ref,
  "G1" = nu_ref, "G2" = nu_ref, "G3" = nu_ref
)
# Gamma priors for reciprocal scale hyperparameters of model parameters
zeta <- c(
  "Q" = 1440/nu_ref,
  "L1" = (0.66 * 0.25 * 0.1)^(-1)/nu_ref, "L2" = (0.75 * 0.1)^(-1)/nu_ref,
  "P1" = (0.66 * 0.25 * 0.2)^(-1)/nu_ref, "P2" = (0.34 * 0.25 * 0.2)^(-1)/nu_ref, "P3" = (0.75 * 0.2)^(-1)/nu_ref,
  "G1" = (0.66 * 0.75 * 1)^(-1)/nu_ref, "G2" = (0.34 * 0.75 * 1)^(-1)/nu_ref, "G3" = (0.25 * 1)^(-1)/nu_ref
  # "Q" = 2880,
  # "L1" = 80, "L2" = 40,
  # "P1" = 40, "P2" = 40, "P3" = 20,
  # "G1" = 8, "G2" = 8, "G3" = 8/3
)

# 3) Fit matrices to the data ---------------------------------------------

# TODO: 
# [] Set up general fitting function that deals with possible errors from phtMCMC2
# [] Keep track of likelihood and AIC of each fit
# [] Decide: how to fix dimensions

# # perform N MCMC iterations to fit data structured phase-type distribution to data "x"
# data_sample = sample(simulated_data$Out_time, sample_size)
# res <- phtMCMC2(data_sample, mech_mat, alpha, nu, zeta, n = 100, method = "DCS")

out_df = tibble(class = as.character(), dimension = as.integer(), data_set = as.integer(),
                # logLik = as.double(), AIC = as.double(), 
                A = list(), alpha = list())
for (i in 1:num_samples) {
  
  data_sample = sample(simulated_data$Out_time, 1000)
  print(paste0("Fitting PH distribution:"))
  print(paste0("Data set: ", i, " out of ", num_samples))
  capture.output(
    res <- phtMCMC2(simulated_data$Out_time, mech_mat, alpha, nu, zeta, n = 1000, 
                    method = "DCS", silent = FALSE), 
    file = nullfile()
  )
  
  # Put matrix on scale of days instead of minutes
  DayToMin = 1440
  out_mat <- mech_mat_func(res$samples, median) * DayToMin
  
  # print(paste0("Likelihood = ", temp_logLik ))
  # temp_param_count = param_count_func(class, dimension)
  # temp_AIC = 2 * temp_param_count - 2 * temp_logLik
  
  out_df = add_row(out_df, class = "Mechanistic", dimension = 4,
                   data_set = i,
                   # logLik = temp_logLik, AIC = temp_AIC,
                   A = list(out_mat), alpha = list(alpha))
}

write_rds(out_df, "data/mech_dists.rds")

# 
# Can we get confidence intervals around each of the parameters?
# Not from this method since the EM algorithm will always give the same results
# Could try boot-strapping or methods like that

# To get initial probability vector alpha
alpha = coef(ph_fit)$alpha

# To get subintensity matrix A
A = coef(ph_fit)$S

# To get logLikelihood
# ph_fit@fit$logLik


# 4) Transform matrix entries to parameters -------------------------------

# Helper function: Put output of phtMCMC2 in matrix form following our bespoke mechanistic matrix form
mech_mat_func <- function(in_mat, stat_func) {
  # in_mat should be of the form of the output of res$samples[i,]
  
  stat_mat <- in_mat %>% 
    as.data.frame() %>% 
    reshape2::melt() %>% 
    group_by(variable) %>% 
    summarise(mean = stat_func(value)) %>% 
    pivot_wider(names_from = variable, values_from = mean)
  
  G_mat = matrix(c(
    0,  stat_mat$L1, stat_mat$P1, stat_mat$G1, 0,
    stat_mat$Q,0,   stat_mat$P2, stat_mat$G2, 0,
    0,  stat_mat$L2,0,    0,    0,
    0,  0,   stat_mat$P3, 0,    0,
    0,  0,   0,    stat_mat$G3,    0
  ),5)
  row_sums = rowSums(G_mat)
  
  temp_mat = G_mat - diag(row_sums)
  
  # remove absorbing state row and column
  mat_dim = dim(temp_mat)[1]
  
  A_mat = temp_mat[1:(mat_dim-1), 1:(mat_dim-1)]
  
  out_mat = t(A_mat)
}



# 5) Create visualizations of model fits ----------------------------------

# QQ plots to assess fit visually


# Figure XXX: Posterior distributions of fitted parameters

plot_parm_dists_func <- function(param_set) {
  dists_plot <- parms %>% 
    ggplot(aes(x = value)) + 
    geom_histogram() +
    geom_density() +
    facet_wrap(~name, scales = "free") +
    theme_cowplot()
  
}

# Parameter set to play with
sample_size = 100
res <- phtMCMC2(sample(simulated_data$Out_time, sample_size),
                mech_mat, alpha, nu, zeta, n = 1000)

# Function: Empirical model set up for input to fitting algorithm
emp_mat_func <- function(dimension) {
  dim <- 5
  dirpi <- rep(0, dim)
  dirpi[1] <-  1
  # Gamma prior: shape hyperparameters (one per matrix element, columnwise)
  nu <- rep(1, dim^2)
  # Gamma prior: reciprocal scale hyperparameters (one per matrix row)
  zeta <- rep(100, dim)
  
}

# Function: Phenomenological model set up for input to fitting algorithm 

# Define dimension of model to fit

res <- phtMCMC(sample(simulated_data$Out_time, sample_size) / 1440, 
               dim, dirpi, nu, zeta, 1000, mhit=10)

# Function: get mechanistic parameters from transient rate matrix
mech_parm_func <- function(param_set) {
  out_df <- param_set %>% 
    mutate(
      lG = G1 + G2 + G3, 
      pG = G3 / lG,
      f = G1 / (G1 + G2),
      lP = P1 + P2 + P3,
      pP = P3 / lP,
      lL = (L1 + f * L2) / f,
      pL = (f * L2) / (L1 + f * L2),
      lQ = Q
    ) %>% 
    select(lG:lQ)
}


# Actual distributions data was sampled from, in form for plot
# KD: This was done incorrectly because the original samples weren't all gamma distributed like this
# all_dist_data <- data.frame(n = 1:sample_size)
# for (variable in actual_dists$name) {
#   temp_nu = filter(actual_dists, name == variable)$nu
#   temp_zeta = filter(actual_dists, name == variable)$zeta
#   dist_data <-
#     rgamma(
#       n = sample_size,
#       shape = temp_nu, 
#       scale = temp_zeta
#     )
#   temp_df <- data.frame(variable = dist_data / 1440) # put on same time scale as parameters
#   names(temp_df)[names(temp_df) == "variable"] <- variable
#   
#   all_dist_data <- cbind(all_dist_data, temp_df)
#   
# }
# all_dist_plot <- pivot_longer(select(all_dist_data, -n), cols = everything())


# Put parameters in workable form
parms <- mech_parm_func(as.data.frame(res$samples)) %>%
  pivot_longer(everything())

# parms <- as.data.frame(res$samples) %>% 
#   pivot_longer(everything())

correct_vals_quote_unquote <- data.frame(
  lG = 1, 
  pG = 0.25,
  f = 0.66,
  lP = 0.2, 
  pP = 0.75, 
  lL = 0.1,
  pL = 0.75,
  lQ = 1/1440
) %>% pivot_longer(everything())

dists_plot <- parms %>% 
  ggplot(aes(x = value)) + 
  geom_histogram() +
  geom_density() +
  geom_vline(data = correct_vals_quote_unquote, aes(xintercept = value),
             color = "blue") +
  # geom_density(data = all_dist_plot, color = "darkgray", linetype = "dashed") +
  facet_wrap(~name, scales = "free")


dists_plot

# AIC_compare_plot <- AIC_df %>% 
#   mutate(actual_dim = ifelse(data_set == "direct", 3, 12)) %>% 
#   group_by(data_set) %>% 
#   mutate(min_AIC = min(AIC)) %>% 
#   mutate(delta_AIC = abs(AIC - min_AIC)) %>% 
#   # mutate(low_AIC = min_AIC - 2) %>% 
#   # mutate(high_AIC = min_AIC + 2) %>% 
#   ggplot() +
#   geom_point(aes(x = dimension, y = delta_AIC, color = class)) +
#   geom_line(aes(x = dimension, y = delta_AIC, color = class)) +
#   geom_ribbon(aes(x = dimension, ymin = 0, ymax = 2), alpha = 0.3) +
#   # geom_hline(aes(yintercept =min_AIC)) +
#   geom_vline(aes(xintercept = actual_dim)) +
#   scale_x_continuous(breaks = dimensions, expand = c(0,0)) +
#   scale_y_continuous(expand = c(0,0)) +
#   ggtitle("Comparing AIC") +
#   facet_wrap( ~ data_set, nrow = 2, scales = "free") +
#   theme_cowplot() +
#   coord_cartesian(xlim = range(dimensions), y = c(0, 100))
# 
# logLik_compare_plot <- AIC_df %>% 
#   mutate(actual_dim = ifelse(data_set == "direct", 3, 12)) %>% 
#   ggplot() +
#   geom_point(aes(x = dimension, y = logLik, color = class)) +
#   geom_line(aes(x = dimension, y = logLik, color = class)) +
#   geom_vline(aes(xintercept = actual_dim)) +
#   scale_x_continuous(breaks = dimensions) +
#   ggtitle("Comparing logLikelihood") +
#   facet_wrap( ~ data_set, nrow = 2, scales = "free") +
#   theme_cowplot()
