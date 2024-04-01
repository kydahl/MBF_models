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
library(Matrix)
library(corrplot)

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

# Function: set up uninformed priors for each parameter of the PH distribution
uninf_priors_func <- function(mean, variance, dimension) {
  shape_k = mean^2 / variance
  scale_theta = variance / mean
  
  nu <- rep(shape_k, dimension^2) #rep(1/mean(test_data), num_states)
  # Gamma prior: reciprocal scale hyperparameters (one per matrix row)
  zeta <- rep(scale_theta, dimension)
  
  out <- data.frame(nu = nu, zeta = zeta)
  return(out)
}

# Function: reshape phtMCMC/2 output to matrix form
reshape_phtMCMC <- function(res) {
  
  # Get median values from the parameter posterior distributions
  medians <- as.data.frame(res$samples) %>% 
    summarise(across(everything(), mean)) %>% 
    as.double()
  
  dimension <- sqrt(length(medians))
  med_seq = seq(1,dimension^2, by = dimension)
  
  out_rates <- medians[med_seq]
  off_diags <- medians[-med_seq]
  positions <- expand.grid(col = seq(1, dimension), row = seq(1, dimension)) %>% 
    filter(row != col)
  
  PH_mat <- Matrix(data = NA, nrow = dimension, ncol = dimension)
  
  # Fill in the off-diagonal entries
  for (index in 1:dim(positions)[1]) {
    row = positions$row[index]
    col = positions$col[index]
    PH_mat[row, col] = off_diags[index]
  }
  
  # Fill in the diagonal entries
  
  
  PH_mat <- t(matrix(as.numeric(PH_mat), nrow = dimension, ncol = dimension, byrow = TRUE))
  
  diags <- out_rates + rowSums(PH_mat, na.rm = TRUE)
  PH_mat[cbind(seq(1:dimension), seq(1:dimension))] <- -diags
  
  return(PH_mat)
  
}

# Function: Calculate expected likelihood of PH distribution parameters to data
PHlikelihood <- function(data, PH_mat) {
  dimension = dim(PH_mat)[1]
  
  alpha = matrix(rep(0, dimension), nrow = 1); alpha[1] = 1
  
  out_rates = -PH_mat %*% matrix(rep(1, dimension), nrow = dimension, ncol = 1)
  
  
  # N = number of observations
  # y_1,...,N = data
  # pi = alpha
  # T = the sub-intensity matrix, out_mat
  # t = -T*e, the row sums of out_mat
  Lik = 1
  for (i in 1:length(test_data)) {
    Lik = alpha %*% expm(PH_mat * data[i]) %*% out_rates
  }
  logLik = log(Lik)
  
  return(logLik[1])
}

# 1) Set up data sets with increasing levels of data ----------------------

# TODO: 
# [] Determine appropriate levels to use
# [] Figure out how to set up submatrices to fit to lower levels of data
# []

#### Get data ----
# Simulated data set
simulated_data = read_csv("data/noQ_data_continuous.csv") #"data/sample_data_continuous.csv"
total_sample_size = length(simulated_data$Out_time)
sample_size = 100

test_data = simulated_data$Out_time[1:100]

# test_data = simulated_data %>% 
#   mutate(not_Q = time_G+time_L+time_P+time_G) %>%  select(not_Q)
# test_data = test_data$not_Q[1:100]

# 2) Set up structured matrices for fitting -------------------------------

# TODO: 
# [] Define each distribution type and its structure in comments
# [] Determine reasonable uninformed priors to use for each type
# [] Write function that creates a structured matrix of the given type
# [] 


### Empirical distribution matrix ----

# alpha is always set to (1,0,...,0)
num_transient_states = 4
alpha = rep(0, num_transient_states); alpha[1] = 1

priors <- uninf_priors_func(mean(test_data) / num_transient_states, var(test_data), num_transient_states)
# Gamma prior: reciprocal scale hyperparameters (one per matrix row)
nu <- priors$nu #rep(1/mean(test_data), num_states)
zeta <- priors$zeta[1:num_transient_states]

num_samples = 100
res <- phtMCMC(test_data, 
               states = num_transient_states, 
               beta = alpha, 
               nu = nu, zeta = zeta, 
               n = num_samples,
               mhit = 100)

# Reshape output to a manageable matrix
out_mat = reshape_phtMCMC(res)

# Visually inspect the parameter distributions
plot.phtMCMC(res) # posterior densities should be roughly unimodal
corrplot(cor(res$samples)) # there should be no strong correlations among parameters
# Check that largest real part of all eigenvalues is negative
if (any(Re(eigen(out_mat)$values)>0)) {print("Eigenvalue problem")}
if (any(inv(-out_mat)<0)) {print("out-rates problem")}

# PH Parameters Profile likelihood analysis ----
matrix_index <- expand.grid(row = 1:num_transient_states, col = 1:num_transient_states)

PLA_df <- tibble(row_index = as.integer(), col_index = as.integer(), 
                 status = as.character(),
                 value = as.double(), likelihood = as.double())
for (iter in 1:nrow(matrix_index)) {
  print(paste0("Iteration #", iter))
  # Isolate one parameter (matrix element)
  row_index = matrix_index$row[iter]
  col_index = matrix_index$col[iter]
  el_val = out_mat[row_index, col_index]
  
  # Determine the maximum amount of variation allowable in the element
  # Must ensure that row sums remain negative
  # Then set range to explore for the element (+/- variation, cutoff if necessary)
  var_mag = 10
  if (row_index == col_index) {
    # Diagonal entries must remain negative
    min_val = -1440 # Set so that durations don't exceed a day
    max_val = -sum(out_mat[row_index, -col_index])
    min_perturb = max(min_val, var_mag * el_val)
    max_perturb = min(max_val, - var_mag * el_val)
  } else {
    # All other entries must be positive
    min_val = 0
    max_val = -sum(out_mat[row_index, -col_index])
    min_perturb = max(min_val, -var_mag * el_val)
    max_perturb = min(max_val, var_mag * el_val)
  }
  
  # Set up sequence of values for perturbations
  perturb_resolution = 100
  perturb_seq = c(el_val, seq(min_perturb, max_perturb, length.out = perturb_resolution))
  
  # Calculate (expected) likelihood of the data for perturbed values of the matrix element
  for (j in 1:perturb_resolution) {
    status = ifelse(j == 1, "base", "perturbed")
    temp_mat = out_mat
    temp_mat[row_index, col_index] = perturb_seq[j]
    
    # Calculate likelihood
    logLik = PHlikelihood(test_data, temp_mat)
    
    # Put together data
    
    PLA_df <- add_row(PLA_df, status = status,
                      row_index = row_index, col_index = col_index,
                      value = perturb_seq[j], likelihood = logLik)
    
  }
  
}

# Visually inspect the likelihood profiles
likelihood_cutoff = log(0.95 * exp(max(PLA_df$likelihood)))

PLA_plot <- PLA_df %>% ggplot(aes(x = value, y = likelihood)) +
  geom_line() +
  geom_point(data = filter(PLA_df, status == "base"), color = "blue") +
  geom_hline(yintercept = likelihood_cutoff, color = "red") +
  facet_wrap(vars(row_index, col_index), nrow = 4, ncol = 4, 
             labeller = label_both, scales = "free_x")
PLA_plot  

ggsave("figures/PLA_plot.pdf", PLA_plot, width = 16, height = 9)

# RowSums Parameters Profile likelihood analysis ----
rowsum_PLA_df <- tibble(row_index = as.integer(), status = as.character(),
                 value = as.double(), likelihood = as.double())
for (iter in 1:nrow(out_mat)) {
  # Isolate one parameter (matrix element)
  row_index = iter
  col_index = row_index
  el_val = out_mat[row_index, col_index]
  
  # Determine the maximum amount of variation allowable in the element
  # Must ensure that row sums remain negative
  # Then set range to explore for the element (+/- variation, cutoff if necessary)
  var_mag = 100
  # row_sum = -sum(out_mat[row_index, ])
  
  min_perturb = -3
  max_perturb = 3
  
  # Set up sequence of values for perturbations
  perturb_resolution = 1000
  perturb_seq = c(1, 10^seq(min_perturb, max_perturb, length.out = perturb_resolution))
  
  # Calculate (expected) likelihood of the data for perturbed values of the matrix element
  for (j in 1:perturb_resolution) {
    status = ifelse(j == 1, "base", "perturbed")
    temp_mat = out_mat
    temp_mat[row_index, ] = perturb_seq[j] * temp_mat[row_index, ]
    
    # Calculate likelihood
    logLik = PHlikelihood(test_data, temp_mat)
    
    # Put together data
    
    rowsum_PLA_df <- add_row(rowsum_PLA_df,
                      row_index = row_index, status = status,
                      value = perturb_seq[j], likelihood = logLik)
    
  }
  
}

# Visually inspect the likelihood profiles
likelihood_cutoff = log(0.95 * exp(max(rowsum_PLA_df$likelihood)))

rowsum_PLA_plot <- rowsum_PLA_df %>% ggplot(aes(x = value, y = likelihood)) +
  geom_line() +
  geom_point(data = filter(rowsum_PLA_df, status == "base"), color = "blue") +
  geom_hline(yintercept = likelihood_cutoff, color = "red") +
  scale_x_log10() +
  facet_wrap(vars(row_index), nrow = 4, ncol = 4, 
             labeller = label_both, scales = "free_x")

ggsave("figures/rowsum_PLA_plot.pdf", rowsum_PLA_plot, width = 16, height = 9)

### Phenomenological distribution matrix ----

### Mechanistic distribution matrix ----

# Define the structure of the phase-type generator
mech_mat = matrix(c(
  0,   "P2", "G2", 0,
  "L2",0,    0,    0,
  0,   "P3", 0,    0,
  0,   0,    "G3", 0
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
