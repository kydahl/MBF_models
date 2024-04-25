############# Fitting phase type distributions to data #########################
# Fit (structured) phase type distributions to simulated mosquito blood-feeding 
# data using the matrixdist, PhaseTypeR, and PhaseType packages
# Initialized: February 2024

# Outline of script: 
# 0) Load libraries, necessary functions, and data
# 1) Set up data sets
# 2) Fit matrices to the data
# 3) Transform matrix entries to model parameters
# 4) Create visualizations of parameter and model fits
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
library(GGally)

# 1) Set up data sets with increasing levels of data ----------------------

# TODO: 
# [] Determine appropriate levels to use
# [] Figure out how to set up submatrices to fit to lower levels of data
# []

#### Get data ----
# Simulated data set
simulated_data = read_csv("data/direct_samples.csv") 
# simulated_data = read_csv("data/noQ_data_continuous.csv")
# total_sample_size = length(simulated_data$Out_time)
sample_size = 100

test_data = simulated_data$direct_samps[1:sample_size]

# test_data = simulated_data$Out_time[1:100]


# 2) Fit matrices to the data ---------------------------------------------

# TODO: 
# [] Define each distribution type and its structure in comments
# [] Determine reasonable uninformed priors to use for each type
# [] Write function that creates a structured matrix of the given type
# [] 

# TODO: 
# [] Set up general fitting function that deals with possible errors from phtMCMC2
# [] Keep track of likelihood and AIC of each fit
# [] Decide: how to fix dimensions

## Empirical ----

# Run through a range of PH orders and select best fit based on AIC
max_order = 6
out_mats = lapply(1:max_order, function(x) matrix(NA, x, x))
AIC_table = tibble(order = as.integer(), logLik = as.double(), AIC = as.integer())

for (PH_order in 1:max_order) {
  print(paste0("PH Order: ", PH_order, " out of ", max_order))
  num_transient_states = PH_order
  
  # alpha is always set to (1,0,...,0)
  alpha = rep(0, num_transient_states); alpha[1] = 1
  
  priors <- uninf_priors_func(mean(test_data) / num_transient_states, var(test_data), num_transient_states)
  # Gamma prior: reciprocal scale hyperparameters (one per matrix row)
  nu <- priors$nu #rep(1/mean(test_data), num_states)
  zeta <- priors$zeta[1:num_transient_states]
  
  num_samples = 1000
  res <- phtMCMC(test_data, 
                 states = num_transient_states, 
                 beta = alpha, 
                 nu = nu, zeta = zeta, 
                 n = num_samples,
                 mhit = 1000)
  
  # Reshape output to a manageable matrix
  medians <- as.data.frame(res$samples) %>% 
    summarise(across(everything(), median)) %>% 
    as.double()
  out_mat = reshape_phtMCMC(medians, "Empirical")
  out_mats[[PH_order]] = out_mat
  
  # # Visually inspect the parameter distributions
  # plot.phtMCMC(res) # posterior densities should be roughly unimodal
  # corrplot(cor(res$samples)) # there should be no strong correlations among parameters
  # Check that largest real part of all eigenvalues is negative
  if (any(Re(eigen(out_mat)$values)>0)) {print("Eigenvalue problem")}
  if (length(out_mat) == 1) {if (-out_mat<0) {print("out-rates problem")}}
  if (length(out_mat) != 1) {if (any(100 * inv(-100 * out_mat)<0)) {print("out-rates problem")}}
  
  
  # Calculate log likelihood
  logLik = PHlikelihood(test_data, out_mat)
  num_parms = PH_order^2
  AIC = 2 * num_parms - 2 * logLik
  AIC_table = add_row(AIC_table, order = PH_order, logLik = logLik, AIC = AIC)
  
}
best_order = filter(AIC_table, AIC == min(AIC))$order
best_PH = out_mats[[best_order]]


## Phenomenological [NYI] ----

# NYI

## Mechanistic  ----

# Define the structure of the phase-type generator
mech_mat = matrix(c(
  0,    "P1", "G1", 0,
  "L1",    0,    0, 0,
  0,    "P2",    0, 0,
  "L2", "P3", "G2", 0
), 4)
# Gamma priors for shape hyperparameters of model parameters
nu_ref = 0.5
nu <- c(
  "L1" = nu_ref, "L2" = nu_ref,
  "P1" = nu_ref, "P2" = nu_ref, "P3" = nu_ref,
  "G1" = nu_ref, "G2" = nu_ref
)
# Gamma priors for reciprocal scale hyperparameters of model parameters
## Use informed priors from known underlying parameters
# Landing
pL = 0.7
lL = 0.1 # 10 minutes
# Probing
pP = 0.8
lP = 0.2 # 5 minutes
# Ingesting
pG = 0.9 #0.75
lG = 1 # 1 minutes
# Fleeing
f = 0.66
zeta <- c(
  "L1" = (pL * lL)^(-1)/nu_ref, "L2" = (f*(1-pL)*lL)^(-1)/nu_ref,
  "P1" = ((1 - f) * (1 - pP) * lP)^(-1)/nu_ref, "P2" = (pP * lP)^(-1)/nu_ref, "P3" = (f * (1 - pP) * lP)^(-1)/nu_ref,
  "G1" = ((1 - f) * (1 - pG) * lG)^(-1)/nu_ref, "G2" = (pG * lG + f * (1 - pG) * lG)^(-1)/nu_ref
)

# alpha is always set to (1,0,0)
alpha = rep(0, 3); alpha[1] = 1

num_samples = 1000

# Get the fitted matrix entry samples
mech_res <- phtMCMC2(test_data, 
                     TT = mech_mat, 
                     beta = alpha, 
                     nu = nu, zeta = zeta, 
                     n = num_samples,
                     method = "MHRS",
                     mhit = 100000)



# 3) Transform matrix entries to parameters -------------------------------

## Empirical ----

## Mechanistic ----

# Make list of actual parameters used to simulate process
actual_mech_parms <- tibble(
  p_L = pL, lambda_L = lL, p_P = pP, lambda_P = lP, p_G = pG, lambda_G = lG,
  f = f) %>% 
  pivot_longer(everything())


# Make matrix of actual values used to simulate process
actual_mech_mat = matrix(c(
  (-1 + (1 - f) * (1- pL)) * lL, pL * lL, 0,
  (1 - f) * (1 - pP) * lP, -lP, pP * lP,
  (1 - f) * (1 - pG) * lG, 0, -lG), ncol = 3, byrow = TRUE) %>% 
  melt() %>% rename(row_index = Var1, col_index = Var2)

# Reshape samples into matrices and mechanistic parameters
mech_mat_df <- tibble(iterate = as.integer(), row_index = as.integer(), 
                      col_index = as.integer(), value = as.double(), 
                      logLik = as.double())

mech_params_df <- tibble(iterate = as.integer(), lambda_L = as.double(), 
                         p_L = as.double(), lambda_P = as.double(), 
                         p_P = as.double(), lambda_G = as.double(), 
                         p_G = as.double(), f = as.double(), logLik = as.double())

for (i in 1:dim(mech_res$samples)[1]) {
  PH_mat = reshape_phtMCMC(mech_res$samples[i,], "Mechanistic")
  logLik = PHlikelihood(test_data, PH_mat)
  
  # Collect entries of PH matrices
  for (row in 1:3) {
    for (col in 1:3) {
      # Replace numeric zeros with actual zeros in mechanistic matrix
      value = case_when(
        row == 1 & col == 3 ~ 0,
        row == 3 & col == 2 ~ 0,
        TRUE ~ PH_mat[row, col]
      )
      mech_mat_df <- add_row(mech_mat_df, 
                             iterate = i, row_index = row, 
                             col_index = col, value = PH_mat[row, col],
                             logLik = logLik)
    }
  }
  
  if (!any(apply(abs(PH_mat[, 1:nrow(PH_mat)]) <= sqrt(.Machine$double.eps), 
                 1, all))) {
    # Collect mechanistic parameters
    lP_fit = -PH_mat[2, 2]
    pP_fit = PH_mat[2,3] / lP_fit
    temp0 = (PH_mat[2,1] / ((1 - pP_fit) * lP_fit))
    f_fit = 1 - (PH_mat[2,1] / ((1 - pP_fit) * lP_fit))
    
    temp1 = - PH_mat[1,1] - PH_mat[1,2]
    temp2 = temp1 / f_fit
    lL_fit = temp2 + PH_mat[1,2]
    pL_fit = PH_mat[1,2] / lL_fit
    
    lG_fit = -PH_mat[3, 3]
    pG_fit = max(0, 1 - (PH_mat[3,1] / ((1 - f_fit) * lG_fit)))
    
    mech_params_df <- add_row(mech_params_df, 
                              iterate = i, lambda_L = lL_fit, p_L = pL_fit, 
                              lambda_P = lP_fit, p_P = pP_fit, lambda_G = lG_fit, 
                              p_G = pG_fit, f = f_fit)
  }
  
}

# Add loglikelihood back in to params dataframe
mech_params_df <- right_join(mech_params_df, select(mech_mat_df, iterate, logLik), 
                             by = "iterate")

# # Visually inspect the parameter distributions
plot.phtMCMC(mech_res) # posterior densities should be roughly unimodal
corrplot(cor(mech_res$samples)) # there should be no strong correlations among parameters
# Check that largest real part of all eigenvalues is negative
if (any(Re(eigen(out_mat)$values)>0)) {print("Eigenvalue problem")}
if (length(out_mat) == 1) {if (-out_mat<0) {print("out-rates problem")}}
if (length(out_mat) != 1) {if (any(inv(-out_mat)<0)) {print("out-rates problem")}}


# 4) Create visualizations of model fits ----------------------------------

# Empirical distribution matrix ----

#### Profile likelihood analyses ----
matrix_index <- expand.grid(row = 1:best_order, col = 1:best_order)

PLA_df <- tibble(row_index = as.integer(), col_index = as.integer(), 
                 status = as.character(),
                 value = as.double(), likelihood = as.double())
for (iter in 1:nrow(matrix_index)) {
  print(paste0("Iteration #", iter, " / ", nrow(matrix_index)))
  # Isolate one parameter (matrix element)
  row_index = matrix_index$row[iter]
  col_index = matrix_index$col[iter]
  el_val = best_PH[row_index, col_index]
  
  # Determine the maximum amount of variation allowable in the element
  # Must ensure that row sums remain negative
  # Then set range to explore for the element (+/- variation, cutoff if necessary)
  var_mag = 10
  if (row_index == col_index) {
    # Diagonal entries must remain negative
    min_val = -1440 # Set so that durations don't exceed a day
    max_val = -sum(best_PH[row_index, -col_index])
    min_perturb = max(min_val, var_mag * el_val)
    max_perturb = min(max_val, - var_mag * el_val)
  } else {
    # All other entries must be positive
    min_val = 0
    max_val = -sum(best_PH[row_index, -col_index])
    min_perturb = max(min_val, -var_mag * el_val)
    max_perturb = min(max_val, var_mag * el_val)
  }
  
  # Set up sequence of values for perturbations
  perturb_resolution = 100
  perturb_seq = c(el_val, seq(min_perturb, max_perturb, length.out = perturb_resolution))
  
  # Calculate (expected) likelihood of the data for perturbed values of the matrix element
  for (j in 1:perturb_resolution) {
    status = ifelse(j == 1, "base", "perturbed")
    temp_mat = best_PH
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
  facet_wrap(vars(row_index, col_index), nrow = best_order, ncol = best_order, 
             labeller = label_both, scales = "free_x")
PLA_plot  

ggsave("figures/full_emp_PLA_plot.pdf", PLA_plot, width = 16, height = 9)

## RowSums Parameters Profile likelihood analysis ##
rowsum_PLA_df <- tibble(row_index = as.integer(), status = as.character(),
                        value = as.double(), likelihood = as.double())
for (iter in 1:best_order) {
  # Isolate one parameter (matrix element)
  row_index = iter
  col_index = row_index
  el_val = best_PH[row_index, col_index]
  
  # Determine the maximum amount of variation allowable in the element
  # Must ensure that row sums remain negative
  # Then set range to explore for the element (+/- variation, cutoff if necessary)
  var_mag = 100
  # row_sum = -sum(best_PH[row_index, ])
  
  min_perturb = -3
  max_perturb = 3
  
  # Set up sequence of values for perturbations
  perturb_resolution = 100
  perturb_seq = c(1, 10^seq(min_perturb, max_perturb, length.out = perturb_resolution))
  
  # Calculate (expected) likelihood of the data for perturbed values of the matrix element
  for (j in 1:perturb_resolution) {
    status = ifelse(j == 1, "base", "perturbed")
    temp_mat = best_PH
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
  facet_wrap(vars(row_index), nrow = 1, ncol = best_order, 
             labeller = label_both, scales = "free_x")

ggsave("figures/full_emp_rowsum_PLA_plot.pdf", rowsum_PLA_plot, width = 16, height = 9)


# Mechanistic distribution ----
mech_params_plot_df <- mech_params_df %>%
  pivot_longer(cols = -c(iterate, logLik))

#### Parameter distributions plot ----
shift_legend(mech_params_plot_df %>% 
               ggplot(aes(x = value)) +
               geom_histogram() +
               # Add mean line
               geom_vline(data = group_by(mech_params_plot_df, name) %>% 
                            summarise(value = mean(value)),
                          aes(xintercept = value, color = "blue"),
                          linewidth = 1) +
               # Add median line
               geom_vline(data = group_by(mech_params_plot_df, name) %>% 
                            summarise(value = median(value)),
                          aes(xintercept = value, color = "orange"),
                          linewidth = 1) +
               # Add actual value line
               geom_vline(data = actual_mech_parms, aes(xintercept = value, color = "green"),
                          linewidth = 1) +
               facet_wrap(~name, scales = "free") +
               scale_color_manual("Type",
                                  breaks = c("green"="green", "blue"="blue", "orange"="orange"), 
                                  values = c("green"="green", "blue"="blue", "orange"="orange"), 
                                  labels = c("actual", "mean", "median")) +
               theme_cowplot(12) +
               theme(legend.direction = "horizontal")
)

#### Parameter correlations plot ----
mech_params_df %>% 
  select(-c(iterate)) %>% 
  cor() %>% 
  corrplot(., type = "lower")

### Parameter scatter plots ----
mech_params_df %>% 
  filter(lambda_L < 300) %>% 
  select(-c(iterate)) %>% 
  ggpairs(aes()) +
  theme_cowplot(11)

#### Parameter PLA plot ----
# # Use the sample with the highest loglikelihood
# chosen_sample = filter(mech_mat_df, logLik == max(logLik)) %>%
#   select(-c(iterate, logLik))

# Use the mean values across samples with the highest loglikelihood
chosen_sample = mech_params_df %>%
  right_join(select(mech_mat_df, iterate, logLik), by = "iterate") %>% 
  filter(logLik == max(logLik)) %>% 
  unique()

# Transform sample into a matrix
out_parms = select(chosen_sample, -iterate)

# Get actual values and their likelihood
actual_logLik = actual_mech_mat %>% 
  pivot_wider(names_from = col_index, values_from = value) %>% 
  select(-c(row_index)) %>% 
  as.matrix() %>% 
  PHlikelihood(test_data, .)
actual_mech_mat_wlik <- mutate(actual_mech_mat, likelihood = actual_logLik)

PLA_df <- tibble(parameter = as.character(), 
                 status = as.character(),
                 value = as.double(), likelihood = as.double())

for (iter in 1:ncol(out_parms)) {
  print(paste0("Iteration #", iter))
  parameter = colnames(out_parms)[iter]
  # Isolate one parameter (matrix element)
  row_index = matrix_index$row[iter]
  col_index = matrix_index$col[iter]
  el_val = out_mat[row_index, col_index]
  
  # Determine the maximum amount of variation allowable in the element
  # Must ensure that row sums remain negative
  # Then set range to explore for the element (+/- variation, cutoff if necessary)
  var_mag = 2
  min_perturb = -abs(var_mag * el_val)
  max_perturb = abs(var_mag * el_val)
  # if (row_index == col_index) {
  #   # Diagonal entries must remain negative
  #   min_val = -1440 # Set so that durations don't exceed a day
  #   max_val = -sum(out_mat[row_index, -col_index])
  #   min_perturb = max(min_val, var_mag * el_val)
  #   max_perturb = min(max_val, - var_mag * el_val)
  # } else {
  #   # All other entries must be positive
  #   min_val = 0
  #   max_val = -sum(out_mat[row_index, -col_index])
  #   min_perturb = max(min_val, -var_mag * el_val)
  #   max_perturb = min(max_val, var_mag * el_val)
  # }
  
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
likelihood_cutoff = log(0.95 * exp(max(PLA_df$likelihood, na.rm = TRUE)))

max_likelihood_set = PLA_df %>% 
  select(-status) %>% 
  group_by(row_index, col_index) %>% 
  # filter(value = value)
  filter(likelihood == max(likelihood, na.rm = TRUE)) %>% 
  unique()

PLA_plot <- PLA_df %>% 
  filter(-likelihood < -1.005 * likelihood_cutoff) %>% 
  ggplot(aes(x = value, y = likelihood)) +
  geom_line() +
  # Add dot showing baseline values
  geom_point(data = filter(PLA_df, status == "base"), color = "blue") +
  # Add dot showing maximum likelihood values
  geom_point(data = max_likelihood_set, color = "orange") +
  # Add dot showing actual underyling value
  geom_point(data = actual_mech_mat_wlik, color = "green") +
  # geom_vline(data = actual_mech_mat_wlik, aes(xintercept = value), color = "green") +
  geom_hline(yintercept = likelihood_cutoff, color = "red") +
  facet_wrap(vars(row_index, col_index), nrow = 3, ncol = 3, 
             labeller = label_both, scales = "free_x")
PLA_plot  

ggsave("figures/full_mech_PLA_plot.pdf", PLA_plot, width = 16, height = 9)

# Summary statistics for mechanistic parameter posterior distributions
mech_params_summary <- group_by(mech_params_plot_df, name) %>% 
  summarise(mean = mean(value),
            median = median(value))

#### Matrix element distribution plots ----
mech_mat_plot_df <- mech_mat_df %>% 
  mutate(row_label = paste0("row = ", row_index),
         col_label = paste0("column = ", col_index))

mech_mat_plot_df %>% ggplot(aes(x = value)) +
  geom_histogram() +
  # Actual values
  geom_vline(actual_mech_mat %>% 
               mutate(row_label = paste0("row = ", row_index),
                      col_label = paste0("column = ", col_index)), 
             mapping = aes(xintercept = value, color = "green")) +
  # Add mean line
  geom_vline(data = group_by(mech_mat_plot_df, row_label, col_label) %>% 
               summarise(value = mean(value)),
             aes(xintercept = value, color = "blue"),
             linewidth = 1) +
  # Add median line
  geom_vline(data = group_by(mech_mat_plot_df, row_label, col_label) %>% 
               summarise(value = median(value)),
             aes(xintercept = value, color = "orange"),
             linewidth = 1) +
  facet_wrap(row_label ~ col_label, scales = "free") +
  scale_color_manual("Type",
                     breaks = c("green"="green", "blue"="blue", "orange"="orange"), 
                     values = c("green"="green", "blue"="blue", "orange"="orange"), 
                     labels = c("actual", "mean", "median")) +
  theme_cowplot(12) 

#### Matrix elements correlations plot ----
mech_mat_df %>% 
  # select(-c(iterate)) %>% 
  unite(entry, row_index:col_index) %>% 
  pivot_wider(names_from = entry, values_from = value) %>% 
  select(-c(`3_2`, `1_3`)) %>% 
  cor() %>% 
  corrplot(., type = "lower")

#### Matrix PLA plot ----
# # Use the sample with the highest loglikelihood
# chosen_sample = filter(mech_mat_df, logLik == max(logLik)) %>%
#   select(-c(iterate, logLik))

# Use the mean values across samples with the highest loglikelihood
chosen_sample = mech_mat_df %>%
  group_by(row_index, col_index) %>%
  summarise(value = mean(value)) %>%
  ungroup()


# Transform sample into a matrix
out_mat = chosen_sample %>% 
  pivot_wider(names_from = col_index, values_from = value) %>% 
  select(-c(row_index)) %>% 
  as.matrix()

# Get actual values and their likelihood
actual_logLik = actual_mech_mat %>% 
  pivot_wider(names_from = col_index, values_from = value) %>% 
  select(-c(row_index)) %>% 
  as.matrix() %>% 
  PHlikelihood(test_data, .)
actual_mech_mat_wlik <- mutate(actual_mech_mat, likelihood = actual_logLik)

matrix_index <- expand.grid(row = 1:3, col = 1:3)

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
  var_mag = 2
  min_perturb = -abs(var_mag * el_val)
  max_perturb = abs(var_mag * el_val)
  # if (row_index == col_index) {
  #   # Diagonal entries must remain negative
  #   min_val = -1440 # Set so that durations don't exceed a day
  #   max_val = -sum(out_mat[row_index, -col_index])
  #   min_perturb = max(min_val, var_mag * el_val)
  #   max_perturb = min(max_val, - var_mag * el_val)
  # } else {
  #   # All other entries must be positive
  #   min_val = 0
  #   max_val = -sum(out_mat[row_index, -col_index])
  #   min_perturb = max(min_val, -var_mag * el_val)
  #   max_perturb = min(max_val, var_mag * el_val)
  # }
  
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
likelihood_cutoff = log(0.95 * exp(max(PLA_df$likelihood, na.rm = TRUE)))

max_likelihood_set = PLA_df %>% 
  select(-status) %>% 
  group_by(row_index, col_index) %>% 
  # filter(value = value)
  filter(likelihood == max(likelihood, na.rm = TRUE)) %>% 
  unique()

PLA_plot <- PLA_df %>% 
  filter(-likelihood < -1.005 * likelihood_cutoff) %>% 
  ggplot(aes(x = value, y = likelihood)) +
  geom_line() +
  # Add dot showing baseline values
  geom_point(data = filter(PLA_df, status == "base"), color = "blue") +
  # Add dot showing maximum likelihood values
  geom_point(data = max_likelihood_set, color = "orange") +
  # Add dot showing actual underyling value
  geom_point(data = actual_mech_mat_wlik, color = "green") +
  # geom_vline(data = actual_mech_mat_wlik, aes(xintercept = value), color = "green") +
  geom_hline(yintercept = likelihood_cutoff, color = "red") +
  facet_wrap(vars(row_index, col_index), nrow = 3, ncol = 3, 
             labeller = label_both, scales = "free_x")
PLA_plot  

ggsave("figures/full_mech_PLA_plot.pdf", PLA_plot, width = 16, height = 9)

#### RowSums Parameters Profile likelihood analysis ----
rowsum_PLA_df <- tibble(row_index = as.integer(), status = as.character(),
                        value = as.double(), likelihood = as.double())
for (iter in 1:best_order) {
  # Isolate one parameter (matrix element)
  row_index = iter
  col_index = row_index
  el_val = out_mat[row_index, col_index]
  
  # Determine the maximum amount of variation allowable in the element
  # Must ensure that row sums remain negative
  # Then set range to explore for the element (+/- variation, cutoff if necessary)
  var_mag = 10
  # row_sum = -sum(out_mat[row_index, ])
  
  min_perturb = -3
  max_perturb = 3
  
  # Set up sequence of values for perturbations
  perturb_resolution = 100
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
  facet_wrap(vars(row_index), nrow = 1, ncol = 3, 
             labeller = label_both, scales = "free_x")

ggsave("figures/full_mech_rowsum_PLA_plot.pdf", rowsum_PLA_plot, width = 16, height = 9)


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
