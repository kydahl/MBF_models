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

# Fixed Parameter values ----
# Oviposition and resting
gammaV = 1 / (5 * 1440) # 1/(5 * 1440) # exit rate from oviposition to resting, including bloodmeal digestion and site search (5 days)
gammaR = 1/ (2 * 1440) # exit rate from resting to return to blood-feeding (2 days)
# Transmission parameters
betaH = betaB = 1
eta = 1/(7 * 1440) # 6 days for infection to develop in vector
mu = 1/(21 * 1440) # 20 day lifespan for vector
gammaH = 1/(2 * 1440) # rate of recovery in hosts (7 days)
muH = 1/(365.25 * 65 * 1440) # host mortality rate (65 years)
KH = 1E3 # host population density

KJ = 0.75 * KH #0.25 * KH
rhoJ = 1 / (12 * 1440) # 12 day larval development period
muJ = 1 / (20 * 1440) # 62.5% probability of larval survival (20 / (20 + 12))
varPhi = 3 / 1440 # on average 3 eggs per female per day

### Standard case functions: ----
get_standard_stable_pop <- function(b) {
  
  check = (mu / (1440 *varPhi * b)) * ((rhoJ + muJ)/ rhoJ)
  
  if (check > 1) {
    J_star = 0
    B_star = 0
  } else {
    J_star = KJ * (1 - check)
    B_star = rhoJ * J_star / mu
  }
  
  return(
    c(
      J_star = J_star,
      B_star = B_star)
  )
}

standard_R0_calc <- function(b) {
  B_star = as.list(get_standard_stable_pop(b))$B_star
  standard_beta = betaH / sqrt(60) #2 * betaH / 400 #369.5324 * betaH / 4000# / (1000*sqrt(3))
  host_infectious_period = (1 / (muH + gammaH))
  standard_beta = 1 * 1440 / sqrt(host_infectious_period * (eta / (mu + eta)) * (1 / mu) * KJ / KH)
  # 
  # b_vec = seq(1 / (20 * 1440), 1 / (2 * 1440), length.out = 10)
  
  R0 = b * sqrt(host_infectious_period * standard_beta * (eta / (mu + eta)) * (1 / mu) * standard_beta * B_star / KH)
}


### Exponential case functions: ----
# Get LambdaH and LambdaB and BetaH and BetaB
exp_EpiTerms <- function(A_matrix, v_alpha) {
  KB = sum(get_stable_pop(A_matrix, v_alpha)$v_B_star)
  A_dim = dim(A_matrix)[1]
  v_ones = matrix(rep(1, A_dim), ncol = 1)
  out_rates = -as.matrix(A_matrix) %*% v_ones
  
  LambdaB = diag(as.vector(out_rates), nrow = length(out_rates)) %>%
    cbind(matrix(rep(0, A_dim))) %>% # need to add in zeros because oviposition class is added in
    rbind(matrix(rep(0, A_dim+1), nrow = 1)) %>%
    list()
  LambdaH = (LambdaB[[1]] * KB / KH) %>% list()
  
  BetaH = (betaH * diag(A_dim)) %>%
    cbind(matrix(rep(0, A_dim))) %>% # need to add in zeros because oviposition class is added in
    rbind(matrix(rep(0, A_dim+1), nrow = 1)) %>%
    list()
  BetaB = (betaB * BetaH[[1]] / betaH) %>% list()
  
  return(
    c(
      LambdaB = LambdaB,
      LambdaH = LambdaH,
      BetaB = BetaB,
      BetaH = BetaH)
  )
}


### Empirical case functions: ----
get_emp_A <- function(theta) {
  ## Generate 100 samples from a lognormal distribution with mean theta, variance 1...
  sample_num = 100
  
  # Correct values of the lognormal parameters to match R's implementation
  lognorm_sigma <- sqrt(log(1 + (1 / theta^2)))  # Corrected sigma
  lognorm_mu <- log(theta) - 0.5 * lognorm_sigma^2  # Corrected mu
  
  samples <- tibble(theta = theta) %>%
    expand_grid(tibble(n = 1:sample_num)) %>%
    rowwise() %>%
    mutate(sample = rlnorm(1, meanlog = lognorm_mu, sdlog = lognorm_sigma))
  
  # Set dimension of PH distribution to be fit
  emp_A_dim = 8
  
  ## Specify initial guess A_matrix and v_alpha...
  ## Should be dimension 3. Use the same for all theta values (for now)
  # # set.seed(90210)
  # v_alpha_in = rgamma(emp_A_dim, 1)
  # v_alpha_in = v_alpha_in/sum(v_alpha_in)
  
  ## Run EM algorithm to obtain best-fit A_matrix and v_alpha...
  library(mapfit)
  sample_in = samples$sample
  # set.seed(90210) # to ensure same initial state
  out = phfit.point(ph = ph(emp_A_dim),
                    x = sample_in)
  v_alpha = as.matrix(out$alpha)
  
  return(c(
    # theta = theta,
    v_alpha = list(v_alpha),
    A_matrix = list(as.matrix(out$Q)))
  )
  #
  #   saveRDS(EMP_EM_df, "data/EM_fits_to_lognormal.rds")
  #
  #   EMP_EM_df = readRDS("data/EM_fits_to_lognormal.rds")
  
}

get_hypererlang_A <- function(theta) {
  ## Generate 100 samples from a lognormal distribution with mean theta, variance 1...
  sample_num = 100
  
  # Correct values of the lognormal parameters to match R's implementation
  lognorm_sigma <- sqrt(log(1 + (1 / theta^2)))  # Corrected sigma
  lognorm_mu <- log(theta) - 0.5 * lognorm_sigma^2  # Corrected mu
  
  samples <- tibble(theta = theta) %>%
    expand_grid(tibble(n = 1:sample_num)) %>%
    rowwise() %>%
    mutate(sample = rlnorm(1, meanlog = lognorm_mu, sdlog = lognorm_sigma))
  
  # Set dimension of PH distribution to be fit
  emp_A_dim = 4
  
  ## Run EM algorithm to obtain best-fit A_matrix and v_alpha...
  library(mapfit)
  sample_in = samples$sample
  # set.seed(90210) # to ensure same initial state
  out = phfit.point(ph = herlang(emp_A_dim),
                    x = sample_in)
  out = as.gph(out$model)
  v_alpha = as.matrix(out$alpha())
  
  return(c(
    v_alpha = list(v_alpha),
    A_matrix = list(as.matrix(out$Q())))
  )
  
}

# Get LambdaH and LambdaB and BetaH and BetaB
emp_EpiTerms <- function(A_matrix, v_alpha) {
  KB = sum(get_stable_pop(A_matrix, v_alpha)$v_B_star)
  A_dim = dim(A_matrix)[1]
  v_ones = matrix(rep(1, A_dim), ncol = 1)
  out_rates = zero_out(-as.matrix(A_matrix) %*% v_ones)
  
  LambdaB = diag(as.vector(out_rates)) %>%
    cbind(matrix(rep(0, A_dim))) %>% # need to add in zeros because oviposition class is added in
    rbind(matrix(rep(0, A_dim+1), nrow = 1)) %>%
    zero_out() %>%  list()
  LambdaH = (LambdaB[[1]] * KB / KH) %>% zero_out() %>%  list()
  
  BetaH = diag(A_dim) %>%
    cbind(matrix(rep(0, A_dim))) %>% # need to add in zeros because oviposition class is added in
    rbind(matrix(rep(0, A_dim+1), nrow = 1)) %>%
    zero_out() %>%  list()
  BetaB = BetaH
  
  return(
    c(
      LambdaB = LambdaB,
      LambdaH = LambdaH,
      BetaB = BetaB,
      BetaH = BetaH)
  )
}

### Phenomenological case functions: ----
get_phenom_A <- function(theta) {
  b = 1 / theta
  
  v_alpha = matrix(c(1/3, 1/3, 0, 1/3, 0, 0), ncol = 1)
  
  temp_matrix = matrix(rep(0, 36), nrow = 6, ncol = 6 )
  diag(temp_matrix) = c(-b, -2*b, -2*b, -3*b, -3*b, -3*b)
  temp_matrix[3,2] = 2*b
  temp_matrix[5,4] = 3*b
  temp_matrix[6,5] = 3*b
  
  A_matrix = t(temp_matrix)
  
  return(c(
    v_alpha = list(v_alpha),
    A_matrix = list(A_matrix)
  ))
  
}

# Get LambdaH and LambdaB and BetaH and BetaB
phenom_EpiTerms <- function(A_matrix, v_alpha) {
  KB = sum(get_stable_pop(A_matrix, v_alpha)$v_B_star)
  A_dim = dim(A_matrix)[1]
  v_ones = matrix(rep(1, A_dim), ncol = 1)
  out_rates = -as.matrix(A_matrix) %*% v_ones
  
  LambdaB = -diag(diag(A_matrix)) %>%
    cbind(matrix(rep(0, A_dim))) %>% # need to add in zeros because oviposition class is added in
    rbind(matrix(rep(0, A_dim+1), nrow = 1)) %>%
    list()
  LambdaH = (LambdaB[[1]] * KB / KH) %>% list()
  
  BetaH = (betaH * diag(A_dim)) %>%
    cbind(matrix(rep(0, A_dim))) %>% # need to add in zeros because oviposition class is added in
    rbind(matrix(rep(0, A_dim+1), nrow = 1)) %>%
    list()
  BetaB = (BetaH[[1]] * betaB / betaH) %>% list()
  
  return(
    c(
      LambdaB = LambdaB,
      LambdaH = LambdaH,
      BetaB = BetaB,
      BetaH = BetaH)
  )
}

### Mechanistic case functions: ----
# Get A_matrix from mechanistic parameters (!!! Not used)
get_mech_A <- function(lQ, pQ, pL, lL, pP, lP, pG, lG, sigma){
  # subintensity matrix
  out = matrix(c(
    -pQ * lQ,                           pQ * lQ,       0,       0,
    (1 - sigma) * (1- pL) * lL, -lL + sigma * (1- pL) * lL, pL * lL,       0,
    (1 - sigma) * (1 - pP) * lP,     sigma * (1 - pP) * lP,     -lP, pP * lP,
    (1 - sigma) * (1 - pG) * lG,     sigma * (1 - pG) * lG,       0,      -lG
  ),
  ncol = 4, byrow = TRUE)
  return(out)
}

# Get LambdaH and LambdaB and BetaH and BetaB
mech_EpiTerms <- function(v_B_star, lP, lG, varied_parameter) {
  
  if (varied_parameter == "pG") {
    correction_for_plot = 0.125
  } else {
    correction_for_plot = 1
  }
  
  KB = sum(unlist(v_B_star))
  v_ones = matrix(rep(1, 4), ncol = 1)
  
  LambdaH = matrix(c(0, 0, 0,            0, 0,
                     0, 0, 0,            0, 0,
                     0, 0, correction_for_plot * lP * KB / KH, 0, 0,
                     0, 0, 0,            0, 0,
                     0, 0, 0,            0, 0
  ),
  nrow = 5) %>% list()
  LambdaB = matrix(c(0, 0, 0, 0,  0,
                     0, 0, 0, 0,  0,
                     0, 0, 0, 0,  0,
                     0, 0, 0, lG, 0,
                     0, 0, 0, 0,  0
  ),
  nrow = 5) %>% list()
  
  BetaH = matrix(c(0, 0, 0,     0, 0,
                   0, 0, 0,     0, 0,
                   0, 0, betaH, 0, 0,
                   0, 0, 0,     0, 0,
                   0, 0, 0,     0, 0
  ),
  nrow = 5) %>% list()
  BetaB = matrix(c(0, 0, 0, 0,     0,
                   0, 0, 0, 0,     0,
                   0, 0, 0, 0,     0,
                   0, 0, 0, correction_for_plot * betaB, 0,
                   0, 0, 0, 0,     0
  ),
  nrow = 5) %>% list()
  return(
    c(
      LambdaB = LambdaB,
      LambdaH = LambdaH,
      BetaB = BetaB,
      BetaH = BetaH)
  )
}


### Generic functions: ----
# Function to replace small numerical values with zeroes
zero_out = function(vector_in) {
  in_rows = dim(vector_in)[1]
  in_cols = dim(vector_in)[2]
  vector_out = vector_in
  vector_out[vector_out < .Machine$double.eps] <- 0
  vector_out = matrix(vector_out, nrow = in_rows, ncol = in_cols)
  return(vector_out)
}


# Modify A matrix to include oviposition and resting
get_tilde_A <- function(A_matrix) {
  
  A_dim = dim(A_matrix)[1]
  v_ones = matrix(rep(1, A_dim), ncol = 1)
  out_rates = zero_out(-as.matrix(A_matrix) %*% v_ones)
  
  matrix(
    rbind(matrix(
      cbind(
        A_matrix,
        out_rates),
      ncol = A_dim + 1),
      c(rep(0, A_dim), -gammaV)),
    ncol = A_dim + 1)
}

# Get stable mosquito population distributions
get_stable_pop <- function(A_matrix, v_alpha) {
  
  A_dim = dim(A_matrix)[1]
  v_ones = matrix(rep(1, A_dim), ncol = 1)
  out_rates = zero_out(-as.matrix(A_matrix) %*% v_ones)
  
  if (A_dim == 1) {
    green_matrix = 1/(mu * diag(A_dim) - t(A_matrix))
  } else {
    green_matrix = inv(mu * diag(A_dim) - t(A_matrix))
  }
  
  tau = t(out_rates) %*% green_matrix %*% v_alpha %>% as.double()
  varrho = 1 - tau * (gammaR / (mu + gammaR)) * (gammaV / (mu + gammaV))
  nG = 1 / varrho
  N_offspring = (varPhi / (mu + gammaV)) * (rhoJ / (muJ + rhoJ)) * tau * nG
  
  if (N_offspring > 1) {
    J_star = KJ * (1 - 1 / N_offspring)
    v_B_star = rhoJ * J_star * green_matrix %*% v_alpha * nG
    V_star = rhoJ * J_star * (1 / (mu + gammaV)) * tau * nG
    R_star = rhoJ * J_star * (1 / (mu + gammaR)) * (gammaV / (mu + gammaV)) * tau * nG
  } else {
    J_star = 0
    v_B_star = matrix(rep(0, A_dim), ncol = 1)
    V_star = 0
    R_star = 0
  }
  
  return(
    c(
      tau = tau,
      nG = nG,
      N_offspring = N_offspring,
      J_star = J_star,
      v_B_star = list(v_B_star),
      V_star = V_star,
      R_star = R_star)
  )
  # to get this to work with mutate, use
  # out_df %>%
  #   rowwise() %>%
  #   mutate(outs = list(get_stable_pop(A_matrix, v_alpha))) %>% unnest_wider(outs)
}

# Calculate theta, biting process duration
theta_calc <- function(A_matrix, v_alpha) {
  A_dim = dim(A_matrix)[1]
  green_A = inv(-(A_matrix)*1440) * 1440
  v_one = matrix(rep(1, A_dim), ncol = 1)
  theta = t(v_alpha) %*% green_A %*% v_one
}

# Calculate GCD, gonotrophic cycle duration
GCD_calc <- function(A_matrix, v_alpha) {
  # Calculate gonotrophic cycle duration as a function of phase-type distribution parameters
  GCD = theta_calc(A_matrix, v_alpha) + (1/gammaV) + (1/gammaR)
}

# Calculate R0
R0_calc <- function(A_matrix, v_alpha, LambdaH, LambdaB, BetaH, BetaB) {
  # Calculates R0 as a function of the phase-type distribution parameters
  # and specific transmission formulation
  if (is.null(dim(A_matrix))) {
    A_matrix = as.matrix(A_matrix)
    v_alpha = as.matrix(v_alpha)
  }
  
  stable_pop = get_stable_pop(as.matrix(A_matrix), as.matrix(v_alpha))
  if (stable_pop$N_offspring < 1) {
    R0 = 0
  } else {
    
  KB = sum(stable_pop$v_B_star)
  C_vec = rbind(stable_pop$v_B_star, stable_pop$V_star)
  A_dim = dim(A_matrix)[1]
  v_ones = matrix(rep(1, A_dim), ncol = 1)
  tilde_v_alpha = v_alpha %>% rbind(0)
  tilde_v_ones = v_ones %>% rbind(1)
  out_rates = zero_out(-as.matrix(A_matrix) %*% v_ones)
  # Independent terms
  host_infectious_period = 1/(gammaH + muH)
  
  # Form tilde_A
  tilde_A = get_tilde_A(A_matrix)
  tilde_out_rates = zero_out(-as.matrix(tilde_A) %*% tilde_v_ones)
  
  # Intermediate terms
  spec_mat = tilde_v_alpha %*% t(tilde_out_rates)
  GammaI = inv(mu * diag(A_dim + 1) - t(tilde_A) - (gammaR / (mu + gammaR)) * spec_mat)
  GammaE = inv((eta + mu) * diag(A_dim + 1) - t(tilde_A) - (gammaR / (mu + eta + gammaR)) * spec_mat)
  
  tauE = (eta * diag(A_dim + 1) + (eta / (mu + eta + gammaR)) * (gammaR / (mu + gammaR)) * spec_mat) %*% GammaE
  
  R02 = host_infectious_period * t(tilde_v_ones) %*% BetaH %*% LambdaH %*% GammaI %*% tauE %*% BetaB %*% LambdaB %*% C_vec  / KB
  R0 = sqrt(as.double(R02))
  
  }
  
  return(R0)
}


# GCD range ----
# Set range of GCD values to consider
# Note: GCD = biting duration + oviposition duration + resting duration
#       Only theta = biting duration will be varied
resolution = 501

theta_min = 0 # minimum biting duration of 3 minutes
theta_max = 1 / mu #21 * 24 * 60 # maximum biting duration of 3 days
theta_vec = seq(theta_min, theta_max, length.out = resolution)
theta_vec = theta_vec[-1]
inv_theta_vec = seq(1/theta_max, 1/theta_vec[1], length.out = resolution)
theta_vec = sort(unique(c(theta_vec, 1/inv_theta_vec, (1/3) * 1440, (1/2) * 1440, 1 * 1440, 2 * 1440)))


GCD_vec = theta_vec + (1/gammaR) + (1/gammaV)
# Standard exponential model ----
Standard_vec = seq(0, min(GCD_vec), length.out = resolution)
Standard_vec = Standard_vec[-1]
Standard_vec = c(Standard_vec, GCD_vec, (1/2) * 1440, 1440, 2 * 1440)

Standard_df <- tibble(theta = Standard_vec) %>%
  mutate(b = 1/theta) %>%
  # stable population terms
  rowwise() %>%
  mutate(stable_pops = list(get_standard_stable_pop(b))) %>% unnest_wider(stable_pops) %>%
  # basic reproduction number
  rowwise() %>%
  mutate(R0 = standard_R0_calc(b))


# Exponential model with resting and ovipositing ----
Exponential_df <- tibble(theta = theta_vec) %>%
  # get A and alpha
  mutate(
    A_matrix = matrix(-1/theta),
    v_alpha = matrix(1)
  ) %>%
  # stable population terms
  rowwise() %>%
  mutate(stable_pops = list(get_stable_pop(A_matrix, v_alpha))) %>% unnest_wider(stable_pops) %>%
  # epidemiological terms
  # epidemiological terms
  rowwise() %>%
  mutate(EpiTerms = list(exp_EpiTerms(A_matrix, v_alpha))) %>% unnest_wider(EpiTerms) %>%
  # basic reproduction number
  rowwise() %>%
  mutate(R0 = R0_calc(A_matrix, v_alpha, LambdaH, LambdaB, BetaH, BetaB))

# Empirical model ----

## For each theta value...
Empirical_df <- tibble(theta = theta_vec) %>%
  # get A and alpha
  rowwise() %>%
  mutate(params = list(get_emp_A(theta))) %>% unnest_wider(params) %>%
  # stable population terms
  rowwise() %>%
  mutate(stable_pops = list(get_stable_pop(A_matrix, v_alpha))) %>% unnest_wider(stable_pops) %>%
  # epidemiological terms
  rowwise() %>%
  mutate(EpiTerms = list(emp_EpiTerms(A_matrix, v_alpha))) %>% unnest_wider(EpiTerms) %>%
  # basic reproduction number
  rowwise() %>%
  mutate(R0 = R0_calc(A_matrix, v_alpha, LambdaH, LambdaB, BetaH, BetaB))

# # Hyper-Erlang model ----
# 
# ## For each theta value...
# HyperErlang_df <- tibble(theta = theta_vec) %>%
#   # get A and alpha
#   rowwise() %>%
#   mutate(params = list(get_hypererlang_A(theta))) %>% unnest_wider(params) %>%
#   # stable population terms
#   rowwise() %>%
#   mutate(stable_pops = list(get_stable_pop(A_matrix, v_alpha))) %>% unnest_wider(stable_pops) %>%
#   # epidemiological terms
#   rowwise() %>%
#   mutate(EpiTerms = list(emp_EpiTerms(A_matrix, v_alpha))) %>% unnest_wider(EpiTerms) %>%
#   # basic reproduction number
#   rowwise() %>%
#   mutate(R0 = R0_calc(A_matrix, v_alpha, LambdaH, LambdaB, BetaH, BetaB))

# Phenomenological model ----
Phenom_df <- tibble(theta = theta_vec) %>%
  rowwise() %>%
  mutate(mats = list(get_phenom_A(theta))) %>% unnest_wider(mats) %>%
  # stable population terms
  rowwise() %>%
  mutate(stable_pops = list(get_stable_pop(A_matrix, v_alpha))) %>% unnest_wider(stable_pops) %>%
  # Epidemiological terms
  rowwise() %>%
  mutate(EpiTerms = list(phenom_EpiTerms(A_matrix, v_alpha))) %>% unnest_wider(EpiTerms) %>%
  # Basic reproduction number
  rowwise() %>%
  mutate(R0 = R0_calc(A_matrix, v_alpha, LambdaH, LambdaB, BetaH, BetaB))

# Mechanistic model ----

# Set up parameter ranges
# Use "flighty" and "persistent" parameterizations
flighty_parameters = tibble(
  mosquito_type = "flighty",
  # Questing
  pQ = 1,
  lQ = 1 / 480, # 8 hours = 480 minutes
  # Landing
  pL =  0.5,
  lL = 0.1, # 10 minutes
  # Probing
  pP = 0.5,
  lP = 0.2, # 5 minutes
  # Ingesting
  pG = 0.5,
  lG = 1, # 1 minutes
  # Fleeing
  sigma = 1 - 0.1
)
persistent_parameters = tibble(
  mosquito_type = "persistent",
  # Questing
  pQ = 1,
  lQ = 1 / 480,  # 8 hours = 480 minutes
  # Landing
  pL =  0.7,
  lL = 0.1, # 10 minutes
  # Probing
  pP = 0.8,
  lP = 0.2, # 5 minutes
  # Ingesting
  pG = 0.9,
  lG = 1, # 1 minutes
  # Fleeing
  sigma = 1 - 0.66
)
base_parameters = rbind(flighty_parameters, persistent_parameters)

# Set up parameter variation
variation_resolution = resolution

vary_parameter_function = function(df_in, varied_parameter) {
  # Collect and keep track of the baseline value of the chosen parameter
  baseline_values = df_in %>%
    select(mosquito_type, any_of(varied_parameter)) %>%
    filter(!is.na(!!as.symbol(varied_parameter))) %>%
    unique() %>%
    mutate(parameter_type = "baseline") %>%
    right_join(df_in %>%
                 select(-c(any_of(varied_parameter))),
               by = "mosquito_type") %>%
    distinct() %>%
    mutate(varied_parameter = varied_parameter)
  
  # Set up variation in probability parameters
  if (varied_parameter %in% c("pQ", "pL", "pP", "pG", "sigma")) {
    # The probability parameters (pL, pP, pG, and f) will vary from 0 to 1
    prob_vec = seq(0, 1, length.out = variation_resolution+1)
    prob_vec = prob_vec[-1] # get rid of the zero entry
    
    vary_df = baseline_values %>%
      cross_join(tibble(prob = prob_vec)) %>%
      mutate(across(any_of(as.name(varied_parameter)), ~  prob))%>%
      select(-prob) %>%
      mutate(parameter_type = "varied") %>%
      rbind(baseline_values)
    # Set up variation in rate parameters
  } else {
    rate_vec = 10^seq(-6, 2, length.out = variation_resolution)
    
    vary_df = baseline_values %>%
      cross_join(tibble(rate = rate_vec)) %>%
      mutate(across(any_of(as.name(varied_parameter)), ~  rate * .x)) %>%
      select(-rate) %>%
      mutate(parameter_type = "varied") %>%
      rbind(baseline_values)
  }
  
  return(vary_df)
}

# Initialize data frame
full_variation_df = as_tibble(matrix(
  nrow = 0,
  ncol = length(colnames(base_parameters))+2,
  dimnames = list(NULL, c("parameter_type", "varied_parameter", colnames(base_parameters)))
))

# Append all the varied parameter sets
parameter_characters = c("lQ","pQ", "pL", "lL", "pP", "lP", "pG", "lG", "sigma")

for (parameter_name in parameter_characters) {
  new_df = vary_parameter_function(base_parameters, parameter_name)
  full_variation_df = rbind(full_variation_df, new_df)
}

Mech_df <- full_variation_df %>%
  relocate(mosquito_type, parameter_type, varied_parameter) %>%
  # Transform parameters into A matrix
  rowwise() %>%
  mutate(A_matrix = list(matrix(c(
    -lQ,                                           pQ * lQ,       0,       0,
    (1 - sigma) * (1- pL) * lL, -lL + sigma * (1- pL) * lL, pL * lL,       0,
    (1 - sigma) * (1 - pP) * lP,     sigma * (1 - pP) * lP,     -lP, pP * lP,
    (1 - sigma) * (1 - pG) * lG,     sigma * (1 - pG) * lG,       0,      -lG
  ),
  ncol = 4, byrow = TRUE) )) %>%
  mutate(v_alpha = list(matrix(c(1, 0, 0, 0), ncol = 1))) %>%
  # Get the theta value for each parameter set
  rowwise() %>%
  mutate(theta = theta_calc(A_matrix, v_alpha)) %>%
  # stable population terms
  rowwise() %>%
  mutate(stable_pops = list(get_stable_pop(A_matrix, v_alpha))) %>% unnest_wider(stable_pops) %>%
  # Epidemiological terms
  rowwise() %>%
  mutate(EpiTerms = list(mech_EpiTerms(v_B_star, lP, lG, varied_parameter))) %>% unnest_wider(EpiTerms) %>%
  # Basic reproduction number
  rowwise() %>%
  mutate(R0 = R0_calc(A_matrix, v_alpha, LambdaH, LambdaB, BetaH, BetaB))

# Save Mechanistic table separately
saveRDS(Mech_df, "data/Mechanistic_results.rds")

# Join all data ----
Full_df <- bind_rows(
  # Standard model
  Standard_df %>% 
    mutate(`Model type` = "Standard",
           A_matrix = as.list(-b),
           v_alpha = as.list(c(1)),
           B_star = as.list(B_star),
           N_offspring = 1/((mu / (1440 *varPhi * b)) * ((rhoJ + muJ)/ rhoJ))) %>% 
    rename(v_B_star = B_star),
  # Exponential model
  Exponential_df %>% 
    mutate(`Model type` = "Exponential",
           v_B_star = as.list(v_B_star),
           A_matrix = as.list(A_matrix),
           v_alpha = as.list(v_alpha)),
  # Empirical model
  Empirical_df %>% 
    mutate(`Model type` = "Empirical"),
  # Phenomenological model
  Phenom_df %>% 
    mutate(`Model type` = "Phenomenological"),
  # Mechanistic model (pG, lQ, and pP)
  Mech_df %>%
    filter(
      mosquito_type == "persistent",
      varied_parameter %in% c("pG", "lQ", "pP")) %>%
    mutate(`Model type` = "Mechanistic")
) %>%
  rowwise() %>%
  mutate(GCD = if_else(
    `Model type` == "Standard",
    theta - (1/gammaR) - (1/gammaV),
    theta))

saveRDS(Full_df, "data/GCD_R0_data.rds")

# Full_df$`Model type` = factor(
#   Full_df$`Model type`,
#   levels = c("Standard", "Exponential", "Empirical", "Phenomenological", "Mech-persistent-lQ", "Mech-persistent-sigma"),
# )

# Plots ----
# library(latex2exp)
# 
# Quick_plot <- Full_df %>%
#   filter(theta > (1/2) * 1440) %>%
#   filter(theta < 21 * 1440) %>%
#   # filter(R0 < 10) %>%
#   # filter(theta > 1440 / 3) %>%
#   ggplot(aes(x = theta/1440, y = R0, color = `Model type`)) +
#   geom_line(lwd = 1) +
#   geom_hline(aes(yintercept = 1), color = "red", linetype = 2) +
#   scale_x_continuous(
#     name = "Biting stage duration [Days]"
#   ) +
#   scale_y_continuous(
#     name = TeX("Basic reproduction number [$R_0$]")
#   ) +
#   theme_minimal(18)
# 
# Inv_Quick_plot <- Full_df %>%
#   # filter(`Model type` %in% c("Empirical", "Exponential", "Phenomenological")) %>% 
#   # filter(GCD > 10 * 1440) %>%
#   filter(theta > (1/2) * 1440) %>%
#   # filter(R0 < 2) %>%
#   # filter(theta > 1440 / 0.25) %>%
#   # filter(GCD > 0.1 * 1440) %>%
#   # filter(GCD > 1440 / 0.25) %>%
#   ggplot(aes(x = 1440 / theta, y = R0, color = `Model type`)) +
#   geom_line(lwd = 1) +
#   geom_hline(aes(yintercept = 1), color = "red", linetype = 2) +
#   scale_x_continuous(
#     name = TeX("Standard biting rate [Days$^{-1}$]")
#   ) +
#   scale_y_continuous(
#     name = TeX("Basic reproduction number [$R_0$]")
#   ) +
#   scale_color_viridis_d() +
#   theme_minimal(18)
# 
# Inv_Quick_plot
# 
# Quick_Mech_plot <- Mech_df %>%
#   filter(mosquito_type == "persistent") %>%
#   filter(theta > 2 * 1440) %>%
#   filter(theta < 31 * 1440) %>%
#   # filter(theta < 8 * 1440) %>%
#   # filter(theta > 1440 / 3) %>%
#   # filter(GCD > 1440 / 3) %>%
#   # filter(GCD > 7 * 1440) %>%
#   ggplot(aes(x = theta/1440, y = R0, color = varied_parameter, linetype = mosquito_type)) +
#   geom_line(lwd = 1) +
#   geom_hline(aes(yintercept = 1), color = "red", linetype = 2) +
#   scale_x_continuous(
#     name = "Biting stage duration [Days]"
#   ) +
#   scale_y_continuous(
#     name = TeX("Basic reproduction number [$R_0$]")
#   ) +
#   theme_minimal(18)
# 
# Inv_Quick_Mech_plot <- Mech_df %>%
#   filter(mosquito_type == "persistent") %>%
#   filter(theta > 1440 / 0.5) %>%
#   # filter(varied_parameter %in% c("lQ", "sigma")) %>%
#   # filter(mosquito_type == "persistent") %>%
#   # filter(varied_parameter %in% c("sigma")) %>%
#   ggplot(aes(x = 1440 / theta, y = R0, color = varied_parameter, linetype = mosquito_type)) +
#   geom_line(lwd = 1) +
#   geom_hline(aes(yintercept = 1), color = "red", linetype = 2) +
#   scale_x_continuous(
#     name = TeX("Standard biting rate [Days$^{-1}$]")
#   ) +
#   scale_y_continuous(
#     name = TeX("Basic reproduction number [$R_0$]")
#   ) +
#   theme_minimal(18)