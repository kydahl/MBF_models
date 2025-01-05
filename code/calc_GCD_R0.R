# Vary questing duration (flee probability) while fixing mean blood-feeding stage duration
library(tidyverse)
library(cowplot)
source("code/utilities.R")
library(actuar)
library(matlib)
library(latex2exp) # to display TeX
library(cols4all)
library(doFuture)
library(progressr)
library(matrixcalc)

# Cases to plot ----

# 1. "Persistent" vs. "Flighty" mosquito parameter sets
#     [x] Set up two parameter sets that can be input to functions, label them
# 2. Contact rates taken as IN vs. OUT
#     [] Make the calculation for R0 take this assumption as an input
# 3. Host availability as in Ross-Macdonald or Chitnis
#     [x] Decide on a value for the Chitnis max. host contact rate (100?)
#     [x] Set up distinct parameter sets
#     [] Make the calculation for R0 take this assumption as an input
# 4. Exponential contact rates vs. mechanistic
#     [x] Set up distinct parameter sets
#     [x] Make this an input for GCD function
#     [x] Make this an input for contact rates function
#     [] Make this an input for R0 function
# 5. Plotting GCD vs contact rates or 1/GCD vs R0
#     [] Set up a plotting function where only the x and y axes need to be specified
#     [] Ensure that it labels the parameters nicely


# Functions ----

#### Function: Put parameters in matrix form (CURRENTLY UNUSED)
params_to_mat_func <- function(parameters){
  with(as.list(parameters), {
    # lQ = base_lQ * KH / B_sum
    # subintensity matrix
    out = matrix(c(
      -lQ,                                        lQ,       0,       0,
      f * (1- pL) * lL, -lL + (1 - f) * (1- pL) * lL, pL * lL,       0,
      f * (1 - pP) * lP,     (1 - f) * (1 - pP) * lP,     -lP, pP * lP,
      f * (1 - pG) * lG,     (1 - f) * (1 - pG) * lG,       0,      -lG
    ), 
    ncol = 4, byrow = TRUE) 
    return(out)
  })
}

#### Function: Calculate GCD, gonotrophic cycle duration ----
mech_GCD_calc <- function(df_in) {
  one_vec = matrix(rep(1, 4), nrow = 1)
  alpha_vec = matrix(rep(0, 4), ncol = 1); alpha_vec[1] = 1
  
  df_out = df_in %>% 
    dplyr::rowwise() %>%
    mutate(A_matrix = list(matrix(c(
      -lQ,                                        lQ,       0,       0,
      f * (1- pL) * lL, -lL + (1 - f) * (1- pL) * lL, pL * lL,       0,
      f * (1 - pP) * lP,     (1 - f) * (1 - pP) * lP,     -lP, pP * lP,
      f * (1 - pG) * lG,     (1 - f) * (1 - pG) * lG,       0,      -lG
    ),
    ncol = 4, byrow = TRUE)
    )) %>%
    mutate(singular_check = is.singular.matrix(A_matrix, tol = sqrt(.Machine$double.eps))) %>% 
    # mutate(rank = pracma::Rank(A_matrix)) %>% 
    # mutate(condition = kappa(A_matrix)) %>% 
    # mutate(theta = (one_vec %*% -inv(t(A_matrix)) %*% alpha_vec)[1]) %>% 
    mutate(GCD = ifelse(singular_check, NA,
                        (one_vec %*% -Inverse(t(A_matrix), tol = sqrt(.Machine$double.eps)) %*% alpha_vec)[1] + (1 / gammaV) + (1 / gammaR))) %>% 
    # select(-A_matrix) %>% 
    ungroup()
}

#### Function: set up mechanistic parameter variation data frames ----
vary_parameter_function = function(df_in, varied_parameter) {
  # Collect and keep track of the baseline value of the chosen parameter
  baseline_values = df_in %>% 
    filter(model_type == "Mechanistic") %>% 
    select(mosquito_type, any_of(varied_parameter)) %>% 
    filter(!is.na(!!as.symbol(varied_parameter))) %>% 
    unique() %>% 
    mutate(parameter_type = "baseline") %>% 
    right_join(df_in %>% 
                 select(-c(any_of(varied_parameter))) %>% 
                 filter(model_type == "Mechanistic"),
               by = "mosquito_type") %>% 
    distinct() %>% 
    mutate(varied_parameter = varied_parameter)
  
  # Set up variation in probability parameters
  if (varied_parameter %in% c("pL", "pP", "pG", "f")) {
    # The probability parameters (pL, pP, pG, and f) will vary from 0 to 1
    prob_vec = seq(0, 1, length.out = variation_resolution+1)
    prob_vec = prob_vec[-1] # get rid of the zero entry
    
    vary_vec = prob_vec
    vary_df = baseline_values %>% 
      cross_join(tibble(prob = prob_vec)) %>% 
      mutate(across(any_of(as.name(varied_parameter)), ~  prob))%>% 
      select(-prob) %>% 
      mutate(parameter_type = "varied") %>% 
      rbind(baseline_values)
    # Set up variation in rate parameters
  } else {
    rate_vec = 10^seq(-5, 1, length.out = variation_resolution)
    
    vary_df = baseline_values %>% 
      cross_join(tibble(rate = rate_vec)) %>% 
      # The rate parameters (lQ, lL, lP, lG) will vary from 10% of their baseline to 10 times their baseline
      # NB: these don't change for flighty/persistent mosquitoes so don't need to worry about that
      mutate(across(any_of(as.name(varied_parameter)), ~  rate * .x)) %>% 
      select(-rate) %>% 
      mutate(parameter_type = "varied") %>% 
      rbind(baseline_values)
  }
  
  return(vary_df)
}

#### Function: calculate equilibrium population densities ONLY FOR MECHANISTIC MODEL ----
mech_population_density_function <- function(df_in) {
  # Accessory vectors
  alpha_vec = matrix(rep(0, 4), ncol = 1); alpha_vec[1] = 1
  one_vec = matrix(rep(1, 4), nrow = 1)
  
  df_out = df_in %>%
    mutate(
      # Pr(survive blood-feeding state)
      tau = as.double(-matrix(rep(1, 4), nrow = 1) %*% t(A_matrix) %*% inv(mu * diag(4) - t(A_matrix)) %*% alpha_vec),
      rho = 1 - (gammaV / (mu + gammaV)) * (gammaR / (mu + gammaR)) * tau,
      # Average number of gonotrophic cycles
      nG = 1 / rho,
      # Basic offspring number
      N_offspring = tau * (varPhi / (mu + gammaV)) * (rhoJ / ( rhoJ + muJ)) * nG,
      # Resting mosquitoes
      R_star = rhoJ * KJ * tau * nG * (1 - (1/N_offspring)) * (gammaV / (mu + gammaV)) / (mu + gammaR),
      r = (N_offspring - 1) * KJ * ((rhoJ + muJ)/ varPhi) * (mu + gammaV),
      # Aquatic stage mosquitoes
      L_star = KJ * (N_offspring - 1) / N_offspring,
      # Ovipositing mosquitoes
      V_star = r / (mu + gammaV),
      # Blood-feeding stages
      B_vec = list(((N_offspring - 1) * KJ * rhoJ / N_offspring) * (1 + (gammaR / (mu + gammaR)) * (gammaV / (mu + gammaV)) * nG * tau) * inv(mu * diag(4) - t(A_matrix)) %*% alpha_vec),
      # Blood-feeding and ovipositing stages
      BV_vec = list(matrix(rbind(B_vec, V_star), ncol = 1)),
      # Total number of blood-feeding mosquitoes
      B_tot = sum(unlist(B_vec))
      
    )
}

#### Function: calculate effective contact rates ----
contact_rate_function <- function(df_in) {
  # Model type can be Exponential or Mechanistic
  # Contact type can be static (RM) or dynamic (Chitnis)
  # Transmission type can be IN or OUT
  
  df_out = df_in %>% 
    # rowwise() %>% 
    
    # Contact rates to the HOST
    # in units of "contacts per host per day"
    mutate(to_host_contact = case_when(
      # Exponential model
      model_type == "Exponential" & contact_type == "RM" ~ (1/GCD) * B_tot / KH,
      model_type == "Exponential" & contact_type == "Chitnis" ~ sigmaH * (1/GCD) * B_tot / ((1/GCD) * B_tot + sigmaH * KH),
      # Mechanistic model
      model_type == "Mechanistic" & contact_type == "RM" & transmission_type == "IN" ~ pL * lL * B_vec[2] / KH,
      model_type == "Mechanistic" & contact_type == "RM" & transmission_type == "OUT" ~ lP * B_vec[3] / KH,
      model_type == "Mechanistic" & contact_type == "Chitnis" & transmission_type == "IN" ~ sigmaH * pL * lL * B_vec[2] / (pL * lL * B_vec[2] + sigmaH * KH),
      model_type == "Mechanistic" & contact_type == "Chitnis" & transmission_type == "OUT" ~ sigmaH * lP * B_vec[3] / (lP *  B_vec[3] + sigmaH * KH)
    )) %>%
    # Contact rates to the VECTOR
    # in units of "contacts per mosquito per day"
    mutate(to_vector_contact = case_when(
      # Exponential model
      model_type == "Exponential" & contact_type == "RM" ~ (1/GCD),
      model_type == "Exponential" & contact_type == "Chitnis" ~ (1/GCD) * sigmaH * KH / ((1/GCD) * B_tot + sigmaH * KH),
      # Mechanistic model
      model_type == "Mechanistic" & contact_type == "RM" & transmission_type == "IN" ~ pP * lP * B_vec[3] / B_tot,
      model_type == "Mechanistic" & contact_type == "RM" & transmission_type == "OUT" ~ lG * B_vec[4] / B_tot,
      model_type == "Mechanistic" & contact_type == "Chitnis" & transmission_type == "IN" ~ pP * lP * (B_vec[3] / B_tot) * sigmaH * KH / (pP * lP * B_vec[3] + sigmaH * KH),
      model_type == "Mechanistic" & contact_type == "Chitnis" & transmission_type == "OUT" ~ lG * (B_vec[4] / B_tot) * sigmaH * KH / (pG * B_vec[4] + sigmaH * KH)
    )) # %>% 
  # ungroup()
}

#### Function: calculate basic reproduction number ----
R0_function <- function(df_in) {
  # Model type can be Exponential or Mechanistic
  # Contact type can be static (RM) or dynamic (Chitnis)
  # Transmission type can be IN or OUT
  
  # Accessory vectors
  alpha_vec = matrix(rep(0, 5), ncol = 1); alpha_vec[1] = 1
  one_vec = matrix(rep(1, 5), nrow = 1)
  
  # Function to replace small numerical values with zeroes
  zero_out = Vectorize(function(vector_in) {
    vector_out = vector_in
    vector_out[vector_out < 1e-16] <- 0
    return(vector_out)
  })
  
  df_out = df_in %>%
    mutate(
      # Modify the matrix to include the oviposition state
      A_tilde = list(matrix(
        rbind(matrix(
          cbind(
            A_matrix,
            c(0, 0, 0, pG * lG)),
          ncol = 5),
          c(0, 0, 0, 0, -gammaV)),
        ncol = 5)),
      # Set up transmission matrices: beta and Lambda
      betaH_mat = list(matrix(c(
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, betaH, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0
      ),ncol = 5, byrow = TRUE)),
      betaV_mat = list(matrix(c(
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, betaB, 0,
        0, 0, 0, 0, 0
      ),ncol = 5, byrow = TRUE)),
      LambdaH_mat = case_when(
        # !!! double-check all these when returning. only RM + OUT is checked
        contact_type == "RM" & transmission_type == "IN" ~ list(t(matrix(c(
          0, 0, 0, 0, 0,
          0, 0, pL * lL * B_vec[2] / KH, 0, 0,
          0, 0, 0, 0, 0,
          0, 0, 0, 0, 0,
          0, 0, 0, 0, 0
        ), ncol = 5, byrow = TRUE))),
        contact_type == "RM" & transmission_type == "OUT" ~ list(t(matrix(c(
          0, 0, 0, 0, 0,
          0, 0, 0, 0, 0,
          0, 0, lP * B_vec[3] / KH, 0, 0,
          0, 0, 0, 0, 0,
          0, 0, 0, 0, 0
        ), ncol = 5, byrow = TRUE))),
        contact_type == "Chitnis" & transmission_type == "IN" ~ list(t(matrix(c(
          0, 0, 0, 0, 0,
          0, 0, sigmaH * pL * lL * B_vec[2] /(pL * lL * B_vec[2] + sigmaH * KH), 0, 0,
          0, 0, 0, 0, 0,
          0, 0, 0, 0, 0,
          0, 0, 0, 0, 0
        ), ncol = 5, byrow = TRUE))),
        contact_type == "Chitnis" & transmission_type == "OUT" ~ list(t(matrix(c(
          0, 0, 0, 0, 0,
          0, 0, 0, 0, 0,
          0, 0,  sigmaH * lP * B_vec[3] /(lP * B_vec[3] + sigmaH * KH), 0, 0,
          0, 0, 0, 0, 0,
          0, 0, 0, 0, 0
        ), ncol = 5, byrow = TRUE)))
      ),
      LambdaV_mat = case_when(
        # !!! double-check all these when returning. only RM + OUT is checked
        contact_type == "RM" & transmission_type == "IN" ~ list(t(matrix(c(
          0, 0, 0, 0, 0,
          0, 0, 0, 0, 0,
          0, 0, 0, pP * lP, 0,
          0, 0, 0, 0, 0,
          0, 0, 0, 0, 0
        ), ncol = 5, byrow = TRUE))),
        contact_type == "RM" & transmission_type == "OUT" ~ list(t(matrix(c(
          0, 0, 0, 0, 0,
          0, 0, 0, 0, 0,
          0, 0, 0, 0, 0,
          0, 0, 0, lG, 0,
          0, 0, 0, 0, 0
        ), ncol = 5, byrow = TRUE))),
        contact_type == "Chitnis" & transmission_type == "IN" ~ list(t(matrix(c(
          0, 0, 0, 0, 0,
          0, 0, 0, 0, 0,
          0, 0, 0, sigmaH * pP * lP * B_vec[3] /(pP * lP * B_vec[3] + sigmaH * KH), 0,
          0, 0, 0, 0, 0,
          0, 0, 0, 0, 0
        ), ncol = 5, byrow = TRUE))),
        contact_type == "Chitnis" & transmission_type == "OUT" ~ list(t(matrix(c(
          0, 0, 0, 0, 0,
          0, 0, 0, 0, 0,
          0, 0, 0, 0, 0,
          0, 0, 0, sigmaH * lG * B_vec[4] /(lG * B_vec[4] + sigmaH * KH), 0,
          0, 0, 0, 0, 0
        ), ncol = 5, byrow = TRUE)))
      ),
      temp_mat = list(matrix((-one_vec) %*% t(A_tilde), ncol = 1)),
      temp_mat = list(matrix(zero_out(temp_mat), nrow = 5)),
      spec_mat = list(alpha_vec %*% t(temp_mat)),
      M1 = list(mu * diag(5) - t(A_tilde) - (gammaR / (mu + gammaR)) * spec_mat),
      M2 = list(eta * diag(5) + (eta / (mu + gammaR + eta)) * (1 / (mu + gammaR)) * spec_mat),
      M3 = list((eta + mu) * diag(5) - t(A_tilde) - (gammaR / (mu + gammaR + eta)) * spec_mat),
      Q = list(inv(M1) %*% M2 %*% inv(M3)),
      R0 = sqrt(one_vec %*% betaH_mat %*% LambdaH_mat %*% Q %*% betaV_mat %*% LambdaV_mat %*% BV_vec / (sum(B_tot) * (muH + gammaH))),
      R0 = ifelse(N_offspring < 1, 0, R0)
    )
}

# Setting up parameter sets ----

# Parameters that won't change throughout
# Oviposition and resting
gammaV = 1/(5 * 1440) # exit rate from oviposition to resting, including bloodmeal digestion and site search (5 days)
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
varPhi = 3 / 1440 # on average 3 eggs per female per day

# In the Chitnis model, maximum per-capita contact rate for hosts
sigmaH = 100 / 1440 # = 100 bites / day = 100 (bites / day) * (1 day / 1440 minutes)

extra_parameters = tibble(betaH, betaB, eta, mu, gammaH, muH, KH, KJ, rhoJ, muJ, varPhi, sigmaH, gammaV, gammaR)

base_params_flighty = tibble(
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
  f =  0.1
)

base_params_persistent = tibble(
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
  f =  0.66 
)

base_parameters = rbind(base_params_flighty, base_params_persistent)

# Set base gonotrophic cycle duration for the exponential model
base_GCD = 3 * 1440 # 3 days

full_parameters = base_parameters %>% 
  cbind(extra_parameters) %>%
  # Calculate GCD for mechanistic parameters
  mech_GCD_calc() %>% 
  # Add in exponential parameters
  cross_join(tibble(model_type = c("Exponential", "Mechanistic"))) %>%
  mutate(GCD = case_when(
    model_type == "Exponential" ~ base_GCD,
    model_type == "Mechanistic" ~ GCD
  )) %>% 
  # Exponential has NAs in unused parameters
  mutate(across(pQ:f, ~ if_else(model_type == "Exponential", NA, .x))) %>% 
  # Add in dynamic vs. static contact rates
  cross_join(tibble(contact_type = c("RM", "Chitnis"))) %>% 
  mutate(sigmaH = case_when(
    contact_type == "RM" ~ Inf,
    contact_type == "Chitnis" ~ sigmaH
  )) %>% 
  # Add in IN vs. OUT for transmission types
  cross_join(tibble(transmission_type = c("IN", "OUT"))) %>% 
  relocate(model_type, mosquito_type, contact_type, transmission_type) %>% 
  as_tibble()

# Set up parameter variation ----

# Set how finely we explore parameter space
variation_resolution = 10001

## Exponential variable ----
# The only relevant parameter for variation here is GCD, the gonotrophic cycle duration
GCD_vec = 1440 * seq(0, 1/mu/1440, length.out = variation_resolution)
GCD_vec = GCD_vec[-1]

exp_df = full_parameters %>% 
  filter(model_type == "Exponential") %>% 
  select(-GCD) %>% 
  cross_join(tibble(GCD = GCD_vec)) %>% # Function: set up parameter variation dataframe
  mutate(varied_parameter = "theta") %>% 
  mutate(parameter_type = "varied") %>% 
  rbind(
    full_parameters %>% 
      filter(model_type == "Exponential") %>% 
      mutate(varied_parameter = "theta") %>% 
      mutate(parameter_type = "baseline")
  ) %>% 
  rowwise() %>% 
  mutate(
    b = 1/(GCD - (1 / gammaV) - (1 / gammaR)),
    J_eq = KJ * (1 - (varPhi * (rhoJ / (rhoJ + muJ) * (b / (b + mu)) * (1 / (gammaV + mu)) * (1 - (b / (b + mu)) * (gammaV / (gammaV + mu)) * (gammaR / (gammaR + mu)))^(-1)))^(-1)),
    V_eq = (1 - (b / (b + mu)) * (gammaV / (gammaV + mu)) * (gammaR / (gammaR + mu)))^(-1) * (b / (b + mu)) * (1 / (gammaV + mu)) * rhoJ * J_eq ,
    B_eq = ((gammaV + mu) / b) * V_eq,
    nG = (1 - (gammaR / (gammaR + mu)) * (gammaV / (gammaV + mu)) * (b / (b + mu)))^(-1),
    nE = (1 - (gammaR / (eta + gammaR + mu)) * (gammaV / (eta + gammaV + mu)) * (b / (eta + b + mu)))^(-1),
    R02 = betaB * b * (B_eq / KH) * (1 / (gammaH + muH)) * betaH * (b / (b + mu)) * (eta / (eta + gammaV + mu)) * (gammaR / (gammaR + mu)) * (gammaV / (gammaV + mu)) * nG + (gammaV / (eta + gammaV + mu)) * (gammaR / (eta + gammaR + mu)) * ( (eta / (eta + b + mu)) + (gammaR / (gammaR + mu)) * (1 + ( (eta / (eta + b + mu)) * (b / (b + mu)) + (b / (eta + b + mu)) * (eta / (eta + gammaV + mu))) * (gammaV / (gammaV + mu))) * nG) * nE,
    R0 = sqrt(R02)
    ) %>% 
  filter(b > 0) %>% 
  # mutate(B_vec = NA) %>% 
  ungroup() %>% 
  # contact_rate_function() %>% 
  # mutate(b = 1/GCD) %>% 
  # mutate(R0 = case_when(
  #   contact_type == "RM" ~ b * sqrt(betaB * B_tot * betaH * (eta / (mu +eta)) / (KH * (muH + gammaH) * mu)),
  #   contact_type == "Chitnis" ~ (sigmaH * b /(sigmaH * KH + b * B_tot)) * sqrt(betaB * b * B_tot * betaH * sigmaH * (eta / (mu +eta)) * KH / ((muH + gammaH) * mu))
  # )) %>% 
  select(-b)



## Mechanistic variables ----
# The relevant parameters here are
# lQ, pL, lL, pP, lP, pG, lG, f
parameter_characters = c("lQ", "pL", "lL", "pP", "lP", "pG", "lG", "f")

# Initialize data frame
full_variation_df = as_tibble(matrix(
  nrow = 0, 
  ncol = length(colnames(full_parameters))+2, 
  dimnames = list(NULL, c("parameter_type", "varied_parameter", colnames(full_parameters)))
))


# Add in variation data frame for each parameter
for (parameter_name in parameter_characters) {
  new_df = vary_parameter_function(full_parameters, parameter_name) %>% 
    # Calculate GCD
    mech_GCD_calc() %>% 
    rowwise() %>%
    mech_population_density_function() %>% 
    ungroup() %>% 
    rowwise() %>%
    contact_rate_function() %>% 
    rowwise() %>%
    R0_function()
  
  
  full_variation_df = rbind(full_variation_df, new_df)
}

# Add back in the exponential model
combined_df <- rbind(
  select(full_variation_df, intersect(colnames(full_variation_df), colnames(exp_df))),
  select(exp_df, intersect(colnames(full_variation_df), colnames(exp_df)))) %>% 
  # Sort columns for easier reading
  relocate(mosquito_type, model_type, contact_type, transmission_type, varied_parameter, parameter_type)

write_rds(combined_df, "data/GCD_R0_vals.rds.gz", compress = "gz")
