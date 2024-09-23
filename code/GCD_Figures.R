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
#     [] Make this an input for GCD function
#     [] Make this an input for contact rates function
#     [] Make this an input for R0 function
# 5. Plotting GCD vs contact rates or 1/GCD vs R0
#     [] Set up a plotting function where only the x and y axes need to be specified
#     [] Ensure that it labels the parameters nicely


# Functions ----

#### Function: Put parameters in matrix form
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
  
  df_out = df_in %>% rowwise() %>%
    mutate(A_matrix = list(matrix(c(
      -lQ,                                        lQ,       0,       0,
      f * (1- pL) * lL, -lL + (1 - f) * (1- pL) * lL, pL * lL,       0,
      f * (1 - pP) * lP,     (1 - f) * (1 - pP) * lP,     -lP, pP * lP,
      f * (1 - pG) * lG,     (1 - f) * (1 - pG) * lG,       0,      -lG
    ),
    ncol = 4, byrow = TRUE)
    )) %>%
    mutate(condition = kappa(A_matrix)) %>% 
    # mutate(theta = (one_vec %*% -inv(t(A_matrix)) %*% alpha_vec)[1]) %>% 
    mutate(GCD = (one_vec %*% -inv(t(A_matrix)) %*% alpha_vec)[1] + (1 / gammaV) + (1 / gammaR)) %>% 
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
    rate_vec = seq(0.1, 10, length.out = variation_resolution)
    
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

#### Function: calculate equilibrium population densities ONLY FOR MECHANISTIC MODEL
mech_population_density_function <- function(df_in) {
  # Accessory vectors
  alpha_vec = matrix(rep(0, 5), ncol = 1); alpha_vec[1] = 1
  one_vec = matrix(rep(1, 5), nrow = 1)
  
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
      # Pr(survive blood-feeding state)
      # type_A_tilde = list(dim(A_tilde)),
      # test = list(mu * diag(5) - t(A_tilde[[1]])),
      tau = as.double(-matrix(rep(1, 5), nrow = 1) %*% t(A_tilde) %*% inv(mu * diag(5) - t(A_tilde)) %*% alpha_vec),
      rho = 1 - (gammaV / (mu + gammaV)) * (gammaR / (mu + gammaR)) * tau,
      # Average number of gonotrophic cycles
      nG = 1 / rho,
      # Basic offspring number
      N_offspring = tau * (varPhi / (mu + gammaV)) * (rhoL / ( rhoL + muL)) * nG,
      # Resting mosquitoes
      R_star = rhoL * KL * tau * nG * (1 - (1/N_offspring)) * (gammaV / (mu + gammaV)) / (mu + gammaR),
      r = (N_offspring - 1) * KL * ((rhoL + muL)/ varPhi) * (mu + gammaV),
      # Aquatic stage mosquitoes
      L_star = KL * (N_offspring - 1) / N_offspring,
      # Blood-feeding stages
      B_vec = list(((N_offspring - 1) * KL * rhoL / N_offspring) * (1 + (gammaR / (mu + gammaR)) * (gammaV / (mu + gammaV)) * nG * tau) * inv(mu * diag(5) - t(A_tilde)) %*% alpha_vec),
      # Total number of blood-feeding mosquitoes
      B_tot = sum(unlist(B_vec)[1:4])) # Don't include the ovipositing mosquitoes
}



#### Function: calculate effective contact rates
contact_rate_function <- function(df_in) {
  df_out = df_in %>% 
    mutate(to_host_contact = case_when(
      model_type == "Exponential" & contact_type == "RM" ~ (1/GCD) * B_tot / H,
      model_type == "Exponential" & contact_type == "Chitnis" ~ sigmaH * (1/GCD) * B_tot / ((1/GCD) * B_tot + sigmaH * H),
      model_type == "Mechanistic" & contact_type == "RM" & transmission_type == "IN" ~ pL * lL * B_tot / H,
      model_type == "Mechanistic" & contact_type == "RM" & transmission_type == "OUT" ~ lP * B_tot / H,
      model_type == "Mechanistic" & contact_type == "Chitnis" & transmission_type == "IN" ~ sigmaH * pL * lL * B_tot / (pL * lL * B_tot + sigmaH * H),
      model_type == "Mechanistic" & contact_type == "Chitnis" & transmission_type == "OUT" ~ sigmaH * lP * B_tot / (lP * B_tot + sigmaH * H)
    )) %>% 
    mutate(to_vector_contact = case_when(
      model_type == "Exponential" & contact_type == "RM" ~ (1/GCD),
      model_type == "Exponential" & contact_type == "Chitnis" ~ (1/GCD) * sigmaH * H / ((1/GCD) * B_tot + sigmaH * H),
      model_type == "Mechanistic" & contact_type == "RM" & transmission_type == "IN" ~ pP * lP,
      model_type == "Mechanistic" & contact_type == "RM" & transmission_type == "OUT" ~ pG,
      model_type == "Mechanistic" & contact_type == "Chitnis" & transmission_type == "IN" ~ pP * lP * sigmaH * H / (pP * lP * B_tot + sigmaH * H),
      model_type == "Mechanistic" & contact_type == "Chitnis" & transmission_type == "OUT" ~ pG * sigmaH * H / (pG * B_tot + sigmaH * H)
    ))
  # When model is Exponential...
  
  # When model is Mechanistic...
  
  # When contact is static (RM)...
  
  # When model is dynamic (Chitnis)...
  
  # When transmission type is IN...
  
  # When transmission type is OUT...
}

# Setting up parameter sets ----

# Parameters that won't change throughout
# Oviposition and resting
gammaV = 1/(5 * 1440) # exit rate from oviposition to resting, including bloodmeal digestion and site search (5 days)
gammaR = 1/(2 * 1440) # exit rate from resting to return to blood-feeding (2 days)
# Transmission parameters
betaH = betaV = 1
eta = 1/(6 * 1440) # 6 days for infection to develop in vector
mu = 1/(20 * 1440) # 20 day lifespan for vector
gammaH = 1/(7 * 1440) # rate of recovery in hosts (7 days)
muH = 1/(365.25 * 65 * 1440) # host mortality rate (65 years)
KH = 1E5 # host population density

KL = 3 * 1E7
rhoL = 1 / (12 * 1440) # 12 day larval development period
muL = 1 / (20 * 1440) # 62.5% probability of larval survival (20 / (20 + 12))
varPhi = 300 / 1440 # on average 3 eggs per female per day

# In the Chitnis model, maximum per-capita contact rate for hosts
sigmaH = 100

extra_parameters = tibble(betaH, betaV, eta, mu, gammaH, muH, KH, KL, rhoL, muL, varPhi, sigmaH, gammaV, gammaR)

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

# Set wow finely we explore parameter space
variation_resolution = 101

## Exponential variable ----
# The only relevant parameter for variation here is GCD, the gonotrophic cycle duration
GCD_vec = seq(1440, 1/mu, length.out = variation_resolution)

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
  mutate(B_tot = (rhoL / mu) * KL * (1 - mu * (muL + rhoL) / (varPhi * rhoL)))


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
    ungroup() #%>% 
    
  
  full_variation_df = rbind(full_variation_df, new_df)
}

# Add back in the exponential model
combined_df = rbind(
  select(full_variation_df, colnames(exp_df)), 
  exp_df)

# Sort columns for easier reading
combined_df = relocate(full_variation_df, mosquito_type, model_type, contact_type, transmission_type, varied_parameter, parameter_type)

#### Exponential variable ----


# Plot Figures ----

## GCD vs Contact rates figures ----


#### IN Contact rates ----


#### OUT Contact rates ----