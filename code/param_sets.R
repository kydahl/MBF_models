################################################################################
# Multiple blood feeding: parameter sets
################################################################################

# Load packages ---
require(tidyverse)
require(matlib)


# Parameters ----

## Fixed

### Demography
#### Host
K_H <- 100
gamma_H <- 1/30
mu_H <- 1/(365 * 50)

#### Mosquito
Lambda_M_base <- 10 # Recruitment rate (individuals / day)
b <- 1/2 # Biting rate (inverse of gonotrophic cycle length)
mu_M_base <- 1/31 # Mortality rate (inverse of lifespan) [update to depend on k]
gamma_W <- 1/1 # Return to biting rate (inverse of resting period length) 
n_G <- 3 # average number of gonotrophic cycles within the lifespan of a mosquito

rho_B <- b / (b + mu_M_base)
rho_W <- gamma_W / (gamma_W + mu_M_base)

# RHS <- 1/rho_B
# LHS <- (rho_W * (1 - (1 / n_G))^(-1))^(1/k) 


### Transmission
beta_HM_base <- 1
beta_MH_base <- 1

fixed_params <- tibble(
  K_H = K_H, gamma_H = gamma_H, mu_H = mu_H,
  Lambda_M = Lambda_M_base, b = b, mu_M = mu_M_base, gamma_W = gamma_W,
  n_G = n_G
)

# Variable parameter functions ----

## Mortality increment calculator (alpha as a function of n_G)

alpha_func <- function(fixed_params, k) {
  with(as.list(fixed_params),{
    # Check to make sure that alpha will be greater than one
    
    rho_W = gamma_W / (gamma_W + mu_M)
    rho_B = b / (b + mu_M)
    
    RHS <- 1/rho_B # decreases with b, increases with mu_M
    LHS <- (rho_W * (1 - (1 / n_G))^(-1))^(1/k)  # increases with gamma_W, decreases with n_G and nu_G
    if(LHS < RHS) {warning("Chosen parameters do not result in alpha > 1")}
    
    alpha_val <- ifelse(k == 1, 1,
                        (rho_B * ((1 - (1 / n_G)) / rho_W)^(-1 / k))^(2/(k-1)))
    return(alpha_val)
  })
}

## Mosquito fecundity
Lambda_func <- function(Lambda_M_base, alpha, k) {
  return((1 + (k - 1) * alpha) * Lambda_M_base)
}

## Mosquito mortality
mu_func <- function(fixed_params, k, jj) {
  with(as.list(fixed_params),{
    
    alpha_val <- alpha_func(fixed_params, k)
    
    k_vec <- seq(1,k)
    
    mu_M_vec <- alpha_val^(k_vec-1) * (mu_M + b) - b
    
    return(mu_M_vec)
  })
  
}

## Transmission
beta_func <- function(beta, k) {
  return(beta / k)
}

# Make big parameter set ----
# 
# ## Set k high for now
# k_vec <- seq(1,20)
# 
# ## Mosquito fecundity
# Lambda_M_increment <- 0.1
# Lambda_M_vec <- Lambda_func(Lambda_M_base, Lambda_M_increment, k_vec)
# 
# ## Mosquito mortality
# mu_M_increment <- 0.1
# mu_M_vec <- mu_func(fixed_params, k_vec)
# 
# ## Host-to-mosquito transmission probability
# beta_HM_vec <- beta_func(beta_HM_base, k_vec)
# 
# ## Mosquito-to-host transmission probability
# beta_MH_vec <- beta_func(beta_MH_base, k_vec)
# 
# ## Put into table form
# params <- tibble(
#   k = k_vec,
#   K_H = K_H, gamma_H = gamma_H, mu_H = mu_H,
#   Lambda_M = Lambda_M, b = b, mu_M = mu_M_vec, gamma_W = gamma_W,
#   # rho_W = rho_W, rho_b = rho_b,
#   beta_HM = beta_HM_vec, beta_MH = beta_MH_vec
# )
# 
# # Case 1: Independent parameters ----
# case1_params <- filter(params, k == 1)
# 
# # Case 2: Mortality depends on k ----
# # !!!! Make sure DFE calculation can deal with mu_M as a function of k!
# case2_params <- params %>% 
#   mutate(beta_HM = beta_HM_base) %>% 
#   mutate(beta_MH = beta_MH_base)
# 
# # Case 3: Transmission prob. depends on k ----
# case3_params <- params %>% 
#   mutate(mu_M_vec = mu_M_base)
# 
# # Case 4: Both mortality and trans. prob. depend on k ----
# case4_params <- params
