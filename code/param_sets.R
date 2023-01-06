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
b <- 1/5 # Biting rate (inverse of gonotrophic cycle length)
K_M <- 200 # mosquito carrying capacity
k_max <- 10 # Maximum number of bites in a single gonotrophic cycle
mu_base <- 1/21 # Mortality rate (inverse of lifespan) [update to depend on k]
gamma_W <- 1/2 # Return to biting rate (inverse of resting period length) 
n_G <- 3 # average number of gonotrophic cycles within the lifespan of a mosquito

# mu_base <- 0.5 * (-(b+gamma_W) + sqrt((b-gamma_W)^2 + 4 * b * gamma_W * (1 - (1 / n_G))^(-1)))
rho_W <- gamma_W / (gamma_W + mu_base)

# Calculate additional feeding-related mortality so that the n_G is correct
mu_1 <- (b * rho_W) * (1 - (1 / n_G))^(-1) - b - mu_base
rho_B <- (b) / (b + mu_1 + mu_base)

Lambda_M_base <- 200 * mu_base # Recruitment rate (individuals / day)

(1 / (1 - rho_B * rho_W))

(RHS <- 1 - (1 / n_G))
(LHS <- ((b * k_max) / (b * k_max + mu_base))^k_max * (gamma_W / (gamma_W + mu_base)))
# Needs LHS > RHS

# Choose large b, small mu_M and small k
# Choose large gamma_W, small n_G and small mu_M

### Transmission
beta_HM_base <- 1
beta_MH_base <- 1

fixed_params <- tibble(
  K_H = K_H, gamma_H = gamma_H, mu_H = mu_H,
  Lambda_M = Lambda_M_base, b = b, mu = mu_base, mu_1 = mu_1, gamma_W = gamma_W,
  n_G = n_G, K_M = K_M
)

# Variable parameter functions ----

## Mortality increment calculator (alpha as a function of n_G)

alpha_func <- function(fixed_params, k) {
  with(as.list(fixed_params),{
    # Check to make sure that alpha will be greater than one
    
    rho_W = gamma_W / (gamma_W + mu)
    rho_B = b / (b + mu + mu_1)
    
    RHS <- 1 - (1 / n_G)
    LHS <- ((b * k_max) / (b * k_max + mu_base))^k_max * (gamma_W / (gamma_W + mu_base))

    if(LHS < RHS) {warning("Chosen parameters will not result in mu_k > 0")}
    
    # alpha_val <- ifelse(k == 1, 1,
    #                     ((rho_B^k) * ((1 - (1 / n_G)) / rho_W)^(-1))^(2/(k * (k + 1)))
    #                     )
    alpha_val <- ifelse(k == 1, 1,
                        (((rho_B^k) * rho_W * (k^k)) / (1 - (1 / n_G)))^(1 / (k * (k-1)))
    )
  
  return(alpha_val)
  })
}

## Mosquito fecundity
Lambda_func <- function(fixed_params,k) {
  with(as.list(fixed_params), {
    
  # Get baseline population size (Lambda_M_base / mu_base)
  # Set Lambda_M_k so that when k = k, population size is equal to baseline
    
    mu_k <- mu_func(fixed_params, k)
    
    rho_B_k <- (b * k) / ((b * k) + mu + mu_k)
    # rho_B_1 <- (b / (b + mu + mu_1))
    
    # prod_Lambda <- (1 - (1 / n_G)) / gamma_W + (1 / b)
    S_B <- vector(mode = "numeric", length = k)
    for (j in seq(1,k)) {
      S_B[j] <- (n_G * rho_B_k^(j-1)) / ((b * k) + mu + mu_k)
    }
    S_W <- (n_G * rho_B_k^k) / (mu + gamma_W)
    
    Lambda_k <-  K_M / (sum(S_B) + S_W)
  
  return(Lambda_k)
  
  })
}

## Mosquito mortality
mu_func <- function(fixed_params, k) {
  with(as.list(fixed_params),{
    
    alpha_val <- alpha_func(fixed_params, k)
    mu_k <- alpha_val^(k-1) * (b + mu + mu_1) - (b * k) - mu
    
    return(mu_k)
  })
  
}

## Transmission
beta_func <- function(beta_in, k) {
  return(beta_in / k)
}

# Make big parameter set ----

## Set k high for now
k_vec <- seq(1, k_max)





## Build parameter table for a given value of k
build_params <- function(fixed_params, max_k) {
  with(as.list(fixed_params), {
    k_vec <- seq(1,max_k)
    
    ## Host-to-mosquito transmission probability
    beta_HM_vec <- beta_func(beta_HM_base, k_vec)
    
    ## Mosquito-to-host transmission probability
    beta_MH_vec <- beta_func(beta_MH_base, k_vec)
    
    ## Mosquito fecundity
    Lambda_k_vec <- vector(mode = "numeric", length = length(k_vec))
    for (i in k_vec) {
      Lambda_k_vec[i] <- Lambda_func(fixed_params, k_vec[i])
    }
    
    ## Mosquito mortality
    mu_k_vec <- vector(mode = "numeric", length = length(k_vec))
    for (i in k_vec) {
      mu_k_vec[i] <- mu_func(fixed_params, k_vec[i])
    }
    
    params <- tibble(
      k = k_vec,
      K_H = K_H, gamma_H = gamma_H, mu_H = mu_H,
      Lambda_k = Lambda_k_vec, b = b, mu = mu_base, mu_k = mu_k_vec, gamma_W = gamma_W,
      # rho_W = rho_W, rho_b = rho_b,
      beta_HM = beta_HM_vec, beta_MH = beta_MH_vec,
      n_G = n_G
    )
    return(params)
  })
}

# Case 1: Independent parameters ----
build_case1 <- function(fixed_params, max_k) {
  temp_params <- build_params(fixed_params, max_k)
  
  case1_params <- temp_params %>%
    mutate(mu_k = temp_params$mu_k[1]) %>% 
    mutate(Lambda_k = temp_params$Lambda_k[1]) %>%
    mutate(beta_HM = beta_HM_base) %>%
    mutate(beta_MH = beta_MH_base)
  
  return(case1_params)
}


# Case 2: Demographic parameters depends on k ----
build_case2 <- function(fixed_params, max_k) {
  temp_params <- build_params(fixed_params, max_k)
  
  case2_params <- temp_params %>%
    mutate(beta_HM = beta_HM_base) %>%
    mutate(beta_MH = beta_MH_base)
  
  return(case2_params)
}


# Case 3: Transmission prob. depends on k ----
build_case3 <- function(fixed_params, max_k) {
  temp_params <- build_params(fixed_params, max_k)
  
  case3_params <- temp_params %>%
    mutate(mu_k = temp_params$mu_k[1]) %>% 
    mutate(Lambda_k = temp_params$Lambda_k[1])
  
  return(case3_params)
}

## Create parameter tables
case1_params <- build_case1(fixed_params, k_max)
case2_params <- build_case2(fixed_params, k_max)
case3_params <- build_case3(fixed_params, k_max)
case4_params <- build_params(fixed_params, k_max)

