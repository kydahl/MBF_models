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
Lambda_M <- 10 # Recruitment rate (individuals / day)
b <- 1/5 # Biting rate (inverse of gonotrophic cycle length)
mu_M <- 1/21 # Mortality rate (inverse of lifespan) [update to depend on k]
gamma_W <- 1/2 # Return to biting rate (inverse of resting period length)

rho_W <- gamma_W / (gamma_W + mu_M)
rho_b <- b / (b + mu_M)


### Transmission
beta_HM <- 1
beta_MH <- 1

## Variable


## Put into table form
params <- tibble(
  K_H = K_H, gamma_H = gamma_H, mu_H = mu_H,
  Lambda_M = Lambda_M, b = b, mu_M = mu_M, gamma_W = gamma_W,
  rho_W = rho_W, rho_b = rho_b, 
  beta_HM = beta_HM, beta_MH = beta_MH
)

# Case 1: Independent parameters ----


# Case 2: Mortality depends on k ----
# !!!! Make sure DFE calculation can deal with mu_M as a function of k!

mu_M_func <- function(mu_M, k) {
  return(k * mu_M)
}

# Case 3: Transmission prob. depends on k ----
beta_func <- function(beta, k) {
  return(beta / k)
}

# Case 4: Both mortality and trans. prob. depend on k ----
