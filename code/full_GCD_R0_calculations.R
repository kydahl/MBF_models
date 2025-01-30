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

# General functions ----

# Calculate R0
R0_calc <- function(A, v_alpha, LambdaB, LambdaH, BetaB, BetaH) {
  # Calculates R0 as a function of the phase-type distribution parameters
  # and specific transmission formulation
}

# Calculate theta, biting process duration
theta_calc <- function(A, v_alpha) {
  A_dim = dim(A)[1]
  green_A = inv(-A)
  v_one = matrix(rep(1, A_dim), 1)
  theta = v_one %*% green_A %*% v_alpha
}

# Calculate GCD, gonotrophic cycle duration
GCD_calc <- function(A, v_alpha) {
  # Calculate gonotrophic cycle duration as a function of phase-type distribution parameters
  GCD = theta_calc(A, v_alpha) + (1/gammaV) + (1/gammaR)
}


# Fixed Parameter values ----
# Oviposition and resting
gammaV = 1 / (2 * 1440) #1/(5 * 1440) # exit rate from oviposition to resting, including bloodmeal digestion and site search (5 days)
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
varPhi = 300 / 1440 # on average 3 eggs per female per day

# GCD range ----
# Set range of GCD values to consider
# Note: GCD = biting duration + oviposition duration + resting duration
#       Only theta = biting duration will be varied
resolution = 100

theta_min = 8 * 60 # minimum biting duration of 8 hours
theta_max = 20 * 24 * 60 # maximum biting duration of 20 days
theta_vec = seq(theta_min, theta_max, resolution)

GCD_vec = theta_vec + (1/gammaR) + (1/gammaV)


# Exponential model ----
b_vec = 1/theta_vec

## Calculate R0...

# Empirical model ----

## For each theta value...

## Generate 100 samples from a lognormal distribution with mean theta, variance 1...

## Specify initial guess A and v_alpha...
## Should be dimension 3. Use the same for all theta values (for now)

## Run EM algorithm to obtain best-fit A and v_alpha...

## [Sanity check] Calculate theta_calc(A, v_alpha) and compare to theta

## Determine LambdaH, LambdaB, BetaH, BetaB from A...

## Calculate R0...

# Phenomenological model ----

## For each theta value...

## Assemble A and v_alpha following structure in paper...

## Determine LambdaH, LambdaB, BetaH, BetaB from A...

## Calculate R0....

# Mechanistic model ----

## !!! Cannibalize most of this from old code

## Define ranges for each of the parameters...
## Might consider flighty vs persistent

## Set up parameter grid(s)...

## Calculate GCD over parameter grid(s)...

## Determine LambdaH, LambdaB, BetaH, BetaB over parameter grid(s)...

## Calculate R0 over parameter grid(s)...