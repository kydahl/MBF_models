# Fitting phase type distributions to data
# Fit (structured) phase type distributions to simulated mosquito blood-feeding 
# data.
# !!! NB: This fitting necessarily assumes that alpha = c(1, 0, ..., 0) !!!
# Initialized: February 2024

# Load libraries and define helper functions ----
library(PhaseType)
library(actuar)
library(tidyverse)

# Function: reshape the output of PhaseType to subintensity matrices
ph_reshape_func <- function(phtMCMC_out) {
  dimension = length(phtMCMC_out$beta)
  # Get median and mean values of parameter posterior samples
  results = data.frame(phtMCMC_out$samples) %>% 
    mutate(id = row_number()) %>% 
    pivot_longer(cols = -id) %>% 
    group_by(name) %>% 
    summarize(mean = mean(value),
              median = median(value)) %>% 
    mutate(row = as.integer(substr(name, 2, 2)),
           col = as.integer(substr(name, 3, 3)))
  
  cens_results <- filter(results, !is.na(col))
  
  A_mat <- matrix(0, ncol = dimension, nrow = dimension)
  
  for(i in 1:nrow(cens_results)){
    A_mat[cens_results[i,]$col, cens_results[i,]$row] <- cens_results[i,]$median
  }
  
  Sii_vals = as.double(-t(t(filter(results, is.na(col))$median)) - rowSums(A_mat))
  
  A_mat = A_mat + diag(Sii_vals)
  
  return(A_mat)
}

# Function: 

# Get data ----

# Define the S matrix (columnwise)
S <- matrix(c(-3.6, 9.5, 9.5, 1.8, -11.3, 0, 1.8, 0, -11.3), 3)

# Define starting state distribution
pi <- c(1, 0, 0)

# Generate 50 random absorption times from the Phase-type with subgenerator S
# and starting distribution pi, which we will try to infer next
x <- rphtype(50, pi, S)

# General fit ----

# One-dimensional (exponential)
dimension = 1
alpha = 1



# Two-dimensional
dimension = 2
alpha = runif(dimension)
alpha = alpha / sum(alpha)
nu = seq(1:dimension^2)
zeta = seq(1:dimension)
res2 <- phtMCMC(x, dimension, alpha, nu, zeta, 200, mhit=1)
print(res2)

plot(res2)

# Three-dimensional
dimension = 3
alpha = rep(0, dimension)
alpha[1] = 1
# alpha = runif(dimension)
# alpha = alpha / sum(alpha)
nu = seq(1:dimension^2)
zeta = seq(1:dimension)
res3 <- phtMCMC(x, dimension, alpha, nu, zeta, 2, mhit=1)
print(res3)

plot(res3)

A_mat = ph_reshape_func(res3)





# ...

# Ten-dimensional (?)

# FIRST: descriptive model fit (Bladt et al. 2003)
# Prior on starting state
dirpi <- c(1, 0, 0)

# Gamma prior: shape hyperparameters (one per matrix element, columnwise)
nu = rep(1, p^2)
nu <- c(24, 24, 1, 180, 1, 24, 180, 1, 24)

# Gamma prior: reciprocal scale hyperparameters (one per matrix row)
zeta = rep(10,)
zeta <- c(16, 16, 16)

# Define dimension of model to fit
n <- 3

# Perform 20 MCMC iterations (fix inner Metropolis-Hastings to one iteration
# since starts in stationarity here).  Do more in practise!!
res1 <- phtMCMC(x, n, dirpi, nu, zeta, 20, mhit=1)
print(res1)

plot(res1)

# Mechanistic fits ----

# Full mechanistic

# Approximated mechanistic (Coxian)

# Disrupt model

# Fate model


# SECOND: mechanistic model fit (Aslett and Wilson 2011)
# Prior on starting state
dirpi <- c(1, 0, 0)

# Define the structure of the Phase-type generator
TT <- matrix(c(0,"R","R",0,"F",0,0,0,"F",0,0,0,0,"F","F",0), 4)

# Gamma prior: shape hyperparameters (one per model parameter)
nu <- list("R"=180, "F"=24)

# Gamma prior: reciprocal scale hyperparameters (one per model parameter)
zeta <- c("R"=16,"F"=16)

# Perform 20 MCMC iterations.  Do more in practise!!
res2 <- phtMCMC2(x, TT, dirpi, nu, zeta, 20)
print(res2)

plot(res2)