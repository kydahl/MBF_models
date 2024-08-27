######### EM algorithm for fitting phase type distributions to data ############
# Initialized: April 2024

# Uses functions from 'fitting.R'

# Outline of script: 
# 0) Load libraries, necessary functions, and data
# 1) Set up data sets
# 2) Set up structured matrices for fitting
# 3) Fit matrices to the data
# 4) Transform matrix entries to model parameters
# 5) Create visualizations of parameter and model fits
# ---------------------------------------------------------------------------- #

# 0) Load libraries, necessary functions, and data ------------------------
# library(matrixdist)
# library(tidyverse)
# library(matlib) # for a bunch of matrix operations
# library(cowplot)
# library(PhaseTypeR) # to create and simulate phase type distributions
# library(PhaseType) # to fit structured phase type distributions
# library(mapfit) # alternative way to fit phase type distributions
# library(reshape2)
# library(Matrix)
# library(corrplot)

# 1) Set up data sets with increasing levels of data ----------------------

#### Get data ----
# Simulated data set
simulated_data = read_csv("data/noQ_data_continuous.csv") # "data/sample_data_continuous.csv") # 
total_sample_size = length(simulated_data$Out_time)
sample_size = 100

Y = simulated_data$Out_time[1:100]


# 2) Initialize with "arbitrary" alpha and A ------------------------------

# Be sure that A is formed correctly
subintensity_init_func <- function(PH_type, PH_order, data = NULL) {
  # alpha is always the same
  alpha = rep(0, PH_order); alpha[1] = 1
  
  if (PH_type == "Empirical") {
    if(is.null(data)) {stop("data must be provided for initializing subintensity matrix in the empirical case")}
    # Use (very) uninformed priors to create an initial "arbitrary" A
    priors <- uninf_priors_func(mean(data) / PH_order, var(data), PH_order)
    nu = priors$nu
    zeta = priors$zeta[1:PH_order]
    entries = rgamma(PH_order^2, shape = nu[1], rate = zeta[1])
    
    
  }
  if (PH_type == "Mechanistic") {
    
  } 
}


