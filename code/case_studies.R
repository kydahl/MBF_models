# Case studies
# Looking at how between-oviposition waiting time distribution changes with K and Gamma_i's 
# Initiated: October 2023


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(cowplot)


# Settings & parameters ----------------------------------------------------

# Maximum number of feedings before reaching repletion
k_max = 5
# Explore range of k values 
k_vec = 1:k_max

# Common parameters ---- 
# Biting rate
sigma = 1/3
# Eggs per ovipositing female per day
f = 100
# Larval carrying capacity
K_L = 300
# Larval intrinsic mortality
mu_L = 1/15
# Larval maturation rate
rho_L = 1 * mu_L
# Rate to leave oviposition stage
gamma = 1/1
# Intrinsic adult mortality
mu = 1/21

# Base lambda
base_lambda = 1

# Case study parameter sets (defined by values of the lambda and rho vectors)
# 2x2 grid: row 1 = equidistributed rho values, column 1 = constant lambda
#           row 2 = decreasing rho values,      column 2 = fixed "total lambda"

# Parameter set 1: simple, 
# Constant lambda_i = lambda, 
# Equidistributed rho = [1/k]

lambdas_case_1 <- function(k) {
  rep(base_lambda, k)
}

rhos_case_1 <- function(k) {
  rep(1/k , k)
}

params_1_func <- function(k) {
  tibble(k_val = rep(k,k),
         lambdas = lambdas_case_1(k),
         rhos = rhos_case_1(k)
  )
}

params_1 <- params_1_func(k_max)

# Parameter set 2: constant lambda, decreasing rho
# rho_i+1 < rho_i
# Choose some "canonical" decreasing distribution for the rho_i's or just make
# it linear (remembering to normalize to a pdf)


# Parameter set 3: fixed time to oviposition
# Differs between the two models!
# For "fate", need to fix 1/lambda = sum 1/lambda_i,j from i to j, for j = 1 to k
# Equidistributed rhos

# Parameter set 4: fixed time to oviposition, decreasing rho
# lambda = same as parameter set 3
# rho = same as parameter set 2

# Set up dynamical models -------------------------------------------------

## The "disruption" model ----

# Define the alpha and A matrices


## The "fate" model ----

# Define the alpha and A matrices



# Run simulations ---------------------------------------------------------

# Vary the maximum number of bites per oviposition cycle, k

## Equilibrium quantities ----

# Equilibrium fecundity
e_func <- function() {
  e <- (rho_L / (rho_L + mu_L)) * sigma * f / (gamma + mu)
}

# Disrupt model: tau
tau_D_func <- function(rhos, lambdas) {
  temp_tau <- 0
  for (i in 1:length(rhos)) {
    temp_lambdas <- lambdas[1:i]
    temp_rhos <- rhos[1:i]
    
    correct_rhos = 1 - temp_rhos
    
    prod_part <- prod(temp_lambdas * (1 / (mu + temp_lambdas)) * correct_rhos)
    
    new_term <- (1 - correct_rhos[i]) * temp_lambdas[i] * prod_part
    
    temp_tau <- temp_tau + new_term
  }
  return(temp_tau)
}

# Fate model: tau
tau_F_func <- function(index, lambdas) {
  if (index > length(lambdas)) {
    warning(paste0("Index is too high. Length of lambda vector is only ", length(lambdas)))
    }
  temp_lambdas <- lambdas[1:index]
  prod(temp_lambdas * (1 / (mu + temp_lambdas)))
}

### Reproductive numbers ----

# Disrupt model:
R_D_func <- function(rhos, lambdas) {
  e <- e_func()
  tau <- tau_D_func(rhos, lambdas)
  ovi_prob <- gamma / (gamma + mu)
  e * tau / (1 - ovi_prob * tau)
}


# Fate model:
R_F_func <- function(rhos, lambdas) {
  e <- e_func()
  ovi_prob <- gamma / (gamma + mu)
  
  temp_R <- 0
  
  for (i in 1:length(rhos)) {
    tau_i <- tau_F_func(i, lambdas)
    new_term <- rhos[i] * tau_i / (1 - ovi_prob * tau_i)
    
    temp_R <- temp_R + new_term
  }
  R <- e * temp_R 
}

# Get table of R values
R_table_func <- function(k_max) {
  R_table <- tibble(k = double(),
                    R = double(),
                    model = character())
  
  k_vec <- 1:k_max
  
  for (k in k_vec) {
    rhos <- rhos_case_1(k)
    lambdas <- lambdas_case_1(k)
    
    mod_rhos <- c(rhos[1:k-1], 1)
    
    R_D <- R_D_func(mod_rhos, lambdas)
    R_F <- R_F_func(rhos, lambdas)
    
    new_row <- tibble(k = k, R = c(R_D, R_F), model = c("disrupt", "fate"))
    
    R_table <- rbind(R_table, new_row)
    
    
  }
  return(R_table)
}


### Equilibria ----




## Transient trajectories ----
# Use deSolve

# Idea: in "fate" model, what is the effect of shifting the mean/mode/median of the rho distribution?
#       Come up with a way to systematically shift the mean from 1, to 2, to 3, etc.
#       How do results change as this shift occurs?


# Summarize results -------------------------------------------------------

## Disease-free case ----

### Equilibrium results ----
# Average lifespan?

# Average number of gonotrophic cycles per lifespan?

# Reproductive number vs # of bites per cycle

### Transient results ----

# Transient population size vs. # of bites per cycle
# Top = disruption model, bottom = fate
# Color time series by # of bites per cycle using an "increasing" color scheme

# Distribution of mosquitoes by # of bites per cycle

## With simple transmission ----


# Create visualizations ---------------------------------------------------

# Reproductive number table
R_table <- R_table_func(100)

case1_R_plot <- R_table %>% 
  ggplot(aes(x = k, y = R, color = model)) +
  geom_line() +
  geom_point()
case1_R_plot