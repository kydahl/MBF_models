# Case studies
# Looking at how between-oviposition waiting time distribution changes with K and Gamma_i's 
# Initiated: October 2023


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(cowplot)


# Settings & parameters ----------------------------------------------------

# Common parameters
# N: Total population size
# f: Eggs per female per oviposition
# mu: Mortality rate
# gamma: Return to blood-feeding rate or 1 / oviposition period length


# Case study parameter sets (defined by values of the lambda and rho vectors)
# 2x2 grid: row 1 = equidistributed rho values, column 1 = constant lambda
#           row 2 = decreasing rho values,      column 2 = fixed "total lambda"

# Parameter set 1: simple, 
# Constant lambda_i = lambda, 
# Equidistributed rho = [1/k]

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

### Reproductive numbers ----

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


