################################################################################
# Multiple blood feeding model and analysis
################################################################################


# 0) Load in necessary packages, functions, settings ----------------------
# Packages
require(tidyverse)
require(reshape2)
require(deSolve)


# 1) Specify parameters ---------------------------------------------------
## Static parameters ----

b <- 0.3 / 4.5 # biting rate # KD: either 1.05 or 4.74
beta_HM <- 0.5 # transmission prob. from M to H
beta_MH <- 0.5 # transmission prob. from H to M
NH <- 1000 # host population size
NM <- 10000 # mosquito population size
gH <- 0.1 # recovery rate of hosts
muM <- 0.1 # mortality of vectors

R0 <- sqrt((b^2 * tHV * tVH * NV) / (NH * gH * mV))

Hstar <- max((NH * (b^2 * NV * tHV * tVH - NH * mV * gH)) / (b * tVH * (b * NV * tHV + NH * gH)), 100)
Vstar <- max((b^2 * NV * tVH * tHV - NH * mV * gH) / (b * tHV * (b * tVH + mV)), 1000)

## Variables ----
k <- c(1, 3, 10) # number of bites per gonotrophic cycle
muM_k <- function(mu, k) {mu * k} # !!! temporary
beta_k <- function(beta, k) {beta * k} # !!! temporary

## Initial conditions ----
y0 <- c(H = 0, I_M = 10) # These are numbers of infected, so better to initialize at smaller values
tMax <- 10 * 365 # 10 years
dt <- 0.01
tspan <- seq(0.0, tMax, dt)

## Parameter vector ----
parms.vec <- c()


# 3) ODE model ------------------------------------------------------------




