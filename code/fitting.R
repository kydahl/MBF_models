# Fitting phase type distributions to data
# Fit (structured) phase type distributions to simulated mosquito blood-feeding 
# data
# Initiated: February 2024

# Packages ----
library(tidyverse)
library(progress)
library(matlib) # for a bunch of matrix operations
library(moments) # to calculate moments of distributions
library(Matrix) # to calculate matrix exponentials
library(expm) # to calculate matrix exponentials more efficiently and accurately
library(PhaseTypeR) # to create and simulate phase type distributions

# Helper functions ----
# Function: sample from a set phase type distribution (alpha, A)
ph_sampler <- function(alpha, A, numSamps) {
  ph = PH(A, alpha)
  samples = rPH(numSamps, ph)
}

# Function: make unit basis vectors
make_basis <- function(basis, dimension) {t(t(replace(numeric(dimension), basis, 1)))}


# Load data ----
simulated_data = read_csv("data/sample_data.csv")

# For testing: set a specific phase type distribution
# subintensity matrix
set_A = matrix(c(
  -1.5, 0, 0,
  1.5, -1, 0,
  0, 1, -0.5), ncol = 3) 
set_alpha = c(0.9, 0.1, 0) # initial probability vector

# Get 100 random samples from the distribution
direct_samps = ph_sampler(set_alpha, set_A, 100)

# Fitting distributions ----
# Fitting is done through the "EM algorithm"

# Choose data for fitting:
# Set up sample vector
Y = direct_samps #sample_data$Out_time/1440 # change units to days
NumSamples = length(Y)


# 0. Initialize estimates alpha and A

# Initial estimates of alpha and A are determined by the model being considered
# because if any entry of alpha or A is set to zero, it will remain zero 
# at all points in the algorithm. Note also that implicit in this definition
# is the dimension 'p' of the phase type distribution

# Define initial estimates for alpha and A corresponding to different sub-models

p = 16
matrix_dimension = p

# Toy values for now:
make_fake_initial <- function(matrix_dimension) {
  # Create a p x p matrix such that t = - A * 1 has all positive entries
  p = matrix_dimension
  if (p == 1) {new_A = -rexp(1, rate = 1)
  } else {
    temp_A = matrix(rexp(p^2, rate = 1), ncol = p)
    offdiag_A = temp_A - diag(diag(temp_A))
    diag_A = -diag(rexp(p, rate = 1))- diag(as.vector(offdiag_A %*% t(t(rep(1, p)))))
    
    new_A = diag_A + offdiag_A
  }
  
  t_vec = - new_A %*% t(t(rep(1, p)))
  
  temp_alpha = runif(p)
  alpha =temp_alpha / sum(temp_alpha)
  
  out <- list(A = new_A, alpha = alpha, t = t_vec)
}

initial_params <- make_fake_initial(p)
A <- t(initial_params$A) # Matrix(rnorm(p^2,0), nrow = p) #set_A #
alpha <- initial_params$alpha # set_alpha + abs(rnorm(1, 0) / 100)
# alpha <- alpha / sum(alpha)
t_vec <- initial_params$t # -rowSums(A)

# EM algorithm parameters
tolerance = 1e-6
Eell = -Inf
Eell_max = -Inf
tolCheck = tolerance + 1
keep_running = TRUE
runCount = 0
tracker <- data.frame(run_num = as.integer(), loglikelihood = as.double(), 
                      best = as.double(), mat_diff = as.double(), alpha_diff = as.double())
# Initialize arrays
Init_1_array <- rep(NA, NumSamples)
Init_2_array <- array(rep(NA, p * NumSamples), c(p, NumSamples))
Init_3_array <- array(rep(NA, p * p * NumSamples), c(p, p, NumSamples))

# Set to true if there's a set phase type distribution for testing (alpha, A)
testing = FALSE

# Run EM algorithm
while (keep_running) {
  runCount = runCount + 1
  # 1. Calculate the expected values of sufficient statistics given observations 
  
  # Calculate the matrices exp(Ay) and J
  
  # Initialize arrays
  # expAy <- Init_3_array
  # J <- Init_3_array
  # denoms <- Init_1_array
  # temp_EB <- Init_2_array
  # temp_EZ <- Init_2_array
  # temp_ENi <- Init_2_array
  # temp_EN <- Init_3_array
  
  EB = EZ = ENi = rep(0, p)
  EN = Matrix(rep(0, p^2), ncol = p)
  
  for (k in 1:length(Y)) {
    y_val = Y[k]
    bigMat = expm(
      rbind(cbind(A * y_val, t_vec %*% t(alpha) * y_val),
            cbind(0 * A, A * y_val)),
      do.sparseMsg = FALSE
    )
    temp_expAy = bigMat[1:p,1:p]
    temp_J = bigMat[1:p,(p+1):(2*p)]
    
    # expAy[1:p,1:p,k] = temp_expAy
    # J[1:p,1:p,k] = temp_J
    
    denom = as.double(t(alpha) %*% temp_expAy %*% t_vec)
    Bi_numerators = alpha * diag(1,p) %*% temp_expAy %*% t_vec
    EB = EB + (Bi_numerators / denom)
    EZ = EZ + (diag(temp_J) / denom)
    ENi = ENi + (t(alpha) %*% temp_expAy %*% diag(1, p) * t(t_vec) / denom)
    EN = EN + (A * t(temp_J) / denom)
  }
  
  # 2. Update initial estimates alpha and A
  
  alpha_hat = EB / NumSamples
  
  tij_hat = EN / EZ
  
  ti_hat = ENi / EZ
  
  tii_hat = - ti_hat - rowSums(tij_hat-diag(diag(tij_hat)))
  
  A_hat = (tij_hat-diag(diag(tij_hat))) + diag(as.vector(tii_hat))
  
  t_hat = - A_hat %*% t(t(rep(1, p)))
  
  # Check to see if likelihood is within the chosen tolerance
  B_part = t(EB) %*% replace(log(alpha_hat), alpha_hat == 0, 0)
  Nij_part = sum((EN - diag(diag(EN))) * log(A_hat - diag(diag(A_hat))), na.rm = TRUE)
  Ni_part = ENi %*% t(log(replace(ti_hat, ti_hat <= 0, .Machine$double.eps)))
  Z_part = sum(diag(EZ) %*% A_hat)
  
  Eell_new = as.double(B_part + Nij_part + Ni_part - Z_part)
  print(Eell_new)
  
  if (Eell_new > Eell_max) {Eell_max = Eell_new}
  tolCheck = abs(Eell_new - Eell) < tolerance | Eell_new > Eell_max
  
  if (tolCheck & runCount > 1000) {
    keep_running = FALSE
  } else {
    if (any(is.nan(A_hat))) {
      print("Getting NaNs...")
      break
    }
    Eell = Eell_new
    A = A_hat
    alpha = alpha_hat
    t_vec = t_hat
    # t_vec = replace(t_hat, t_hat < 0, .Machine$double.eps) # in case we get numerical errors making these values negative
  }
  
  print(runCount)
  # # Sys.sleep(.01)
  if (testing) {
    # Compare actual matrix to one set for testing
    mat_norm = norm(A - set_A)
    alpha_norm = norm(Matrix(alpha_hat - set_alpha))
    tracker <- add_row(tracker, run_num = runCount, loglikelihood = Eell, best = Eell_max, 
                       mat_diff = mat_norm, alpha_diff = alpha_norm)
  } else {
    tracker <- add_row(tracker, run_num = runCount, loglikelihood = Eell, best = Eell_max, 
                       mat_diff = NA, alpha_diff = NA)
  }
}

# Compare fit and set PH ----
fit_PH = PH(matrix(as.numeric(A_hat), ncol = p), matrix(alpha_hat))
set_PH = PH(set_A, set_alpha)

summary(fit_PH)

summary(set_PH)

all.moments(rPH(1000, fit_PH), order.max = 5)

all.moments(rPH(1000, set_PH), order.max = 5)

# Compare moments

# 3. Return to 1 with updated alpha and A

