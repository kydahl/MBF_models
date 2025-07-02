# Fitting phase-type distribution using EM algorithm



# Load data and libraries -------------------------------------------------
library(matrixdist)
library(tidyverse)
library(matlib) # for a bunch of matrix operations
library(cowplot)
library(PhaseTypeR) # to create and simulate phase type distributions
library(PhaseType) # to fit structured phase type distributions
library(mapfit) # alternative way to fit phase type distributions
library(reshape2)
library(Matrix)
library(lhs) # easy Latin Hypercube sampling
library(actuar) # to create and simulate phase type distributions
library(data.table)
library(progressr) # Progress bars
library(doFuture) # Faster parallel processing
library(philentropy)


# Define parameters -------------------------------------------------------

# Define phase type distribution
# Questing
pQ = 1
lQ = 1 / (8 * 60) # 8 hours = 480 minutes

# Landing
pL = 0.7
lL = 0.1 # 10 minutes

# Probing
pP = 0.8
lP = 0.2 # 5 minutes

# Ingesting
pG = 0.9 #0.75
lG = 1 # 1 minutes

# Fleeing
f = 0.66

param_table = tibble(
  name = c("lambda_Q", "p_L", "lambda_L", "p_P", "lambda_P", "p_G", "lambda_G", "f"),
  Symbol = c("$\\lambda_Q$", "$p_L$", "$\\lambda_L$", "$p_P$", "$\\lambda_P$", "$p_G$", "$\\lambda_G$", "$f$"),
  Description = c("Exit rate from questing stage (per minute)", 
                  "Probability of progressing from landing to probing ",
                  "Exit rate from landing stage (per minute)", 
                  "Probability of progressing from probing to ingesting",
                  "Exit rate from probing stage (per minute)",
                  "Probability of progressing from ingesting to ovipositing",
                  "Exit rate from ingestion stage (per minute)",
                  "Probability of seeking a new vertebrate host given feeding failure"
  ),
  value = c(lQ, pL, lL, pP, lP, pG, lG, f)
)

eta = 1/(6 * 1440) # 6 days for infection to develop in vector
mu = 1/(21 * 1440) # 20 day lifespan for vector
gamma = 1/(3 * 1440) # rate of return from oviposition to blood-feeding (3 days)
gammaH = 1/(7 * 1440) # rate of recovery in hosts (7 days)
muH = 1/(365.25 * 65 * 1440) # host mortality rate (65 years)
KH = 10^7 # host population density

KL = 300
rhoL = 1/(12 * 1440)
muL = 1/(20 * 1440)
varPhi = 3 / 1440

# Epidemiological parameters
epi_parms = tibble(betaP = 1, betaG = 1, betaH = 1, betaV = 1)

parameters = tibble(
  eta = eta, mu = mu, gamma = gamma, gammaH = gammaH, muH = muH, KH = KH, KL = KL, rhoL = rhoL, muL = muL, varPhi = varPhi, epi_parms
)


# Simulate data -----------------------------------------------------------
set_A = matrix(c(
  -lQ,               lQ,                           0,       0,
  f * (1- pL) * lL,  -lL + (1 - f) * (1- pL) * lL, pL * lL, 0,
  f * (1 - pP) * lP, (1 - f) * (1 - pP) * lP,      -lP,     pP * lP,
  f * (1 - pG) * lG, (1 - f) * (1 - pG) * lG,      0,       -lG), ncol = 4, byrow = TRUE) 
set_alpha = Matrix(c(1, 0, 0, 0), nrow = 1) # initial probability vector
ph = PH(set_A, as.matrix(set_alpha))

# Get 100 random samples from the distribution
direct_samps = rphtype(10000, as.matrix(set_alpha), set_A)

test_data = sample(direct_samps, 1000)


# Helper: Calculate J matrix ----------------------------------------------
J_mat_func <- function(data_point, subint_mat, alpha, t_vec, PH_dim) {
  temp_subint = Matrix(rep(0, 4 * PH_dim^2), ncol = 2 * PH_dim, nrow = 2 * PH_dim, sparse = TRUE)
  temp_subint[1:PH_dim,1:PH_dim] = subint_mat
  temp_subint[(PH_dim+1):(2 * PH_dim),(PH_dim+1):(2 * PH_dim)] = subint_mat
  temp_subint[1:PH_dim,(PH_dim+1):(2 * PH_dim)] = t_vec %*% alpha
  
  temp_matexp = Matrix(expm(temp_subint * data_point), sparse = TRUE)
  
  J_k = temp_matexp[1:PH_dim,(PH_dim+1):(2 * PH_dim)]
  
  return(J_k)
}


# Expectation Step --------------------------------------------------------
PH_expect_func <- function(data_in, subint_mat, alpha) {
  # Intermediate quantities
  data_length = length(data_in)
  
  A_mat = subint_mat
  PH_dim = dim(A_mat)[1]
  e_vec = matrix(rep(1, PH_dim), ncol = 1)
  t_vec = - A_mat %*% e_vec
  y_vec = data_in
  
  EB = EZ = ENi = rep(0, PH_dim)
  ENij = Matrix(rep(0, PH_dim^2), ncol = PH_dim, nrow = PH_dim)
  
  plan(sequential)
  # plan(multisession)
  # plan(list(tweak(multisession, workers = PH_dim),
  #           tweak(multisession, workers = availableCores()/PH_dim)))
  # progressr::handlers(global = TRUE)
  # progressr::handlers("cli")
  
  expectations = foreach(
    i = 1:PH_dim,
    .combine = 'rbind'
  ) %do% {
    e_i_vec = matrix(diag(PH_dim)[,i], ncol = 1)
    EB_val = EZ_val = ENi_val = 0
    ENij_row = rep(0, PH_dim)
    
    
    inter_table = foreach(
      k = 1:data_length,
      .combine = 'rbind'
    ) %do% {
      # for (k in 1:data_length) {
      # Calculate common denominator
      common_denom = alpha %*% expm(A_mat * y_vec[k]) %*% t_vec
      
      # Calculate J matrix
      J_k = J_mat_func(y_vec[k], A_mat, alpha, t_vec, PH_dim)
      
      # E[B_i]: Expected number of processes starting in state i
      #         sum_k=1^N {(alpha_i * e_i' * exp(A * y_k) * t) / (alpha * exp(A * y_k) * t)}
      temp_EB = (alpha[i] * t(e_i_vec) %*% expm(A_mat * y_vec[k]) %*% t_vec) / common_denom
      # EB_val = EB_val + temp_EB
      
      # E[Z_i]: Expected total time spent in state i in all processes
      temp_EZ = J_k[i,i] / common_denom
      # EZ_val = EZ_val + temp_EZ
      
      # E[N_i]: Expected number of processes that exit to the absorbing state from state i
      temp_ENi = alpha %*% expm(A_mat * y_vec[k]) %*% e_i_vec * t_vec[i] / common_denom
      # ENi_val = ENi_val + temp_ENi
      
      temp = tibble(EB = temp_EB[1,1], EZ = temp_EZ[1,1], ENi = temp_ENi[1,1])
      # E[N_ij]: Expected number of observed jumps from state i to state j
      for (j in 1: PH_dim) {
        # print(paste0("Column number # ", j, " out of ", PH_dim))
        temp_ENij_row = A_mat[i,j] * J_k[j,i] / common_denom
        temp = mutate(temp, "ENi{j}" := temp_ENij_row[1,1])
        
      }
      temp
    }
    
    final = inter_table %>% summarise(across(everything(), ~ sum(., na.rm = TRUE)))
    EB[i] = final$EB
    EZ[i] = final$EZ
    ENi[i] = final$ENi
    ENij[i,] = Matrix(as.double(final[1,4:(3+PH_dim)]), nrow = 1)
    
    tibble(EB = EB[i], EZ = EZ[i], ENi = ENi[i], as_tibble(t(ENij[i,]), .name_repair = "unique"))
    
  }
  
  EB = Matrix(expectations$EB, nrow = 1)
  EZ = Matrix(expectations$EZ, nrow = 1)
  ENi = Matrix(expectations$ENi, nrow = 1)
  ENij = Matrix(as.matrix(expectations[1:PH_dim, 4:(3+PH_dim)]))
  # 
  # EZ = foreach(
  #   i = 1:PH_dim,
  #   .combine = rbind
  # ) %do% {
  # }
  # 
  # ENi = foreach(
  #   i = 1:PH_dim,
  #   .combine = rbind
  # ) %do% {
  # }
  # 
  # EN = foreach(
  #   i = 1:PH_dim,
  #   j = 1:PH_dim,
  #   .combine = rbind
  # ) %do% {
  # }
  # 
  # 
  # 
  expectations = foreach(
    i = 1:PH_dim,
    .combine = cbind
  ) %do% {
  # for (i in 1: PH_dim) {
    # print(paste0("Row number # ", i, " out of ", PH_dim))
    # Set up standard basis vector
    e_i_vec = matrix(diag(PH_dim)[,i], ncol = 1)

    for (k in 1:data_length) {
      # print(paste0("Iteration # ", k, " out of ", data_length))

      # Calculate J matrix
      J_k = J_mat_func(y_vec[k], A_mat, alpha, t_vec, PH_dim)

      # Calculate common denominator
      common_denom = alpha %*% expm(A_mat * y_vec[k]) %*% t_vec

      # E[B_i]: Expected number of processes starting in state i
      temp_EB = (alpha[i] * t(e_i_vec) %*% expm(A_mat * y_vec[k]) %*% t_vec) / common_denom
      EB[i] = EB[i] + temp_EB

      # E[Z_i]: Expected total time spent in state i in all processes
      temp_EZ = J_k[i,i] / common_denom
      EZ[i] = EZ[i] + temp_EZ

      # E[N_ij]: Expected number of observed jumps from state i to state j
      for (j in 1: PH_dim) {
        # print(paste0("Column number # ", j, " out of ", PH_dim))
        temp_ENij = A_mat[i,j] * J_k[j,i] / common_denom
        EN[i,j] = EN[i,j] + temp_ENij
      }

      # E[N_i]: Expected number of processes that exit to the absorbing state from state i
      temp_ENi = alpha %*% expm(A_mat * y_vec[k]) %*% e_i_vec * t_vec[i] / common_denom
      ENi[i] = ENi[i] + temp_ENi

    }
    return(list(EB = EB, EZ = EZ, EN = EN, ENi = ENi))
  }
  
  # plan(sequential)
  return(list(EB = EB, EZ = EZ, EN = ENij, ENi = ENi))
}


# Update Step -------------------------------------------------------------
PH_update_func <- function(expect_out, data_length) {
  with(as.list(expect_out), {
    PH_dim = length(EB)
    
    alpha_out = t_vec_out = Matrix(rep(0, PH_dim), nrow = 1)
    temp_mat_out = Matrix(rep(0, PH_dim^2), nrow = PH_dim, ncol = PH_dim)
    
    for (i in 1:PH_dim) {
      alpha_out[i] = EB[i] / data_length
      t_vec_out[i] = ENi[i] / EZ[i]
      for (j in 1:PH_dim) {
        temp_mat_out[i,j] = ifelse(i==j, 0, EN[i,j] / EZ[i])
      }
    }
    temp_rowsums = Matrix(rowSums(temp_mat_out), ncol = 1)
    
    A_mat_out = temp_mat_out
    for (i in 1:PH_dim) {
      A_mat_out[i,i] = -temp_rowsums[i] - t_vec_out[i]
    }
    return(list(alpha = alpha_out, A = A_mat_out, t = t_vec_out))
  })
}


# Likelihood Estimation ---------------------------------------------------

PH_likelihood <- function(data_in, A_mat, alpha_vec) {
  
  out_rates = -Matrix(rowSums(A_mat), ncol = 1)
  likelihood = tibble(x = test_data) %>% 
    rowwise() %>% 
    transmute(temp = log(alpha_vec %*% expm(A_mat * x) %*% out_rates)[1]) %>% 
    sum()
  
  return(likelihood)
  
}

out_rates = -Matrix(rowSums(set_A), ncol = 1)
test_likelihood = tibble(x = test_data) %>% 
  rowwise() %>% 
  mutate(temp = log(set_alpha %*% expm(set_A * x) %*% out_rates)[1]) %>% 
  select(temp) %>% 
  sum()


# Iterate the EM algorithm ------------------------------------------------
num_iterations = 100

# A_mat = set_A
# alpha_vec = set_alpha
# 
# likelihood = rep(-Inf, num_iterations)
# 
# for (iter in 1:num_iterations) {
#   print(paste0("Iteration number # ", iter, " out of ", num_iterations))
#   expect = PH_expect_func(test_data, A_mat, alpha_vec)
#   update = PH_update_func(expect, length(test_data))
#   
#   A_mat = update$A
#   alpha_vec = update$alpha
#   
#   likelihood[iter] = PH_likelihood(direct_samps, A_mat, alpha_vec)
# }

system.time(test <- PH_update_func(PH_expect_func(test_data[1:100], A_mat, alpha_vec), length(test_data)))

