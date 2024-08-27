# Calculating equilibria and R0 for model with generalized blood-feeding dynamics
# Initiated: February 2024

# Load libraries and define helper functions ----
library(matrixdist)
library(tidyverse)
library(matlib) # for a bunch of matrix operations
library(cowplot)

# Helper function: extract matrices from a tibble
mat_extract <- function(in_df, col_name, index) {
  if (!is.character(col_name)) {stop("Make sure col_name is a character string!")}
  
  if(!col_name %in% colnames(in_df)) {stop("Incorrect column name")}
  
  if (index > dim(in_df)[1]) {stop("Index out of range of data frame")}
  
  out = in_df[[col_name]][index][[1]]
}

# Helper function: turn a list into a matrix
make_mat <- function(list, nrow, ncol) {
  matrix(list[[1]], nrow = nrow, ncol = ncol)
}

# Helper function: temporarily assign variables using column names from dataframe
cheap_assign <- function(df, envir) {
  for (i in 1:ncol(df)) {
    var_str = paste0(names(df)[i])
    
    assign(paste0(var_str), df[[var_str]][[1]], envir = envir)
  }
}

# Set up data frames ----

# Aquatic stage parameter set
rhoL_vals = 1/12
muL_vals = 1/20
varPhi_vals = 3
KL_vals = 300
types = c("logistic", "constant")

aquatic_params = tibble(rhoL = rhoL_vals, muL = muL_vals, varPhi = varPhi_vals, 
                        KL = KL_vals, type = types)

# Adult stage parameter set
mu_vals = 1/21 #seq(1/31, 1/21, length.out = 10)
gamma_vals = c(1/2)

adult_params = tibble(
  mu = mu_vals,
  gamma = gamma_vals)

# Matrix parameters
# # Define some A and alpha sets for testing for now. These will be loaded in separately later
# test_A2 = matrix(c(-2,1,1,-2), 2, 2)
# test_alpha2 = matrix(c(1,0), 2, 1)
# 
# A_dim = 3
# test_A3 = matrix(c(-3,2,0,1,-3,1,0,1,-3), 3, 3)
# test_alpha3 = matrix(c(1,0,0), 3, 1)
# 
# 
# mats = list(test_A2, test_A3)
# alphas = list(test_alpha2, test_alpha3)
# dims = c(2,3)
# 
# bite_params <- tibble(A = mats, alpha = alphas, dims = dims, class = NA, data_set = NA)


# Epidemiological parameters
# !!! For testing for now. Need to set them later
# betaH = list(diag(1,2),diag(1,3))
# eta = 1
# betaV = list(diag(1,2),diag(1,3))
# 
# Lambda_2 = diag(c(-t(rep(1,2)) %*% test_A2), 2, 2)
# Lambda_3 = diag(c(-t(rep(1,3)) %*% test_A3), 3, 3)
# 
# LambdaH = list(Lambda_2, Lambda_3)
# LambdaV = list(Lambda_2, Lambda_3)

epi_setup_func <- function(dimension, A_mat) {
  
  eta = 1/3 # !!! Temporary. Will be varied in sensitivity analysis later
  if (dimension == 1) {
    betaH = 1
    betaV = 1
    Lambda = -A_mat
  } else {
    
    Lambda_vec = -diag(A_mat)
    betaH_vec = rep(0, dimension)
    # !!! Temporary: assume transmission to host in "fastest" stage
    # betaH_vec[Lambda_vec[1:(dimension-1)] == max(Lambda_vec[1:(dimension-1)])] = 1
    
    
    betaH_vec[1:(dimension - 1)] = 1 
    betaH = diag(betaH_vec)
    betaV_vec = rep(0, dimension)
    # !!! Temporary: assume transmission occurs only when mosquitoes are going to oviposition at a substantial rate
    betaV_vec[-colSums(A_mat)>1E-3] = 1 
    # betaV_vec[dimension] = 1
    betaV = diag(betaV_vec)
  
  # !!! Assume that contact rates are the same as "successful" transition rates through the blood feeding stages
  Lambda = diag(-diag(A_mat))
  
  # Lambda= diag(c(-t(matrix(rep(1, dimension), ncol = 1)) %*% t(A_mat)), dimension, dimension)
  }
  LambdaH = Lambda
  LambdaV = Lambda
  
  out <- tibble(betaH = list(betaH), betaV = list(betaV), eta, 
                LambdaH = list(LambdaH), LambdaV = list(LambdaV))
}

temp_df <- select(bite_params, c(dimension, A, class, data_set)) 
epi_params = tibble(dimension = as.integer(), A = list(), class = as.character(), data_set = as.integer(),
                    betaH = list(), betaV = list(), eta = as.double(), LambdaH = list(), LambdaV = list())

for (i in 1:dim(temp_df)[1]) {
  cur_dim = temp_df$dimension[i]
  cur_A = temp_df$A[i][[1]]
  cur_class = temp_df$class[i]
  cur_dataset = temp_df$data_set[i]
  epi_row = epi_setup_func(cur_dim, cur_A)
  
  epi_params <- add_row(epi_params, class = cur_class, data_set = cur_dataset, dimension = cur_dim, A = list(cur_A), epi_row)
}


# Host parameters
gammaH = 1/7
muH = 1/(365 * 65)
KH = 10^4

# epi_params = tibble(betaH, eta, betaV, LambdaH, LambdaV, gammaH, muH, KH, dims)

mat_params = right_join(bite_params, epi_params, by = c("dimension", "class", "A", "data_set"))

# Combine all data
full_df <- expand_grid(aquatic_params, adult_params, mat_params)


# Functions for calculating outputs ----

# # Function: recruitment rate
# phi_func <- function(params) {
#   with(as.list(params), {
#     # type = params$type
#     # V = params$V
#     # L = params$L
#     # varPhi = params$phi
#     # 
#     out = case_when(
#       type == "logistic" ~ varPhi * V * (1 - (L / KL)),
#       type == "constant" ~ varPhi * V
#     )
#     
#     # if (type == "logistic") {
#     #   KL = params$KL
#     #   out = 
#     # }
#     # if (type == "constant") {
#     #   out = varPhi * V
#     # }
#   })
# }
# 
# # Function: calculate tau from alpha, A, and mu
# tau_func <- function(alpha, A, mu) {
#   p = dim(A)[1]
#   
#   ones_vec = matrix(rep(1, p), ncol = 1)
#   first_prod = -1 %*% t(ones_vec) %*% A
#   
#   print(length(mu))
#   
#   # mu %*% diag(p)
#   
#   # print(t(A))
#   
#   second_prod = inv(mu * diag(p) - t(A))
#   
#   out = as.double(first_prod %*% second_prod %*% alpha)
#   
# }
# 
# # Function: find the mosquito population reproductive number R
# R_func <- function(in_df, A, alpha) {
#   with(as.list(in_df), {
#   
#   diag_mu = make_mat(list_diag_mu, dims, dims)
#   
#   print(diag_mu)
#   
#   p = dim(A)[1]
#   
#   ones_vec = matrix(rep(1, p), ncol = 1)
#   
#   first_prod = -1 %*% t(ones_vec) %*% A
# 
#   second_prod = inv(diag_mu - t(A))
# 
#   tau = as.double(first_prod %*% second_prod %*% alpha_vec)
#   # 
#   # mu
#   
#   # tau = tau_func(alpha, A, mu)
#   # print(tau)
#   
#   # 
#   # phi_V = case_when(
#   #   type == "logistic" ~ varPhi,
#   #   type == "constant" ~ varPhi
#   # )
#   # phi_L = case_when(
#   #   type == "logistic" ~ 0,
#   #   type == "constant" ~ 0
#   # )
#   # 
#   # R = tau * (phi_V / (gamma + mu)) * (rhoL / (rhoL + muL - phi_L)) * ((1 - tau * (gamma / (gamma + mu)))^(-1)) 
#   # 
#   # 
#   })
# }
# 
# 
# 
# 
# # Function: find equilibrium value of L (based on choice of recruitment function)
# L_func <- function(in_df) {
#   if (is.null(in_df$R)) {in_df <- mutate(in_df, R = R_func(.))}
#   
#   index = 1
#   A = mat_extract(in_df, "list_A", index)
#   alpha = mat_extract(in_df, "list_alpha", index)
#   
#   # with(as.list(in_df), {
#   
#   tau = tau_func(alpha, A, mu)
#   phi = phi_func(params)
#   
#   L_star = case_when(
#     type == "logistic" ~ KL * (1 - 1/R),
#     type == "constant" ~ 0
#   )
#   
#   firstRHS = tau * (1 / (gamma + mu)) * (rhoL / (rhoL + muL))
#   invRHS = firstRHS * ((1 - tau * (gamma / (gamma + mu)))^(-1))
#   RHS = 1/invRHS
#   
#   bar_phi = phi / V
#   
#   # })
#   
# }

# Calculations ----

# Function: calculate equilibrium values for a given parameter set
eq_func <- function(in_df) {
  # Get dimensions
  num_rows = dim(in_df)[1]
  max_A_dims = max(in_df$dimension)
  
  # Initialize arrays
  init_vec = rep(NA, num_rows)
  tau = R = omega = nu = theta = init_vec
  L_star = V_star = init_vec
  B_star = matrix(rep(NA, num_rows * max_A_dims), nrow = num_rows)
  # Give informative names to B_star columns
  B_star_colnames = purrr::map_chr(1:max_A_dims, ~paste0("B", .x,"_star"))
  colnames(B_star) = B_star_colnames
  
  calc_df = cbind(in_df, tau, R, L_star, V_star, B_star)
  
  # Calculate equilibrium quantities
  for (i in 1:num_rows) {
    # Load in parameter values from the dataframe and assign them in the 
    # environment for convenient reference
    parms = in_df[i,]
    cheap_assign(parms, environment())
    A = t(A)
    alpha = matrix(alpha, ncol = 1)
    p = dim(A)[1]
    
    # Calculate tau, nu, and theta
    diag_mu = diag(mu, p)
    ones_vec = matrix(rep(1, p), ncol = 1)
    first_prod = -1 %*% t(ones_vec) %*% t(A)
    second_prod = if (p == 1){(diag_mu- t(A))^(-1)} else {inv(diag_mu - t(A))}
    tau[i] = as.double(first_prod %*% second_prod %*% alpha)
    nu[i] = (1 - tau[i] * (gamma / (gamma + mu)))^(-1)
    theta[i] = tau[i] * nu[i]
    
    # Calculate varPhi so that the positive equilibrium exists
    varPhi = case_when(
      type == "logistic" ~ varPhi,
      type == "constant" ~ (theta[i] * (rhoL / (rhoL + muL)) * (1 / (gamma + mu)))^-1
    )
    phi_V = case_when(
      type == "logistic" ~ varPhi,
      type == "constant" ~ varPhi
    )
    phi_L = case_when(
      type == "logistic" ~ 0,
      type == "constant" ~ 0
    )
    
    omega[i] = (rhoL / (rhoL + muL - phi_L)) * phi_V * (1 / (gamma + mu))
    
    R[i] = tau[i] * (phi_V / (gamma + mu)) * (rhoL / (rhoL + muL - phi_L)) * ((1 - tau[i] * (gamma / (gamma + mu)))^(-1))
    
    if (R[i] < (1 - .Machine$double.eps)) {
      L_star[i] = 0
      V_star[i] = 0
      B_star[i,1:dimension] = rep(0, dimension)
    } else {
      L_star[i] = case_when(
        type == "logistic" ~ KL * (R[i] - 1) / (R[i]),
        type == "constant" ~ KL
      )
      
      
      bar_phi = case_when(
        type == "logistic" ~ varPhi * (1 - (L_star[i]) / KL),
        type == "constant" ~ varPhi
      )
      
      V_star[i] = L_star[i] * (rhoL + muL) / bar_phi 
      
      B_star[i,1:dimension] = (rhoL * L_star[i] + gamma * V_star[i]) * second_prod %*% alpha
    }
    
  }
  
  B_star_colnames = purrr::map_chr(1:max(in_df$dimension), ~paste0("B", .x,"_star"))
  colnames(B_star) = B_star_colnames
  
  # Add calculated column values in
  calc_df = cbind(in_df, tau, R, L_star, V_star, B_star) %>% 
    mutate(B_star = rowSums(across(B1_star:B3_star), na.rm = TRUE))
}

# Function: vectorial capacity
veccap_func <- function(in_df) {
  # Get dimensions
  num_rows = dim(in_df)[1]
  max_A_dims = max(in_df$dimension)
  
  # Initialize arrays
  init_vec = rep(NA, num_rows)
  veccap = init_vec
  
  
  out_df = cbind(in_df, veccap)
  
  # Calculate equilibrium quantities
  for (i in 1:num_rows) {
    # Load in parameter values from the dataframe and assign them in the 
    # environment for convenient reference
    parms = in_df[i,]
    cheap_assign(parms, environment())
    # A = t(A)
    alpha = matrix(alpha, ncol = 1)
    p = dim(A)[1]
    
    # Get equilibrium B values
    B_star_colnames = purrr::map_chr(1:max(in_df$dimension), ~paste0("B", .x,"_star"))
    B_star = matrix(unlist(parms[B_star_colnames[1:p]], p, 1))
    
    veccap[i] = sum(betaH %*%LambdaH %*% B_star)
    
  }
  
  # Add calculated column values in
  calc_df = cbind(in_df, veccap)
}

# Function: reproduction numbers and endemic equilibria
repnum_func <- function(in_df) {
  # Get dimensions
  num_rows = dim(in_df)[1]
  max_A_dims = max(in_df$dimension)
  
  # Initialize arrays
  init_vec = rep(NA, num_rows)
  R0 = init_vec
  
  
  out_df = cbind(in_df, R0)
  
  # Calculate equilibrium quantities
  for (i in 1:num_rows) {
    # Load in parameter values from the dataframe and assign them in the 
    # environment for convenient reference
    parms = in_df[i,]
    cheap_assign(parms, environment())
    # A = t(A)
    alpha = matrix(alpha, ncol = 1)
    p = dim(A)[1]
    
    # Get equilibrium B values
    B_star_colnames = purrr::map_chr(1:max(in_df$dimension), ~paste0("B", .x,"_star"))
    B_star = matrix(unlist(parms[B_star_colnames[1:p]], p, 1))
    
    # Calculate intermediate quantities
    ones_vec = matrix(rep(1, p), ncol = 1)
    
    r_mat = alpha %*% t(ones_vec) %*% t(A)
    
    first_prod = t(ones_vec) %*% betaH %*% LambdaH
    
    second_prod = eta * diag(1, p) - (eta / (eta + gamma + mu)) * (gamma / (gamma + mu)) * r_mat
    
    third_quot =  (eta + mu) * diag(1, p) - t(A) + (gamma / (eta + gamma + mu)) * r_mat
    third_prod = if (p == 1) { 1 / third_quot} else {inv(third_quot)}
    
    fourth_quot =  mu * diag(1, p) - t(A) + (gamma / (gamma + mu)) * r_mat
    fourth_prod = if (p == 1) {1 / fourth_quot} else {inv(fourth_quot)}
    
    last_prod = betaV %*% LambdaV %*% B_star / (KH * (gammaH + muH))
    
    
    R0[i] = sqrt(first_prod %*% second_prod %*% third_prod %*% fourth_prod %*% last_prod)
    
  }
  
  # Add calculated column values in
  calc_df = cbind(in_df, R0)
}

# Calculations ----
# calc_df = eq_func(full_df)
# 
# in_df = calc_df
# 
# out_df = repnum_func(in_df)

# Visualizations ----

# # Look at the phi functions
# phi_vals <- aquatic_params %>% 
#   expand_grid(expand_grid(V = seq(1,300), L = seq(1,3000))) %>% 
#   mutate(phi = phi_func(.))
# 
# phi_plot <- phi_vals %>% 
#   filter(phi > 0) %>% 
#   ggplot(aes(x = L, y = V, z = phi)) +
#   geom_contour_filled() +
#   facet_wrap(~ type, scales = "free") +
#   theme_cowplot()
