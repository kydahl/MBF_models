################################################################################
# Multiple blood feeding: matrix calculations
################################################################################

# Load packages ---
require(tidyverse)
require(matlib)

# from Jthorpe on stackoverflow :
# https://stackoverflow.com/questions/28974507/efficient-creation-of-tridiagonal-matrices
tridiag <- function(upper, lower, main){
  out <- matrix(0,length(main),length(main))
  diag(out) <- main
  indx <- seq.int(length(upper))
  out[cbind(indx+1,indx)] <- lower
  out[cbind(indx,indx+1)] <- upper
  return(out)
}

# Equilibria calculators ----
DFE_func <- function(params, k) {
  with(as.list(params), {

    S_H <- K_H
    I_H <- 0
    R_H <- 0
    S_B <- vector(mode = "numeric", length = k)
    I_B <- S_B
    # S_B[1] <- Lambda_k * n_G / ((b * k) + mu_k + mu)
    
    
    rho_B_k <- (b * k) / ((b * k) + mu + mu_k)
    rho_W <- gamma_W / (mu + gamma_W)
    
    for (j in 1:k) {
      S_B[j] <- (Lambda_k * n_G * rho_B_k^(j-1)) / ((b * k) + mu + mu_k)
      I_B[j] <- 0
    }
    S_W <- (Lambda_k * n_G * rho_B_k^k) / (mu + gamma_W)
    I_W <- 0
    return(tibble(
      S_H = S_H,
      I_H = I_H,
      R_H = R_H,
      S_B = S_B,
      I_B = I_B,
      S_W = S_W,
      I_W = I_W
    ))
  })
}

M_total_func <- function(params) {
  with(as.list(params), {
    k_vec <- params$k
    M_total <- vector(mode = "numeric", length = max(k_vec))
    for (j in k_vec) {
      k_val <- k_vec[j]
      param_set <- filter(params, k == k_val)
      DFE <- DFE_func(param_set, k_val)
      
      M_total[j] <- sum(DFE$S_B) + DFE$S_W[1]
    }
    return(M_total)
  })
}

# Matrix builders ----

## F, the new infections matrix
F_func <- function(params, k_in) {
  param_set <- filter(params, k == k_in)
  with(as.list(c(param_set)), {
    # Get necessary DFE values
    DFE <- DFE_func(param_set, k_in)
    S_H <- DFE$S_H[1]
    S_B <- DFE$S_B
    
    row_1 <- vector(mode = "numeric", length = k+2)
    column_1 <- vector(mode = "numeric", length = k+2)
    F_mat <- matrix(0, nrow = k+2, ncol = k+2)
    
    # Build k+2 by k+2 matrix 
    for (j in 1:k+2) {
      row_1[j] <- ifelse(j>2, beta_MH * b * k, 0)
      column_1[j] <- ifelse(j>3, beta_HM * (b * k) * S_B[j-3] / K_H, 0)
    }
    column_1[2] <- beta_HM * (b * k) * S_B[k] / K_H
    F_mat[1,] <- row_1
    F_mat[,1] <- column_1
    return(F_mat)
  })
}

## V, the net-rate out matrix
V_func <- function(params, k_in) {
  param_set <- filter(params, k == k_in)
  with(as.list(c(param_set)), {
    V_mat <- matrix(0, nrow = k+2, ncol = k+2)
    
    V_mat[1,1] <- gamma_H + mu_H
    V_mat[2,2] <- gamma_W + mu
    V_mat[2,k+2] <- -b * k
    V_mat[3,2] <- -gamma_W
    V_mat[3,3] <- (b * k) + mu + mu_k
    
    if (k>1) {
    main <- as.vector(rep((b * k) + mu + mu_k, k))
    lower <- as.vector(rep(-(b * k), k-1))
    upper <- as.vector(rep(0, k-1))
    lower_block <- tridiag(upper, lower, main)
    V_mat[3:(k+2), 3:(k+2)] <- lower_block
    }
    return(V_mat)
  })
}


## K, the next-generation matrix
K_func <- function(params, k_in) {
  F_mat <- F_func(params, k_in)
  V_mat <- V_func(params, k_in)
  inv_V <- inv(V_mat)
  
  K <- F_mat %*% inv_V
  return(K)
}

# Calculators ----

## R0, the spectral radius of K
R0_func <- function(params, k_in) {
  K_mat <- K_func(params, k_in)
  K_eigen <- eigen(K_mat, only.values = TRUE)
  R0 <- max(Re(K_eigen$values))
}