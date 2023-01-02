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
    for (j in 1:k) {
      S_B[j] <- Lambda_M * (1/b) * (rho_b^j) * ( 1 - (rho_W)*(rho_b)^k)^-1
      I_B[j] <- 0
    }
    S_W <- Lambda_M * (rho_b^k) * ( 1 - (rho_W)*(rho_b)^k)^-1 / (gamma_W + mu_M)
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

# Matrix builders ----

## F, the new infections matrix
F_func <- function(params, k) {
  with(as.list(c(params)), {
    # Get necessary DFE values
    DFE <- DFE_func(params,k)
    S_H <- DFE$S_H[1]
    S_B <- DFE$S_B
    
    row_1 <- vector(mode = "numeric", length = k+2)
    column_1 <- vector(mode = "numeric", length = k+2)
    F_mat <- matrix(0, nrow = k+2, ncol = k+2)
    
    # Build k+2 by k+2 matrix 
    for (j in 1:k+2) {
      row_1[j] <- ifelse(j>2, beta_MH*b*S_H / K_H, 0)
      column_1[j] <- ifelse(j>3, beta_HM * b * S_B[j -3] / K_H, 0)
    }
    column_1[2] <- beta_HM * b * S_B[k] / K_H
    F_mat[1,] <- row_1
    F_mat[,1] <- column_1
    return(F_mat)
  })
}

## V, the net-rate out matrix
V_func <- function(params, k) {
  with(as.list(c(params)), {
    V_mat <- matrix(0, nrow = k+2, ncol = k+2)
    
    V_mat[1,1] <- gamma_H + mu_H
    V_mat[2,2] <- gamma_W + mu_M
    V_mat[2,k+2] <- -b
    V_mat[3,2] <- -gamma_W
    V_mat[3,3] <- b + mu_M
    
    if (k>1) {
    main <- as.vector(rep(b + mu_M, k))
    lower <- as.vector(rep(-b, k-1))
    upper <- as.vector(rep(0, k-1))
    lower_block <- tridiag(upper, lower, main)
    V_mat[3:(k+2), 3:(k+2)] <- lower_block
    }
    return(V_mat)
  })
}

## V^-1
Vinv_func <- function(params, k) {inv(V)}

## K, the next-generation matrix
K_func <- function(params, k) {
  F_mat <- F_func(params,k)
  V_mat <- V_func(params,k)
  inv_V <- inv(V_mat)
  
  K <- F_mat %*% inv_V
  return(K)
}

# Calculators ----

## R0, the spectral radius of K
R0_func <- function(params, k) {
  K_mat <- K_func(params, k)
  K_eigen <- eigen(K_mat, only.values = TRUE)
  R0 <- max(Re(K_eigen$values))
}