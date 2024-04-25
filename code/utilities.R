#########################  Utility functions  ##################################

# Initialized: April 2024

# Maintained by Kyle Dahlin (kydahlin@gmail.com)

# ---------------------------------------------------------------------------- #

# Misc. helper functions ----
# Helper function: place legends in empty facets of plot grids
# Code by: Artem Sokolov, found here: https://stackoverflow.com/questions/54438495/shift-legend-into-empty-facets-of-a-faceted-plot-in-ggplot2
shift_legend <- function(p) {
  pnls <- cowplot::plot_to_gtable(p) %>%
    gtable::gtable_filter("panel") %>%
    with(setNames(grobs, layout$name)) %>%
    purrr::keep(~ identical(.x, zeroGrob()))
  
  if (length(pnls) == 0) stop("No empty facets in the plot")
  
  lemon::reposition_legend(p, "center", panel = names(pnls))
}

# Phase-type distribution fitting functions ----

# Function: Get the number of parameters for a phase-type distribution of
#           a given class and order
param_count_func <- function(class, order) {
  p = order
  out <- case_when(
    class == "general" ~ p^2 + p,
    class == "coxian" ~ max(1, 2*p - 2),
    class == "hyperexponential" ~ 2*p,
    class == "gcoxian" ~ 3*p - 2,
    class == "gerlang" ~ p,
    class == "mechanistic" ~ 9 * p / 4 # 9 parameters per 4x4 motif. p must be a multiple of 4
  )
} 

# Function: Calculate expected likelihood of PH distribution parameters to data
PHlikelihood <- function(data, PH_mat) {
  
  dimension = dim(PH_mat)[1]
  
  alpha = matrix(rep(0, dimension), nrow = 1); alpha[1] = 1
  
  out_rates = -PH_mat %*% matrix(rep(1, dimension), nrow = dimension, ncol = 1)
  
  # N = number of observations
  # y_1,...,N = data
  # pi = alpha
  # T = the sub-intensity matrix, out_mat
  # t = -T*e, the row sums of out_mat
  Lik = 0
  for (i in 1:length(data)) {
    sub_prod = log(alpha %*% expm(PH_mat * data[i]) %*% out_rates)
    # print(paste0(i, ", ", sub_prod))
    Lik = Lik + sub_prod
  }
  
  return(Lik[1])
}

# Function: put together likelihood and AIC for each model and data set
model_selection_func <- function(data, PH_mat, class, order) {
  logLik <- PHlikelihood(data, PH_mat)
  temp_param_count = param_count_func(class, order)
  temp_AIC = 2 * temp_param_count - 2 * temp_logLik
}

# Function: Plot phtMCMC outputs nicely (taken directly from PhaseType package)
plot.phtMCMC <- function(x, ...) {
  # Get chain info
  StartEndThin <- attr(x$samples, "mcpar")
  
  # Trace plots
  print(
    ggplot(melt(as.data.frame(x$samples[seq(StartEndThin[1], StartEndThin[2], StartEndThin[3]),]), id.vars = NULL)) +
      geom_line(aes_string(x = "1:length(value)", y = "value")) +
      geom_smooth(aes_string(x = "1:length(value)", y = "value"), method = "glm") +
      #geom_hline(aes(yintercept = value), data = truth, colour = "red", linetype = "dashed") +
      facet_wrap(~variable, scales = "free") +
      #theme_grey(base_family = "serif", base_size = 11) +
      ggtitle("Parameter Traces") + #, plot.title = theme_text(size = 14, face = "bold", family = "serif")) +
      xlab("Iteration") + ylab("Parameter Value")
  )
  
  # Marginal posterior densities
  print(
    ggplot() +
      geom_density(aes_string(x = "value"), melt(as.data.frame(x$samples[seq(StartEndThin[1], StartEndThin[2], StartEndThin[3]),]), id.vars = NULL)) +
      #geom_vline(aes(xintercept = value), data = truth, colour = "red") +
      facet_wrap(~variable, scales = "free") +
      #theme_grey(base_family = "serif", base_size = 11) +
      ggtitle("Marginal Posterior Densities") + #, plot.title = theme_text(size = 14, face = "bold", family = "serif")) +
      xlab("Parameter Value") + ylab("Density")
  )
}

# Function: set up uninformed priors for each parameter of the PH distribution
uninf_priors_func <- function(mean, variance, order) {
  shape_k = mean^2 / variance
  scale_theta = variance / mean
  # Gamma prior: shape hyperparameters (one per matrix entry)
  nu <- rep(shape_k, order^2)
  # Gamma prior: reciprocal scale hyperparameters (one per matrix row)
  zeta <- rep(scale_theta, order)
  
  out <- data.frame(nu = nu, zeta = zeta)
  return(out)
}

# Function: reshape phtMCMC/2 output to matrix form
# Function: Transform samples from phtMCMC2 into Phase-type matrix form
reshape_phtMCMC <- function(res_sample, class) {
  # Get median values from the parameter posterior distributions
  order <- case_when(
    class == "Empirical" ~ sqrt(length(res_sample)),
    # class == "Phenomenological" ~ length(res_sample)/2, # NYI
    class == "Mechanistic" ~ 3
  )
  
  if (abs(order - round(order))>0) {stop("Incorrect number of entries in input")}
  
  if (order > 1) {
    PH_mat <- Matrix(data = 0, nrow = order, ncol = order)
    med_seq = seq(1,order^2, by = order)
    
    if (class == "Mechanistic") {
      PH_mat[1,1] = -(res_sample[3] + res_sample[4]) # -(L1 + L2)
      PH_mat[1,2] = res_sample[3]
      
      PH_mat[2,1] = res_sample[5]
      PH_mat[2,2] = -(res_sample[5] + res_sample[6] + res_sample[7]) # P1 + P2 + P3
      PH_mat[2,3] = res_sample[6]
      
      PH_mat[3,1] = res_sample[1]
      PH_mat[3,3] = -(res_sample[1] + res_sample[2])# G1 + G2
      PH_mat = as.matrix(PH_mat)
    } else {
      out_rates <- res_sample[med_seq]
      off_diags <- res_sample[-med_seq]
      positions <- expand.grid(col = seq(1, order), row = seq(1, order)) %>% 
        filter(row != col)
      
      # Fill in the off-diagonal entries
      for (index in 1:dim(positions)[1]) {
        row = positions$row[index]
        col = positions$col[index]
        PH_mat[row, col] = off_diags[index]
      }
      
      # Fill in the diagonal entries
      PH_mat <- t(matrix(as.numeric(PH_mat), nrow = order, ncol = order, byrow = TRUE))
      
      diags <- -out_rates - rowSums(PH_mat, na.rm = TRUE)
      PH_mat[cbind(seq(1:order), seq(1:order))] <- diags
    }
    
  } else {
    PH_mat <- as.matrix(-res_sample,1,1)
  }
  
  return(PH_mat)
}


# Helper function: Put output of phtMCMC2 in matrix form following our bespoke mechanistic matrix form
mech_mat_func <- function(in_mat, stat_func) {
  # in_mat should be of the form of the output of res$samples[i,]
  
  stat_mat <- in_mat %>% 
    as.data.frame() %>% 
    reshape2::melt() %>% 
    group_by(variable) %>% 
    summarise(mean = stat_func(value)) %>% 
    pivot_wider(names_from = variable, values_from = mean)
  
  G_mat = matrix(c(
    0,  stat_mat$L1, stat_mat$P1, stat_mat$G1, 0,
    stat_mat$Q,0,   stat_mat$P2, stat_mat$G2, 0,
    0,  stat_mat$L2,0,    0,    0,
    0,  0,   stat_mat$P3, 0,    0,
    0,  0,   0,    stat_mat$G3,    0
  ),5)
  row_sums = rowSums(G_mat)
  
  temp_mat = G_mat - diag(row_sums)
  
  # remove absorbing state row and column
  mat_dim = dim(temp_mat)[1]
  
  A_mat = temp_mat[1:(mat_dim-1), 1:(mat_dim-1)]
  
  out_mat = t(A_mat)
}

# Parameter manipulation functions ----
get_beta_func <- function(PH_mat, PH_type, epi_parms) {
  
  order = dim(PH_mat)[1]
  
  if (PH_type == "Mechanistic") {
    # Mechanistic
    # epi_parms = betaP, betaG
    zero_mat = matrix(rep(0, (order + 1)^2), nrow = order + 1, ncol = order + 1)
    betaH = zero_mat; betaV = zero_mat; LambdaH = zero_mat; LambdaV = zero_mat
    betaH[3,3] = epi_parms$betaP
    betaV[4,4] = epi_parms$betaG
    LambdaH[3,3] = -PH_mat[2,2]
    LambdaV[4,4] = -PH_mat[3,3]  
  } else if (PH_type == "Empirical") {
    # Empirical
    # epi_parms = betaH = betaV = 1
    betaH = diag(order)
    betaV = betaH
    LambdaH = diag(-rowSums(PH_mat))
    LambdaV = LambdaH
  } else if (PH_type == "Phenomenological") {
    # Phenomenological 
    # epi_parms = beta1V, ..., betapV, beta1H, ..., betapH,
    #             sigma1, ..., sigmap, c1, d1, ..., cp, dp
    betaV = diag(epi_parms$betaV)
    betaH = diag(epi_parms$betaH)
    LambdaH = diag(-epi_parms$d * epi_parms$sigma * diag(PH_mat))
    LambdaV = diag(-epi_parms$c * epi_parms$sigma * diag(PH_mat))
  }
  
  out <- tibble(betaH = betaH, betaV = betaV, LambdaH = LambdaH, LambdaV = LambdaV)
}


# Quantity of interest calculators ----
R0_calc <- function(params, PH_mat, PH_type) {
  with(as.list(params), {
    
    if (PH_type == "Mechanistic") {
      
      
      # !!! temporary values for testing
      eta = 1/(3 * 1440) #1/1
      mu = 1/(14 * 1440)
      gamma = 1/(3 * 1440)
      gammaH = 1/(4 * 1440)
      muH = 1/(365 * 65 * 1440)
      KH = 10^4
      
      KL = 300
      rhoL = 1/(12 * 1440)
      muL = 1/(20 * 1440)
      varPhi = 3 / 1440
      
      iterate_choice = 884
      test_mat = filter(mech_mat_df, iterate == iterate_choice) %>% 
        pivot_wider(names_from = col_index, values_from = value) %>% 
        select(-c(row_index, logLik, iterate)) %>% 
        as.matrix()
      test_parms = mech_params_df[1,]
      PH_mat = test_mat
      
      lambdaL = mech_params_df[iterate_choice,]$lambda_L
      pL = mech_params_df[iterate_choice,]$p_L
      lambdaP = mech_params_df[iterate_choice,]$lambda_P
      pP = mech_params_df[iterate_choice,]$p_P
      lambdaG = mech_params_df[iterate_choice,]$lambda_G
      pG = mech_params_df[iterate_choice,]$p_G
      f = mech_params_df[iterate_choice,]$f
      
      order = dim(PH_mat)[1]
      A = matrix(0, ncol = order+1, nrow = order+1)
      # !!! do this step separately later
      A[2:(order+1), 2:(order+1)] = PH_mat
      A[1,1] = -lambdaQ; A[1,2] = lambdaQ
      A[2,1] = -rowSums(PH_mat)[1] # f * (1 - pL) * lambdaL # 
      A[3,1] = -rowSums(PH_mat)[2] # f * (1 - pP) * lambdaP # -rowSums(PH_mat)[2]
      A[4,1] = f * (1 - pG) * lambdaG # -rowSums(PH_mat)[3]
      
      alpha = matrix(rep(0, order + 1), ncol = 1); alpha[1] = 1
      one_vec = matrix(rep(1, order + 1), nrow = 1)
      
      # Set up transmission matrices: beta and Lambda
      betas = get_beta_func(PH_mat, PH_type, epi_parms)
      betaH = betas$betaH
      betaV = betas$betaV
      LambdaH = betas$LambdaH
      LambdaV = betas$LambdaV
      
      spec_mat = alpha %*% one_vec %*% t(A)
      M1 = (eta + mu) * diag(order + 1) - t(A) + (gamma / (eta + gamma + mu)) * spec_mat
      M2 = eta * diag(order + 1) - (eta / (eta + gamma + mu)) * (gamma / (gamma + mu)) * spec_mat
      M3 = mu * diag(order + 1) - t(A) + (gamma / (gamma + mu)) * spec_mat
      
      # Pr(complete a gonotrophic cycle)
      tau = as.double(-one_vec %*% t(A) %*% inv(mu * diag(order + 1) - t(A)) %*% alpha)
      # Reproductive number
      R = tau * (varPhi / (gamma + mu)) * (rhoL / ( rhoL + muL)) * (1 - tau * gamma / (gamma + mu))^-1
      # Stable distribution of feeding classes
      B_star = KL * ((R-1)/R) * (rhoL + muL) * (1/varPhi) * (R + varPhi * (rhoL/(rhoL + muL))) * inv(mu * diag(order + 1) - t(A)) %*% alpha
      # B_star = 1 * KH * B_star / (sum(B_star))
      
      Q = inv(M3) %*% M2 %*% inv(M1)
      
      R0 = as.double(one_vec %*% betaH %*% LambdaH %*% Q %*% betaV %*% LambdaV %*% B_star / ((gammaH + muH) * KH))
      (R0)
      
    }
  })
}
