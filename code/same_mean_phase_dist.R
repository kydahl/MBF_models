# Vary questing duration (flee probability) while fixing mean blood-feeding stage duration
library(tidyverse)
library(cowplot)
source("code/utilities.R")
library(actuar)
library(matlib)
library(latex2exp) # to display TeX
library(cols4all)

# Set up parameters ----

# Baseline parameters

# Define phase type distribution
# Questing
pQ = 1
lQ = 1 / (8 * 60) # 1 / (3 * 1440) # 8 hours = 480 minutes

# Landing
pL = 0.5 # 0.7 # 0.5 #0.3
lL = 0.1 # 10 minutes

# Probing
pP = 0.5 # 0.8 # 0.5 #0.2
lP = 0.2 # 5 minutes

# Ingesting
pG = 0.5 # 0.9 # 0.5 #0.75
lG = 1 # 1 minutes

# Fleeing
f = 0.5 # 0.66 # 0.5 #0.8

# Average total duration of blood-feeding stage
base_theta = 1440 * 3
theta = base_theta # three days for feeding, five days for ovipositing, and two days of resting

gammaV = 1/(0.5 * 1440) # exit rate from oviposition to resting, including bloodmeal digestion and site search (5 days)
gammaR = 1/(0.5 * 1440) # exit rate from resting to return to blood-feeding (2 days)

parameters = tibble(lQ, pL, lL, pP, lP, pG, lG, f, gammaV, gammaR, theta)

# Transmission parameters
betaH = betaV = betaP = betaG = 1
eta = 1/(6 * 1440) # 6 days for infection to develop in vector
mu = 1/(20 * 1440) # 20 day lifespan for vector
# gamma = 1/(2 * 1440) # rate of return from oviposition to blood-feeding (3 days)
gammaH = 1/(7 * 1440) # rate of recovery in hosts (7 days)
muH = 1/(365.25 * 65 * 1440) # host mortality rate (65 years)
KH = 1E7 # host population density

KL = 3 * 1E9
rhoL = 1/(12 * 1440) # 12 day larval development period
muL = 1/ (20 * 1440) # 62.5% probability of larval survival (20 / (20 + 12))
varPhi = 3 / 1440 # on average 3 eggs per female per day

extra_parameters = tibble(betaH, betaV, betaP, betaG, eta, mu, gammaH, muH, KH, KL, rhoL, muL, varPhi)

(GCD_func(cbind(parameters, extra_parameters))/1440)

1/(GCD_func(cbind(parameters, extra_parameters))/1440)


R_calc_function_mech(cbind(parameters, extra_parameters))
R_calc_function_exp(cbind(parameters, extra_parameters))

# Compare waiting time distributions ----

# Function: transform mechanistic parameters to mechanistic matrix
params_to_mat_func <- function(parameters){
  with(as.list(parameters), {
    # subintensity matrix
    out = matrix(c(
      -lQ,                                        lQ,       0,       0,
      f * (1- pL) * lL, -lL + (1 - f) * (1- pL) * lL, pL * lL,       0,
      f * (1 - pP) * lP,     (1 - f) * (1 - pP) * lP,     -lP, pP * lP,
      f * (1 - pG) * lG,     (1 - f) * (1 - pG) * lG,       0,      -lG
    ), 
    ncol = 4, byrow = TRUE) 
  })
}

params_to_tilde_mat_func <- function(parameters){
  with(as.list(parameters), {
    out = diag(5)
    out[1:4, 1:4] = params_to_mat_func(parameters)
    
    out[4,5] = pG * lG
    out[5,] = c(0, 0, 0, 0, -gammaV)
    return(out)
  })
}

params_to_prime_mat_func <- function(parameters){
  with(as.list(parameters), {
    out = diag(6)
    out[1:5, 1:5] = params_to_tilde_mat_func(parameters)
    out[5,6] = gammaV
    out[6,] = c(0, 0, 0, 0, 0, -gammaR)
    
    return(out)
  })
}

# Function to calculate mean duration (to check that other calculations are correct)
GCD_func <- function(parameters) {
  with(as.list(parameters), {
    GCD = ((1 - (1 - f) * (1 - pL * pP * pG))/lQ + 1 / lL + pL / lP + pL * pP / lG) / (pL * pP * pG) + (1 / gammaV) + (1 / gammaR)
  })
}

# Basic offspring number and reproduction number comparisons ----

# Function: calculate basic offspring and reproductive numbers
R_calc_function_exp <- function(params) {
  with(as.list(params), {
    order = 1
    A = -1 / theta
    b = -A
    
    # Set up transmission matrices: beta and Lambda
    betaH = 1
    betaV = 1
    LambdaH = b
    LambdaV = b
    
    # Basic offspring number
    # N_offspring = tau * (varPhi / (mu)) * (rhoL / (muL + rhoL)) * nG
    N_offspring = (rhoL / (muL + rhoL)) * (varPhi / mu)
    B_star = ((N_offspring - 1) * KL * rhoL / N_offspring) / mu
    # B_star = (N_offspring - 1) * KL * ((rhoL / N_offspring) + (muL + rhoL) / varPhi) * (1/(mu+b))
    
    R0 = b * sqrt(betaV * B_star * (1 / (KH * (muH + gammaH))) * betaH * (eta / (mu + eta)) * (1 / mu))
    
    return(tibble(N_offspring = N_offspring, R0 = R0))
  })
}

R_calc_function_mech <- function(params) {
  with(as.list(params), {
    # # !!!
    # list2env(as.list(params), envir = .GlobalEnv)
    # # !!!
    
    order = 4
    A = params_to_mat_func(params)
    A_tilde = params_to_tilde_mat_func(params)
    # A_prime = params_to_prime_mat_func(params)
    
    alpha = matrix(rep(0, order), ncol = 1); alpha[1] = 1
    one_vec = matrix(rep(1, order), nrow = 1)
    
    # Pr(complete a gonotrophic cycle)
    tau = as.double(-one_vec %*% t(A) %*% inv(mu * diag(order) - t(A)) %*% alpha)
    rho = 1 - (gammaV / (mu + gammaV)) * (gammaR / (mu + gammaR)) * tau
    nG = 1 / rho
    # Basic offspring number
    N_offspring = tau * (varPhi / (mu + gammaV)) * (rhoL / ( rhoL + muL)) * nG
    r = (N_offspring - 1) * KL * ((rhoL + muL)/ varPhi) * (mu + gammaV)
    # r = (N_offspring - 1) * KL * nG * tau * rhoL
    
    # Stable distribution of feeding classes
    L_star = KL * (N_offspring - 1) / N_offspring
    # L_star2 = KL * (1 - (tau * varPhi * nG * (rhoL / (muL + rhoL)) / (mu + gammaV))^(-1))
    # L_star2 = KL * (1 - (1 / varPhi) * (tau * (rhoL / (muL + rhoL)) * (1 - tau * (gammaR / (mu + gammaR)) * (gammaV / (mu + gammaV)))^(-1) / (mu + gammaV))^(-1))
    
    # B_star = (N_offspring - 1) * KL * ((rhoL / N_offspring) + (gammaR / (mu + gammaR)) * ((rhoL + muL) / varPhi)) * inv(mu * diag(order) - t(A)) %*% alpha
    B_star = ((N_offspring - 1) * KL * rhoL / N_offspring) * (1 + (gammaR / (mu + gammaR)) * (gammaV / (mu + gammaV)) * nG * tau) * (inv(mu * diag(order) - t(A))) %*% alpha
    
    # V_star = (N_offspring - 1) * KL * (rhoL + muL) / varPhi
    V_star = r / (mu + gammaV)
    
    R_star = r * (gammaV / (mu + gammaV)) / (mu + gammaR)
    
    B_vec = rbind(B_star, V_star)
    
    if (N_offspring < 1) {R0 = 0
    } else {
      
      # Adjust matrix to include oviposition (where new infections can occur)
      alpha = rbind(alpha, 0)
      one_vec = cbind(one_vec, 1)
      
      # Set up transmission matrices: beta and Lambda
      zero_mat = matrix(rep(0, (order+1)^2), nrow = order+1, ncol = order+1)
      betaH = zero_mat; betaV = zero_mat; LambdaH = zero_mat; LambdaV = zero_mat
      betaH[3,3] = 1
      betaV[4,4] = 1
      LambdaH[3,3] = -A[3,3]
      LambdaV[4,4] = -A[4,4]
      
      temp_mat = (-one_vec) %*% t(A_tilde)
      temp_mat[temp_mat < 1e-16] = 0 # remove small numerical errors to ease computation
      
      spec_mat = alpha %*% temp_mat
      M1 = mu * diag(order+1) - (A_tilde) - (gammaR / (mu + gammaR)) * spec_mat
      M2 = eta * diag(order+1) + (eta / (mu + gammaR + eta)) * (1 / (mu + gammaR)) * spec_mat
      M3 = (eta + mu) * diag(order+1) - (A_tilde) - (gammaR / (mu + gammaR + eta)) * spec_mat
      Q = inv(M1) %*% M2 %*% inv(M3)
      
      R0 = sqrt(one_vec %*% betaH %*% LambdaH %*% Q %*% betaV %*% LambdaV %*% B_vec / (KH * (muH + gammaH)))
    }
    
    return(tibble(N_offspring = N_offspring, R0 = R0))
  })
}

# Set up data tables to compute R and R0 across ranges of...

f_vec = seq(0, 1, length.out = 101)

# Vary f
f_vec2 = seq(0, 1, length.out = 1001)
vary_f_parameters = cbind(extra_parameters, parameters) %>%
  dplyr::select(-c(theta, f)) %>%
  cbind(tibble(f = f_vec2)) %>%
  mutate(theta = GCD_func(.)) %>% 
  mutate(type = "Mechanistic") %>%
  mutate(vary_parm = "Seek-new-host probability")

# Vary pL
pL_vec = seq(0, 1, length.out = 1001)
pL_vec = pL_vec[-1]
vary_pL_parameters = cbind(extra_parameters, parameters) %>%
  dplyr::select(-c(theta, pL)) %>%
  cbind(tibble(pL = pL_vec)) %>%
  mutate(theta = GCD_func(.)) %>% 
  mutate(type = "Mechanistic") %>%
  mutate(vary_parm = "Landing success probability")

# Vary pP
pP_vec = seq(0, 1, length.out = 1001)
pP_vec = pP_vec[-1]
vary_pP_parameters = cbind(extra_parameters, parameters) %>%
  dplyr::select(-c(theta, pP)) %>%
  cbind(tibble(pP = pP_vec)) %>%
  mutate(theta = GCD_func(.)) %>% 
  mutate(type = "Mechanistic") %>%
  mutate(vary_parm = "Probing success probability")

# Vary pG
pG_vec = seq(0, 1, length.out = 1001)
pG_vec = pG_vec[-1]
vary_pG_parameters = cbind(extra_parameters, parameters) %>%
  dplyr::select(-c(theta, pG)) %>%
  cbind(tibble(pG = pG_vec)) %>%
  mutate(theta = GCD_func(.)) %>% 
  mutate(type = "Mechanistic") %>%
  mutate(vary_parm = "Ingestion success probability")

# Vary lQ
inv_lQ_vec = seq(0, 100 * 1440, length.out = 10002)
inv_lQ_vec = inv_lQ_vec[-1]
vary_lQ_parameters = cbind(extra_parameters, parameters) %>%
  dplyr::select(-c(theta, lQ)) %>%
  cbind(tibble(lQ = 1/inv_lQ_vec)) %>%
  mutate(theta = GCD_func(.)) %>% 
  mutate(type = "Mechanistic") %>%
  mutate(vary_parm = "Questing rate")

# Vary lL
inv_lL_vec = seq(0, 100 * 1440, length.out = 10002)
inv_lL_vec = inv_lL_vec[-1]
vary_lL_parameters = cbind(extra_parameters, parameters) %>%
  dplyr::select(-c(theta, lL)) %>%
  cbind(tibble(lL = 1/inv_lL_vec)) %>%
  mutate(theta = GCD_func(.)) %>% 
  mutate(type = "Mechanistic") %>%
  mutate(vary_parm = "Landing rate")

# Vary lP
inv_lP_vec = seq(0, 100 * 1440, length.out = 10002)
inv_lP_vec = inv_lP_vec[-1]
vary_lP_parameters = cbind(extra_parameters, parameters) %>%
  dplyr::select(-c(theta, lP)) %>%
  cbind(tibble(lP = 1/inv_lP_vec)) %>%
  mutate(theta = GCD_func(.)) %>% 
  mutate(type = "Mechanistic") %>%
  mutate(vary_parm = "Probing rate")

# Vary lG
inv_lG_vec = seq(0, 100 * 1440, length.out = 10002)
inv_lG_vec = inv_lG_vec[-1]
vary_lG_parameters = cbind(extra_parameters, parameters) %>%
  dplyr::select(-c(theta, lG)) %>%
  cbind(tibble(lG = 1/inv_lG_vec)) %>%
  mutate(theta = GCD_func(.)) %>% 
  mutate(type = "Mechanistic") %>%
  mutate(vary_parm = "Ingesting rate")

all_vary_parameters = rbind(
  vary_lQ_parameters, vary_f_parameters, vary_pL_parameters, vary_pP_parameters, vary_pG_parameters,
  vary_lL_parameters, vary_lP_parameters, vary_lG_parameters
)

all_vary_R_vals = tibble(R0 = as.double(), N_offspring = as.double(), theta = as.double(), type = as.character(), vary_parm = as.character())

for (i in 1:dim(all_vary_parameters)[1]) {
  params = all_vary_parameters[i,]
  Rs = R_calc_function_mech(params)
  N_val = Rs$N_offspring
  R0_val = Rs$R0
  temp_tibble = tibble(
    R0 = R0_val,
    N_offspring = N_val,
    type = params$type,
    vary_parm = params$vary_parm,
    theta = params$theta)
  all_vary_R_vals <- add_row(all_vary_R_vals, temp_tibble)
}


# Set up theta variation for the exponential case
exp_vary_R_vals = tibble(R0 = as.double(), N_offspring = as.double(), theta = as.double(), type = as.character(), vary_parm = as.character())

theta_vec = 1440 * (seq(0, 100, length.out = 101))
theta_vec = theta_vec[-1]

vary_exp_parameters = tibble(
  cbind(extra_parameters,
        dplyr::select(parameters, -theta),
        tibble(theta = theta_vec)),
  type = "Exponential")

for (i in 1:dim(vary_exp_parameters)[1]) {
  params = vary_exp_parameters[i,]
  Rs = R_calc_function_exp(params)
  N_val = Rs$N_offspring
  R0_val = Rs$R0
  temp_tibble = tibble(
    R0 = R0_val,
    N_offspring = N_val,
    type = params$type,
    vary_parm = "Exponential",
    theta = params$theta)
  exp_vary_R_vals <- add_row(exp_vary_R_vals, temp_tibble)
}

param_table = tibble(
  name = c("lambda_Q", "p_L", "lambda_L", "p_P", "lambda_P", "p_G", "lambda_G", "f"),
  Symbol = c("$\\lambda_Q$", "$p_L$", "$\\lambda_L$", "$p_P$", "$\\lambda_P$", "$p_G$", "$\\lambda_G$", "$f$"),
  Description = c("Exit rate from questing stage (per minute)", 
                  "Probability of progressing from landing to probing",
                  "Exit rate from landing stage (per minute)", 
                  "Probability of progressing from probing to ingesting",
                  "Exit rate from probing stage (per minute)",
                  "Probability of progressing from ingesting to ovipositing",
                  "Exit rate from ingestion stage (per minute)",
                  "Probability of seeking a new vertebrate host given feeding failure"
  ),
  Label = c("Questing rate", "Landing success probability", "Landing rate", "Probing success probability", "Probing rate", "Ingestion success probability", "Ingesting rate", "Seek-new-host probability"), 
  Type = c("Rate", "Probability", "Rate", "Probability", "Rate", "Probability", "Rate", "Probability"),
  Prefix = c("Questing", "Landing", "Landing", "Probing", "Probing", "Ingesting",  "Ingesting", "Seeking"), 
  value = c(lQ, pL, lL, pP, lP, pG, lG, f)
)

all_vary_plot_df = all_vary_R_vals %>% 
  rbind(exp_vary_R_vals) %>% 
  filter(!is.na(R0)) %>% 
  filter(!is.na(vary_parm)) %>% 
  left_join(rename(param_table, vary_parm = Label)) %>% 
  rename(Label = vary_parm) %>% 
  mutate(Label = ifelse(is.na(Label), "Exponential", Label)) %>% 
  mutate(Prefix = ifelse(is.na(Prefix), "Exponential", Prefix)) %>% 
  mutate(Type = ifelse(is.na(Type), "Rate", Type)) %>%
  filter(is.finite(R0))

all_vary_plot_df$Label <- factor(all_vary_plot_df$Label, levels = c("Exponential", "Questing rate", "Seek-new-host probability", "Landing success probability", "Landing rate", "Probing success probability", "Probing rate", "Ingestion success probability", "Ingesting rate"))

# Labels

all_vary_plot_df$Prefix = factor(all_vary_plot_df$Prefix, levels = c("Exponential", "Questing", "Seeking", "Landing", "Probing", "Ingesting"))
all_vary_plot_df$Type = factor(all_vary_plot_df$Type, levels = c("Exponential", "Rate", "Probability"))


# Plot R0 vs "GCD" ----
# !!! Make distinction between rates and probabilities as in main RMD file
R0_vary_plot_inverted_allparams = ggplot() +
  geom_line(data = all_vary_plot_df %>% 
              # filter(type != "Exponential") %>%
              # filter(Type == "Probability") %>%
              # mutate(b = 1440 / theta) %>% 
              mutate(b = case_when(
                type == "Exponential" ~ 1440 / (theta),# - (1/gammaV) - (1/gammaR)), 
                type == "Mechanistic" ~ 1440 / (theta) # - (1/gammaV) - (1/gammaR))
              )) %>% 
              # restrict to a feasible domain
              # filter(between(b, 0, 0.5)) %>% 
              # filter(between(R0, 0, 1.1)) %>% 
              as.tibble(),
            aes(x = b, y = R0, 
                group = Label,
                color = Prefix,
                linetype = Type
            ),
            lwd = 1, alpha = 0.9) +
  geom_hline(aes(yintercept = 1), color = "red") +
  scale_linetype_manual("Type:", values = c(1, 2), breaks = c("Rate", "Probability")) +
  scale_x_continuous("``Bites per mosquito per day`` = 1/Gonotrophic cycle duration", 
                     # limits = c(0, 0.33),
                     expand = c(0,0)) +
  scale_y_continuous(name = TeX("Basic reproduction number ($R_0$)"),
                     # limits = c(0, NA),
                     # breaks = seq(0,4),
                     # trans = 'log10',
                     expand = c(0,0)
  ) +
  scale_color_manual("Parameter", values = c("Black", c4a("poly.dark24", 5))) +
  # coord_cartesian(xlim = c(0,0.5)) +
  # coord_cartesian(ylim = c(0, 4.1)) +
  # ggtitle("The response of R0 to 1/OCL depends on the mechanism causing variations in OCL") +
  theme_minimal_hgrid(16)

R0_vary_plot_inverted_allparams

ggsave("figures/R0_model_invert_all_params.png", plot = R0_vary_plot_inverted_allparams, 
       width = 13.333, height = 6.5, units = "in")

# Rel'ships among parameters with fixed GCD ----

# Function to calculate values of questing rate to ensure flee probability is feasible (between zero and one)
lQ_range_func <- function(parameters) {
  with(as.list(parameters), {
    denom = ((theta - (1 / gammaV) - (1 / gammaR)) * pL * pP * pG) - (1 / lL) - (pL / lP) - (pL * pP / lG)
    min_lQ = pL * pP * pG / denom
    max_lQ = 1 / denom
    return(tibble(min_lQ, max_lQ))
  })
}

# Function to calculate flee probability (f) for a fixed value of theta
f_theta_func <- function(parameters) {
  with(as.list(parameters), {
    denom = ((theta - (1 / gammaV) - (1 / gammaR)) * pL * pP * pG) - (1 / lL) - (pL / lP) - (pL * pP / lG)
    f = 1 - ( (1 - lQ * denom) / (1 - pL * pP * pG))
  })
}

# Function to calculate questing rate (lQ) as a function of flee probability (f) when mean duration is fixed
lQ_theta_func <- function(parameters) {
  with(as.list(parameters), {
    denom = ((theta - (1 / gammaV) - (1 / gammaR)) * pL * pP * pG) - (1 / lL) - (pL / lP) - (pL * pP / lG)
    lQ = (1 - (1-f)*(1-pL*pP*pG))/denom
  })
}

# Set feasible range of questing rate so that flee probability is between zero and one
f_vec = seq(0, 1, length.out = 101)
lQ_range = lQ_range_func(parameters)
lQ_vec = seq((1+1E-10) * lQ_range$min_lQ, lQ_range$max_lQ, length.out = 10)

theta_vec = c(10 * 1440, 1440* seq(1, 21, length.out = 101))

varied_parameters = parameters %>% 
  select(-f, -lQ, -theta) %>% 
  cross_join(tibble(f = f_vec)) %>%
  cross_join(tibble(theta = theta_vec)) %>% 
  mutate(lQ = lQ_theta_func(.)) %>% 
  filter(lQ > 0)

# Check that the mean duration is being maintained
check_mean = varied_parameters %>% 
  mutate(mean_duration = GCD_func(.)) %>% 
  mutate(test = abs(theta - mean_duration)) %>% 
  filter(test > 1E-10)
if (dim(check_mean)[1] > 0) {warning("Mean duration is not being fixed!")}

# Plot f as a function of lambda_Q when mean duration is fixed
varied_parameters %>% 
  filter(theta == base_theta) %>%
  mutate(Q_duration = 1 / lQ) %>% 
  ggplot(aes(x = Q_duration/1440, y = f)) +
  geom_line() +
  labs(x = "Questing duration (days)", y = "Fleeing probability") +
  coord_cartesian(xlim = c(0.3, 0.4)) +
  theme_minimal()


# Function: Get waiting time distributions from the mechanistic phase-type model with above parameters
wait_time_func <- function(parameters, f_val, resolution, resting_bool = TRUE) {
  with(as.list(parameters), {
    new_params = parameters %>% 
      dplyr::select(-f, -lQ) %>% 
      mutate(f = f_val) %>% 
      mutate(lQ = lQ_theta_func(.))
    
    if (resting_bool) {
      out_mat = params_to_prime_mat_func(new_params)
      alpha_vec = c(1,0,0,0,0,0)
    } else {
      out_mat = params_to_mat_func(new_params)
      alpha_vec = c(1,0,0,0)
    }
    
    out_times = rphtype(resolution, alpha_vec, out_mat)
    
    out = tibble(x = out_times)
  })
}

# Put together list of waiting time distributions
# Set resolution
compare_res = 1E7

# exp_mat = matrix(c(-1/theta))

big_list = tibble(type = "exp", duration = rphtype(compare_res, c(1), as.matrix(-1/theta))) %>% 
  rbind(tibble(type = "f_10", duration = wait_time_func(parameters, 0.1, compare_res)$x)) %>% 
  rbind(tibble(type = "f_90", duration = wait_time_func(parameters, 0.9, compare_res)$x))

big_list_no_resting = tibble(
  type = "exp", 
  duration = rphtype(compare_res, c(1), as.matrix(-1/1/(theta - 1/gammaV - 1/gammaR)))) %>% 
  rbind(
    tibble(
      type = "f_10", 
      duration = wait_time_func(parameters, 0.1, compare_res, FALSE)$x)
  ) %>% 
  rbind(
    tibble(
      type = "f_90", 
      duration = wait_time_func(parameters, 0.9, compare_res, FALSE)$x))

# Plot the comparison of the distributions
color_labs = c("Mechanistic (f = 0.10)", "Mechanistic (f = 0.90)", "Exponential")


# Approx. PDFs of waiting time distributions
GCD_density_compare_plot = big_list %>% ggplot(aes(x = duration/1440, fill = type)) +
  # Exponential durations
  geom_density(alpha = 0.5) +
  scale_x_continuous("Gonotrophic cycle duration (days)",
                     trans = 'log10'
  ) +
  scale_y_continuous("") +
  scale_fill_manual("Model",
                    values = c("blue","orange", "lightgreen"),
                    breaks = c("f_10", "f_90", "exp"), 
                    labels = color_labs) +
  # coord_cartesian(xlim = c(1E-1, 10^(5))) +
  theme_minimal_grid()

biting_density_compare_plot = big_list_no_resting %>% ggplot(aes(x = duration/1440, fill = type)) +
  # Exponential durations
  geom_density(alpha = 0.5) +
  scale_x_continuous("Blood-feeding stage duration (days)",
                     trans = 'log10'
  ) +
  scale_y_continuous("") +
  scale_fill_manual("Model",
                    values = c("blue","orange", "lightgreen"),
                    breaks = c("f_10", "f_90", "exp"), 
                    labels = color_labs) +
  # coord_cartesian(xlim = c(1E-1, 10^(5))) +
  theme_minimal_grid()

# Save density plots
ggsave("figures/GCD_density_comparison.png", plot = GCD_density_compare_plot, 
       width = 13.333, height = 6.5, units = "in")
ggsave("figures/biting_density_comparison.png", plot = biting_density_compare_plot, 
       width = 13.333, height = 6.5, units = "in")

# Approx. PDFs of waiting time distributions
density_compare_zoom = big_list %>% ggplot(aes(x = duration/1440, fill = type)) +
  # Exponential durations
  geom_density(alpha = 0.5) +
  scale_x_continuous("Gonotrophic cycle duration (days)",
                     trans = 'log10') +
  scale_y_continuous("") +
  scale_fill_manual("Model",
                    values = c("blue","orange", "lightgreen"),
                    breaks = c("f_10", "f_90", "exp"), 
                    labels = color_labs) +
  # coord_cartesian(xlim = c(1E-1, 1E2),
  # ylim = c(0, 0.06)) +
  theme_minimal_grid()

# Alternate visualization: violin or box plots

violin_compare_plot = big_list %>% ggplot(aes(x = as.factor(type), y = duration/1440)) +
  geom_violin(fill = "turquoise") +
  scale_x_discrete("Model", 
                   labels = c(
                     "exp" = "Exponential", 
                     "f_75" = "Mechanistic (f = 0.75)", 
                     "f_25" = "Mechanistic (f = 0.25)")
  ) +
  scale_y_continuous("Gonotrophic cycle duration (days)", trans = 'log10') +
  # coord_cartesian(ylim = c(1E-1, 1E4)) +
  theme_minimal_grid()

ggsave("figures/violin_comparison.png", plot = violin_compare_plot, 
       width = 13.333, height = 6.5, units = "in")



# Old stuff ----
vary_mech_parameters = cbind(extra_parameters, parameters) %>%
  dplyr::select(-c(theta, lQ, f)) %>%
  cbind(tibble(theta = theta_vec)) %>% 
  cross_join(tibble(f = f_vec)) %>%
  # cross_join(tibble(lQ = lQ_vec)) %>% 
  # mutate(f = f_theta_func(.)) %>%  
  mutate(lQ = lQ_theta_func(.)) %>%
  mutate(type = "mech") %>%
  # filter(f > 0 & f<1)
  filter(lQ > 0)

mech_R_vals_list = tibble(R = as.double(), R0 = as.double(), type = as.character(), f = as.double(), lQ = as.double(), theta = as.double())
exp_R_vals_list = tibble(R = as.double(), R0 = as.double(), type = as.character(), f = as.double(), lQ = as.double(), theta = as.double())

for (i in 1:dim(vary_mech_parameters)[1]) {
  params = vary_mech_parameters[i,]
  cur_type = params$type
  R_val = R_calc_function_mech(params)$R
  R0_val = R_calc_function_mech(params)$R0
  temp_tibble = tibble(R = R_val,
                       R0 = R0_val,
                       type = cur_type,
                       f = params$f,
                       lQ = params$lQ,
                       theta = params$theta)
  mech_R_vals_list <- add_row(mech_R_vals_list, temp_tibble)
}
for (i in 1:dim(vary_exp_parameters)[1]) {
  params = vary_exp_parameters[i,]
  cur_type = params$type
  R_val = R_calc_function_exp(params)$R
  
  R0_val = R_calc_function_exp(params)$R0
  temp_tibble = tibble(R = R_val,
                       R0 = R0_val,
                       type = cur_type,
                       f = params$f,
                       lQ = params$lQ,
                       theta = params$theta)
  exp_R_vals_list <- add_row(exp_R_vals_list, temp_tibble)
}

# Plot changes in R or R0 as parameters are varied
R_vals_list = rbind(exp_R_vals_list, mech_R_vals_list)


R0_vary_plot = ggplot() +
  geom_line(data = R_vals_list %>% filter(type == "mech") %>% as.tibble() %>%
              filter(f == 0.5),
            # filter(lQ == min(lQ)),
            aes(x = theta / 1440, y = R0, color = type),
            lwd = 1) +
  geom_line(data = R_vals_list %>% filter(type == "Exponential", f == 0.66) %>% as.tibble(),
            aes(x = theta / 1440, y = R0, color = type),
            lwd = 1) +
  geom_hline(aes(yintercept = 1), color = "red", lwd = 1) +
  scale_x_continuous("Blood-feeding stage duration (days)", expand = c(0,0)) +
  scale_y_continuous(name = TeX("Basic reproduction number, $R_0$"),
                     breaks = seq(1,6), expand = c(0,0)) +
  scale_color_manual("Model", 
                     breaks = c("Exponential", "mech"),
                     values = c("blue", "green"),
                     labels = c("Exponential", "Mechanistic")) +
  coord_cartesian(xlim = c(0,10)) +
  theme_minimal_grid()

R0_vary_plot
ggsave("figures/R0_model_diff.png", plot = R0_vary_plot, 
       width = 13.333, height = 6.5, units = "in")


R0_vary_plot_inverted = ggplot() +
  geom_line(data = R_vals_list %>% 
              filter(type == "mech")  %>%
              mutate(b = 1440 / (theta + 1 / gamma)) %>% 
              as.tibble() %>%
              filter(f == 0.5),
            # filter(lQ == min(lQ)),
            aes(x = b, y = R0, color = type),
            lwd = 1) +
  geom_line(data = R_vals_list %>%
              filter(type == "Exponential", f == 0.66) %>%
              mutate(b = 1440 / (theta + 1 / gamma)) %>% # 1440 /(theta + 2)
              as.tibble(),
            aes(x = b, y = R0, color = type),
            lwd = 1) +
  geom_hline(aes(yintercept = 1), color = "red", lwd = 1) +
  scale_x_continuous("``Bites per mosquito per day``", expand = c(0,0)) +
  scale_y_continuous(name = TeX("Basic reproduction number, $R_0$"),
                     breaks = seq(1,6), expand = c(0,0)) +
  scale_color_manual("Model", 
                     breaks = c("Exponential", "mech"),
                     values = c("blue", "green"),
                     labels = c("Exponential", "Mechanistic")) +
  # coord_cartesian(xlim = c(0, 20)) +
  theme_minimal_grid()

R0_vary_plot_inverted

ggsave("figures/R0_model_invert.png", plot = R0_vary_plot_inverted, 
       width = 13.333, height = 6.5, units = "in")

R0_explain_plot = ggplot() +
  geom_line(data = R_vals_list %>% 
              filter(type == "mech")  %>%
              mutate(b = 1440 / (theta + 1 / gamma)) %>% 
              mutate(dQ = 1/lQ/1440) %>% 
              as.tibble() %>%
              filter(f == 0.5),
            # filter(lQ == min(lQ)),
            aes(x = dQ, y = R0, color = type),
            lwd = 1) +
  # geom_line(data = R_vals_list %>%
  #             filter(type == "exp", f == 0.66) %>%
  #             mutate(b = 1440 / (theta + 1 / gamma)) %>% # 1440 /(theta + 2)
  #             as.tibble(),
  #           aes(x = b, y = R0, color = type),
  #           lwd = 1) +
  geom_hline(aes(yintercept = 1), color = "red", lwd = 1) +
  scale_x_continuous("Questing duration", expand = c(0,0)) +
  scale_y_continuous(name = TeX("Basic reproduction number, $R_0$"),
                     breaks = seq(1,6), expand = c(0,0)) +
  scale_color_manual("Model", 
                     breaks = c("Exponential", "mech"),
                     values = c("blue", "green"),
                     labels = c("Exponential", "Mechanistic")) +
  # coord_cartesian(xlim = c(0, 20)) +
  theme_minimal_grid()

R0_explain_plot

ggsave("figures/R0_explain.png", plot = R0_explain_plot, 
       width = 13.333, height = 6.5, units = "in")


test = R_vals_list %>%
  filter(type == "mech") %>% 
  rbind(R_vals_list %>% filter(type == "Exponential") %>% dplyr::select(-f) %>% cross_join(tibble(f = f_vec))) %>% 
  dplyr::select(-lQ) %>% 
  group_by(theta) %>% 
  pivot_wider(names_from = type, values_from = c(R, R0)) %>% 
  mutate(R_ratio = R0_Exponential / R0_mech)


R0_ratio_plot = R_vals_list %>%
  filter(type == "mech") %>% 
  rbind(R_vals_list %>% filter(type == "Exponential") %>% dplyr::select(-f) %>% cross_join(tibble(f = f_vec))) %>% 
  dplyr::select(-lQ) %>% 
  group_by(theta) %>% 
  pivot_wider(names_from = type, values_from = c(R, R0)) %>% 
  mutate(R_ratio = R0_Exponential / R0_mech) %>% 
  ggplot(aes(x = theta / 1440, y = 1/R_ratio)) +
  geom_line(lwd = 1) +
  scale_x_continuous("Blood-feeding stage duration (days)") +
  coord_cartesian(xlim = c(0,10)) +
  theme_minimal_grid()

R0_ratio_plot


