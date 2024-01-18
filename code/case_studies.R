# Case studies
# Looking at how between-oviposition waiting time distribution changes with
# numbers and timings between blood meals
# Initiated: October 2023


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(cowplot)
library(deSolve)
library(latex2exp)

source("code/parameterization.R")

# Settings & parameters ----------------------------------------------------

k_max = 5
mu = 1/21
alpha = 0.25
lambda = 3

all_params <- params_func(mu, k_max, alpha, lambda)

# Functions ---------------------------------------------------------------

# Calculate eggs per female per oviposition cycle
EFOC_func <- function(parms) {
  with(as.list(parms), {
    EFOC = (rho_L / (rho_L + mu_L)) * sigma_f / (gamma + mu)
  })
}

# Calculate total prob. of surviving to reach oviposition
tau_func_disrupt <- function(parms) {
  out <- parms %>% 
    mutate(part_prod = ifelse(type == 1, 
                              1,
                              lag((lambdas / (mu + lambdas)) * (1 - probs)))) %>%
    mutate(part_prod2 = cumprod(part_prod)) %>% 
    mutate(part_prod3 = (lambdas / (mu + lambdas)) * probs * part_prod2) %>%  
    mutate(tau = sum(part_prod3)) %>% 
    select(-c(part_prod, part_prod2, part_prod3))
}
tau_func_fate <- function(parms) {
  out <- parms %>% 
    mutate(tau = cumprod(lambdas / (mu + lambdas)))
}

# Calculate avg. num. of oviposition cycles
nu_func <- function(parms) {
  out <- parms %>% 
    mutate(nu = 1 / (1 - (gamma / (mu + gamma)) * tau))
}

# Vector reproductive numbers
ReproNum_func <- function(parms) {
  ReproNum <- parms %>% 
    group_by(case, model) %>% 
    mutate(big_Nu = sum(probs * tau * nu)) %>% 
    mutate(ReproNum = case_when(
      model == "Disrupt" ~ EFOC * tau * nu,
      model == "Fate" ~ EFOC * big_Nu
    ))
}

# Get population sizes at equilibrium
EquiNums_func <- function(parms) {
  with(as.list(parms), {
    k_max = unique(k_max)
    if (unique(model) == "Disrupt") {
      # Immature mosquitoes
      L = (EFOC * tau * nu - 1) * K_L / (EFOC * tau * nu)
      # Ovipositing mosquitoes
      V = rho_L * (EFOC * tau * nu - 1) * K_L / (EFOC * (gamma + mu))
      # Feeding mosquitoes
      B = rep(0,k_max)
      B[1] = unique(rho_L * K_L * (EFOC * tau * nu - 1) / ((EFOC * tau) * (mu + lambdas[1])))
      B_vec = cumprod(lag((1 - probs) * (lambdas / (mu + lambdas)))[-1]) * (rho_L * K_L * (EFOC * tau * nu - 1) / ((EFOC * tau) * (mu + lambdas)))[-1]
      if (k_max > 1) {B[2:k_max] = B_vec[-k_max]}
      B_df <- tibble(type = 1:k_max, feed_count = 1:k_max, B = B)
      out <- cbind(parms, tibble(L = L, V = V)) %>% 
        right_join(B_df, by = "type")
    }
    # Fate model
    else if (unique(model) == "Fate") {
      # Immature mosquitoes
      L = (EFOC * big_Nu - 1) * K_L / (EFOC * big_Nu)
      
      # Ovipositing mosquitoes
      V = probs * nu * tau * rho_L * K_L * (EFOC * big_Nu - 1) / (EFOC * big_Nu * (gamma + mu))
      
      # Feeding mosquitoes
      B_mat = matrix(nrow = k_max, ncol = k_max)
      mus = rep(mu, k_max)
      B_mat[,1] = probs * nu * rho_L * K_L * (EFOC * big_Nu - 1) / ((EFOC * big_Nu) * (mu + lambdas[1]))
      
      for (i in 1:k_max) {
        if (k_max > 1) {
          for (j in 2:k_max) {
            temp_lambda = lambdas[1:j]
            temp_mu = mu[1:j]
            B_mat[i,j] = unique(prod((temp_lambda / (temp_mu + temp_lambda))) * probs[i] * nu[i] * rho_L * K_L * (EFOC * big_Nu - 1) / ((EFOC * big_Nu) * (lambdas[i])))
          }
        }
      }
      B_df <- B_mat %>% as_tibble() %>%
        mutate(type = row_number()) %>% 
        pivot_longer(cols = !"type", names_to = "feed_count", names_prefix = "V",
                     values_to = "B") %>% 
        mutate(feed_count = as.integer(feed_count)) %>% 
        filter(feed_count < type + 1)
      
      out <- cbind(parms, tibble(L = L, V = V)) %>% 
        right_join(B_df, by = "type")
    }
  })
}

# Get biting pressure
BitingPressure_func <- function(parms) {
  with(as.list(parms), {
    if (unique(model) == "Disrupt") {
      biting_pressure = sum(lambdas * B)
    }
    # Fate model
    else if (unique(model) == "Fate") {
      
      # Immature mosquitoes
      L = (EFOC * big_Nu - 1) * K_L / (EFOC * big_Nu)
      
      # Ovipositing mosquitoes
      V = probs * nu * tau * rho_L * K_L * (EFOC * big_Nu - 1) / (EFOC * big_Nu * (gamma + mu))
      
      # Feeding mosquitoes
      B_mat = matrix(nrow = k_max, ncol = k_max)
      mus = rep(mu, k_max)
      B_mat[,1] = probs * nu * rho_L * K_L * (EFOC * big_Nu - 1) / ((EFOC * big_Nu) * (mu + lambdas[1]))
      
      for (i in 1:k_max) {
        for (j in 2:k_max) {
          temp_lambda = lambdas[1:j]
          temp_mu = mu[1:j]
          B_mat[i,j] = unique(prod((temp_lambda / (temp_mu + temp_lambda))) * probs[i] * nu[i] * rho_L * K_L * (EFOC * big_Nu - 1) / ((EFOC * big_Nu) * (lambdas[i])))
        }
      }
      B_df <- B_mat %>% as_tibble() %>%
        mutate(type = row_number()) %>% 
        pivot_longer(cols = V1:V10, names_to = "feed_count", names_prefix = "V",
                     values_to = "B") %>% 
        mutate(feed_count = as.integer(feed_count)) %>% 
        filter(feed_count < type + 1)
      
      out <- cbind(parms, tibble(L = L, V = V)) %>% 
        right_join(B_df, by = "type")
    }
  })
}


# Do equilibrium analyses -------------------------------------------------

# Variables to explore
# Set length of variable vectors
base_resolution = 10

# Maximum number of feedings before reaching repletion
k_max_vec = c(1, 5, 10)

# Intrinsic adult mortality
base_mu = 1/21
mu_res = 100
lf_vec <- 7:28
mu_vec = 1/lf_vec#, length.out = base_resolution)

# Base lambda
hbb_vec = seq(0,72, by = 3)[-1] # hours between blood meals
# hbb = 24*1/base_lambda
base_lambda = 1/3
mu_res = 100 
lambda_vec = 24/hbb_vec

# Base reduction
base_alpha = 0.25
alpha_res = 100
alpha_vec = seq(0, 0.5, by = 0.05)

# Determine maximum allowable values of alpha as a function of k_max
root_func <- function(x, k_max) {1 - 2 * x + x^k_max}

max_alpha_vec2 <- tibble(k_max = k_max_vec) %>%
  rowwise() %>% 
  filter(k_max > 2) %>% 
  do(max_alpha = uniroot(root_func, c(0,1-1E-6), k_max = .$k_max, extendInt = "yes")$root) %>% 
  unlist() %>% 
  as.double()
max_alpha_vec = tibble(k_max = k_max_vec, max_alpha = c(1, max_alpha_vec2))

# Set up dataframe of variables
var_df <- expand.grid(k_max = k_max_vec, alpha = alpha_vec) %>%
  right_join(max_alpha_vec, by = "k_max") %>% 
  filter(alpha < max_alpha) %>%
  select(-max_alpha) %>% 
  expand_grid(mu = mu_vec, lambda = lambda_vec)


###* Reproductive numbers ----

ReproNum_df <- tibble(model = character(), case = character(), k_max = integer(),
                      mu = double(), alpha = double(), lambda = double(),
                      ReproNum = double())

BitingPressure_df <- tibble(model = character(), case = character(), k_max = integer(),
                            mu = double(), alpha = double(), lambda = double(),
                            BitingPressure = double())

options(dplyr.summarise.inform = FALSE)
library(foreach)
library(doParallel) 
cl<-makeCluster(13) 
registerDoParallel(cl)

ReproNum_df <- foreach(i = 1:dim(var_df)[1],
                       .combine = 'rbind',
                       .packages = "tidyverse") %dopar% {
                         print(i)
                         vars = var_df[i,]
                         parms = params_func(vars$mu, vars$k_max, vars$alpha, vars$lambda)
                         
                         
                         temp_params <- rbind(
                           filter(parms, model == "Disrupt", case == 1) %>% tau_func_disrupt(.),
                           filter(parms, model == "Fate", case == 1) %>% tau_func_fate(.),
                           filter(parms, model == "Disrupt", case == 2) %>% tau_func_disrupt(.),
                           filter(parms, model == "Fate", case == 2) %>% tau_func_fate(.),
                           filter(parms, model == "Disrupt", case == 3) %>% tau_func_disrupt(.),
                           filter(parms, model == "Fate", case == 3) %>% tau_func_fate(.)
                         ) %>%
                           nu_func(.) %>% 
                           mutate(EFOC = EFOC_func(.)) %>% 
                           ReproNum_func(.)
                         
                         final_params <- rbind(
                           filter(temp_params, model == "Disrupt", case == 1) %>% EquiNums_func(.),
                           filter(temp_params, model == "Fate", case == 1) %>% EquiNums_func(.),
                           filter(temp_params, model == "Disrupt", case == 2) %>% EquiNums_func(.),
                           filter(temp_params, model == "Fate", case == 2) %>% EquiNums_func(.),
                           filter(temp_params, model == "Disrupt", case == 3) %>% EquiNums_func(.),
                           filter(temp_params, model == "Fate", case == 3) %>% EquiNums_func(.)
                         ) 
                         
                         # temp_BitingPressure_df <- final_params %>% 
                         #   select(model, case, lambdas, B) %>% 
                         #   group_by(model, case) %>% 
                         #   summarise(biting_pressure = sum(lambdas * B))%>% 
                         #   ungroup() %>%
                         #   mutate(k_max = vars$k_max, base_lambda = vars$lambda, mu = vars$mu,
                         #          alpha = vars$alpha) %>%
                         #   select(model, case, k_max, mu, alpha, base_lambda,
                         #          BitingPressure = biting_pressure)
                         #  
                         # temp_BitingPressure_df
                         # BitingPressure_df <- rbind(BitingPressure_df, temp_BitingPressure_df)
                         
                         temp_ReproNum_df <- final_params %>%
                           mutate(k_max = vars$k_max, base_lambda = vars$lambda, mu = vars$mu) %>%
                           select(model, case, k_max, mu, alpha, base_lambda, ReproNum) %>%
                           distinct()
                         temp_ReproNum_df
                         # ReproNum_df <- rbind(ReproNum_df, temp_ReproNum_df)
                         
                       }

saveRDS(ReproNum_df, file = "data/ReproNumbers.rds")
# saveRDS(BitingPressure_df, file = "data/BitingPressure.rds")

###* Equilibrium proportions ----


###* Biting pressure ----


# Early visualization -----------------------------------------------------

###* Biting pressure
BitingPressure_df <- readRDS("data/BitingPressure.rds")


# Comparing biting pressure for different k_max
kmax_compare_plot <- BitingPressure_df %>% 
  group_by(case, model, mu, alpha, base_lambda) %>% 
  mutate(BitePres_diff = (BitingPressure[k_max == 5] / BitingPressure[k_max == 1])) %>% 
  mutate(lf = 1/mu) %>% 
  mutate(hbb = 24*1/base_lambda) %>% 
  mutate(dbb = hbb/24) %>% 
  filter(#k_max %in% c(1,5),
         # hbb == 9,
         # lf == 21,
         alpha %in% c(0.5)) %>% 
  distinct() %>% 
  ggplot(aes(x = lf, y = dbb, z = BitePres_diff)) +
  geom_contour_filled() +
  scale_x_continuous("Average lifespan (days)",
                     expand = c(0,0)) +
  scale_y_continuous("Days between blood meals",
                     expand = c(0,0)) +
  scale_fill_viridis_d("-fold increase in biting pressure\nfrom k_max = 1 to 5", 
                       guide = guide_colorsteps(barheight = 20),
                       # labels = scales::percent
                       ) +
  facet_grid(rows = vars(case),
             cols = vars(model),
             labeller = labeller(
               case = as_labeller(function(x) {unname(TeX(paste("$Case = $", x)))},
                                   default = label_parsed),
               model =  as_labeller(label_parsed),
             ))



BitePres_plot_df <- BitingPressure_df %>% 
  group_by(case, k_max, mu, alpha, base_lambda) %>% 
  mutate(BitePres_diff = (BitingPressure[model == "Disrupt"] - BitingPressure[model == "Fate"])/BitingPressure[model == "Disrupt"]) %>%
  # mutate(BitePres_diff = scales::percent(BitePres_diff,pct)) %>% 
  filter(k_max %in% c(5,10),
         # alpha %in% c(NA, 0.25),
         # base_lambda == 1
  ) %>% 
  mutate(lf = 1/mu) %>% 
  mutate(hbb = 24*1/base_lambda) %>% 
  mutate(dbb = hbb/24) %>% 
  group_by(model, case, base_lambda) %>% 
  arrange(lf, BitingPressure)

Case1_plot <- BitePres_plot_df %>% 
  filter(case == 1) %>% 
  select(-alpha) %>% 
  distinct() %>% 
  ggplot(aes(x = lf, y = dbb, z = BitePres_diff)) +
  geom_contour_filled() +
  scale_x_continuous("Average lifespan (days)",
                     expand = c(0,0)) +
  scale_y_continuous("Days between blood meals",
                     expand = c(0,0)) +
  facet_grid(cols = vars(k_max),
             # rows = vars(alpha), 
             labeller = labeller(
               alpha = as_labeller(function(x) {unname(TeX(paste("$\\alpha = $", x)))},
                                   default = label_parsed),
               k_max =  as_labeller(function(x) {unname(TeX(paste("$k = $", x)))},
                                    default = label_parsed),
             )) +
  scale_fill_viridis_d("Change in biting pressure\nfrom Disrupt to Fate", 
                       guide = guide_colorsteps(barheight = 20),
                       labels = scales::percent)
Case1_plot

Case2_plot <- BitePres_plot_df %>%
  filter(case == 2,
         alpha %in% c(0.25, 0.5),
         # hbb < 24
  ) %>% 
  ggplot(aes(x = lf, y = dbb, z = BitePres_diff)) +
  geom_contour_filled() +
  facet_grid(cols = vars(k_max),
             rows = vars(alpha), 
             labeller = labeller(
               alpha = as_labeller(function(x) {unname(TeX(paste("$\\alpha = $", x)))},
                                   default = label_parsed),
               k_max =  as_labeller(function(x) {unname(TeX(paste("$k = $", x)))},
                                    default = label_parsed),
             )) +
  scale_x_continuous("Average lifespan (days)",
                     expand = c(0,0)) +
  scale_y_continuous("Days between blood meals",
                     expand = c(0,0)) +
  scale_fill_viridis_d("Change in biting pressure\nfrom Disrupt to Fate", 
                       guide = guide_colorsteps(barheight = 20),
                       labels = scales::percent)
Case2_plot


Case3_plot <- BitePres_plot_df %>%
  filter(case == 3,
         alpha %in% c(0.25, 0.5),
         # hbb < 24
         ) %>% 
  ggplot(aes(x = lf, y = dbb, z = BitePres_diff)) +
  geom_contour_filled() +
  facet_grid(cols = vars(k_max),
             rows = vars(alpha), 
             labeller = labeller(
               alpha = as_labeller(function(x) {unname(TeX(paste("$\\alpha = $", x)))},
                                   default = label_parsed),
               k_max =  as_labeller(function(x) {unname(TeX(paste("$k = $", x)))},
                                    default = label_parsed),
             )) +
  scale_x_continuous("Average lifespan (days)",
                     expand = c(0,0)) +
  scale_y_continuous("Days between blood meals",
                     expand = c(0,0)) +
  scale_fill_viridis_d("Change in biting pressure\nfrom Disrupt to Fate", 
                       guide = guide_colorsteps(barheight = 20),
                       labels = scales::percent)
Case3_plot

BitePres_plot <- BitePres_plot_df %>% 
  ggplot(aes(x = lf, y = BitingPressure, color = model,
             fill = model, alpha = base_lambda)) +
  geom_path() +
  facet_wrap(~ case)

BitePres_plot

###* Reproductive numbers
ReproNum_plot_df <- ReproNum_df %>% 
  filter(k_max == 5,
         # base_lambda == 1/3
  ) %>% 
  mutate(lf = 1/mu)

ReproNum_plot <- ReproNum_plot_df %>% 
  ggplot(aes(x = lf, y = ReproNum, color = model,
             fill = model, alpha = base_lambda)) +
  geom_line() +
  facet_wrap(~ case)

ReproNum_plot


# Set up full data frame --------------------------------------------------

temp_params <- rbind(
  filter(all_params, model == "Disrupt", case == 1) %>% tau_func_disrupt(.),
  filter(all_params, model == "Fate", case == 1) %>% tau_func_fate(.),
  filter(all_params, model == "Disrupt", case == 2) %>% tau_func_disrupt(.),
  filter(all_params, model == "Fate", case == 2) %>% tau_func_fate(.),
  filter(all_params, model == "Disrupt", case == 3) %>% tau_func_disrupt(.),
  filter(all_params, model == "Fate", case == 3) %>% tau_func_fate(.)
) %>%
  nu_func(.) %>% 
  mutate(EFOC = EFOC_func(.)) %>% 
  ReproNum_func(.)

final_params <- rbind(
  filter(temp_params, model == "Disrupt", case == 1) %>% EquiNums_func(.),
  filter(temp_params, model == "Fate", case == 1) %>% EquiNums_func(.),
  filter(temp_params, model == "Disrupt", case == 2) %>% EquiNums_func(.),
  filter(temp_params, model == "Fate", case == 2) %>% EquiNums_func(.),
  filter(temp_params, model == "Disrupt", case == 3) %>% EquiNums_func(.),
  filter(temp_params, model == "Fate", case == 3) %>% EquiNums_func(.)
) 

ReproNum_df <- final_params %>% 
  select(model, case, ReproNum) %>% 
  mutate(k_max = k_max) %>% 
  distinct()

BitingPressure_df <- final_params %>% 
  select(model, case, lambdas, B) %>% 
  group_by(model, case) %>% 
  summarise(biting_pressure = sum(lambdas * B))

# Set up dynamical models -------------------------------------------------

model.ODE <- function(max_time, vars, parms) {}

# !!! Sample code from Paul Hurtado
# SEIR with Coxian dwell times in a GLCT framework
SEIRpt.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE =ps[['kE']]
  kI =ps[['kI']]
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  dS = -b*z[1]*Itot
  dE =  b*z[1]*Itot                        *  aE + AEt %*% z[1+(1:kE)]  # t(AE) %*% Evec
  dI =  aI %*% (-OnesE%*%AEt %*% z[1+(1:kE)]) + AIt %*% z[1+kE+(1:kI)] #  t(AI) %*% Ivec
  dR = as.numeric(-OnesI%*%AIt %*% z[1+kE+(1:kI)]) 
  
  return(list(c(dS, as.numeric(dE), as.numeric(dI), dR)))
}

SEIRpt.init <- function(ps) {
  
  #print("Defining ICs, AEt, AIt, OnesE, and OnesI ... ")
  
  # Unpack some parameter values...
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = ps[['kE']] # Number of substates in E, and...
  kI = ps[['kI']] # ... I.
  rhoE=ps[['rhoE']] # fraction advancing to the next "phase" (for Coxian distribution)
  rhoI=ps[['rhoI']]
  
  # These are Coxian distributions framed in a Phase-type distribution context,
  # where vector a = (1 0 ... 0) and matrix A is as follows...
  
  aE = matrix(0,nrow = kE, ncol=1); aE[1] = 1;
  AE = sum(rhoE^(0:(kE-1)))/muE*(diag(rep(-1,kE),kE)); 
  if(kE>1) for(i in 1:(kE-1)) {AE[i,i+1] = rhoE*sum(rhoE^(0:(kE-1)))/muE}
  
  aI = matrix(0,nrow = kI, ncol=1); aI[1] = 1;
  AI = sum(rhoI^(0:(kI-1)))/muI*(diag(rep(-1,kI),kI)); 
  if(kI>1) for(i in 1:(kI-1)) {AI[i,i+1] = rhoI*sum(rhoI^(0:(kI-1)))/muI}
  
  
  # Initial conditions
  PopSize=10000
  z0=numeric(kE+kI+2) # initialize the 1+kE+kI+1 state variables
  z0[2] <- 1/PopSize  # 1 in the initial exposed class 
  z0[1] <- 1-z0[2]    # susceptibles = (PopSize - 1)/PopSize
  
  # Set some global variables...
  aE  <<- aE
  AEt <<- t(AE)
  aI  <<- aI
  AIt <<- t(AI)
  OnesE <<- matrix(1,ncol=kE,nrow=1)
  OnesI <<- matrix(1,ncol=kI,nrow=1)
  ICs <<- z0
}


## The "disruption" model ----

# Define the alpha and A matrices


## The "fate" model ----

# Define the alpha and A matrices



# Run simulations ---------------------------------------------------------

# Vary the maximum number of bites per oviposition cycle, k

## Equilibrium quantities ----
# 1. Reproductive numbers
# 2. Distribution of mosquito types
# 3. "Biting pressure"
# 4. Number of mosquitoes taking a second bite within a "window"
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