# Parameterization
# Looking at how between-oviposition waiting time distribution changes with
# numbers and timings between blood meals
# Initiated: October 2023

# Common parameters ---- 
# Eggs per ovipositing female per cycle
sigma_f = 30
# Larval carrying capacity
K_L = 300
# Larval intrinsic mortality
mu_L = 1/15
# Larval maturation rate
rho_L = 1/7
# Rate to leave oviposition stage
gamma = 1/3

## Variable parameters ----
# Set length of variable vectors
base_resolution = 10

# Intrinsic adult mortality
base_mu = 1/21
mu_res = 100
mu_vec = 1/seq(1 * 7, 4 * 7, length.out = base_resolution)

# Base lambda
base_lambda = 1/3
mu_res = 100 
lambda_vec = 1/seq(1, 3, length.out = base_resolution)

# Base reduction
base_alpha = 0.25
alpha_res = 100
# alpha_vec = seq(0, max_alpha, length.out = base_resolution)[-1]



# Combine all parameters into a single dataframe ----
params_func <- function(mu, k_max, alpha, lambda) {
  
  root_func <- function(x, k_max) {1 - 2 * x + x^k_max}
  
  max_alpha = ifelse(k_max < 3, 1,
                     uniroot(root_func, c(0,1-1E-6), k_max = k_max, extendInt = "yes")$root
                     )
  
  if (alpha > max_alpha) {
    warning(paste0(
      "alpha should be less than ", round(max_alpha, digits = 5))
    )}
  
  # Explore range of k values 
  k_vec = 1:k_max
  
  # Case study parameter sets (defined by values of the lambda and rho vectors)
  # 2x2 grid: row 1 = equidistributed rho values, column 1 = constant lambda
  #           row 2 = decreasing rho values,      column 2 = fixed "total lambda"
  
  #* Parameter set 1: ----
  # All rates equal
  # Equidistributed prob`abilities
  
  # Disruption model
  params_1_disrupt = tibble(lambdas = rep(lambda, k_max), 
                            probs = 1/(1+k_max-k_vec),
                            type = k_vec,
                            model = "Disrupt")
  # Fate model
  params_1_fate = tibble(lambdas = rep(lambda, k_max), 
                         probs = 1/k_max,
                         type = k_vec,
                         model = "Fate")
  
  
  params_1 <- rbind(params_1_disrupt, params_1_fate) %>% 
    mutate(case = 1, alpha = NA)
  
  #* Parameter set 2: ----
  # All rates equal
  # Increased type probabilities reduced by alpha
  
  # Disruption model
  sigma_vec <- 1/(rev(cumsum(alpha^(k_vec-1))))
  
  params_2_disrupt = tibble(lambdas = rep(lambda, k_max), 
                            probs = sigma_vec,
                            type = k_vec,
                            model = "Disrupt")
  # Fate model
  kappa_vec <- alpha^(k_vec - 1) / sum(alpha^(k_vec - 1))
  params_2_fate = tibble(lambdas = rep(lambda, k_max), 
                         probs = kappa_vec,
                         type = k_vec,
                         model = "Fate")
  
  
  params_2 <- rbind(params_2_disrupt, params_2_fate) %>% 
    mutate(case = 2, alpha = alpha)
  
  #* Parameter set 3: ---- 
  # Cycle durations are equal among mosquito types
  # Increased type probabilities reduced by alpha
  
  # Disruption model
  sigma_vec <- 1/(rev(cumsum(alpha^(k_vec-1))))
  temp_sigma <- ifelse(k_max == 1, 1, sigma_vec[-k_max])
  prod_vec <- as.double(cumprod(1 - temp_sigma))
  sigma_check <- c(sigma_vec[1] > 1/2,
                   sigma_vec[-1] > prod_vec / (1 + prod_vec)
  )
  
  lambda_vec  = rep(0, k_max)
  lambda_vec[1] = sigma_vec[1] * lambda
  if (k_max > 1) {
    for (i in 2:k_max) {
      lambda_vec[i] = lambda_vec[i-1] * sigma_vec[i] / (2 * sigma_vec[i-1] - 1)
    }
  }
  
  params_3_disrupt = tibble(lambdas = lambda_vec, 
                            probs = sigma_vec,
                            type = k_vec,
                            model = "Disrupt")
  # Fate model
  kappa_vec <- alpha^(k_vec - 1) / sum(alpha^(k_vec - 1))
  params_3_fate = tibble(lambdas = k_vec * lambda, 
                         probs = kappa_vec,
                         type = k_vec,
                         model = "Fate")
  
  params_3 <- rbind(params_3_disrupt, params_3_fate) %>% 
    mutate(case = 3, alpha = alpha)
  
  # Add in fixed parameters
  out <- rbind(params_1, params_2, params_3) %>% 
    mutate(sigma_f = sigma_f, K_L = K_L, mu_L = mu_L, rho_L = rho_L, 
           gamma = gamma, mu = mu, k_max = k_max)
}