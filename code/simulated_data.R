# Simulated blood-feeding data
# Create simulated data for the timing between oviposition events
# numbers and timings between blood meals
# Initiated: February 2024

# Packages ----
library(tidyverse)
library(progress)
library(moments) # to calculate moments of distributions
library(Matrix) # to calculate matrix exponentials
library(expm) # to calculate matrix exponentials more efficiently and accurately

# Parameters ----

# Time steps = minutes
DayToMin = 1440
MaxTime = 30 * DayToMin # maximum time spent blood feeding is maximum lifespan

# Questing
pQ = 1
lQ = 1/(2*DayToMin) # 2 days average questing duration

# Landing
pL = 0.75
lL = 1/10 # 10 minutes

# Probing
pP = 0.75
lP = 1/5 # 5 minutes

# Ingesting
pG = 0.75
lG = 1/1 # 1 minutes

# Fleeing
f = 1/2

# Maximum number of partial blood meals
MaxPartial = 5

# Total number of iterates
NumIter = 1E2

# Algorithm ----

# Data to keep track of:
# B_i = number of observed processes initiating in state i
# N_ij = number of observed jumps from state i to state j in ALL the processes
# Z_i = total time spent in state i in ALL the processes
# N_i = the number of processes that exit the absorbing state from state i

# Initialize data frame
data_df <- data.frame(Iterate = as.integer(), Exit_code = as.character(),
                      Out_time = as.integer(), time_Q = as.integer(),
                      time_L = as.integer(), time_P = as.integer(),
                      time_G = as.integer())

# Set up progress bar
iterations <- NumIter
pb <- progress_bar$new(
  format = ":spin :system progress = :percent [:bar] :elapsed | eta: :eta",
  total = iterations,
  width = 100)   
progress <- function(n){
  pb$tick()
}
opts <- list(progress = progress)

for (Iter in 1:NumIter) {
  pb$tick()
  finish_bool = 0
  t = 0
  t_Q = 0
  t_L = 0
  t_P = 0
  t_G = 0
  
  NumPartial = 0
  to_Q = TRUE
  to_L = FALSE
  to_P = FALSE
  to_G = FALSE
  to_V = FALSE
  
  while (!finish_bool) {
    # Questing
    leave_Q = FALSE
    # if (to_Q) {print("Entering Q!")}
    while (!leave_Q & to_Q) {
      t = t + 1
      t_Q = t_Q + 1
      if (runif(1) < lQ) {
        leave_Q = TRUE
        to_L = TRUE}
    }
    
    # if (to_L) {print("Entering L!")}
    # Landing
    leave_L = FALSE
    while (!leave_L & to_L) {
      t = t + 1
      t_L = t_L + 1
      # Check if finished landing
      if (runif(1) < lL) {
        leave_L = TRUE
        # Check if succesful landing
        if (runif(1) < pL) {
          to_P = TRUE
        } else {
          # Check if flee to new host
          if (runif(1) < f) {
            to_L = FALSE
            to_Q = TRUE
          } else {
            to_L = TRUE
          }
        }
      }
    }     
    
    # if (to_P) {print("Entering P!")}
    # Probing
    leave_P = FALSE
    while (!leave_P & to_P) {
      t = t + 1
      t_P = t_P + 1
      # Check if finished probing
      if (runif(1) < lP) {
        leave_P = TRUE
        # Check if successful probing
        if (runif(1) < pP) {
          leave_P = TRUE
          to_G = TRUE
        } else {
          # Check if flee to new host
          if (runif(1) < f) {
            to_L = FALSE
            to_Q = TRUE
          } else {to_L = TRUE}
        }
      }
    }  
    
    # if (to_G) {print("Entering G!")}
    # Ingesting
    leave_G = FALSE
    while (!leave_G & to_G) {
      NumPartial = NumPartial + 1
      t = t + 1
      t_G = t_G + 1
      # Check if finished ingesting
      if (runif(1) < lG) {
        leave_G = TRUE
        # Check if successful probing
        if (runif(1) < pG) {
          leave_G = TRUE
          to_V = TRUE
        } else {
          # Check if flee to new host
          if (runif(1) < f) {
            to_L = FALSE
            to_Q = TRUE
          } else {to_L = TRUE}
        }
      }
    }   
    
    # Finish
    if (t > MaxTime) {
      finish_bool = 1
      # print("Feeding exceeds maximum lifespan")
      code = "fail"
    }
    if (to_V) {
      finish_bool = 1
      # print("Successful blood feeding")
      code = "success"
    }
    if (NumPartial > MaxPartial) {
      finish_bool = 1
      # print("Hit max partial blood meals")
      code = "partial"
    }
  }
  
  new_row = list(Iterate = Iter, Exit_code = code, Out_time = t, partial_feeds = NumPartial,
                 time_Q = t_Q, time_L = t_L, time_P = t_P, time_G = t_G)
  
  data_df <- rbind(data_df, new_row)
  
}

# Calculating moments ----

# Remove unrealistic data where blood feeding exceeded maximum lifespan
sample_data <- filter(data_df, Exit_code != "failed")

sample_moments = all.moments(sample_data$Out_time, order.max = 7)
sample_moments_day = all.moments(sample_data$Out_time/(24*60), order.max = 7)

write_csv(sample_data, "data/sample_data.csv")

# Fitting distributions ----
# Fitting is done through the "EM algorithm"

# Function: make unit basis vectors
make_basis <- function(basis, dimension) {t(t(replace(numeric(dimension), basis, 1)))}

# Define data used for fitting

# 0. Initialize estimates alpha and A

# Initial estimates of alpha and A are determined by the model being considered
# because if any entry of alpha or A is set to zero, it will remain zero 
# at all points in the algorithm. Note also that implicit in this definition
# is the dimension 'p' of the phase type distribution

# Define initial estimates for alpha and A corresponding to different sub-models

p = 5
matrix_dimension = p

# Toy values for now:
(A = matrix(rexp(p^2, rate = 0.1), ncol = p) * 10^(-ceiling(log10(max(Y)))-2))
# A = t(Matrix(c(-17,1,2,3,4,-26,5,6,7,8,-35,9,10,11,12,-43)/10000, nrow = p, ncol = p))
temp_alpha = runif(p)
alpha =temp_alpha / sum(temp_alpha)
t_vec = - A %*% t(t(rep(1, p)))

# Set up sample vector
Y = sample_data$Out_time
NumSamples = length(Y)

# Start EM algorithm
tolerance = 1e-6
tolCheck = tolerance + 1
runCount = 0


Init_1_array <- rep(NA, NumSamples)
Init_2_array <- array(rep(NA, p * NumSamples), c(p, NumSamples))
Init_3_array <- array(rep(NA, p * p * NumSamples), c(p, p, NumSamples))

while (tolCheck > tolerance) {
  runCount = runCount + 1
  # 1. Calculate the expected values of sufficient statistics given observations 
  
  # Calculate the matrices exp(Ay) and J
  
  # Initialize arrays
  expAy <- Init_3_array
  J <- Init_3_array
  denoms <- Init_1_array
  temp_EB <- Init_2_array
  temp_EZ <- Init_2_array
  temp_ENi <- Init_2_array
  temp_EN <- Init_3_array
  
  for (k in 1:length(Y)) {
    y_val = Y[k]
    bigMat = expm(
      rbind(cbind(A*y_val, t_vec%*%t(alpha)*y_val),
            cbind(0 * A, A*y_val))
    )
    temp_expAy = bigMat[1:p,1:p]
    temp_J = bigMat[1:p,(p+1):(2*p)]
    
    # expAy[1:p,1:p,k] = temp_expAy
    # J[1:p,1:p,k] = temp_J
    
    denom = t(alpha) %*% temp_expAy %*% t_vec
    
    for (i in 1:p) {
      ei = make_basis(i, p)
      Bi_numerator = alpha[i] * t(ei) %*% temp_expAy %*% t_vec
      temp_EB[i,k] = as.double(Bi_numerator / denom)
      
      Zi_numerator = temp_J[i,i]
      temp_EZ[i,k] = as.double(Zi_numerator / denom)
      
      Ni_numerator = t(alpha) %*% temp_expAy %*% ei * t_vec[i]
      temp_ENi[i,k] = as.double(Ni_numerator / denom)
      
      for (j in 1:p) {
        Nij_numerator = A[i,j] * temp_J[j,i]
        
        temp_EN[i,j,k] = as.double(Nij_numerator / denom)
      }
    }  
  }
  
  EB = rowSums(temp_EB)
  EZ = rowSums(temp_EZ)
  EN = rowSums(temp_EN, dims = 2)
  ENi = rowSums(temp_ENi)
  
  # 2. Update initial estimates alpha and A
  
  alpha_hat = EB / NumSamples
  
  tij_hat = diag(1/EZ) %*% EN
  
  ti_hat = diag(1/EZ) %*% ENi
  
  tii_hat = -ti_hat - colSums(tij_hat-diag(diag(tij_hat)))
  
  A_hat = (tij_hat-diag(diag(tij_hat))) + diag(as.vector(tii_hat))
  
  t_hat = - A_hat %*% rep(1, p)
  
  # Check to see if we're within the chosen tolerance
  tolCheck = max(
    max(abs(t_hat-t_vec)), 
    max(abs(A_hat-A)), 
    max(abs(alpha_hat - alpha))
  )
  print(tolCheck)
  print(runCount)
  
  if (any(is.nan(A_hat))) {break}
  A = A_hat
  alpha = alpha_hat
  t_vec = t_hat
  
}


# 3. Return to 1 with updated alpha and A

# Visualizations ----

# Histograms
hist_plots <- sample_data %>% 
  pivot_longer(cols = Out_time:time_G) %>% 
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_wrap(~name, scales = "free", nrow = 1)
