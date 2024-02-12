# Simulated blood-feeding data
# Create simulated data for the timing between oviposition events
# numbers and timings between blood meals
# Initiated: February 2024

# Packages ----
library(tidyverse)
library(progress)
library(matlib) # for a bunch of matrix operations
library(moments) # to calculate moments of distributions
library(Matrix) # to calculate matrix exponentials
library(expm) # to calculate matrix exponentials more efficiently and accurately
library(PhaseTypeR) # to create and simulate phase type distributions

# Parameters ----

# Time steps = minutes
DayToMin = 1440
MaxTime = 30 * DayToMin # maximum time spent blood feeding is maximum lifespan

# Questing
pQ = 1
lQ = 1/(1*DayToMin) # 2 days average questing duration

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
f = 0.8

# Maximum number of partial blood meals
MaxPartial = 5

# Total number of iterates
NumIter = 1E4

# Simulated data Algorithm ----

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

data_df <- data_df %>% mutate(out_days = Out_time / 1440)
write_csv(data_df, "data/sample_data.csv")

# Simulate data directly from phase type distribution ----

# Define phase type distribution

# subintensity matrix
set_A = matrix(c(
  -1.5, 0, 0,
  1.5, -1, 0,
  0, 1, -0.5), ncol = 3) 
set_alpha = c(0.9, 0.1, 0) # initial probability vector
ph = PH(set_A, set_alpha)

# Get 100 random samples from the distribution
direct_samps = rPH(1000, ph)