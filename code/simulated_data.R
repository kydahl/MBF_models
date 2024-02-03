# Simulated blood-feeding data
# Create simulated data for the timing between oviposition events
# numbers and timings between blood meals
# Initiated: February 2024

# Packages ----
library(tidyverse)
library(progress)
library(moments)

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
NumIter = 1E4

# Algorithm ----

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
  
  new_row = list(Iterate = Iter, Exit_code = code, Out_time = t, 
                 time_Q = t_Q, time_L = t_L, time_P = t_P, time_G = t_G)
  
  data_df <- rbind(data_df, new_row)
  
}

# Calculating moments ----

# Remove unrealistic data where blood feeding exceeded maximum lifespan
sample_data <- filter(data_df, Exit_code != "failed")

sample_moments = all.moments(sample_data$Out_time, order.max = 7)
sample_moments_day = all.moments(sample_data$Out_time/(24*60), order.max = 7)

# Visualizations ----

# Histograms
hist_plots <- sample_data %>% 
  pivot_longer(cols = Out_time:time_G) %>% 
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_wrap(~name, scales = "free", nrow = 1)
