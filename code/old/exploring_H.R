# H exploration
# Looking at how between-oviposition waiting time distribution changes with K and Gamma_i's 
# Initiated: October 2023


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(cowplot)



# Settings & parameters ----------------------------------------------------

resolution = 10000

# Define K and Gamma_i's
K_prob = 0.75
K = rgeom(resolution, K_prob) + 1

maxK = max(K)

# rates = rep(1/maxK, maxK)
rate_length = 50
rate_min = 1 /(2/24) # 2 hours
rate_max = 1 / (48/24) # 2 days
rates = seq(rate_min, rate_max, length.out = rate_length)
rates[rate_length:maxK] = rate_max
rates = rates[1:maxK]

# Build data frame of sampled values of the Gamma_is
Gammas = tibble(i = 1:maxK, rate = rates) %>% 
  cross_join(tibble(rep = seq(1, resolution))) %>% 
  rowwise() %>% 
  mutate(Gamma_i = rexp(1, rate))


# Calculate H -------------------------------------------------------------

# H = sum_i^K Gamma_i

H_df <- tibble(K = K) %>% 
  mutate(rep = row_number()) %>% 
  right_join(Gammas) %>% 
  group_by(rep) %>% 
  # Get rid of irrelevant Gamma_i values
  filter(i <= K) %>% 
  # Calculate H
  summarise(H = sum(Gamma_i))

mean(H_df$H)
var(H_df$H)
# Visualize distribution of H ---------------------------------------------


H_df %>% 
  ggplot(aes(x = H)) +
  geom_histogram() +
  theme_cowplot()
