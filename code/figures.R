################################################################################
# Multiple blood feeding: Plots
################################################################################

# Load libraries & packages ----
library(ggplot2)
require(deSolve)
require(reshape2)
require(tidyverse)
require(cowplot)
require(pracma)
require(gridExtra)
require(RColorBrewer)
require(MetBrewer)
require(scales)
source("~/GitHub/MBF_models/code/matrix_funcs.R")
source("~/GitHub/MBF_models/code/param_sets.R")

# Distribution across biting stages plot ----

build_DFE_table <- function(params) {
  with(as.list(params), {
    k_max <- max(params$k)
    DFE_table <- tibble(k = vector("numeric", k_max),
                        SW = vector("numeric", k_max))
    for (i in 1:k_max){DFE_table <- mutate(DFE_table, "SB.{i}" := k * i)}
    DFE_table$k <- 1:k_max
    
    for (i in seq(1,k_max)) {
      k_val <- i
      param_set <- filter(params, k == k_val)
      DFE <- DFE_func(param_set, k_val)
      SB_out <- vector("numeric",k_max)
      SB_out[1:k_val] <- as.list(DFE$S_B)
      
      DFE_table[i, 3:(k_max+2)] <- SB_out
      DFE_table[i, 2] <- as.list(DFE$S_W[1])
    }
    return(DFE_table)
  })
}

SB_dist_plot_function <- function(DFE_table) {
SB_dist_plot <- DFE_table %>%
  select(-"SW") %>% 
  melt(id = c("k"),
       value.name = "SB_vals") %>%
  mutate(variable = as.integer(substr(variable, 4,5))) %>% 
  group_by(variable) %>%
  ggplot(aes(x = variable-1, y = SB_vals, color = as.factor(k))) +
  geom_point() +
  geom_line(lwd = 1) +
  # color:
  scale_color_manual(
    name = "# bites necessary\nfor repletion (k)",
    values = met.brewer("Nizami", 10, direction = -1)
  ) +
  theme_minimal_grid(16) +
  scale_x_continuous(limits = c(0,9),
                   breaks = seq(0,9),
                   name = "# of previous feeding attempts") +
  scale_y_continuous(name = "Population size")
  
return(SB_dist_plot)
}

## Plots by case
SB_dist_case1 <- case1_params %>% 
  build_DFE_table() %>% 
  SB_dist_plot_function() +
  ggtitle("Fitness independent of k")

SB_dist_case2 <- case2_params %>% 
  build_DFE_table() %>% 
  SB_dist_plot_function()

SB_dist_case3 <- case3_params %>% 
  build_DFE_table() %>% 
  SB_dist_plot_function()

SB_dist_case4 <- case4_params %>% 
  build_DFE_table() %>% 
  SB_dist_plot_function() +
  ggtitle("Fitness depends on k")

#* SB distribution plot ----

SB_dist_plot <- grid.arrange(SB_dist_case1, SB_dist_case4, ncol = 1)

ggsave("SB_dist.png",
       plot = SB_dist_plot,
       width = 16,
       height = 9)

# R0 as a function of k plots ----

k_vec <- seq(1,k_max)
R0_table <- tibble(k = k_vec,
                   case1 = numeric(length = k_max),
                   case2 = numeric(length = k_max),
                   case3 = numeric(length = k_max),
                   case4 = numeric(length = k_max)
                   )
for (i in k_vec) {
  case1_params <- build_case1(fixed_params, k_vec[i])
  case2_params <- build_case2(fixed_params, k_vec[i])
  case3_params <- build_case3(fixed_params, k_vec[i])
  case4_params <- build_params(fixed_params, k_vec[i])
  
  R0_table$case1[i] <- R0_func(case1_params, k_vec[i])
  R0_table$case2[i] <- R0_func(case2_params, k_vec[i])
  R0_table$case3[i] <- R0_func(case3_params, k_vec[i])
  R0_table$case4[i] <- R0_func(case4_params, k_vec[i])
}


R0_plot <- R0_table %>% 
  melt(id = c("k")) %>% 
  ggplot(aes(x = k, y = value, color = variable)) +
  geom_line(lwd = 1.5) +
  geom_point() +
  # color:
  scale_color_manual(
    name = "k dependence",
    values = met.brewer("Egypt", 4, direction = -1),
    labels = function(x){
      case_when(
        x == "case1" ~ "independent",
        x == "case2" ~ "only mortality & fecundity",
        x == "case3" ~ "only transmission",
        x == "case4" ~ "full"
        )
    }
  ) +
  # x axis
  scale_x_continuous(name = "# bites necessary for repletion (k)", 
                     limits = c(1,k_max),
                     breaks = seq(1,k_max)
                     ) +
  # y axis
  scale_y_continuous(name = "Basic reproduction number") +
  theme_minimal_grid(16)


ggsave("R0_plot.png",
       plot = R0_plot,
       width = 16,
       height = 9)

R0_case12_plot <- R0_table %>% 
  select(-c(case3,case4)) %>% 
  melt(id = c("k")) %>% 
  ggplot(aes(x = k, y = value, color = variable)) +
  geom_line(lwd = 1.5) +
  geom_point() +
  # color:
  scale_color_manual(
    name = "k dependence",
    values = met.brewer("Egypt", 4, direction = -1),
    labels = function(x){
      case_when(
        x == "case1" ~ "independent",
        x == "case2" ~ "only mortality & fecundity",
        x == "case3" ~ "only transmission",
        x == "case4" ~ "full"
      )
    }
  ) +
  # x axis
  scale_x_continuous(name = "# bites necessary for repletion (k)", 
                     limits = c(1,k_max),
                     breaks = seq(1,k_max)
  ) +
  # y axis
  scale_y_continuous(name = "Basic reproduction number") +
  theme_minimal_grid(16)

ggsave("R0_case12_plot.png",
       plot = R0_case12_plot,
       width = 16,
       height = 9)

R0_case34_plot <- R0_table %>% 
  select(-c(case1,case2)) %>% 
  melt(id = c("k")) %>% 
  ggplot(aes(x = k, y = value, color = variable)) +
  geom_line(lwd = 1.5) +
  geom_point() +
  # color:
  scale_color_manual(
    name = "k dependence",
    values = met.brewer("Egypt", 4, direction = -1)[3:4],
    labels = function(x){
      case_when(
        x == "case1" ~ "independent",
        x == "case2" ~ "only mortality & fecundity",
        x == "case3" ~ "only transmission",
        x == "case4" ~ "full"
      )
    }
  ) +
  # x axis
  scale_x_continuous(name = "# bites necessary for repletion (k)", 
                     limits = c(1,k_max),
                     breaks = seq(1,k_max)
  ) +
  # y axis
  scale_y_continuous(name = "Basic reproduction number") +
  theme_minimal_grid(16)

ggsave("R0_case34_plot.png",
       plot = R0_case34_plot,
       width = 16,
       height = 9)

#* Case 1: R0 as function of k ----
#         parameters independent of k

#* Case 2: R0 as function of k ----
#         mortality & fecundity depend on k

#* Case 3: R0 as function of k ----
#         transmission depends on k

#* Case 4: R0 as function of k ----
#         mortality & fecundity & transmisison depend on k