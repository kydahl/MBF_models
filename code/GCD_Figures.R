# Vary questing duration (flee probability) while fixing mean blood-feeding stage duration
library(tidyverse)
library(cowplot)
source("code/utilities.R")
library(actuar)
library(matlib)
library(latex2exp) # to display TeX
library(cols4all)
library(doFuture)
library(progressr)
library(matrixcalc)

# Cases to plot ----

# 1. "Persistent" vs. "Flighty" mosquito parameter sets
#     [x] Set up two parameter sets that can be input to functions, label them
# 2. Contact rates taken as IN vs. OUT
#     [] Make the calculation for R0 take this assumption as an input
# 3. Host availability as in Ross-Macdonald or Chitnis
#     [x] Decide on a value for the Chitnis max. host contact rate (100?)
#     [x] Set up distinct parameter sets
#     [] Make the calculation for R0 take this assumption as an input
# 4. Exponential contact rates vs. mechanistic
#     [x] Set up distinct parameter sets
#     [x] Make this an input for GCD function
#     [x] Make this an input for contact rates function
#     [] Make this an input for R0 function
# 5. Plotting GCD vs contact rates or 1/GCD vs R0
#     [] Set up a plotting function where only the x and y axes need to be specified
#     [] Ensure that it labels the parameters nicely


# Load in data ----
combined_df = readRDS(file = "data/GCD_R0_vals.rds")

# Plot Figures ----

chosen_mosquito_type = "persistent"

folder_name = paste0("figures/GCD_figures/",chosen_mosquito_type)

# For now, just look at persistent mosquitoes
plot_df <- combined_df %>% 
  filter(mosquito_type == chosen_mosquito_type) %>% 
  arrange(across(lQ:f)) %>% 
  select(c(mosquito_type:f,GCD, to_host_contact, to_vector_contact, R0)) %>% 
  mutate(mod_GCD = case_when(
    model_type == "Exponential" ~ GCD,
    model_type == "Mechanistic" ~ GCD - (1/gammaV) - (1/gammaR)
  )) %>% 
  mutate(GCD_day = mod_GCD / 1440) %>% 
  filter(GCD_day < 1/mu/1440)

# Dictionary for labeling
param_table = tibble(
  Symbol = c("$\\lambda_Q$", "$p_L$", "$\\lambda_L$", "$p_P$", "$\\lambda_P$", "$p_G$", "$\\lambda_G$", "$f$", "$b$"),
  Description = c("Exit rate from questing stage (per minute)", 
                  "Probability of progressing from landing to probing",
                  "Exit rate from landing stage (per minute)", 
                  "Probability of progressing from probing to ingesting",
                  "Exit rate from probing stage (per minute)",
                  "Probability of progressing from ingesting to ovipositing",
                  "Exit rate from ingestion stage (per minute)",
                  "Probability of seeking a new vertebrate host given feeding failure",
                  "Biting rate (exponential model)"
                  
  ),
  Label = c("Questing rate", "Landing success probability", "Landing rate", "Probing success probability", "Probing rate", "Ingestion success probability", "Ingesting rate", "Seek-new-host probability", "Exp. biting rate"), 
  Type = c("A Rates", "B Probabilities", "A Rates", "B Probabilities", "A Rates", "B Probabilities", "A Rates", "B Probabilities", "A Rates"),
  Prefix = c("Host-seeking", "Landing", "Landing", "Probing", "Probing", "Ingesting",  "Ingesting", "Stick-around", "Exponential"), 
  short_label = c("lQ", "pL", "lL", "pP", "lP", "pG", "lG", "f", "theta")
)

# Add in nice labels for the parameters
plot_df = plot_df %>% 
  left_join(rename(param_table, varied_parameter = short_label))

plot_df$Label <- factor(plot_df$Label, levels = c("Exp. biting rate", "Host detection rate", "Seek-new-host probability", "Landing success probability", "Landing rate", "Probing success probability", "Probing rate", "Ingestion success probability", "Ingesting rate"))

# Labels
plot_df$Prefix = factor(plot_df$Prefix, levels = c("Exponential", "Host-seeking", "Stick-around", "Landing", "Probing", "Ingesting"))
plot_df$Type = factor(plot_df$Type, levels = c("A Rates", "B Probabilities"))


## GCD vs Contact rates figures ----

#### IN Contact rates ----
# 
# GCD_contact_IN_df <- plot_df %>% 
#   filter(transmission_type == "IN") %>%
#   filter(contact_type == "RM",
#          GCD_day > 1/3,
#          GCD_day < 7) %>% 
#   distinct()
# 
# # Calculate midpoints of curves to place an arrow showing the direction of increase of the parameter
# 
# midpoints <- GCD_contact_IN_df %>% 
#   group_by(varied_parameter) %>% 
#   mutate(median_GCD_day = median(GCD_day)) %>% 
#   filter(GCD_day < median_GCD_day)
# 
# ###### To host
# GCD_contact_IN_host_plot <- GCD_contact_IN_df %>% 
#   ggplot(aes(x = GCD_day, y = to_host_contact, group = interaction(Label, Prefix, contact_type), color = Prefix)) +
#   # Curves
#   geom_line(lwd = 1) +
#   # Arrows showing direction of increasing parameter values
#   geom_line(
#     data = GCD_contact_IN_df %>% 
#       group_by(varied_parameter) %>% 
#       mutate(median_GCD_day = median(GCD_day)) %>% 
#       filter(GCD_day > 0 + rnorm(1,0.5)),
#     lwd = 1,
#     arrow = arrow(ends = "first", type = "closed")
#   ) +
#   geom_line(
#     data = GCD_contact_IN_df %>% 
#       group_by(varied_parameter) %>% 
#       mutate(median_GCD_day = median(GCD_day)) %>% 
#       filter(GCD_day > 2 + rnorm(1,0.5)),
#     lwd = 1,
#     arrow = arrow(ends = "first", type = "closed")
#   ) +
#   # geom_line(
#   #   data = GCD_contact_IN_df %>% 
#   #     group_by(varied_parameter) %>% 
#   #     mutate(median_GCD_day = median(GCD_day)) %>% 
#   #     filter(GCD_day > 4 + rnorm(1,0.5)),
#   #   lwd = 1,
#   #   arrow = arrow(ends = "first", type = "closed")
#   # ) +
#   geom_line(
#     data = GCD_contact_IN_df %>% 
#       group_by(varied_parameter) %>% 
#       mutate(median_GCD_day = median(GCD_day)) %>% 
#       filter(GCD_day > 6 + rnorm(1,0.5)),
#     lwd = 1,
#     arrow = arrow(ends = "first", type = "closed")
#   ) +
#   facet_wrap(vars(Type), ncol = 1, scales = "free") +
#   scale_x_continuous("Gonotrophic cycle duration (days)",
#                      # trans = 'log10'
#   ) +
#   scale_y_continuous("Effective contact rate to hosts (bites per human per day)",
#                      # trans = 'log10'
#   ) +
#   scale_color_manual("Parameter", values = c("Black", c4a("poly.dark24", 5))) +
#   ggtitle("To-HOST contact rate as a function of GCD (using IN-rates)") +
#   theme_minimal_grid()
# 
# ggsave(paste0(folder_name, "/GCD_contact_IN_host.png"), GCD_contact_IN_host_plot, 
#        width = 13.333, height = 7.5, units = "in")
# 
# ###### To vector
# GCD_contact_IN_vector_plot <- GCD_contact_IN_df %>% 
#   ggplot(aes(x = GCD_day, y = to_vector_contact, group = interaction(Label, Prefix, contact_type), color = Prefix)) +
#   geom_line(lwd = 1,
#             arrow = arrow(ends = "first", type = "closed"),
#             alpha = 0.5
#   ) +
#   facet_wrap(vars(Type), ncol = 1, scales = "free") +
#   scale_x_continuous("Gonotrophic cycle duration (days)",
#                      # trans = 'log10'
#   ) +
#   scale_y_continuous("Effective contact rate to vectors (bites per mosquito per day)",
#                      # trans = 'log10'
#   ) +
#   scale_color_manual("Parameter", values = c("Black", c4a("poly.dark24", 5))) +
#   ggtitle("To-VECTOR contact rate as a function of GCD (using IN-rates)") +
#   theme_minimal_grid()
# 
# ggsave(paste0(folder_name, "/GCD_contact_IN_vector.png"), GCD_contact_IN_vector_plot, 
#        width = 13.333, height = 7.5, units = "in")

#### OUT Contact rates ----

GCD_contact_OUT_df <- plot_df %>% 
  filter(transmission_type == "OUT") %>% 
  filter(contact_type == "RM",
         GCD_day > 1/3,
         GCD_day < 7
  ) %>% 
  distinct() %>% 
  group_by(varied_parameter) %>% 
  mutate(nudge = rnorm(1, 0.1)) %>% 
  mutate(GCD_nudge = GCD_day + nudge) %>% 
  ungroup()

small_nudge_func = function(sd_in) {1+(rnorm(1, mean = 0, sd = sd_in))}

GCD_contact_OUT_df <- plot_df %>% 
  filter(transmission_type == "OUT",
         contact_type == "RM",
         GCD_day > 1/3,
         GCD_day < 4
  ) %>% 
  distinct() %>% 
  group_by(varied_parameter) %>% 
  mutate(nudge0 = min(GCD_day)) %>% 
  mutate(nudge1 = small_nudge_func(0.1) * (1 * max(GCD_day)/4)) %>% 
  mutate(nudge2 = small_nudge_func(0.1) * (2 * max(GCD_day)/4)) %>% 
  mutate(nudge3 = small_nudge_func(0.1) * (3 * max(GCD_day)/4)) %>% 
  mutate(nudge4 = 4 * max(GCD_day)/4) %>% 
  ungroup()


###### To host
alt_GCD_contact_OUT_host_plot <- GCD_contact_OUT_df %>% 
  ggplot(aes(x = GCD_day, y = to_host_contact, group = interaction(Label, Prefix, contact_type), color = Prefix)) +
  # Curves
  geom_line(lwd = 1) +
  # Arrows showing direction of increasing parameter values
  geom_line(
    data = GCD_contact_OUT_df %>%
      group_by(varied_parameter) %>%
      filter(between(GCD_day, nudge0, nudge0 + 0.05)),
    lwd = 1,
    arrow = arrow(ends = "first", type = "closed")
  ) +
  # geom_line(
  #   data = R0_invGCD_OUT_df %>%
  #     group_by(varied_parameter) %>%
  #     filter(between(invGCD, nudge1, nudge2)),
  #   lwd = 1,
  #   arrow = arrow(ends = "first", type = "closed")
  # ) +
  geom_line(
    data = GCD_contact_OUT_df %>%
      group_by(varied_parameter) %>%
      filter(between(GCD_day, nudge2 - 0.05, nudge2)),
    lwd = 1,
    arrow = arrow(ends = "first", type = "closed")
  ) +
  geom_line(
    data = GCD_contact_OUT_df %>%
      group_by(varied_parameter) %>%
      filter(between(GCD_day, nudge4 - 0.05, nudge4)),
    lwd = 1,
    arrow = arrow(ends = "first", type = "closed")
  ) +
  facet_wrap(vars(Type), ncol = 1, scales = "free") +
  scale_x_continuous("Gonotrophic cycle duration (days)",
                     # trans = 'log10'
  ) +
  scale_y_continuous("Effective contact rate to hosts (bites per human per day)",
                     # trans = 'log10'
  ) +
  scale_color_manual("Parameter", values = c("Black", rev(c4a("poly.dark24", 5)))) +
  # ggtitle("R0 as a function of 1/GCD (using OUT-rates)") +
  theme_minimal_grid() +
  theme(strip.text.x = element_text(hjust = 0))

ggsave(paste0(folder_name, "/alt_GCD_contact_OUT_host.png"), alt_GCD_contact_OUT_host_plot, 
       width = 13.333, height = 7.5, units = "in")


# GCD_contact_OUT_host_plot <- GCD_contact_OUT_df %>% 
#   ggplot(aes(x = GCD_day, y = to_host_contact, group = interaction(Label, Prefix, contact_type), color = Prefix)) +
#   # Curves
#   geom_line(lwd = 1) +
#   # Arrows showing direction of increasing parameter values
#   geom_line(
#     data = GCD_contact_OUT_df %>%
#       group_by(varied_parameter) %>%
#       filter(between(GCD_nudge, 0, 1.1)),
#     lwd = 1,
#     arrow = arrow(ends = "first", type = "closed")
#   ) +
#   geom_line(
#     data = GCD_contact_OUT_df %>% 
#       group_by(varied_parameter) %>% 
#       filter(between(GCD_nudge, 3, 3.1)),
#     lwd = 1,
#     arrow = arrow(ends = "first", type = "closed")
#   ) +
#   # geom_line(
#   #   data = GCD_contact_OUT_df %>% 
#   #     group_by(varied_parameter) %>% 
#   #     filter(GCD_day > 4 + rnorm(1,0.5)),
#   #   lwd = 1,
#   #   arrow = arrow(ends = "first", type = "closed")
#   # ) +
#   geom_line(
#     data = GCD_contact_OUT_df %>% 
#       group_by(varied_parameter) %>% 
#       filter(between(GCD_nudge, 5, 5.1)),
#     lwd = 1,
#     arrow = arrow(ends = "first", type = "closed")
#   ) + 
#   facet_wrap(vars(Type), ncol = 1, scales = "free") +
#   scale_x_continuous("Gonotrophic cycle duration (days)",
#                      # trans = 'log10'
#   ) +
#   scale_y_continuous("Effective contact rate to hosts (bites per human per day)",
#                      # trans = 'log10'
#   ) +
#   scale_color_manual("Parameter", values = c("Black", c4a("poly.dark24", 5))) +
#   ggtitle("To-HOST contact rate as a function of GCD (using OUT-rates)") +
#   theme_minimal_grid()
# 
# ggsave(paste0(folder_name, "/GCD_contact_OUT_host.png"), GCD_contact_OUT_host_plot, 
#        width = 13.333, height = 7.5, units = "in")

###### To vector
GCD_contact_OUT_vector_plot <- GCD_contact_OUT_df %>% 
  ggplot(aes(x = GCD_day, y = to_vector_contact, group = interaction(Label, Prefix, contact_type), color = Prefix)) +
  # Arrows showing direction of increasing parameter values
  geom_line(
    data = GCD_contact_OUT_df %>% 
      group_by(varied_parameter) %>% 
      mutate(median_GCD_day = median(GCD_day)) %>% 
      filter(GCD_day > 0 + rnorm(1,0.5)),
    lwd = 1,
    arrow = arrow(ends = "first", type = "closed")
  ) +
  geom_line(
    data = GCD_contact_OUT_df %>% 
      group_by(varied_parameter) %>% 
      mutate(median_GCD_day = median(GCD_day)) %>% 
      filter(GCD_day > 2 + rnorm(1,0.5)),
    lwd = 1,
    arrow = arrow(ends = "first", type = "closed")
  ) +
  # geom_line(
  #   data = GCD_contact_IN_df %>% 
  #     group_by(varied_parameter) %>% 
  #     mutate(median_GCD_day = median(GCD_day)) %>% 
  #     filter(GCD_day > 4 + rnorm(1,0.5)),
  #   lwd = 1,
  #   arrow = arrow(ends = "first", type = "closed")
  # ) +
  geom_line(
    data = GCD_contact_OUT_df %>% 
      group_by(varied_parameter) %>% 
      mutate(median_GCD_day = median(GCD_day)) %>% 
      filter(GCD_day > 6 + rnorm(1,0.5)),
    lwd = 1,
    arrow = arrow(ends = "first", type = "closed")
  ) + 
  # Curves
  geom_line(lwd = 1) +
  facet_wrap(vars(Type), ncol = 1, scales = "free") +
  scale_x_continuous("Gonotrophic cycle duration (days)",
                     # trans = 'log10'
  ) +
  scale_y_continuous("Effective contact rate to vectors (bites per mosquito per day)",
                     # trans = 'log10'
  ) +
  scale_color_manual("Parameter", values = c("Black", c4a("poly.dark24", 5))) +
  ggtitle("To-VECTOR contact rate as a function of GCD (using OUT-rates)") +
  theme_minimal_grid()

ggsave(paste0(folder_name, "/GCD_contact_OUT_vector.png"), GCD_contact_OUT_vector_plot, 
       width = 13.333, height = 7.5, units = "in")


alt_GCD_contact_OUT_vector_plot <- GCD_contact_OUT_df %>% 
  ggplot(aes(x = GCD_day, y = to_vector_contact, group = interaction(Label, Prefix, contact_type), color = Prefix)) +
  # Curves
  geom_line(lwd = 1) +
  # Arrows showing direction of increasing parameter values
  geom_line(
    data = GCD_contact_OUT_df %>%
      group_by(varied_parameter) %>%
      filter(between(GCD_day, nudge0, nudge0 + 0.05)),
    lwd = 1,
    arrow = arrow(ends = "first", type = "closed")
  ) +
  # geom_line(
  #   data = R0_invGCD_OUT_df %>%
  #     group_by(varied_parameter) %>%
  #     filter(between(invGCD, nudge1, nudge2)),
  #   lwd = 1,
  #   arrow = arrow(ends = "first", type = "closed")
  # ) +
  geom_line(
    data = GCD_contact_OUT_df %>%
      group_by(varied_parameter) %>%
      filter(between(GCD_day, nudge2 - 0.05, nudge2)),
    lwd = 1,
    arrow = arrow(ends = "first", type = "closed")
  ) +
  geom_line(
    data = GCD_contact_OUT_df %>%
      group_by(varied_parameter) %>%
      filter(between(GCD_day, nudge4 - 0.05, nudge4)),
    lwd = 1,
    arrow = arrow(ends = "first", type = "closed")
  ) +
  facet_wrap(vars(Type), ncol = 1, scales = "free") +
  scale_x_continuous("Gonotrophic cycle duration (days)",
                     # trans = 'log10'
  ) +
  scale_y_continuous("Effective contact rate to vectors (bites per mosquito per day)",
                     # trans = 'log10'
  ) +
  scale_color_manual("Parameter", values = c("Black", rev(c4a("poly.dark24", 5)))) +
  # ggtitle("R0 as a function of 1/GCD (using OUT-rates)") +
  theme_minimal_grid() +
  theme(strip.text.x = element_text(hjust = 0))

ggsave(paste0(folder_name, "/alt_GCD_contact_OUT_vector.png"), alt_GCD_contact_OUT_vector_plot, 
       width = 13.333, height = 7.5, units = "in")

## GCD vs R0 figures ----

#### IN Contact rates ----
# 
# R0_GCD_IN_df <- plot_df %>% 
#   filter(transmission_type == "IN") %>% 
#   distinct()
# 
# R0_GCD_IN_plot <- R0_GCD_IN_df %>% 
#   filter(contact_type == "RM") %>% 
#   ggplot(aes(x = GCD_day, y = R0, group = interaction(Label, Prefix, contact_type), color = Prefix)) +
#   geom_line(lwd = 1,
#             arrow = arrow(ends = "first", type = "closed"),
#             alpha = 0.5
#   ) +
#   facet_wrap(vars(Type), ncol = 1, scales = "free") +
#   scale_x_continuous("Gonotrophic cycle duration (days)",
#                      breaks = seq(1,21),
#                      # trans = 'log10'
#   ) +
#   scale_y_continuous("Basic reproduction number",
#                      # trans = 'log10'
#   ) +
#   scale_color_manual("Parameter", values = c("Black", c4a("poly.dark24", 5))) +
#   ggtitle("R0 as a function of GCD (using IN-rates)") +
#   theme_minimal_grid()
# 
# ggsave(paste0(folder_name, "/R0_GCD_IN.png"), R0_GCD_IN_plot, 
#        width = 13.333, height = 7.5, units = "in")

#### OUT Contact rates ----

R0_GCD_OUT_df <- plot_df %>% 
  filter(transmission_type == "OUT") %>% 
  distinct()

R0_GCD_OUT_plot <- R0_GCD_OUT_df %>% 
  filter(contact_type == "RM") %>% 
  ggplot(aes(x = GCD_day, y = R0, group = interaction(Label, Prefix, contact_type), color = Prefix)) +
  geom_line(lwd = 1,
            arrow = arrow(ends = "first", type = "closed"),
            alpha = 0.5
  ) +
  facet_wrap(vars(Type), ncol = 1, scales = "free") +
  scale_x_continuous("Gonotrophic cycle duration (days)",
                     breaks = seq(1,21),
                     # trans = 'log10'
  ) +
  scale_y_continuous("Basic reproduction number",
                     # trans = 'log10'
  ) +
  scale_color_manual("Parameter", values = c("Black", c4a("poly.dark24", 5))) +
  ggtitle("R0 as a function of GCD (using OUT-rates)") +
  theme_minimal_grid()

ggsave(paste0(folder_name, "/R0_GCD_OUT.png"), R0_GCD_OUT_plot, 
       width = 13.333, height = 7.5, units = "in")


## 1/GCD vs R0 figures ----

#### IN Contact rates ----

# R0_invGCD_IN_df <- plot_df %>% 
#   mutate(invGCD = 1/GCD_day) %>% 
#   filter(transmission_type == "IN",
#          # invGCD < 0.4
#   ) %>% 
#   distinct()
# 
# R0_invGCD_IN_plot <- R0_invGCD_IN_df %>% 
#   filter(contact_type == "RM") %>% 
#   ggplot(aes(x = 1/GCD_day, y = R0, group = interaction(Label, Prefix, contact_type), color = Prefix)) +
#   geom_line(lwd = 1,
#             arrow = arrow(ends = "last", type = "closed"),
#             alpha = 0.5
#   ) +
#   facet_wrap(vars(Type), ncol = 1, scales = "free") +
#   scale_x_continuous("1/Gonotrophic cycle duration = ``Biting rate''",
#                      # breaks = seq(1,21),
#                      # trans = 'log10'
#   ) +
#   scale_y_continuous("Basic reproduction number",
#                      # trans = 'log10'
#   ) +
#   scale_color_manual("Parameter", values = c("Black", c4a("poly.dark24", 5))) +
#   ggtitle("R0 as a function of 1/GCD (using IN-rates)") +
#   theme_minimal_grid()
# 
# ggsave(paste0(folder_name, "/R0_invGCD_IN.png"), R0_invGCD_IN_plot, 
#        width = 13.333, height = 7.5, units = "in")

#### OUT Contact rates ----

R0_invGCD_OUT_df <- plot_df %>% 
  mutate(invGCD = 1/GCD_day) %>% 
  filter(transmission_type == "OUT",
         contact_type == "RM",
         invGCD < 2
  ) %>% 
  distinct() %>% 
  group_by(varied_parameter) %>% 
  mutate(nudge = rnorm(1, 0.1)) %>% 
  mutate(invGCD_nudge = invGCD + nudge) %>% 
  ungroup()


R0_invGCD_OUT_plot <- R0_invGCD_OUT_df %>% 
  ggplot(aes(x = invGCD, y = R0, group = interaction(Label, Prefix, contact_type), color = Prefix)) +
  geom_line(lwd = 1,
            arrow = arrow(ends = "last", type = "closed"),
            alpha = 1
  ) +
  facet_wrap(vars(Type), ncol = 1, scales = "free") +
  scale_x_continuous("1/Gonotrophic cycle duration = ``Biting rate''",
                     # breaks = seq(0, 0.4, by = 0.1),
                     # trans = 'log10'
  ) +
  scale_y_continuous("Basic reproduction number",
                     # trans = 'log10'
  ) +
  scale_color_manual("Parameter", values = c("Black", c4a("poly.dark24", 5))) +
  ggtitle("R0 as a function of 1/GCD (using OUT-rates)") +
  theme_minimal_grid() +
  theme(strip.text.x = element_text(hjust = 0))

ggsave(paste0(folder_name, "/R0_invGCD_OUT.png"), R0_invGCD_OUT_plot, 
       width = 13.333, height = 7.5, units = "in")

##### alternative R0 plot

R0_invGCD_OUT_df <- plot_df %>% 
  mutate(invGCD = 1/GCD_day) %>% 
  filter(transmission_type == "OUT",
         contact_type == "RM",
         invGCD < 2.75
  ) %>% 
  distinct() %>% 
  group_by(varied_parameter) %>% 
  mutate(nudge0 = min(invGCD)) %>% 
  mutate(nudge1 = small_nudge_func(0.1) * (1 * max(invGCD)/4)) %>% 
  mutate(nudge2 = small_nudge_func(0.1) * (2 * max(invGCD)/4)) %>% 
  mutate(nudge3 = small_nudge_func(0.1) * (3 * max(invGCD)/4)) %>% 
  mutate(nudge4 = 4 * max(invGCD)/4) %>% 
  ungroup()


# nudge_points = seq(0, range(R0_invGCD_OUT_df$invGCD_nudge)[2], length.out = 4)

alt_R0_invGCD_OUT_plot <- R0_invGCD_OUT_df %>% 
  ggplot(aes(x = invGCD, y = R0, group = interaction(Label, Prefix, contact_type), color = Prefix)) +
  # Curves
  geom_line(lwd = 1) +
  # Arrows showing direction of increasing parameter values
  geom_line(
    data = R0_invGCD_OUT_df %>%
      group_by(varied_parameter) %>%
      filter(between(invGCD, nudge0, nudge0 + 0.01)),
    lwd = 1,
    arrow = arrow(ends = "first", type = "closed")
  ) +
  # geom_line(
  #   data = R0_invGCD_OUT_df %>%
  #     group_by(varied_parameter) %>%
  #     filter(between(invGCD, nudge1, nudge2)),
  #   lwd = 1,
  #   arrow = arrow(ends = "first", type = "closed")
  # ) +
  geom_line(
    data = R0_invGCD_OUT_df %>%
      group_by(varied_parameter) %>%
      filter(between(invGCD, nudge2 - 0.01, nudge2)),
    lwd = 1,
    arrow = arrow(ends = "first", type = "closed")
  ) +
  geom_line(
    data = R0_invGCD_OUT_df %>%
      group_by(varied_parameter) %>%
      filter(between(invGCD, nudge4 - 0.01, nudge4)),
    lwd = 1,
    arrow = arrow(ends = "first", type = "closed")
  ) +
  facet_wrap(vars(Type), ncol = 1, scales = "free") +
  scale_x_continuous("1/Gonotrophic cycle duration = ``Biting rate''",
                     # breaks = seq(0, 3)
                     # breaks = seq(0, 0.4, by = 0.1),
                     # trans = 'log10'
  ) +
  scale_y_continuous("Basic reproduction number",
                     # trans = 'log10'
  ) +
  scale_color_manual("Parameter", values = c("Black", rev(c4a("poly.dark24", 5)))) +
  # ggtitle("R0 as a function of 1/GCD (using OUT-rates)") +
  theme_minimal_grid() +
  theme(strip.text.x = element_text(hjust = 0))

ggsave(paste0(folder_name, "/alt_R0_invGCD_OUT.png"), alt_R0_invGCD_OUT_plot, 
       width = 13.333, height = 7.5, units = "in")


### Plot GCD against each parameter individually ----
var_GCD_plot_rates <- plot_df %>% 
  pivot_longer(cols = pQ:f) %>% 
  filter(name == varied_parameter) %>% 
  filter(Type == unique(plot_df$Type)[1]) %>% 
  ggplot(aes(x = value, y = GCD_day, color = Prefix)) +
  geom_line(lwd = 1) +
  scale_x_log10("") +
  scale_y_continuous("Gonotrophic cycle duration (days)") +
  scale_color_manual("Parameter", values = c("Black", rev(c4a("poly.dark24", 5)))) +
  facet_wrap(vars(Label), scales = 'free', nrow = 1) +
  theme_minimal_grid()

var_GCD_plot_probs <- plot_df %>% 
  pivot_longer(cols = pQ:f) %>% 
  filter(name == varied_parameter) %>% 
  filter(Type == unique(plot_df$Type)[2]) %>% 
  ggplot(aes(x = value, y = GCD_day, color = Prefix)) +
  geom_line(lwd = 1) +
  scale_x_log10("Parameter value") +
  scale_y_continuous("Gonotrophic cycle duration (days)") +
  scale_color_manual("Parameter", values = c("Black", rev(c4a("poly.dark24", 5))), guide = "none") +
  facet_wrap(vars(Label), scales = 'free', nrow = 1) +
  theme_minimal_grid()

plot_grid(var_GCD_plot_rates, var_GCD_plot_probs, nrow = 2)

### Plot R0 against each parameter individually ----
coeff = 10

var_R0_plot_df = plot_df %>%
  pivot_longer(cols = pQ:f, names_to = "variable", values_to = "variable_value") %>% 
  filter(variable == varied_parameter,
         transmission_type == "OUT") %>% 
  select(contact_type, R0, Label, Prefix, Type, variable, variable_value) %>% 
  pivot_wider(names_from = contact_type, values_from = R0)

var_R0_plot_rates <- var_R0_plot_df %>%
  filter(Type == unique(plot_df$Type)[1]) %>%
  ggplot(aes(x = variable_value, color = Prefix)) +
  # Ross-Macdonald
  geom_line(aes(y = RM), lwd = 1, lty = 1) +
  # Chitnis
  geom_line(aes(y = Chitnis * coeff), lwd = 1, lty = 2) +
  scale_x_log10("Parameter value") +
  scale_y_continuous("Basic reproduction number\n(Ross-Macdonald, solid)",
                     sec.axis = sec_axis(~./coeff, name = "Basic reproduction number\n(Chitnis, dashed)")) +
  scale_linetype("Contact function") +
  scale_color_manual("", values = c(rev(c4a("poly.dark24", 5))[-2]), guide = "none") +
  facet_wrap(vars(Label), scales = 'free', nrow = 1) +
  theme_minimal_grid()

var_R0_plot_probs <- var_R0_plot_df %>%
  filter(Type == unique(plot_df$Type)[2]) %>%
  ggplot(aes(x = variable_value, color = Prefix)) +
  # Ross-Macdonald
  geom_line(aes(y = RM), lwd = 1, lty = 1) +
  # Chitnis
  geom_line(aes(y = Chitnis * coeff), lwd = 1, lty = 2) +
  scale_x_continuous("Parameter value") +
  scale_y_continuous("Basic reproduction number\n(Ross-Macdonald, solid)",
                     sec.axis = sec_axis(~./coeff, name = "Basic reproduction number\n(Chitnis, dashed)")) +
  scale_linetype("Contact function") +
  scale_color_manual("", values = c(rev(c4a("poly.dark24", 5))[-1]), guide = "none") +
  facet_wrap(vars(Label), scales = 'free', nrow = 1) +
  theme_minimal_grid()

vars_R0_out_plot = plot_grid(var_R0_plot_rates, var_R0_plot_probs, nrow = 2)

ggsave(paste0(folder_name, "/vars_R0_out.png"), vars_R0_out_plot, 
       width = 13.333, height = 7.5, units = "in")



















