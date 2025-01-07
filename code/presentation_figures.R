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

# TODO
# [] Reduce to only show key parameters
# [] PRCC plots


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
# 6. PRCC plots


# Load in data ----
combined_df = readRDS(file = "data/GCD_R0_vals.rds.gz")

# Plot Figures ----

chosen_mosquito_type = "persistent"

folder_name = paste0("figures/GCD_figures/",chosen_mosquito_type)

# For now, just look at persistent mosquitoes
plot_df <- combined_df %>% 
  filter(mosquito_type == chosen_mosquito_type) %>% 
  mutate(f = 1-f, .keep = "unused") %>% 
  mutate(varied_parameter = ifelse(varied_parameter == "f", "sigma", varied_parameter)) %>% 
  rename(sigma = f) %>% 
  arrange(across(pQ:sigma)) %>% 
  # select(c(mosquito_type:gammaR,GCD, to_host_contact, to_vector_contact, R0)) %>% 
  select(c(mosquito_type:gammaR, GCD, R0)) %>% 
  mutate(GCD_day = GCD / 1440) %>%   
  mutate(b = 1/GCD,
         varied_parameter = if_else(
           varied_parameter == "theta",
           "b",
           varied_parameter
         )
  ) %>% 
  relocate(b, .after = "lQ") #%>% 
# mutate(mod_GCD = case_when(
#   model_type == "Exponential" ~ GCD,
#   model_type == "Mechanistic" ~ GCD - (1/gammaV) - (1/gammaR)
# )) %>% 
# mutate(GCD_day = mod_GCD / 1440) %>% 
# filter(GCD_day < 1/mu/1440)


# Dictionary for labeling
param_table = tibble(
  Symbol = c("$\\lambda_Q$", "$p_L$", "$\\lambda_L$", "$p_P$", "$\\lambda_P$", "$p_G$", "$\\lambda_G$", "$\\sigma$", "$b$"),
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
  Label = c("Host-seeking rate", "Landing success probability", "Landing rate", "Probing success probability", "Probing rate", "Ingestion success probability", "Ingesting rate", "Persistence probability", "Biting rate (exponential)"), 
  Type = c("A Rates", "B Probabilities", "A Rates", "B Probabilities", "A Rates", "B Probabilities", "A Rates", "B Probabilities", "A Rates"),
  Prefix = c("Host-seeking", "Landing", "Landing", "Probing", "Probing", "Ingesting",  "Ingesting", "Persistence", "Exponential"), 
  short_label = c("lQ", "pL", "lL", "pP", "lP", "pG", "lG", "sigma", "b")
)

# Add in nice labels for the parameters
plot_df = plot_df %>% 
  left_join(rename(param_table, varied_parameter = short_label))

plot_df$Label <- factor(plot_df$Label, levels = c("Biting rate (exponential)", "Host-seeking rate", "Landing success probability", "Landing rate", "Probing success probability", "Probing rate", "Ingestion success probability", "Ingesting rate", "Persistence probability"))

# Labels
plot_df$Prefix = factor(plot_df$Prefix, levels = c("Exponential", "Host-seeking", "Persistence", "Landing", "Probing", "Ingesting"))
plot_df$Type = factor(plot_df$Type, levels = c("A Rates", "B Probabilities"))

## Settings for plots ----

selected_parameters <- c(
  # "theta", "pP", "sigma", "lP"
  # try the following
  # "pP", "pG", "lQ", "theta" # maybe not pG. definitely not pG or pP
  "b", "lQ", "pP" # persistent
)

parameter_palette = c("Black", (c4a("tol.muted", 2))) # !!! update to use consistent colors for each parameter across all the plots

## GCD vs R0 figures ----

#### OUT Contact rates ----

R0_GCD_OUT_df <- plot_df %>%
  filter(varied_parameter %in% selected_parameters) %>%
  filter(transmission_type == "OUT") %>% 
  distinct()

R0_GCD_OUT_plot <- R0_GCD_OUT_df %>% 
  filter(contact_type == "RM") %>% 
  ggplot(aes(x = GCD_day, y = R0, group = Label, color = Label)) +
  # Curves
  geom_line(lwd = 1,
            # arrow = arrow(ends = "first", type = "closed")
  ) +
  scale_x_continuous("Gonotrophic cycle duration [days]",
                     breaks = seq(1,21),
                     limits = c(4, 21)
                     # limits = c(NA, max(filter(R0_GCD_OUT_df, varied_parameter == "lP")$GCD_day))
  ) +
  scale_y_continuous(TeX("Basic reproduction number, $R_0$"),
                     # limits = c(NA, 20)
                     # limits = c(NA, max(filter(R0_GCD_OUT_df, varied_parameter != "theta")$R0))
  ) +
  scale_color_manual("Parameter", values = parameter_palette) +
  # ggtitle("R0 as a function of GCD (using OUT-rates)") +
  theme_minimal_grid(16)

ggsave(paste0(folder_name, "/R0_GCD_OUT.png"), R0_GCD_OUT_plot, 
       width = 13.333, height = 7.5, units = "in")


## 1/GCD vs R0 figures ----

R0_invGCD_OUT_df <- plot_df %>% 
  filter(varied_parameter %in% selected_parameters) %>%
  mutate(invGCD = 1/GCD_day) %>% 
  filter(transmission_type == "OUT",
         contact_type == "RM",
         # invGCD < 2
  ) %>% 
  distinct() %>% 
  group_by(varied_parameter) %>% 
  mutate(nudge = rnorm(1, 0.1)) %>% 
  mutate(invGCD_nudge = invGCD + nudge) %>% 
  ungroup()


R0_invGCD_OUT_plot <- R0_invGCD_OUT_df %>% 
  ggplot(aes(x = invGCD, y = R0, group = Label, color = Label)) +
  geom_line(lwd = 1,
            # arrow = arrow(ends = "last", type = "closed"),
            # alpha = 1
  ) +
  scale_x_continuous("1/Gonotrophic cycle duration = \"Biting rate\" [per day]",
                     limits = c(0.05, 1/4)
  ) +
  scale_y_continuous(TeX("Basic reproduction number, $R_0$"),
                     breaks = seq(0,7),
                     limits = c(NA, 7)
  ) +
  scale_color_manual(
    "Parameter", 
    values = parameter_palette
    ) +
  theme_minimal_grid(16) +
  theme(strip.text.x = element_text(hjust = 0))

ggsave(paste0(folder_name, "/R0_invGCD_OUT.png"), R0_invGCD_OUT_plot, 
       width = 13.333, height = 7.5, units = "in")

### Plot GCD against each parameter individually ----
var_GCD_plot_probs <- plot_df %>% 
  pivot_longer(cols = lQ:sigma) %>% 
  filter(name == varied_parameter) %>% 
  filter(Type == unique(plot_df$Type)[1]) %>% 
  ggplot(aes(x = value, y = GCD_day, color = Prefix)) +
  geom_line(lwd = 1) +
  scale_x_continuous("",
                     limits = c(10^(-2), NA)
  ) +
  scale_y_continuous("Gonotrophic cycle duration (days)",
                     limits = c(NA, 21)
  ) +
  # scale_color_manual("Parameter", values = parameter_palette, guide = "none") +
  facet_wrap(vars(Label), scales = 'free', nrow = 1) +
  theme_minimal_grid()

var_GCD_plot_rates <- plot_df %>% 
  pivot_longer(cols = lQ:sigma) %>% 
  filter(name == varied_parameter) %>% 
  filter(Type == unique(plot_df$Type)[2]) %>% 
  ggplot(aes(x = value, y = GCD_day, color = Prefix)) +
  geom_line(lwd = 1) +
  scale_x_log10("Parameter value") +
  scale_y_continuous("Gonotrophic cycle duration (days)",
                     limits = c(NA, 21)
  ) +
  # scale_color_manual("Parameter", values = parameter_palette, guide = "none") +
  facet_wrap(vars(Label), scales = 'free', nrow = 1) +
  theme_minimal_grid()


var_GCD_plot <- plot_df %>% 
  filter(varied_parameter %in% selected_parameters) %>%
  pivot_longer(cols = lQ:sigma) %>% 
  filter(name == varied_parameter) %>%
  filter(GCD_day < 21) %>% 
  # filter(Type == unique(plot_df$Type)[2]) %>% 
  ggplot(aes(x = value, y = GCD_day, color = Label)) +
  geom_line(lwd = 1) +
  scale_x_log10("Parameter value") +
  scale_y_continuous("Gonotrophic cycle duration (days)"
                     # limits = c(NA, 21)
  ) +
  # scale_color_manual("Parameter", values = parameter_palette, guide = "none") +
  # facet_wrap(vars(Prefix), nrow = 2) +
  theme_minimal_grid()

plot_grid(var_GCD_plot_rates, var_GCD_plot_probs, nrow = 2)

### Plot R0 against each parameter individually ----
coeff = 10

var_R0_plot_df = plot_df %>%
  pivot_longer(cols = pQ:sigma, names_to = "variable", values_to = "variable_value") %>% 
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

# PRCC Plots ----

# Load in data
LHS_data = read_csv("data/julia_outputs.csv.gz")

rank_data = LHS_data %>% 
  group_by(type) %>% 
  mutate(across(lQ:R0, ~ rank(.x)))

PRCC_data <- rank_data %>%
  pivot_longer(cols = lQ:pG, names_to = "input", values_to = "input_value") %>%
  pivot_longer(cols = GCD:R0, names_to = "output", values_to = "output_value") %>%
  group_by(type, input, output) %>%
  summarise(
    PRCC = cor(input_value, output_value),
    .groups = "drop"
  )


plot_data <- PRCC_data %>% 
  right_join(
    param_table %>% 
      rename(input = short_label) %>% 
      rbind(list(
        Symbol = "$\\sigma$",
        Description = "Persistence probability",
        Label = "Persistence probability",
        Type = "B Probabilities",
        Prefix = "Persistence",
        input = "sigma"
      ))
  ) %>% 
  mutate(
    output_label = case_when(
      output == "R0" ~ "Basic reproduction number",
      output == "GCD" ~ "Gonotrophic cycle duration",
      output == "N_offspring" ~ "Basic offspring number",
    )
  ) %>% 
  mutate(
    type_label = case_when(
      type == "flighty" ~ "Flighty baseline",
      type == "persistent" ~ "Persistent baseline",
      type == "max" ~ "Maximum variation",
    )
  ) 

plot_data$type = factor(plot_data$type, levels = c("flighty", "persistent", "max"))
plot_data$input = factor(plot_data$input, levels = c(
  "lQ", "lL", "lP", "lG", "sigma", "pL", "pP", "pG"
))
plot_data$output_label = factor(plot_data$output_label, levels = c(
  "Gonotrophic cycle duration", "Basic offspring number","Basic reproduction number"
))

plot_data$type_label = factor(plot_data$type_label, levels = c(
  "Flighty baseline","Persistent baseline","Maximum variation"
))


plot_data$Label = factor(plot_data$Label, levels = c(
  c("Questing rate", "Landing rate","Probing rate", "Ingesting rate", "Persistence probability", "Landing success probability",  "Probing success probability", "Ingestion success probability")
))


PRCC_plots <- plot_data %>% 
  arrange(input) %>% 
  filter(type %in% c("flighty", "persistent", "max")) %>% 
  ggplot(aes(x = Label, y = PRCC, fill = type_label)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, color = "black", size =1 ) +
  facet_wrap(~output_label, ncol = 1, scales = "free_y") +
  scale_fill_discrete(name = "Parameter set:") +
  scale_x_discrete(name = "") +
  theme_minimal(16)+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/GCD_figures/PRCCs_draft.png", PRCC_plots, 
       width = 13.333, height = 7.5, units = "in")


