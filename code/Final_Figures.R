# New figures

# Load libraries ----
library(tidyverse)
library(latex2exp)
library(cols4all)
library(cowplot)
library(expm)
library(matlib)
library(scales)

# Load data ----
Full_df = read_rds("data/GCD_R0_data.rds") %>% 
  mutate(varied_parameter = if_else(
    is.na(varied_parameter),
    "none",
    varied_parameter)) %>% 
  # Set up labels
  mutate(Type = if_else(
    `Model type` == "Mechanistic",
    paste0(`Model type`, " (", varied_parameter,")"),
    `Model type`))

# Set varied_parameter levels
Full_df$varied_parameter = factor(Full_df$varied_parameter, 
                                  levels = c("none", "lQ", "pP", "pG"))

# Labels for the model types
type_order = c("Standard", "Exponential", "Empirical", "Phenomenological", "Mechanistic")

Full_df$`Model type` = factor(Full_df$`Model type`, 
                              levels = type_order)

# Define functions ----

# Define pdf for Phase-type distribution
PH_pdf <- function(x, A_matrix, v_alpha) {
  A_dim = dim(A_matrix)[1]
  if (is.null(A_dim)) {
    v_ones = 1
    pdf_val = exp(A_matrix * x) * (-A_matrix)
  } else {
    v_ones = matrix(rep(1, A_dim), ncol = 1)
    
    pdf_val = t(v_alpha) %*% expm(A_matrix * x) %*% (-A_matrix %*% v_ones)
  }
  return(pdf_val)
}

# Define mean for Phase-type distribution
PH_mean <- function(A_matrix, v_alpha) {
  A_dim = dim(A_matrix)[1]
  if (is.null(A_dim)) {
    v_ones = 1
    mean_val = (-1 / A_matrix[[1]])
  } else {
    v_ones = matrix(rep(1, A_dim), ncol = 1)
    
    mean_val = t(v_ones) %*% inv(-t(A_matrix)) %*% v_alpha
  }
  return(mean_val)
}

# Helper function: place legends in empty facets of plot grids
# Code by: Artem Sokolov, found here: https://stackoverflow.com/questions/54438495/shift-legend-into-empty-facets-of-a-faceted-plot-in-ggplot2
shift_legend <- function(p) {
  pnls <- cowplot::plot_to_gtable(p) %>%
    gtable::gtable_filter("panel") %>%
    with(setNames(grobs, layout$name)) %>%
    purrr::keep(~ identical(.x, zeroGrob()))
  
  if (length(pnls) == 0) stop("No empty facets in the plot")
  
  lemon::reposition_legend(p, "center", panel = names(pnls))
}

# Plot parameters ----

# Colors for the model types
type_colors = c("black", c4a("brewer.dark2",4))

# 1. Pdfs of distributions ----

# Plot the pdf (using the actual function!) of each model type for a couple 
# values of theta: 1/3 days, 1/2 days, 1 days, 2 days
theta_vals = c((1/2) * 1440, 1440, 2 * 1440)
tolerance = 4.5

pdf_max = 4 * 1440
resolution = 10000

Fig1_labeller <- function(value) TeX(paste0("Mean of ", as.double(value)/1440, " days"))

# For mechanistic model, we can't guarantee these values so choose the closest ones
Mech_df_filtered <- Full_df %>% 
  filter(`Model type` %in% c("Mechanistic")) %>% 
  filter(varied_parameter == "lQ") %>% # the curves look almost identical, so just look at one
  # mutate(Type = "Mechanistic (lQ, pP, pG)") %>% 
  group_by(`Model type`, varied_parameter) %>% 
  filter(theta %in% sapply(theta_vals, function(x) theta[which.min(abs(theta - x))])) %>%
  ungroup()

# Arrange in a grid or row
Figure1_df <- Full_df %>% 
  # Remove the Mechanistic rows, which are handled separately
  filter(!(`Model type` %in% c("Mechanistic"))) %>% 
  # Combine Standard and Exponential since they are equivalent in this case
  mutate(Type = if_else(
    Type %in% c("Standard", "Exponential"),
    "Standard / Exponential",
    Type)
  ) %>% 
  # Just keep chosen theta values
  filter(theta %in% theta_vals) %>% 
  # Add the filtered Mechanistic rows
  bind_rows(Mech_df_filtered) %>% 
  # Assign the closest theta to each row (to relate Mechanistic with the rest)
  group_by(Type) %>% 
  mutate(closest_theta = theta_vals[sapply(theta, function(x) which.min(abs(x - theta_vals)))]) %>%
  ungroup() %>% 
  # Select only relevant rows
  dplyr::select(Type, theta, closest_theta, v_alpha, A_matrix) %>% 
  # Add in values to plug into the pdf
  cross_join(tibble(x = seq(0, pdf_max, length.out = resolution))) %>% 
  # Calculate values of the pdf
  rowwise() %>% 
  mutate(pdf_val = PH_pdf(x, A_matrix, v_alpha))
# Cut off when pdf_val is really small
# group_by(Type, closest_theta)# %>% 
# filter(pdf_val > 1e-5) %>% 

Figure1_df$Type = factor(
  Figure1_df$Type,
  levels = c("Empirical", "Phenomenological", "Mechanistic (lQ, pP, pG)","Standard / Exponential"))

Figure1_labels = c(expression("Standard / Exponential"), expression("Empirical"), expression("Phenomenological"), 
                   expression("Mechanistic " (lambda[Q]))
)

Fig1_color_vals = c("Standard / Exponential" = "black",
                    # "Exponential" = c4a("brewer.dark2", 4)[2],
                    "Empirical" = c4a("brewer.dark2", 3)[1],
                    "Phenomenological" = c4a("brewer.dark2", 3)[2],
                    "Mechanistic (lQ, pP, pG)" = c4a("brewer.dark2", 3)[3]#,
                    # "Mechanistic (lQ)" = c4a("brewer.dark2", 4)[3],
                    # "Mechanistic (pP)" = c4a("brewer.dark2", 4)[3],
                    # "Mechanistic (pG)" = c4a("brewer.dark2", 4)[3]
)

Fig1_lty_vals = c("Standard / Exponential" = 1,
                  # "Exponential" = 1,
                  "Empirical" = 1,
                  "Phenomenological" = 1,
                  "Mechanistic (lQ)" = 2,
                  "Mechanistic (pP)" = 3,
                  "Mechanistic (pG)" = 4)


Figure1 = Figure1_df %>%
  ggplot(aes(x = x / 1440, y = pdf_val, color = Type)) +
  geom_path(
    lwd = 0.75,
    alpha = 0.9
  ) + 
  # # Add line showing the means of the distributions (which are pretty much equal)
  # geom_vline(
  #   data = Figure1_df %>%
  #     select(-c(x, pdf_val)) %>% unique() %>%
  #     # group_by(`Model type`, closest_theta) %>%
  #     rowwise() %>%
  #     mutate(mean = PH_mean(A_matrix, v_alpha)),
  #   aes(xintercept = mean / 1440, color = Type),
  #   lty = 2, alpha = 0.1, lwd = 2) +
  facet_wrap( ~ closest_theta,
              labeller = as_labeller(Fig1_labeller),
              ncol = 1,
              scales = "free_y") +
  scale_x_continuous(
    name = "Gonotrophic cycle duration [Days]",
    limits =  c(0, 3),
    expand = c(0.01,0),
    breaks = seq(0,10)
  ) +
  scale_y_continuous(
    name = "Density",
    limits = c(0, NA),
    expand = c(0.01,0)
  ) +
  # scale_linetype_manual(
  #   name = "Model type",
  #   values = Fig1_lty_vals,
  #   breaks = unique(Figure1_df$Type)
  # ) +
  scale_color_manual(
    name = "Model type:",
    values = Fig1_color_vals,
    breaks = unique(Figure1_df$Type),
    labels = Figure1_labels
    # labels = scales::label_parse()
  ) +
  guides(
    color = guide_legend(
      position = "top",
      direction = "horizontal",
      nrow = 2
    )
  ) +
  theme_half_open(11) +
  theme(
    strip.background = element_rect(color = "white", fill = "white"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

Figure1

# shift_legend(Figure1)

ggsave("figures/Figure1.pdf", Figure1, width = 4.25, height = 3.25 * 9/6.5, units = "in")

# 2. R0 vs. standard biting rates ----

# Plot R0 as a function of the standard biting rate for each model type
# Distinguish types by color and linetype (make a backup with linetype)
# Make sure the labels are nice

Fig2_color_vals = c("Standard" = "black",
                    "Exponential" = "black",
                    "Empirical" = c4a("brewer.dark2", 3)[1],
                    "Phenomenological" = c4a("brewer.dark2", 3)[2],
                    "Mechanistic~(lambda[Q])" = c4a("brewer.dark2", 3)[3],
                    "Mechanistic~(p[P])" = c4a("brewer.dark2", 3)[3],
                    "Mechanistic~(p[G])" = c4a("brewer.dark2", 3)[3])

Fig2_lty_vals = c("Standard" = 1,
                  "Exponential" = 2,
                  "Empirical" = 1,
                  "Phenomenological" = 1,
                  "Mechanistic~(lambda[Q])" = 2,
                  "Mechanistic~(p[P])" = 3,
                  "Mechanistic~(p[G])" = 4)

Figure2_df = Full_df %>%
  filter(theta > (1/2) * 1440) %>%
  mutate(Type = case_when(
    Type == "Mechanistic (lQ)" ~ "Mechanistic~(lambda[Q])",
    Type == "Mechanistic (pP)" ~ "Mechanistic~(p[P])",
    Type == "Mechanistic (pG)" ~ "Mechanistic~(p[G])",
    TRUE ~ Type
  ))

Figure2_labels = c(expression("Standard"), expression("Exponential"), expression("Empirical"), expression("Phenomenological"), 
                   expression("Mechanistic " (lambda[Q])), expression("Mechanistic " (p[P])),expression("Mechanistic " (p[G]))
)



Figure2 <- Figure2_df %>% 
  ggplot(aes(x = 1440 / theta, y = R0, color = Type, lty = Type)) +
  geom_hline(aes(yintercept = 1), color = "grey", lwd = 2) +
  geom_line(lwd = 0.75) +
  scale_x_continuous(
    name = TeX("Standard biting rate [Days$^{-1}$]"),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = TeX("Basic reproduction number \\, [$R_0$]"),
    expand = c(0,0)
  ) +
  scale_linetype_manual(
    name = "Model type:",
    values = Fig2_lty_vals,
    breaks = unique(Figure2_df$Type),
    labels = Figure2_labels
  ) +
  scale_color_manual(
    name = "Model type:",
    values = Fig2_color_vals,
    breaks = unique(Figure2_df$Type),
    labels = Figure2_labels
  ) +
  theme_half_open(11) + 
  theme(
    legend.key.width = unit(0.4, "in")
  )

Figure2

ggsave("figures/Figure2.pdf", Figure2, width = 7.5, height = 3.25 * 9/6.5, units = "in")

Figure2_alt <- Figure2 +
  guides(
    color = guide_legend(
      position = "top",
      direction = "horizontal",
      nrow = 2,
      byrow = T
    ),
    linetype = guide_legend(
      position = "top",
      direction = "horizontal",
      nrow = 2,
      byrow = T
    )
  )

ggsave("figures/Figure2_alt.pdf", Figure2_alt, width = 7.5, height = 3.25 * 9/6.5, units = "in")

# 3. R0 vs. mechanistic parameters ----

Mech_df <- read_rds("data/Mechanistic_results.rds")


# Distinguish parameters by color -- different from 2.
# Distinguish sets by linetype (flighty vs. persistent)

nice_mech_labels = data.frame(
  pQ = "Questing~success~probability~(p[Q])",
  pL = "Landing~success~probability~(p[L])",
  pP = "Probing~success~probability~(p[P])",
  pG = "Ingesting~success~probability~(p[G])",
  sigma = "Persistence~probability~(sigma)",
  lQ = "Questing~rate~(lambda[Q])",
  lL = "Landing~rate~(lambda[P])",
  lP = "Probing~rate~(lambda[P])",
  lG = "Ingesting~rate~(lambda[G])"
) %>% pivot_longer(everything(), values_to = "nice_labels")

# Define a labeller function that applies TeX()
TeX_labeller <- function(x) {
  sapply(x, function(label) {
    # Escape underscores for LaTeX
    label <- gsub("_", "\\\\_", label)  
    
    # Replace spaces with LaTeX-friendly spacing
    label <- gsub(" ", "~", label)      
    
    # Convert label to expression string
    paste0("expression(", label, ")")
  })
}


Figure3_df <- Mech_df %>% 
  filter(parameter_type == "varied") %>% 
  select(-parameter_type) %>% 
  pivot_longer(lQ:sigma) %>% 
  filter(name == varied_parameter) %>% 
  mutate(parameter_type = case_when(
    name %in% c("lQ", "lL", "lP", "lG") ~ "rate",
    name %in% c("pQ", "pL", "pP", "pG", "sigma") ~ "probability"
  )) %>% 
  # Reduce the range for rates
  filter(!(parameter_type == "rate" & value >= (1/60))) %>% 
  right_join(nice_mech_labels)

Figure3_df$name <- factor(
  Figure3_df$name,
  levels = c("pQ", "pL", "pP", "pG", "sigma", "lQ", "lL", "lP", "lG"))

Figure3_df$nice_labels <- factor(
  Figure3_df$nice_labels,
  levels = nice_mech_labels$nice_labels)

Figure3 <- Figure3_df %>% 
  ggplot(aes(x = value, y = R0, linetype = mosquito_type, color = parameter_type)) +
  geom_line(lwd = 1) +
  facet_wrap(
    ~ nice_labels,
    ncol = 5,
    labeller = labeller(nice_labels = label_parsed),
    scales = "free_x"
  ) +
  scale_x_continuous(
    name = "Parameter value",
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = TeX("Basic reproduction number \\, [$R_0$]"),
    expand = c(0,0)
  ) +
  scale_color_discrete(
    name = "Parameter type:"
  ) +
  scale_linetype_discrete(
    name = "Mosquito type:"
  ) +
  theme_half_open(18)

shift_legend(Figure3)
ggsave("figures/Figure3.pdf", shift_legend(Figure3), width = 20, height = 8, units = "in")

# Supplementary: theta vs. mechanistic parameters ----


Figure4 <- Figure3_df %>%
  filter(theta < 21 * 1440) %>% 
  ggplot(aes(x = value, y = theta / 1440, linetype = mosquito_type, color = parameter_type)) +
  geom_line(lwd = 1) +
  facet_wrap(
    ~ nice_labels,
    ncol = 5,
    labeller = labeller(nice_labels = label_parsed),
    scales = "free"
  ) +
  scale_x_continuous(
    name = "Parameter value",
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = TeX("Gonotrophic cylce duration \\, [Days]"),
    expand = c(0,0)
  ) +
  scale_color_discrete(
    name = "Parameter type:"
  ) +
  scale_linetype_discrete(
    name = "Mosquito type:"
  ) +
  theme_half_open(18)

shift_legend(Figure4)
ggsave("figures/SuppFigure.pdf", shift_legend(Figure4), width = 20, height = 8, units = "in")
# Distinguish parameters by color
# Distinguish sets by linetype (flighty vs. persistent)