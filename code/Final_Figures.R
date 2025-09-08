# New figures

# Load libraries ----
library(tidyverse)
library(latex2exp)
library(cols4all)
library(cowplot)
library(expm)
library(matlib)
library(scales)
library(ggh4x)
library(varhandle) # to deal with factor problems in plotting
library(egg) # for arranging plots

# Load data ----
Full_df = read_rds("data/GCD_R0_data.rds") |> 
  mutate(varied_parameter = if_else(
    is.na(varied_parameter),
    "none",
    varied_parameter)) |> 
  # Set up labels
  mutate(Type = if_else(
    `Model type` == "Mechanistic",
    paste0(`Model type`, " (", varied_parameter,")"),
    `Model type`))

# Dictionary for labeling
param_table = tibble(
  Symbol = c("$p_Q$", "$\\lambda_Q$", "$p_L$", "$\\lambda_L$", "$p_P$", "$\\lambda_P$", "$p_G$", "$\\lambda_G$", "$\\sigma$", "$b$"),
  Description = c(
    "Probability of progressing from seeking to landing",
    "Exit rate from seeking stage (per minute)", 
    "Probability of progressing from landing to probing",
    "Exit rate from landing stage (per minute)", 
    "Probability of progressing from probing to ingesting",
    "Exit rate from probing stage (per minute)",
    "Probability of progressing from ingesting to ovipositing",
    "Exit rate from ingestion stage (per minute)",
    "Probability of seeking a new host given feeding failure",
    "Biting rate (exponential model)"
    
  ),
  Label = c("Seeking success", "Seeking rate", "Landing success", "Landing rate", "Probing success", "Probing rate", "Ingesting success", "Ingesting rate", "Persistence probability", "Biting rate (exponential)"), 
  Type = c("B Probabilities", "A Rates", "B Probabilities", "A Rates", "B Probabilities", "A Rates", "B Probabilities", "A Rates", "B Probabilities", "A Rates"),
  Prefix = c("Seeking", "Seeking", "Landing", "Landing", "Probing", "Probing", "Ingesting",  "Ingesting", "Persistence", "Exponential"), 
  short_label = c("pQ", "lQ", "pL", "lL", "pP", "lP", "pG", "lG", "sigma", "b")
)


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
  pnls <- cowplot::plot_to_gtable(p) |>
    gtable::gtable_filter("panel") |>
    with(setNames(grobs, layout$name)) |>
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
theta_vals = c((1/4) * 1440, (1/2) * 1440, 1440, 2 * 1440)
tolerance = 4.5

pdf_max = 8 * 1440
resolution = 10000

Fig1_labeller <- function(value) TeX(paste0("Mean = ", as.double(value)/1440, " days"))

# For mechanistic model, we can't guarantee these values so choose the closest ones
Mech_df_filtered <- Full_df |> 
  filter(`Model type` %in% c("Mechanistic")) |> 
  filter(varied_parameter == "lQ") |> # the curves look almost identical, so just look at one
  # mutate(Type = "Mechanistic (lQ, pP, pG)") |> 
  group_by(`Model type`, varied_parameter) |> 
  filter(theta %in% sapply(theta_vals, function(x) theta[which.min(abs(theta - x))])) |>
  ungroup()

# Arrange in a grid or row
Figure1_df <- Full_df |> 
  # Remove the Mechanistic rows, which are handled separately
  filter(!(`Model type` %in% c("Mechanistic"))) |> 
  # Combine Standard and Exponential since they are equivalent in this case
  mutate(Type = if_else(
    Type %in% c("Standard", "Exponential"),
    "Standard / Exponential",
    Type)
  ) |> 
  # Just keep chosen theta values
  filter(theta %in% theta_vals) |> 
  # Add the filtered Mechanistic rows
  bind_rows(Mech_df_filtered) |> 
  # Assign the closest theta to each row (to relate Mechanistic with the rest)
  group_by(Type) |> 
  mutate(closest_theta = theta_vals[sapply(theta, function(x) which.min(abs(x - theta_vals)))]) |> 
  ungroup() |> 
  # Select only relevant rows
  dplyr::select(Type, theta, closest_theta, v_alpha, A_matrix) |> 
  mutate(theta_label = case_when(
    closest_theta == 360 ~ "A",
    closest_theta == 720 ~ "B",
    closest_theta == 1440 ~ "C",
    closest_theta == 2880 ~ "D"
  )) |>
  # Add in values to plug into the pdf
  cross_join(tibble(x = seq(0, pdf_max, length.out = resolution)))  |> 
  filter(case_when(
    closest_theta == 360 ~ x < 0.75 * 1440,
    closest_theta == 720 ~ x < 1.5 * 1440,
    closest_theta == 1440 ~ x < 3 * 1440,
    closest_theta == 2880 ~ x < 6 * 1440
  )) |>
  # Calculate values of the pdf
  rowwise() |> 
  mutate(pdf_val = PH_pdf(x, A_matrix, v_alpha))

Figure1_df$Type = factor(
  Figure1_df$Type,
  levels = c("Empirical", "Phenomenological", "Mechanistic (lQ)","Standard / Exponential"))

Figure1_labels = c(expression("Standard / Exponential"), expression("Empirical"), expression("Phenomenological"), 
                   expression("Mechanistic " (lambda[Q]))
)

Fig1_color_vals = c("Standard / Exponential" = "black",
                    "Empirical" = c4a("brewer.dark2", 3)[1],
                    "Phenomenological" = c4a("brewer.dark2", 3)[2],
                    "Mechanistic (lQ)" = c4a("brewer.dark2", 3)[3]#
)

facet_labels <- Figure1_df |>
  group_by(closest_theta) |> 
  filter(x > 0) |> 
  mutate(label = LETTERS[1:n()],  # Assign "A", "B", "C", ...
         x_pos = 360*min(x, na.rm = TRUE)/1440,  # Align left
         y_pos = max(pdf_val, na.rm = TRUE) * .95) |>  # Slightly above max y
  distinct(closest_theta, x_pos, y_pos, theta_label)

# pdf figure
Figure1 = Figure1_df |>
  filter(closest_theta < 2 * 1440) |> 
  ggplot(aes(x = 24 * x / (1440), y = pdf_val, color = Type)) +
  # Plot grey dotted line showing the mean
  geom_vline(
    aes(xintercept = 24 * closest_theta/1440),
    lwd = 0.5, color = "grey", lty = 2
  ) +
  # Plot pdfs
  geom_line(
    lwd = 0.5,
    alpha = 0.75
  ) + 
  # Label subplots A, B, C, ...
  geom_text(data = facet_labels |> filter(closest_theta < 2 * 1440),
            aes(x = 0.01, y = y_pos, label = theta_label, group = theta_label),
            color = "black", size = 4,
            # vjust = "inward",
            hjust ="inward",
            # nudge_x = 360
  ) +
  # Wrap by theta value
  facet_wrap( ~ closest_theta,
              # labeller = as_labeller(Fig1_labeller),
              nrow = 1,
              scales = "free") +
  scale_x_continuous(
    name = "Gonotrophic cycle duration [Hours]",
    limits =  c(0, NA),
    expand = c(0.0,0),
    breaks = seq(0, 72, by = 6)
  ) +
  scale_y_continuous(
    name = "Density",
    limits = c(0, NA),
    expand = c(0.01,0),
  ) +
  scale_color_manual(
    name = "Model type:",
    values = Fig1_color_vals,
    breaks = unique(Figure1_df$Type),
    labels = Figure1_labels
  ) +
  guides(
    # Put the color legend at the top of the plot
    color = guide_legend(
      position = "top",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  theme_half_open(10) +
  theme(
    strip.background = element_rect(color = "white", fill = "white"),
    strip.text = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

Figure1

# shift_legend(Figure1)

ggsave("figures/Figure1.pdf", Figure1, width = 6.5, height = 1.75 * 9/6.5, units = "in")

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

Figure2_df = Full_df |>
  mutate(
    Type = case_when(
      Type == "Mechanistic (lQ)" ~ "Mechanistic~(lambda[Q])",
      Type == "Mechanistic (pP)" ~ "Mechanistic~(p[P])",
      Type == "Mechanistic (pG)" ~ "Mechanistic~(p[G])",
      TRUE ~ Type
    ),
    sbr = 1440/theta  # standard biting rate
  )

Figure2_ticks <- Figure2_df |> 
  filter(sbr < 2.01) |> # only keep values up to 2 bites per day
  # Get x-coordinate where R0 first exceeds/is less than one
  group_by(Type) |> 
  arrange(sbr) |> 
  summarise(
    first_R0_greater_1 = sbr[R0 > 1][1],
    first_R0_less_1 = max(sbr[R0 > 1])
  ) |> 
  # select(-sbr) |> 
  unique()

Figure2_labels = c(expression("Standard"), expression("Exponential"), expression("Empirical"), expression("Phenomenological"), 
                   expression("Mechanistic " (lambda[Q])), expression("Mechanistic " (p[P])),expression("Mechanistic " (p[G]))
)

Figure2_arrows = Figure2_df |>
  filter(`Model type` == "Mechanistic") |> 
  pivot_longer(cols = lQ:sigma) |> 
  filter(varied_parameter == name) |> 
  group_by(varied_parameter) |>
  arrange(value) |> 
  mutate(sbr = 1440 / theta) |> 
  filter(between(sbr, 0.9, 1)) |> 
  group_by(varied_parameter) |> 
  mutate(
    x_first = min(sbr),
    x_last = max(sbr)
  ) |> 
  group_by(varied_parameter) |> 
  filter(sbr %in% c(x_first, x_last)) |> 
  mutate(
    R0_first = R0[sbr == min(sbr)],
    R0_last = R0[sbr == max(sbr)]
  ) |> ungroup() |> 
  select(Type, x_first, x_last, R0_first, R0_last) |> distinct()


Figure2 <- Figure2_df |> 
  filter(sbr < 2.01) |> # only keep values up to 2 bites per day
  ggplot(aes(color = Type, lty = Type)) +
  # Grey line for R0 = 1
  geom_hline(aes(yintercept = 1), color = "grey", lwd = 1) +
  # R0-GCD curves
  geom_line(aes(x = 1440 / theta, y = R0),
            lwd = 0.75) +
  # Add ticks for R0 = 1 crossing below the x-axis
  geom_rug(
    data = Figure2_ticks,
    aes(x = first_R0_greater_1),
    sides = "b", size = 0.75, outside = TRUE,
    length = unit(0.3, "in"),
    show.legend = F
  ) +
  # Arrows showing direction of increasing parameter values
  geom_segment(
    data = Figure2_arrows,
    aes(x = x_first, y = R0_first, xend = x_last, yend = R0_last, color = Type),
    lwd = 0, 
    # alpha = 0,
    arrow = arrow(length = unit(0.2, "inches"), ends = "last", type = "closed"),
    show.legend = F,
    inherit.aes = F
  ) +
  scale_x_continuous(
    name = TeX("Standard biting rate [Days$^{-1}$]"),
    # limits = c(0,1.9),
    expand = c(0, 0, 0, 0.01),
    # breaks = seq(0, 2, by = 0.25)
  ) +
  scale_y_continuous(
    name = TeX("Basic reproduction number \\, [$R_0$]"),
    # breaks = seq(0, 2, by = 0.25),
    limits = c(0, NA),
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
  coord_cartesian(
    xlim = c(0, 1.8),
    clip = "off" # needed to let ticks get plotted out of the axes
  ) + 
  theme_half_open(11) + 
  theme(
    axis.title.x = element_text(margin = margin(t = 10)),
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

Figure2_alt

ggsave("figures/Figure2_alt.pdf", Figure2_alt, width = 7.5, height = 3.25 * 9/6.5, units = "in")

# Alternate figure with non-mechanistic and mechanistic model curves separated for easier viewing

Figure2_separated_nonmech <- Figure2_df |> 
  filter(`Model type` != "Mechanistic") |> 
  filter(sbr < 1.6) |> # only keep values up to 2 bites per day
  ggplot(aes(color = Type, lty = Type)) +
  # Grey line for R0 = 1
  geom_hline(aes(yintercept = 1), color = "grey", lwd = 1) +
  # R0-GCD curves
  geom_line(aes(x = 1440 / theta, y = R0),
            lwd = 0.75) +
  # Add ticks for R0 = 1 crossing below the x-axis
  geom_rug(
    data = filter(Figure2_ticks, !(Type %in% c("Mechanistic~(lambda[Q])", "Mechanistic~(p[G])", "Mechanistic~(p[P])"))),
    aes(x = first_R0_greater_1),
    sides = "b", size = 0.75, outside = TRUE,
    length = unit(0.3, "in"),
    show.legend = F
  ) +
  scale_x_continuous(
    name = "",
    # limits = c(0,1.9),
    expand = c(0, 0, 0, 0.01),
    # breaks = seq(0, 2, by = 0.25)
  ) +
  scale_y_continuous(
    name = TeX("Basic reproduction number \\, [$R_0$]"),
    breaks = seq(0, 1.75, by = 0.25),
    limits = c(0, NA),
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
  coord_cartesian(
    xlim = c(0, 1.51),
    clip = "off" # needed to let ticks get plotted out of the axes
  ) + 
  theme_half_open(8) + 
  theme(
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.key.width = unit(0.35, "in")
    ) +
  guides(
    color = guide_legend(
      position = "top",
      direction = "horizontal",
      nrow = 1,
      byrow = T
    ),
    linetype = guide_legend(
      position = "top",
      direction = "horizontal",
      nrow = 1,
      byrow = T
    )
  )

Figure2_separated_mech <- Figure2_df |> 
  filter(`Model type` == "Mechanistic") |> 
  filter(sbr < 1.6) |> # only keep values up to 2 bites per day
  ggplot(aes(color = Type, lty = Type)) +
  # Grey line for R0 = 1
  geom_hline(aes(yintercept = 1), color = "grey", lwd = 1) +
  # R0-GCD curves
  geom_line(aes(x = 1440 / theta, y = R0),
            lwd = 0.75) +
  # Add ticks for R0 = 1 crossing below the x-axis
  geom_rug(
    data = filter(Figure2_ticks, Type %in% c("Mechanistic~(lambda[Q])", "Mechanistic~(p[G])", "Mechanistic~(p[P])")),
    aes(x = first_R0_greater_1),
    sides = "b", size = 0.75, outside = TRUE,
    length = unit(0.3, "in"),
    show.legend = F
  ) +
  # Arrows showing direction of increasing parameter values
  geom_segment(
    data = Figure2_arrows,
    aes(x = x_first, y = R0_first, xend = x_last, yend = R0_last, color = Type, lty = Type),
    lwd = 0, 
    # alpha = 0,
    arrow = arrow(length = unit(0.15, "inches"), ends = "last", type = "closed"),
    show.legend = F,
    inherit.aes = F
  ) +
  scale_x_continuous(
    name = TeX("Standard biting rate [Days$^{-1}$]"),
    # limits = c(0,1.9),
    expand = c(0, 0, 0, 0.01),
    # breaks = seq(0, 2, by = 0.25)
  ) +
  scale_y_continuous(
    name = TeX("Basic reproduction number \\, [$R_0$]"),
    breaks = seq(0, 1.75, by = 0.25),
    limits = c(0, NA),
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
  coord_cartesian(
    xlim = c(0, 1.51),
    clip = "off" # needed to let ticks get plotted out of the axes
  ) + 
  theme_half_open(8) + 
  theme(
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.key.width = unit(0.35, "in")
  ) +
  guides(
    color = guide_legend(
      position = "top",
      direction = "horizontal",
      nrow = 1,
      byrow = T
    ),
    linetype = guide_legend(
      position = "top",
      direction = "horizontal",
      nrow = 1,
      byrow = T
    )
  )

Figure2_separated_mech

Figure2_separated = egg::ggarrange(Figure2_separated_nonmech, Figure2_separated_mech, ncol = 1)

ggsave("figures/Figure2_separated.pdf", Figure2_separated, width = 7.5, height = 3.25 * 9/6.5, units = "in")

# Table: Characteristics of R0 curves ----
R0_characteristics_table <- Full_df |>
  filter(between(theta, 1, 20*1440)) |> # remove unrealistically short GCD (less than 1 second)
  mutate(sbr = 1440/theta) |>  # standard biting rate
  arrange(sbr) |> 
  mutate(
    Type = case_when(
      Type == "Mechanistic (lQ)" ~ "Mechanistic~(lambda[Q])",
      Type == "Mechanistic (pP)" ~ "Mechanistic~(p[P])",
      Type == "Mechanistic (pG)" ~ "Mechanistic~(p[G])",
      TRUE ~ Type
    )
  ) |> 
  # Get x-coordinate where R0 first exceeds/is less than one
  group_by(Type) |> 
  summarise(
    max_sbr = max(sbr),
    crit_min_sbr = sbr[R0 > 1][1],
    max_test = length(sbr[R0>1]),
    temp_max = max(sbr[R0 > 1], na.rm = T),
    crit_max_sbr = ifelse(temp_max > 0.99*max_sbr | max_test == 0, NA, temp_max),
    max_R0 = max(R0)
  ) |> 
  select(-c(max_sbr, max_test, temp_max)) |> 
  unique()

# 3. R0 vs. mechanistic parameters ----

Mech_df <- read_rds("data/Mechanistic_results.rds")


# Distinguish parameters by color -- different from 2.

nice_mech_labels = data.frame(
  pQ = "Seeking~success*','~p[Q]",
  pL = "Landing~success*','~p[L]",
  pP = "Probing~success*','~p[P]",
  pG = "Ingesting~success*','~p[G]",
  sigma = "Persistence~probability*','~sigma",
  lQ = "Seeking~rate*','~lambda[Q]",
  lL = "Landing~rate*','~lambda[L]",
  lP = "Probing~rate*','~lambda[P]",
  lG = "Ingesting~rate*','~lambda[G]",
  dummy = "Dummy~variable"
) %>% pivot_longer(everything(), values_to = "nice_labels")


Figure3_df <- Mech_df |> 
  filter(parameter_type == "varied") |> 
  select(-parameter_type) |> 
  pivot_longer(lQ:sigma) |> 
  filter(name == varied_parameter) |> 
  mutate(parameter_type = case_when(
    name %in% c("lQ", "lL", "lP", "lG") ~ "rate",
    name %in% c("pQ", "pL", "pP", "pG", "sigma") ~ "probability"
  )) |> 
  # Reduce the range for rates
  filter(!(parameter_type == "rate" & value >= (1/60))) |> 
  right_join(filter(nice_mech_labels, name != "dummy"))

Figure3_df$name <- factor(
  Figure3_df$name,
  levels = c("pQ", "pL", "pP", "pG", "sigma", "lQ", "lL", "lP", "lG"))

Figure3_df$nice_labels <- factor(
  Figure3_df$nice_labels,
  levels = filter(nice_mech_labels, name != "dummy")$nice_labels)

Figure3 <- Figure3_df |> 
  ggplot(aes(x = value, y = R0, linetype = mosquito_type, color = parameter_type)) +
  geom_hline(aes(yintercept = 1), color = "grey") +
  geom_line(lwd = 0.5) +
  facet_wrap(
    ~ nice_labels,
    ncol = 5,
    labeller = labeller(nice_labels = label_parsed),
    scales = "free_x"
  ) +
  scale_x_continuous(
    name = "",
    breaks = waiver(),
    n.breaks = 5,
    labels = function(x) sub("\\.?0+$", "", format(x, nsmall = 2)),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = TeX("Basic reproduction number \\, [$R_0$]"),
    labels = function(x) sub("\\.?0+$", "", format(x, nsmall = 2)),
    expand = c(0,0)
  ) +
  scale_color_manual(
    name = "Parameter type:",
    values = c4a("parks.saguaro",2)
    # values = c("#1B9E77", "#D952")
  ) +
  scale_linetype_discrete(
    name = "Mosquito type:",
    labels = c("Flighty", "Persistent")
  ) +
  theme_half_open(font_size = 8) +
  guides(
    color = guide_none()
  ) +
  theme(
    strip.background = element_rect(color = "white", fill = "white"),
    legend.key.width = unit(0.4, "in")
  )

shift_legend(Figure3)
ggsave("figures/Figure3.pdf", shift_legend(Figure3), width = 6.5, height = 2.25 * 9/6.5, units = "in")

# 4. PRCCs of R0 against mechanistic parameters ----

# Load in final PRCC results
PRCC_data = read_csv("data/julia_PRCC.csv") |>
  pivot_longer(cols = lQ:dummy, names_to = "input", values_to = "PRCC")


PRCC_plot_data <- PRCC_data |> 
  group_by(type, output) |> 
  # Add in nice labels
  left_join(
    param_table |> 
      rename(input = short_label)
  ) |>
  left_join(rename(nice_mech_labels, input = name), by = "input")  |> 
  mutate(
    output_label = case_when(
      output == "R0" ~ "Basic~reproduction~number~'['*R[0]*']'",
      output == "RVH" ~ "Probing~'mosquito-to-host'*~reproduction~number~'['*R[PH]*']'",
      output == "RHV" ~ "'Host-to-ingesting'*~mosquito~reproduction~number~'['*R[HG]*']'",
      output == "N_offspring" ~ "Basic~offspring~number~'['*N[0]*']'",
    ),
    type_label = case_when(
      type == "flighty" ~ "Flighty",
      type == "persistent" ~ "Persistent",
      type == "max" ~ "Maximum variation",
      type == "inv_flighty" ~ "Flighty (inverse)",
      type == "inv_persistent" ~ "Persistent (inverse)",
      type == "inv_max" ~ "Maximum variation (inverse)",
    ),
    dummy_min = -abs(PRCC[input == "dummy"]),
    dummy_max = abs(PRCC[input == "dummy"])
  ) |> 
  filter(!is.na(output), input != "dummy") |> 
  mutate(input_num = as.numeric(factor(input)))


PRCC_plot_data$type = factor(PRCC_plot_data$type, levels = c("flighty", "inv_flighty", "persistent", "inv_persistent", "max", "inv_max"))
PRCC_plot_data$input = factor(PRCC_plot_data$input, levels = c(
  "pQ", "pL", "pP", "pG", "sigma", "lQ", "lL", "lP", "lG", "dummy"
))


PRCC_plot_data$type_label = factor(PRCC_plot_data$type_label, levels = c(
  "Flighty", "Flighty (inverse)","Persistent (inverse)", "Persistent","Maximum variation", "Maximum variation (inverse)"
))

PRCC_plot_data$output_label = factor(PRCC_plot_data$output_label, levels = c(
  "Basic~offspring~number~'['*N[0]*']'",
  "Basic~reproduction~number~'['*R[0]*']'",
  "'Host-to-ingesting'*~mosquito~reproduction~number~'['*R[HG]*']'",
  "Probing~'mosquito-to-host'*~reproduction~number~'['*R[PH]*']'"
))

PRCC_plot_data$Label = factor(PRCC_plot_data$Label, levels = rev(c(
  c("Seeking success",
    "Landing success",  "Probing success", "Ingesting success"),
  "Persistence probability", 
  "Seeking rate", "Landing rate", "Probing rate", "Ingesting rate", "Dummy variable"
)))

PRCC_plot_data$nice_labels <- factor(
  PRCC_plot_data$nice_labels,
  levels = c(rev(nice_mech_labels$nice_labels)))

plot_ribbon <- PRCC_plot_data |>
  group_by(type, output_label) |> 
  group_modify( ~ bind_rows(
    tibble(
      input_num = min(.x$input_num) - 0.5,
      dummy_min = .x$dummy_min[1],
      dummy_max = .x$dummy_max[1],
      type = .x$type[1],
      output_label = .x$output_label[1],
    ),
    .x,
    tibble(
      input_num = max(.x$input_num) + 0.5,
      dummy_min = .x$dummy_min[nrow(.x)],
      dummy_max = .x$dummy_max[nrow(.x)],
      type = .x$type[nrow(.x)],
      output_label = .x$output_label[nrow(.x)],
    )
  )
  ) |> ungroup() #|> 
  # select(input_num, 
  #        dummy_min, dummy_max,
  #        type, output_label) |> distinct()

PRCC_plots <- PRCC_plot_data |> 
  # filter(type_label %in% c("Maximum variation", "Maximum variation (inverse)")) |> 
  arrange(input) |> 
  # filter((input %in% c("lP"))) |>
  filter(output %in% c("N_offspring","R0")) |>
  # filter(type %in% c("flighty", "persistent")) |> 
  ggplot() +
  geom_col(
    aes(x = nice_labels, y = PRCC, fill = type_label),
    position = "dodge", alpha = 0.75
  ) +
  # Add light grey lines to divide up categories
  geom_vline(xintercept = seq(1.5, length(levels(PRCC_plot_data$input)) - 0.5, by = 1),
             color = "grey60", linetype = "dashed", linewidth = 0.25) +
  # Add zero line
  geom_hline(yintercept = 0, color = "black", linewidth =0.5 ) +
  # Add grey ribbon to show dummy variable values
  geom_ribbon(
    data = plot_ribbon |> filter(output_label %in% c(TeX("Basic offspring number [$N_{0}$]"),TeX("Basic reproduction number [$R_{0}$]"))),
    aes(
      x = input_num,
      ymin = dummy_min, ymax = dummy_max,
      group = output_label),
    color = NA, fill = "grey60",
    alpha = 0.5
  ) +
  facet_wrap( ~ output_label, ncol = 1, scales = "free_y", labeller = label_parsed) +
  scale_fill_manual(
    name = "Parameter set:",
    values = c(c4a("met.juarez", 6))#, "black")
  ) +
  scale_x_discrete(
    name = "",
    labels = function(x) parse(text = x),
  ) +
  scale_y_continuous(
    TeX("Partial Rank Correlation Coefficient"),
    limits = c(-1,1)
  ) +
  theme_half_open(11) +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1),
    strip.background = element_rect(color = "white", fill = "white"),
    legend.key.width = unit(0.4, "in"),
    legend.position = "top",
    legend.direction = "horizontal"
  ) +
  coord_cartesian(xlim = c(1.1, 8.9))

PRCC_plots

ggsave("figures/Figure4.pdf", PRCC_plots, width = 7, height = 3.25 * 9/6.5, units = "in")

PRCC_plots_row <- PRCC_plot_data |> 
  arrange(input) |> 
  filter(output %in% c("N_offspring","R0")) |>
  ggplot() +
  geom_col(aes(x = nice_labels, y = PRCC, fill = type_label), position = "dodge") +
  # Add light grey lines to divide up categories
  geom_vline(xintercept = seq(1.5, length(levels(PRCC_plot_data$input)) - 0.5, by = 1), 
             color = "grey60", linetype = "dashed", linewidth = 0.25) +
  # Add zero line
  geom_hline(yintercept = 0, color = "black", linewidth =1 ) +
  # Add grey ribbon to show dummy variable values
  geom_ribbon(
    data = plot_ribbon |> filter(output_label %in% c(TeX("Basic offspring number [$N_{0}$]"),TeX("Basic reproduction number [$R_{0}$]"))),
    aes(x = input_num,
        ymin = dummy_min, ymax = dummy_max,
        group = output_label),
    color = NA, fill = "grey60",
    alpha = 0.5
  ) +
  facet_wrap(~output_label, nrow = 1, scales = "free_x") +
  # scale_fill_manual(
  #   name = "Parameter set:",
  #   values = c(c4a("met.juarez",3))#, "black")
  # ) +
  scale_x_discrete(
    name = "",
    labels = function(x) parse(text = x)
  ) +
  scale_y_continuous(
    TeX("Partial Rank Correlation Coefficient"),
    limits = c(-1,1)
  ) +
  theme_half_open(11) +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1),
    strip.background = element_rect(color = "white", fill = "white"),
    legend.key.width = unit(0.4, "in"),
    legend.position = "top",
    legend.direction = "horizontal",
    # legend.justification = "center"
  )

PRCC_plots_row

ggsave("figures/Figure4_row.pdf", PRCC_plots_row, width = 7, height = 3.25 * 9/6.5, units = "in")

PRCC_plots_max_only <- PRCC_plot_data |> 
  # Just keep maximum variation for now
  filter(type == "max", input != "dummy") |>
  arrange(input) |> 
  filter(output %in% c("N_offspring","R0")) |>
  group_by(type, output) |> 
  mutate(star_flag = abs(PRCC) < abs(dummy_min),
         star_xpos = PRCC + sign(PRCC) * 0.025) |>
  # Plot
  ggplot() +
  geom_col(aes(y = nice_labels, x = PRCC, fill = output_label), position = "dodge", width = 0.75) +
  # Add light grey lines to divide up categories
  geom_hline(yintercept = seq(1.5, length(levels(PRCC_plot_data$input)) - 0.5, by = 1), 
             color = "grey60", linetype = "dashed", linewidth = 0.125) +  # Add grey ribbon to show dummy variable values
  geom_ribbon(
    data = plot_ribbon |> 
      filter(output_label %in% c(TeX("Basic offspring number [$N_{0}$]"),TeX("Basic reproduction number [$R_{0}$]")), type == "max", input !="dummy"),
    aes(y = nice_labels,
        xmin = dummy_min, xmax = dummy_max,
        group = output_label),
    color = NA, fill = "grey60",
    alpha = 0.5
  ) +
  # Add zero line
  geom_vline(xintercept = 0, color = "black", lwd = 0.25) +
  scale_fill_manual(
    name = "",
    values = c(c4a("met.juarez",3))
  ) +
  # Add stars to flagged bars
  geom_text(
    data = . %>% filter(star_flag),
    aes(y = nice_labels, x = star_xpos, label = "n.s."),
    position = position_dodge(width = 0.75),
    size = 2, vjust = -0.55
  ) +
  scale_y_discrete(
    name = "",
    labels = function(x) parse(text = x)
  ) +
  scale_x_continuous(
    TeX("Partial Rank Correlation Coefficient"),
    breaks = seq(-1,1,by = 0.25),
    expand = c(0.05,0)
    # limits = c(-1,1)
  ) +
  theme_half_open(11) +
  theme(
    strip.background = element_rect(color = "white", fill = "white"),
    legend.key.width = unit(0.25, "in"),
    legend.key.height = unit(0.03125, "in"),
    legend.position = "top",
    legend.direction = "horizontal"
  )

PRCC_plots_max_only

ggsave("figures/Figure4_max_only.pdf", PRCC_plots_max_only, width = 6.5, height = 2.25 * 9/6.5, units = "in")

PRCC_plots_max_only <- PRCC_plot_data |> 
  # filter(index_type == "S1") |>
  # Just keep maximum variation for now
  filter(type == "max") |>
  # filter(input != "dummy") |> 
  arrange(input) |> 
  # filter(output %in% c("N_offspring","R0")) |>
  group_by(type, output) |> 
  mutate(
    star_flag = abs(PRCC) < min(abs(dummy_min), abs(dummy_max)),
    star_label = ifelse(star_flag, "*", ""),
    star_xpos = PRCC + 2 * dummy_min + 0.1
  ) %>% 
  # Plot
  ggplot() +
  geom_col(aes(y = nice_labels, x = PRCC, fill = output_label), color = "black", position = "identity", width = 0.75, linewidth = 0.25) +
  # Add light grey lines to divide up categories
  geom_hline(yintercept = seq(1.5, length(levels(PRCC_plot_data$input)) - 0.5, by = 1), 
             color = "grey60", linetype = "dashed", linewidth = 0.125) +
  # Add zero line
  geom_vline(xintercept = 0, color = "black", lwd = 0.25) +  
  # Add n.s. to flagged bars
  geom_text(
    data = . %>% filter(star_flag),
    aes(y = nice_labels, x = star_xpos, label = "n.s."),
    position = position_dodge(width = 0.75),
    size = 2, 
    # vjust = -0.55
  ) +
  scale_fill_manual(
    name = "",
    values = c(c4a("met.juarez"))
  ) +
  scale_y_discrete(
    name = "",
    labels = function(x) parse(text = x)
  ) +
  scale_x_continuous(
    "",
    labels = label_number(accuracy = 0.01, trim = TRUE, drop0trailing = TRUE),
    expand = expansion(mult = c(0.01, 0))
  ) +
  theme_half_open(8) +
  facet_wrap(~ output_label, ncol = 2, labeller = label_parsed) +
  guides(
    # alpha = guide_none(),
    # linetype = guide_none(),
    fill = guide_none()
  ) + 
  theme(
    strip.background = element_rect(color = "white", fill = "white"),
    legend.key.width = unit(0.25, "in"),
    legend.key.height = unit(0.03125, "in"),
    legend.position = "top",
    legend.direction = "horizontal"
  ) #+
  # ggtitle("Partial Rank Correlation Coefficients")

PRCC_plots_max_only

ggsave("figures/PRCCFigure4_max_only.pdf", PRCC_plots_max_only, width = 6.5, height = 2.25 * 9/6.5, units = "in")

# 5. eFAST indices of N0 and R0 wrt parameters ----
# Load in data

eFAST_data = read_csv("data/new_eFAST_test.csv") |> 
  mutate(type = factor(type, levels = c("max", "flighty", "persistent", "inv_persistent", "inv_flighty")),
         param_type = ifelse(input %in% c("sigma", "pQ", "pL", "pP", "pG"), "probability", "rate"))

# eFAST_stability_plot = eFAST_data |> 
#   filter(sample_size > 10000) |> 
#   group_by(input, output) |> 
#   arrange(sample_size) |> 
#   ggplot(aes(x = sample_size, y = value, color = input)) +
#   geom_path(lwd = 1, alpha = 0.75) +
#   facet_wrap(output~type+index_type, scales = "free", ncol = 6) +
#   scale_x_continuous(
#     "Sample size for sensitivity analysis",
#     labels = scales::label_scientific(),
#     n.breaks = 3,
#     trans = 'log2'
#   ) +
#   scale_y_continuous("") + 
#   scale_color_discrete("Parameter") +
#   theme_minimal()
# eFAST_stability_plot

# eFAST sensitivity plots
plot_data <- eFAST_data |>
  # filter(sample_size == max(sample_size)) |> 
  group_by(type, output, index_type) |> 
  # Add in nice labels
  left_join(
    param_table |> 
      rename(input = short_label)
  ) |>
  left_join(rename(nice_mech_labels, input = name), by = "input")  |> 
  mutate(
    output_label = case_when(
      output == "R0" ~ "Basic~reproduction~number~'['*R[0]*']'",
      output == "RVH" ~ "Probing~'mosquito-to-host'*~reproduction~number~'['*R[PH]*']'",
      output == "RHV" ~ "'Host-to-ingesting'*~mosquito~reproduction~number~'['*R[HG]*']'",
      output == "N_offspring" ~ "Basic~offspring~number~'['*N[0]*']'",
    ),
    type_label = case_when(
      type == "flighty" ~ "Flighty",
      type == "persistent" ~ "Persistent",
      type == "max" ~ "Maximum variation",
      type == "inv_persistent" ~ "Persistent (durations)",
      type == "inv_flighty" ~ "Flighty (durations)",
    ),
    dummy_min = -abs(mean_value[input == "dummy"]),
    dummy_max = abs(mean_value[input == "dummy"])
  ) %>% 
  # filter(!is.na(output), input != "dummy") %>% 
  mutate(input_num = as.numeric(factor(input))) |> 
  mutate(
    index_type_label = case_when(
      index_type == "S1" ~ "First order eFAST sensitivity index",
      index_type == "ST" ~ "Total order eFAST sensitivity index"
    ),
    index_type_shortlabel = case_when(
      index_type == "S1" ~ "First order",
      index_type == "ST" ~ "Total order"
    )
  )


plot_data$type = factor(plot_data$type, levels = c("max", "flighty", "persistent", "inv_persistent", "inv_flighty"))
plot_data$input = factor(plot_data$input, levels = c(
  "pQ", "pL", "pP", "pG", "sigma", "lQ", "lL", "lP", "lG", "dummy"
))
plot_data$output_label = factor(plot_data$output_label, levels = c(
  "Basic~offspring~number~'['*N[0]*']'",
  "Basic~reproduction~number~'['*R[0]*']'",
  "'Host-to-ingesting'*~mosquito~reproduction~number~'['*R[HG]*']'",
  "Probing~'mosquito-to-host'*~reproduction~number~'['*R[PH]*']'"
))

plot_data$type_label = factor(plot_data$type_label, levels = c(
  "Maximum variation", "Flighty", "Flighty (durations)", "Persistent", "Persistent (durations)"
))

plot_data$index_type_label = factor(plot_data$index_type_label, levels = c("First order eFAST sensitivity index", "Total order eFAST sensitivity index"))

plot_data$index_type_shortlabel = factor(plot_data$index_type_shortlabel, levels = c("First order", "Total order"))

plot_data$Label = factor(plot_data$Label, levels = (c(
  c("Seeking success",
    "Landing success",  "Probing success", "Ingesting success"),
  "Persistence probability", 
  "Seeking rate", "Landing rate", "Probing rate", "Ingesting rate", "Dummy variable"
)))

plot_data$nice_labels <- factor(
  plot_data$nice_labels,
  levels = c(rev(nice_mech_labels$nice_labels)))

# plot_ribbon <- plot_data |>
#   group_by(type, output_label) |> 
#   group_modify( ~ bind_rows(
#     tibble(
#       input_num = min(.x$input_num) - 0.5,
#       dummy_min = .x$dummy_min[1],
#       dummy_max = .x$dummy_max[1],
#       type = .x$type[1],
#       output_label = .x$output_label[1],
#     ),
#     .x,
#     tibble(
#       input_num = max(.x$input_num) + 0.5,
#       dummy_min = .x$dummy_min[nrow(.x)],
#       dummy_max = .x$dummy_max[nrow(.x)],
#       type = .x$type[nrow(.x)],
#       output_label = .x$output_label[nrow(.x)],
#     )
#   )
#   ) %>% ungroup() %>% 
#   # select(input_num, type, output_label) %>% 
#   distinct()

FirsteFAST_plots <- plot_data |> 
  filter(index_type == "S1") |>
  arrange(input) |> 
  # filter((input %in% c("lP"))) |>
  # filter(output %in% c("N_offspring","R0")) |>
  # filter(type %in% c("flighty", "persistent")) |> 
  ggplot() +
  geom_col(
    aes(x = nice_labels, y = mean_value, fill = type_label),
    position = "dodge", alpha = 0.75
  ) +
  # Add light grey lines to divide up categories
  geom_vline(xintercept = seq(1.5, length(levels(plot_data$input)) - 0.5, by = 1),
             color = "grey60", linetype = "dashed", linewidth = 0.25) +
  # Add zero line
  geom_hline(yintercept = 0, color = "black", linewidth =0.5 ) +
  facet_wrap( ~ output_label, ncol = 1, scales = "free_y") +
  scale_fill_manual(
    name = "Parameter set:",
    values = c(c4a("met.juarez"))#, "black")
  ) +
  scale_x_discrete(
    name = "",
    labels = function(x) parse(text = x),
  ) +
  scale_y_continuous(
    TeX("First Order eFAST Sensitivity Index"),
    expand = c(0.0,0)
  ) +
  theme_half_open(11) +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1),
    strip.background = element_rect(color = "white", fill = "white"),
    legend.key.width = unit(0.4, "in"),
    legend.position = "top",
    legend.direction = "horizontal"
  ) +
  coord_cartesian(xlim = c(1.1, 8.9))

FirsteFAST_plots

ggsave("figures/FirsteFASTFigure4.pdf", FirsteFAST_plots, width = 7, height = 3.25 * 9/6.5, units = "in")

TotaleFAST_plots <- plot_data |> 
  filter(index_type == "ST") |>
  arrange(input) |> 
  # filter((input %in% c("lP"))) |>
  # filter(output %in% c("N_offspring","R0")) |>
  # filter(type %in% c("flighty", "persistent")) |> 
  ggplot() +
  geom_col(
    aes(x = nice_labels, y = mean_value, fill = type_label),
    position = "dodge", alpha = 0.75
  ) +
  # Add light grey lines to divide up categories
  geom_vline(xintercept = seq(1.5, length(levels(plot_data$input)) - 0.5, by = 1),
             color = "grey60", linetype = "dashed", linewidth = 0.25) +
  # Add zero line
  geom_hline(yintercept = 0, color = "black", linewidth =0.5 ) +
  facet_wrap( ~ output_label, ncol = 1, scales = "free_y") +
  scale_fill_manual(
    name = "Parameter set:",
    values = c(c4a("met.juarez"))#, "black")
  ) +
  scale_x_discrete(
    name = "",
    labels = function(x) parse(text = x),
  ) +
  scale_y_continuous(
    TeX("Total eFAST Sensitivity Index"),
    expand = c(0.0,0)
  ) +
  theme_half_open(11) +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1),
    strip.background = element_rect(color = "white", fill = "white"),
    legend.key.width = unit(0.4, "in"),
    legend.position = "top",
    legend.direction = "horizontal"
  ) +
  coord_cartesian(xlim = c(1.1, 8.9))

TotaleFAST_plots

ggsave("figures/TotaleFASTFigure4.pdf", TotaleFAST_plots, width = 7, height = 3.25 * 9/6.5, units = "in")

FirsteFAST_plots_row <- plot_data |> 
  filter(index_type == "S1") |>
  arrange(input) |> 
  # filter(output %in% c("N_offspring","R0")) |>
  ggplot() +
  geom_col(aes(x = nice_labels, y = mean_value, fill = type_label), position = "dodge") +
  # Add light grey lines to divide up categories
  geom_vline(xintercept = seq(1.5, length(levels(plot_data$input)) - 0.5, by = 1), 
             color = "grey60", linetype = "dashed", linewidth = 0.25) +
  # Add zero line
  geom_hline(yintercept = 0, color = "black", linewidth =1 ) +
  # # Add grey ribbon to show dummy variable values
  # geom_ribbon(
  #   data = plot_ribbon |> filter(output_label %in% c(TeX("Basic offspring number [$N_{0}$]"),TeX("Basic reproduction number [$R_{0}$]")),
  #   aes(x = input_num, ymin = dummy_min, ymax = dummy_max, group = output_label),
  #   color = NA, fill = "grey60",
  #   alpha = 0.5
  # ) +
  facet_wrap(~output_label, nrow = 1, scales = "free_x") +
  scale_fill_manual(
    name = "Parameter set:",
    values = c(c4a("met.juarez"))#, "black")
  ) +
  scale_x_discrete(
    name = "",
    labels = function(x) parse(text = x)
  ) +
  scale_y_continuous(
    TeX("First Order eFAST Sensitivity Index"),
    expand = c(0.0,0)
  ) +
  theme_half_open(11) +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1),
    strip.background = element_rect(color = "white", fill = "white"),
    legend.key.width = unit(0.4, "in"),
    legend.position = "top",
    legend.direction = "horizontal",
    # legend.justification = "center"
  )

FirsteFAST_plots_row

ggsave("figures/FirsteFASTFigure4_row.pdf", FirsteFAST_plots_row, width = 7, height = 3.25 * 9/6.5, units = "in")

TotaleFAST_plots_row <- plot_data |> 
  filter(index_type == "ST") |>
  arrange(input) |> 
  # filter(output %in% c("N_offspring","R0")) |>
  ggplot() +
  geom_col(aes(x = nice_labels, y = mean_value, fill = type_label), position = "dodge") +
  # Add light grey lines to divide up categories
  geom_vline(xintercept = seq(1.5, length(levels(plot_data$input)) - 0.5, by = 1), 
             color = "grey60", linetype = "dashed", linewidth = 0.25) +
  # Add zero line
  geom_hline(yintercept = 0, color = "black", linewidth =1 ) +
  # # Add grey ribbon to show dummy variable values
  # geom_ribbon(
  #   data = plot_ribbon |> filter(output_label %in% c(TeX("Basic offspring number [$N_{0}$]"),TeX("Basic reproduction number [$R_{0}$]")),
  #   aes(x = input_num, ymin = dummy_min, ymax = dummy_max, group = output_label),
  #   color = NA, fill = "grey60",
  #   alpha = 0.5
  # ) +
  facet_wrap(~output_label, nrow = 1, scales = "free_x") +
  scale_fill_manual(
    name = "Parameter set:",
    values = c(c4a("met.juarez"))#, "black")
  ) +
  scale_x_discrete(
    name = "",
    labels = function(x) parse(text = x)
  ) +
  scale_y_continuous(
    TeX("Total eFAST Sensitivity Index"),
    expand = c(0.0,0)
  ) +
  theme_half_open(11) +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1),
    strip.background = element_rect(color = "white", fill = "white"),
    legend.key.width = unit(0.4, "in"),
    legend.position = "top",
    legend.direction = "horizontal",
    # legend.justification = "center"
  )

TotaleFAST_plots_row

ggsave("figures/TotaleFASTFigure4_row.pdf", TotaleFAST_plots_row, width = 7, height = 3.25 * 9/6.5, units = "in")

FirsteFAST_plots_max_only <- plot_data |> 
  # filter(index_type == "S1") |>
  # Just keep maximum variation for now
  filter(type == "max") |>
  filter(input != "dummy") |>
  arrange(input) |> 
  # filter(output %in% c("N_offspring","R0")) |>
  group_by(type, output) |> 
  mutate(
    star_flag = p_value > 0.01,
    star_label = ifelse(star_flag, "n.s.", ""),
    star_xpos = mean_value + 2 * std_value + 0.05
  ) %>% 
  # Plot
  ggplot() +
  geom_col(aes(y = nice_labels, x = mean_value, linetype = index_type_shortlabel, alpha = index_type_shortlabel, fill = output_label), color = "black", position = "identity", width = 0.75, linewidth = 0.25) +
  geom_errorbar(
    aes(y = nice_labels, xmin = mean_value - 2 * std_value, xmax = mean_value + 2 * std_value, group = output_label),
    width = 0.25, color = "black", position = position_dodge(width = 0.9), linewidth = 0.25
  ) +
  # Add light grey lines to divide up categories
  geom_hline(yintercept = seq(1.5, length(levels(plot_data$input)) - 0.5, by = 1), 
             color = "grey60", linetype = "dashed", linewidth = 0.125) +
  # Add zero line
  geom_vline(xintercept = 0, color = "black", lwd = 0.25) +
  scale_fill_manual(
    name = "",
    values = c(c4a("met.juarez"))
  ) +
  # Add stars to flagged bars
  geom_text(
    aes(y = nice_labels, x = star_xpos, label = star_label, group = output_label),
    position = position_dodge(width = 0.9),
    size = 2, vjust = 0.75
  ) +
  scale_y_discrete(
    name = "",
    labels = function(x) parse(text = x)
  ) +
  scale_x_continuous(
    "",
    # TeX("First Order eFAST Sensitivity Index"),
    labels = label_number(accuracy = 0.01, trim = TRUE, drop0trailing = TRUE),
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_alpha_manual(
    "",
    values = c(1, 0.5)
  ) +
  scale_linetype_manual(
    "",
    values = c(1, 2)
  ) +
  theme_half_open(8) +
  facet_wrap(~ output_label, label = label_parsed) +
  guides(
    alpha = guide_legend(position = "bottom"),
    linetype = guide_legend(position = "bottom"),
    fill = guide_none()
  ) + 
  theme(
    strip.background = element_rect(color = "white", fill = "white"),
    legend.key.width = unit(0.25, "in"),
    legend.key.height = unit(0.03125, "in"),
    legend.position = "top",
    legend.direction = "horizontal"
  ) #+
  # ggtitle("eFAST Sensitivity Indices")

FirsteFAST_plots_max_only

ggsave("figures/eFASTFigure4_max_only.pdf", FirsteFAST_plots_max_only, width = 6.5, height = 2.25 * 9/6.5, units = "in")

TotaleFAST_plots_max_only <- plot_data |> 
  filter(index_type == "ST") |>
  # Just keep maximum variation for now
  filter(type == "max") |>
  arrange(input) |> 
  # filter(output %in% c("N_offspring","R0")) |>
  group_by(type, output) |> 
  mutate(
    star_flag = p_value < 0.01,
    star_label = ifelse(star_flag, "*", ""),
    star_xpos = mean_value + 2 * std_value + 0.01
  ) %>% 
  # Plot
  ggplot() +
  geom_col(aes(y = nice_labels, x = mean_value, fill = output_label), position = position_dodge(width = 1), width = 0.75) +
  geom_errorbar(
    aes(y = nice_labels, xmin = mean_value - 2 * std_value, xmax = mean_value + 2 * std_value, group = output_label),
    width = 0.5, color = "black", position = position_dodge(width = 1), linewidth = 0.25
  ) +
  # Add light grey lines to divide up categories
  geom_hline(yintercept = seq(1.5, length(levels(plot_data$input)) - 0.5, by = 1), 
             color = "grey60", linetype = "dashed", linewidth = 0.125) +
  # Add zero line
  geom_vline(xintercept = 0, color = "black", lwd = 0.25) +
  scale_fill_manual(
    name = "",
    values = c(c4a("met.juarez",4))
  ) +
  # Add stars to flagged bars
  geom_text(
    aes(y = nice_labels, x = star_xpos, label = star_label, group = output_label),
    position = position_dodge(width = 1),
    size = 4, vjust = 0.75
  ) +
  scale_y_discrete(
    name = "",
    labels = function(x) parse(text = x)
  ) +
  scale_x_continuous(
    TeX("Total eFAST Sensitivity Index"),
    labels = label_number(accuracy = 0.01, trim = TRUE),
    expand = expansion(mult = c(0, 0.01))
  ) +
  theme_half_open(11) +
  theme(
    strip.background = element_rect(color = "white", fill = "white"),
    legend.key.width = unit(0.25, "in"),
    legend.key.height = unit(0.03125, "in"),
    legend.position = "top",
    legend.direction = "horizontal"
  )

TotaleFAST_plots_max_only

ggsave("figures/TotaleFASTFigure4_max_only.pdf", TotaleFAST_plots_max_only, width = 6.5, height = 2.25 * 9/6.5, units = "in")

# Supplementary: theta vs. mechanistic parameters ----
SuppFigure1 <- Figure3_df |>
  filter(theta < 21 * 1440) |> 
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
    name = TeX("Gonotrophic cycle duration \\, [Days]"),
    expand = c(0,0)
  ) +
  scale_color_discrete(
    name = "Parameter type:"
  ) +
  scale_linetype_discrete(
    name = "Mosquito type:"
  ) +
  theme_half_open(18)

shift_legend(SuppFigure1)
ggsave("figures/SuppFigure.pdf", shift_legend(SuppFigure1), width = 20, height = 8, units = "in")
# Distinguish parameters by color


# Supplementary: Parameter set illustrations ------------------------------

# Create plots illustrating the ranges for the three parameter sets

max_lbs = c(1/((1/2)*1440.0), 1/(30), 1/(30), 1/(30), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
max_ubs = c(1/160, 2.0, 2.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 100.0)
flighty_baseline = c(1/480, 1/10, 1/5, 1/1, 1 - 0.9, 1.0, 0.5, 0.5, 0.5, 100.0)
persistent_baseline = c(1/480, 1/10, 1/5, 1/1, 1 - 0.66, 1.0, 0.7, 0.8, 0.9, 100.0)
stretch_val = 0.2
flighty_lbs = (1 - stretch_val) * flighty_baseline
flighty_ubs = (1 + stretch_val) * flighty_baseline
persistent_lbs = (1 - stretch_val) * persistent_baseline
persistent_ubs = (1 + stretch_val) * persistent_baseline
paramset_ranges <- expand_grid(
  paramset = c("max", "flighty", "persistent"),
  params = factor(unique(eFAST_data$input), levels = rev(unique(eFAST_data$input)))
)
paramset_ranges$lbs = c(max_lbs, flighty_lbs, persistent_lbs)
paramset_ranges$ubs = c(max_ubs, flighty_ubs, persistent_ubs)

paramset_ranges =  paramset_ranges |> 
  rowwise() |> 
  mutate(
    paramset = factor(paramset, levels = c("max", "flighty", "persistent")),
    param_type = ifelse(params %in% c("sigma", "pQ", "pL", "pP", "pG"),
                        "probability",
                        "rate"),
    lbs = ifelse(param_type == "probability",
                 max(0, min(lbs, 1)),
                 lbs),
    ubs = ifelse(param_type == "probability",
                 max(0, min(ubs, 1)),
                 ubs)
  ) |> 
  group_by(params, paramset, param_type) |> 
  mutate(
    row_num = cur_group_id(),
    param_name = ifelse(row_num %% 3 == 2, unfactor(params), " ")
  ) |> 
  ungroup() |> 
  arrange(params, paramset)

paramset_prob_plot <- paramset_ranges |> 
  filter(param_type == "probability") |> 
  ggplot(aes(y = factor(row_num), color = paramset)) +
  geom_linerange(
    aes(xmin = lbs, xmax = ubs),
    lwd = 8
  ) +  
  # Add light grey lines to divide up categories
  geom_hline(yintercept = seq(3.5, 14.5, by = 3), 
             color = "grey60", linetype = "dashed", linewidth = 0.5) +
  scale_x_continuous(
    "Probability",
    labels = scales::label_percent(),
    limits = c(0,NA),
    expand = c(0,0)
  ) +
  scale_y_discrete(
    "",
    labels = paramset_ranges$param_name[1:15]
  ) +
  guides(
    color = guide_none()
  ) +
  ggtitle("Probabilities") +
  theme_minimal_vgrid()

paramset_rate_plot <- paramset_ranges |> 
  filter(param_type == "rate") |> 
  ggplot(aes(y = (factor(row_num)), color = paramset)) +
  geom_linerange(
    aes(xmin = 1/ubs, xmax = 1/lbs),
    lwd = 8
  ) +  
  # Add light grey lines to divide up categories
  geom_hline(yintercept = seq(3.5, 14.5, by = 3), 
             color = "grey60", linetype = "dashed", linewidth = 0.5) +
  scale_x_log10(
    "Minutes",
    expand = c(0,0)
  ) +
  scale_y_discrete(
    "",
    labels = (paramset_ranges$param_name[16:27])
  ) +
  ggtitle("Inverted rates") +
  theme_minimal_vgrid()

paramset_plot = egg::ggarrange(paramset_prob_plot, paramset_rate_plot, nrow = 1)

