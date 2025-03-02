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
  Label = c("Seeking success probability", "Seeking rate", "Landing success probability", "Landing rate", "Probing success probability", "Probing rate", "Ingestion success probability", "Ingesting rate", "Persistence probability", "Biting rate (exponential)"), 
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
  pnls <- cowplot::plot_to_gtable(p) %>%
    gtable::gtable_filter("panel") %>%
    with(setNames(grobs, layout$name)) %>%
    purrr::keep(~ identical(.x, zeroGrob()))
  
  if (length(pnls) == 0) stop("No empty facets in the plot")
  
  lemon::reposition_legend(p, "center", panel = names(pnls))
}

# Plot parameters ----

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
  mutate(theta_label = case_when(
    closest_theta == 360 ~ "A",
    closest_theta == 720 ~ "B",
    closest_theta == 1440 ~ "C",
    closest_theta == 2880 ~ "D"
  )) %>%
  # Add in values to plug into the pdf
  cross_join(tibble(x = seq(0, pdf_max, length.out = resolution)))  %>% 
  filter(case_when(
    closest_theta == 360 ~ x < 0.75 * 1440,
    closest_theta == 720 ~ x < 1.5 * 1440,
    closest_theta == 1440 ~ x < 3 * 1440,
    closest_theta == 2880 ~ x < 6 * 1440
  )) %>%
  # Calculate values of the pdf
  rowwise() %>% 
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

facet_labels <- Figure1_df %>%
  group_by(closest_theta) %>% 
  filter(x > 0) %>% 
  mutate(label = LETTERS[1:n()],  # Assign "A", "B", "C", ...
         x_pos = 360*min(x, na.rm = TRUE)/1440,  # Align left
         y_pos = max(pdf_val, na.rm = TRUE) * .95) %>%  # Slightly above max y
  distinct(closest_theta, x_pos, y_pos, theta_label)


quick_fig_1 <- function(keep_names) {
  Fig_out = Figure1_df %>%
    filter(Type %in% keep_names) %>% 
    filter(closest_theta == 360) %>% 
    # filter(closest_theta < 2 * 1440) %>% 
    ggplot(aes(x = 24 * x / (1440), y = pdf_val, color = Type)) +
    # Plot grey dotted line showing the mean
    geom_vline(
      aes(xintercept = 24 * closest_theta/1440),
      lwd = 0.5, color = "grey", lty = 2
    ) +
    # Plot pdfs
    geom_line(
      lwd = 1
    ) + 
    scale_x_continuous(
      name = "Gonotrophic cycle duration [Hours]",
      limits =  c(0, NA),
      expand = c(0.0,0),
      breaks = seq(0, 72, by = 6)
    ) +
    scale_y_continuous(
      name = "Density",
      limits = c(0, max(Figure1_df$pdf_val)),
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
    theme_half_open(18) +
    theme(
      strip.background = element_rect(color = "white", fill = "white"),
      strip.text = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  Fig_out
}

Figure_quarter_1 = quick_fig_1(unique(Figure1_df$Type)[1])
Figure_quarter_2 = quick_fig_1(unique(Figure1_df$Type)[1:2])
Figure_quarter_3 = quick_fig_1(unique(Figure1_df$Type)[1:3])
Figure_quarter_4 = quick_fig_1(unique(Figure1_df$Type)[1:4])

ggsave("presentations/figures/Distribution_quarter_1.png", Figure_quarter_1, 
       width = 13.333, height = 7.5, units = "in",
       dpi = 1200)
ggsave("presentations/figures/Distribution_quarter_2.png", Figure_quarter_2, 
       width = 13.333, height = 7.5, units = "in",
       dpi = 1200)
ggsave("presentations/figures/Distribution_quarter_3.png", Figure_quarter_3, 
       width = 13.333, height = 7.5, units = "in",
       dpi = 1200)
ggsave("presentations/figures/Distribution_quarter_4.png", Figure_quarter_4, 
       width = 13.333, height = 7.5, units = "in",
       dpi = 1200)

# pdf figure
Figure_quarter = Figure1_df %>%
  filter(closest_theta == 360) %>% 
  # filter(closest_theta < 2 * 1440) %>% 
  ggplot(aes(x = 24 * x / (1440), y = pdf_val, color = Type)) +
  # Plot grey dotted line showing the mean
  geom_vline(
    aes(xintercept = 24 * closest_theta/1440),
    lwd = 0.5, color = "grey", lty = 2
  ) +
  # Plot pdfs
  geom_line(
    lwd = 1
  ) + 
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
  theme_half_open(18) +
  theme(
    strip.background = element_rect(color = "white", fill = "white"),
    strip.text = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

Figure_quarter

ggsave("presentations/figures/Distribution_quarter.png", Figure_quarter, 
       width = 13.333, height = 7.5, units = "in",
       dpi = 1200)

# Distribution for half a day
Figure_half = Figure1_df %>%
  filter(closest_theta == 720) %>% 
  # filter(closest_theta < 2 * 1440) %>% 
  ggplot(aes(x = 24 * x / (1440), y = pdf_val, color = Type)) +
  # Plot grey dotted line showing the mean
  geom_vline(
    aes(xintercept = 24 * closest_theta/1440),
    lwd = 0.5, color = "grey", lty = 2
  ) +
  # Plot pdfs
  geom_line(
    lwd = 1
  ) + 
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
  theme_half_open(18) +
  theme(
    strip.background = element_rect(color = "white", fill = "white"),
    strip.text = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

Figure_half

ggsave("presentations/figures/Distribution_half.png", Figure_half, 
       width = 13.333, height = 7.5, units = "in",
       dpi = 1200)

# Distribution for one day
Figure_one = Figure1_df %>%
  filter(closest_theta == 1440) %>% 
  # filter(closest_theta < 2 * 1440) %>% 
  ggplot(aes(x = 24 * x / (1440), y = pdf_val, color = Type)) +
  # Plot grey dotted line showing the mean
  geom_vline(
    aes(xintercept = 24 * closest_theta/1440),
    lwd = 0.5, color = "grey", lty = 2
  ) +
  # Plot pdfs
  geom_line(
    lwd = 1
  ) + 
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
  theme_half_open(18) +
  theme(
    strip.background = element_rect(color = "white", fill = "white"),
    strip.text = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

Figure_one

ggsave("presentations/figures/Distribution_one.png", 
       width = 13.333, height = 7.5, units = "in",
       dpi = 1200)




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
  filter(theta > (1/2.01) * 1440) %>%
  mutate(Type = case_when(
    Type == "Mechanistic (lQ)" ~ "Mechanistic~(lambda[Q])",
    Type == "Mechanistic (pP)" ~ "Mechanistic~(p[P])",
    Type == "Mechanistic (pG)" ~ "Mechanistic~(p[G])",
    TRUE ~ Type
  ))

Figure2_ticks <- Figure2_df %>% 
  # Get x-coordinate where R0 first exceeds/is less than one
  group_by(Type) %>% 
  mutate(sbr = 1440/theta) %>%  # standard biting rate
  arrange(sbr) %>% 
  summarise(
    first_R0_greater_1 = sbr[R0 > 1][1],
    # first_R0_less_1 = sbr[R0 < 1][1]
  ) %>% 
  # select(-sbr) %>% 
  unique()

Figure2_labels = c(expression("Standard"), expression("Exponential"), expression("Empirical"), expression("Phenomenological"), 
                   expression("Mechanistic " (lambda[Q])), expression("Mechanistic " (p[P])),expression("Mechanistic " (p[G]))
)

quick_fig_two <- function(keep_names) {
  Figure2_df %>% 
    filter(Type %in% keep_names) %>% 
    ggplot(aes(color = Type, lty = Type)) +
    geom_hline(aes(yintercept = 1), color = "grey", lwd = 2) +
    geom_line(aes(x = 1440 / theta, y = R0),
              lwd = 1) +
    geom_rug(
      data = Figure2_ticks%>% 
        filter(Type %in% keep_names),
      aes(x = first_R0_greater_1),
      sides = "b", size = 0.75, outside = TRUE,
      length = unit(0.3, "in"),
      show.legend = F
    ) +
    scale_x_continuous(
      name = TeX("Standard biting rate [Days$^{-1}$]"),
      expand = c(0,0),
      breaks = seq(0,2, by = 0.25)
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
    coord_cartesian(clip = "off") +
    theme_half_open(18) + 
    theme(
      axis.title.x = element_text(margin = margin(t = 10)),
      legend.key.width = unit(0.4, "in")
    )  +
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
}

R0theta_one <- quick_fig_two("Standard")
R0theta_two <- quick_fig_two(c("Standard", "Exponential"))
R0theta_three <- quick_fig_two(c("Standard", "Exponential", "Empirical"))
R0theta_four <- quick_fig_two(c("Standard", "Exponential", "Empirical",
                                "Phenomenological"))
R0theta_five <- quick_fig_two(unique(Figure2_df$Type))
R0theta_mech_only <- quick_fig_two(c("Mechanistic~(lambda[Q])", "Mechanistic~(p[P])","Mechanistic~(p[G])"))

ggsave("presentations/figures/R0theta_1.png", R0theta_one,
       width = 13.333, height = 7.5, units = "in")
ggsave("presentations/figures/R0theta_2.png", R0theta_two,
       width = 13.333, height = 7.5, units = "in")
ggsave("presentations/figures/R0theta_3.png", R0theta_three,
       width = 13.333, height = 7.5, units = "in")
ggsave("presentations/figures/R0theta_4.png", R0theta_four,
       width = 13.333, height = 7.5, units = "in")
ggsave("presentations/figures/R0theta_5.png", R0theta_five,
       width = 13.333, height = 7.5, units = "in")
ggsave("presentations/figures/R0theta_mech.png", R0theta_mech_only,
       width = 13.333, height = 7.5, units = "in")

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
  lG = "Ingesting~rate*','~lambda[G]"
) %>% pivot_longer(everything(), values_to = "nice_labels")

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
    breaks = waiver(),
    n.breaks = 4,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = TeX("Basic reproduction number \\, [$R_0$]"),
    expand = c(0,0)
  ) +
  scale_color_manual(
    name = "Parameter type:",
    values = c4a("met.egypt",4)[1:2]
    # values = c("#1B9E77", "#D95F02")
  ) +
  scale_linetype_discrete(
    name = "Mosquito type:",
    labels = c("Flighty", "Persistent")
  ) +
  theme_half_open(18) +
  guides(
    color = guide_none()
  ) +
  theme(
    strip.background = element_rect(color = "white", fill = "white"),
    legend.key.width = unit(0.4, "in")
  )

shift_legend(Figure3)

ggsave("presentations/figures/R0_vs_mech.png", shift_legend(Figure3),
       width = 13.333, height = 7.5, units = "in")
# 4. PRCCs of R0 against mechanistic parameters ----


# Load in final PRCC results
PRCC_data <- read_csv("data/PRCC_data.csv")

plot_data <- PRCC_data %>% 
  group_by(type, output) %>% 
  # Add in nice labels
  left_join(
    param_table %>% 
      rename(input = short_label)
  ) %>%
  left_join(rename(nice_mech_labels, input = name), by = "input")  %>% 
  mutate(
    output_label = case_when(
      output == "R0" ~ "Basic reproduction number",
      output == "GCD" ~ "Gonotrophic cycle duration",
      output == "N_offspring" ~ "Basic offspring number",
    ),
    type_label = case_when(
      type == "flighty" ~ "Flighty",
      type == "persistent" ~ "Persistent",
      type == "max" ~ "Maximum variation",
    ),
    dummy_min = -abs(PRCC[input == "dummy"]),
    dummy_max = abs(PRCC[input == "dummy"])
  ) %>% 
  filter(!is.na(output), input != "dummy") %>% 
  mutate(input_num = as.numeric(factor(input)))


plot_data$type = factor(plot_data$type, levels = c("flighty", "persistent", "max"))
plot_data$input = factor(plot_data$input, levels = c(
  "pQ", "pL", "pP", "pG", "sigma", "lQ", "lL", "lP", "lG", "dummy"
))
plot_data$output_label = factor(plot_data$output_label, levels = c(
  "Gonotrophic cycle duration", "Basic offspring number","Basic reproduction number"
))

plot_data$type_label = factor(plot_data$type_label, levels = c(
  "Flighty","Persistent","Maximum variation"
))


plot_data$Label = factor(plot_data$Label, levels = rev(c(
  c("Seeking success",
    "Landing success",  "Probing success", "Ingesting success"),
  "Persistence probability", 
  "Seeking rate", "Landing rate", "Probing rate", "Ingesting rate", "Dummy variable"
)))

plot_data$nice_labels <- factor(
  plot_data$nice_labels,
  levels = c(rev(nice_mech_labels$nice_labels), "Dummy~variable"))


PRCC_plots_max_only <- plot_data %>% 
  # Just keep maximum variation for now
  filter(type == "max") %>%
  arrange(input) %>% 
  filter(output %in% c("N_offspring","R0")) %>%
  group_by(type, output) %>% 
  mutate(star_flag = abs(PRCC) < abs(dummy_min),
         star_xpos = PRCC + sign(PRCC) * 0.025) %>% 
  # Plot
  ggplot() +
  geom_col(aes(y = nice_labels, x = PRCC, fill = output_label), position = "dodge", width = 0.75) +
  # Add light grey lines to divide up categories
  geom_hline(yintercept = seq(1.5, length(levels(plot_data$input)) - 0.5, by = 1), 
             color = "grey60", linetype = "dashed", linewidth = 0.125) +
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
    size = 4, vjust = -0.55
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
  theme_half_open(18) +
  theme(
    strip.background = element_rect(color = "white", fill = "white"),
    legend.key.width = unit(0.25, "in"),
    legend.key.height = unit(0.03125, "in"),
    legend.position = "top",
    legend.direction = "horizontal"
  )

PRCC_plots_max_only

ggsave("presentations/figures/PRCCs.png", PRCC_plots_max_only,
       width = 13.333, height = 7.5, units = "in")