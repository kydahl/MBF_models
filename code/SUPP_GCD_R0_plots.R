### Supplemental Stuff ####

## [Sanity check] Calculate theta_calc(A_matrix, v_alpha) and compare to theta

# For some chosen values of theta, plot pdfs of lognormal and fitted PH
library(PhaseTypeR)
library(gridExtra)
library(cowplot)
library(latex2exp)

test_theta_vals = c(min(samples$theta), median(samples$theta), max(samples$theta))

test_df = tibble(test_theta = test_theta_vals) %>%
  rowwise() %>%
  mutate(ph_model = list(phfit.point(ph = ph(3),
                                     x = filter(samples, theta == test_theta)$sample)$model),
         row_num = row_number()) %>%
  group_by(row_num)

x_maxs = c(2750, 75000, 200000)
scale_factors = c(1, 10, 10)
legend_list = c("none", "none", "legend")

plot_list = list()

for (i in 1:length(test_df)) {
  test_theta_val = test_df$test_theta[i]
  x_max = x_maxs[i]
  temp_plot = ggplot(data.frame(x_vals = seq(0, x_max, length.out = 1000)), aes(x = x_vals)) +
    stat_function(fun = dphase, args = list(ph = test_df$ph_model[[i]]),
                  aes(color = "PH"),
                  lwd = 1
    ) +
    stat_function(fun = dlnorm, args = list(mean = log(test_theta_val), sd = 1),
                  aes(color = "lognormal"),
                  lwd = 1
    ) +
    scale_x_continuous(name = "Days", breaks = seq(0, x_max, by = scale_factors[i]*1440), labels = seq(0, x_max / (scale_factors[i]*1440), by = 1)) +
    scale_y_continuous(name = "") +
    guides(
      y = "none"
    ) +
    scale_color_discrete(name = "Distribution") +
    ggtitle(TeX(paste0("$\\theta \\approx ", round(test_theta_val/1440, 2)," days"))) +
    theme_minimal(18) +
    theme(
      legend.position = "right",  # Legend positioned to the right
      legend.key.size = unit(1, "cm")  # Adjust legend key size
    )

  plotgrob = ggplotGrob(temp_plot)
  legend_index <- which(sapply(plotgrob$grobs, function(x) x$name) == "guide-box")
  legend <- plotgrob$grobs[[legend_index]]

  # panel_index <- which(plotgrob$layout$name == "panel")
  # plot_panel <- gtable::gtable_filter(plotgrob, "panel")

  plotgrob$grobs[[legend_index]] <- nullGrob()

  plot_list[[i]] = temp_plot + theme(legend.position = "none")
}
all_plots <- grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], legend, nrow = 1, ncol = 4, widths = c(3, 3, 3, 1))

