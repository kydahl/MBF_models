# A re-analysis of the data from:
#     Christofferson, Rebecca C., Helen J. Wearing, Erik A. Turner, 
#     Christine S. Walsh, Henrik Salje, Cécile Tran-Kiem, and Simon Cauchemez. 
#     2022. “How Do i Bite Thee? Let Me Count the Ways: Exploring the 
#     Implications of Individual Biting Habits of Aedes Aegypti for Dengue 
#     Transmission.” PLoS Neglected Tropical Diseases 16 (10): e0010818. 
#     https://doi.org/10.1371/journal.pntd.0010818.

# The goal is to look at correlation between the time between bites and the
# number of bites taken in a single cycle

library(tidyverse)
# library(readr)
library(readxl)
library(cowplot)
library(cols4all)

# Read in data

# Biting data
bites_data = read_excel("data/raw/Christofferson_2022.xlsx", sheet = 1) %>% 
  pivot_longer(cols = "1":"25", names_to = "day", values_to = "bite_bool") %>% 
  group_by(MID, ID, Rep) %>% 
  mutate(RealID = cur_group_id()) %>% 
  ungroup() %>% 
  select(-c(MID, ID, Rep))

bites_test = bites_data %>% 
  group_by(RealID) %>% 
  summarise(
    TimeToFirstBite = min(day[bite_bool == 1]),
    TotalBites = sum(bite_bool),
    # TimeBetweenFirstAndLastBite
  )


# Summary data
summary_data = read_excel("data/raw/Christofferson_2022.xlsx", sheet = 2) %>% 
  select(-c(Rep, `Mosquito ID`, ID)) %>% 
  rename(
    Temperature = Temp,
    TimeToFirstBite = TTFB,
    TotalBites = Sum,
    TimeBetweenFirstAndLastBite = TB1Last,
    TimeBetweenFirstAndSecondBite = TB1Second
  ) %>% 
  mutate(RealID = row_number())

# Three temperature treatments: 24, 28, 32 C
# 16 individual females per temperature group (3 removed)
# Two replicates at each temperature were performed - then collapsed (i.e. can drop the "replicate" label)
# Post-emergence, mosquitoes were offered a 20-minute blood meal daily for 23 days
# Blood feeding was recorded at each offering

temperature_summary = summary_data %>% 
  group_by(Temperature) %>% 
  summarise(
    across(TimeToFirstBite:TimeBetweenFirstAndSecondBite, ~ mean(.x, na.rm = T)),
    # Num_Time_Cor = cor(TimeToFirstBite, TotalBites, use = "everything")
    )


# Correlation between number of bites per 25 days and time to first bite

Cor_df = summary_data %>% 
  group_by(Temperature) %>% 
  summarise(
    cor = cor(TotalBites, TimeToFirstBite, method = "spearman", use = "complete.obs"),
    p_value = cor.test(TotalBites, TimeToFirstBite, method = "spearman", use = "complete.obs")$p.value,
    .groups = "drop"
    )

Cor_plot = summary_data %>% 
  ggplot(aes(x = TimeToFirstBite, y = TotalBites, color = factor(Temperature))) +
  geom_point(
    size = 5,
    alpha  = 0.25
  ) +
  geom_smooth(
    data = summary_data,
    se = F,
    lwd = 2
  ) +
  geom_text(
    data = Cor_df, 
    aes(
      x = 2, y = 11,  # Adjust these values if necessary per Temperature group
      label = paste("Spearman's correlation =", round(cor, 3), "\np-value =", round(p_value, 4))
    ),
    color = "black", hjust = 0
  ) +
  
  facet_wrap(
    ~ Temperature, 
    nrow = 1,
    labeller = as_labeller(function(str) paste0(str," °C"))
  ) +
  scale_color_manual(
    values = c4a("met.moreau", 3),
    labels = c("24°C", "28°C", "32°C")
  ) +
  scale_shape(
    "Temperature"
  ) + 
  scale_x_continuous(
    "Time to first bite [days]",
    expand = c(0.01,0),
    # limits = c(2,20), 
    breaks = seq(0, 20, 2)
    # trans = 'log10'
  ) +
  scale_y_continuous(
    "Total number of bites over 25 days",
    limits = c(0, 12),
    expand = c(0,0),
    # trans = 'log10'
  ) +
  theme_half_open(18) + 
  guides(
    fill = "none",
    color = "none",
    linetype = guide_legend(
      key_glyph = "path",
      keywidth = 5,
      override.aes = list(size = 1, shape = NA, fill = NA)
    )
  ) +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.position = "top",
    legend.direction = "horizontal"
  )

Cor_plot


ggsave("figures/presentations/Christofferson_bite_time_cor.png", Cor_plot,
       width = 16, height = 7, units = "in")


# Total bites density
TotalBiteHist_plot = summary_data %>% 
  ggplot(aes(x = TotalBites, y = ..density.., color = as.factor(Temperature), fill = as.factor(Temperature))) +
  # geom_histogram(
  #   aes(y = ..density..),
  #   position = "dodge"
  # ) +
  geom_density(
    alpha = 0,
    lwd = 2
  ) +
  scale_fill_manual(
    "Temperature",
    values = c4a("met.moreau", 3)
  ) +  
  scale_color_manual(
    "Temperature",
    values = c4a("met.moreau", 3)
  ) +
  scale_shape(
    "Temperature"
  ) + 
  scale_y_continuous(
    "",
    expand = c(0,0)
  ) +
  scale_x_continuous(
    breaks = seq(0, 13),
    expand = c(0,0)
  ) +
  theme_half_open(12) 

# Time to first bite density

TimeToFirstBite_Exp_estimates = summary_data %>% 
  select(Temperature, TimeToFirstBite) %>% 
  group_by(Temperature) %>% 
  summarise(lambda = 1 / mean(TimeToFirstBite, na.rm = T)) %>%
  group_by(Temperature) %>% 
  summarise(
    TimeToFirstBite = seq(0, max(summary_data$TimeToFirstBite, na.rm = T), length.out = 100),
    Density = lambda * exp(-lambda * TimeToFirstBite),
    lambda = lambda
  ) %>%
  mutate(DistributionType = "Fitted Exponential")

# Adjusted Chi-square goodness of fit test function for exponential distribution
chi_square_test_adjusted <- function(empirical_values, exp_lambda, bins = 10) {
  # Bin the empirical data (discretize the continuous data into bins)
  breaks <- seq(min(empirical_values, na.rm = TRUE), max(empirical_values, na.rm = TRUE), length.out = bins + 1)
  observed_counts <- hist(empirical_values, breaks = breaks, plot = FALSE)$counts
  
  # Calculate the expected counts under the exponential distribution
  expected_counts <- diff(pgamma(breaks, shape = 1, rate = exp_lambda)) * length(empirical_values)
  
  # Ensure no NA values in expected_counts (replace NAs with 0 if needed)
  expected_counts[is.na(expected_counts)] <- 0
  
  # If any expected count is below 5, merge bins until all expected counts are >= 5
  while (any(expected_counts < 5)) {
    small_bin <- which.min(expected_counts)  # Find the smallest expected count
    if (small_bin > 1) {
      expected_counts[small_bin - 1] <- expected_counts[small_bin - 1] + expected_counts[small_bin]
      observed_counts[small_bin - 1] <- observed_counts[small_bin - 1] + observed_counts[small_bin]
    } else {
      expected_counts[small_bin + 1] <- expected_counts[small_bin + 1] + expected_counts[small_bin]
      observed_counts[small_bin + 1] <- observed_counts[small_bin + 1] + observed_counts[small_bin]
    }
    # Remove the merged bin from expected_counts and observed_counts
    expected_counts <- expected_counts[-small_bin]
    observed_counts <- observed_counts[-small_bin]
    breaks <- breaks[-small_bin]  # Adjust the bin breaks
    
    # Ensure no NA values after bin merging
    expected_counts[is.na(expected_counts)] <- 0
  }
  
  # Perform Chi-square test (observed vs. expected)
  chi_square_result <- chisq.test(observed_counts, p = expected_counts / sum(expected_counts), rescale.p = TRUE)
  
  return(list(chi_square_statistic = chi_square_result$statistic,
              p_value = chi_square_result$p.value))
}

# Apply to each temperature group
chi_square_results <- summary_data %>%
  group_by(Temperature) %>%
  do({
    # Get the corresponding lambda for this Temperature from TimeToFirstBite_Exp_estimates
    lambda <- unique(filter(TimeToFirstBite_Exp_estimates, Temperature == unique(.$Temperature))$lambda)
    
    # Perform the Chi-square test for this temperature group
    result <- chi_square_test_adjusted(.$TimeToFirstBite, lambda)
    
    # Return the Chi-square statistic and p-value
    data.frame(Temperature = unique(.$Temperature), 
               chi_square_statistic = result$chi_square_statistic,
               p_value = result$p_value)
  })


# Under usual assumptions, this should be distributed exponentially
TimeToFirstBiteHist_plot = summary_data %>%
  mutate(DistributionType = "Empirical Density") %>% 
  ggplot() +
  geom_density(
    aes(x = TimeToFirstBite, y = ..density.., color = factor(Temperature), linetype = DistributionType),
    lwd = 2,
    show.legend = T,
    key_glyph = "path"
  ) +
  geom_line(
    data = TimeToFirstBite_Exp_estimates,
    aes(x = TimeToFirstBite, y = Density, color = factor(Temperature), linetype = DistributionType),
    lwd = 2,
    key_glyph = "path"
  ) +
  facet_wrap(
    ~ Temperature, 
    nrow = 1,
    labeller = as_labeller(function(str) paste0(str," °C")),
    scales = "free_y"
    ) +
  scale_linetype_manual(
    "",
    values = c("Empirical Density" = "solid",
               "Fitted Exponential" = "dashed")
  ) +
  scale_color_manual(
    "Temperature",
    values = c4a("met.moreau", 3)
  ) +
  scale_y_continuous(
    "",
    expand = c(0,0)
  ) +
  scale_x_continuous(
    "Time to first bite [days]",
    breaks = seq(0, 20, 2),
    expand = c(0.01,0)
  ) +
  theme_half_open(18) + 
  guides(
    fill = "none",
    color = "none",
    linetype = guide_legend(
      key_glyph = "path",
      keywidth = 5,
      override.aes = list(size = 1, shape = NA, fill = NA)
    )
  ) +
  theme(
    strip.background = element_rect(fill = "white"),
    axis.text.y = element_blank(),
    legend.position = "top",
    legend.direction = "horizontal"
    )

ggsave("figures/presentations/Christofferson_not_exp.png", TimeToFirstBiteHist_plot,
       width = 16, height = 7, units = "in")