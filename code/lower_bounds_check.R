# Test for feasible lower bounds for the probability and rate parameters

# Load libraries
library(tidyverse)

# Import data from Julia
bounds_data = read_csv("data/lower_bounds_test.csv") |> 
  # translate durations back to rates
  mutate(R = 1/r, G = 1/g)

# Coarse tests
bounds_test = filter(bounds_data, N_offspring >= 1)

min(bounds_test$p)
min(bounds_test$R)
min(bounds_test$G)


# Plot regions
temp_plot <- bounds_data |> 
  filter(R == median(R)) |> 
  ggplot(aes(x = p, y = G, z = N_offspring)) + 
  geom_contour(breaks = c(1)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))

temp_plot

# Create slider plot
library(plotly)
library(akima)

df = bounds_data |> 
  mutate(x1 = p, x2 = r, x3 = g, y = N_offspring) |> 
  filter(x1 > 0, x2 > 0, x3 > 0, y > 0) |> 
  na.omit()

# create function that extracts contour data for given x3
make_contour <- function(x3_val, tol = 0.05) {
  df_sub <- df[abs(df$x3 - x3_val) < tol, ]
  interp_res <- with(df_sub, akima::interp(x1, x2, y, duplicate = "mean"))
  
  contour_df <- data.frame(
    x = rep(interp_res$x, each = length(interp_res$y)),
    y = rep(interp_res$y, times = length(interp_res$x)),
    z = as.vector(interp_res$z)
  )
  contour_df
}

# build plotly widget with slider
x3_vals <- pretty(range(df$x3), 20)
plots <- lapply(x3_vals, function(val) {
  dat <- make_contour(val)
  plot_ly(dat, x = ~x, y = ~y, z = ~z, type="contour", contours=list(showlabels=TRUE, start=1, end=1, size=1)) %>%
    layout(title = paste("x3 =", round(val,2)))
})

# combine into one slider-controlled widget
subplot(plots, nrows = 1, margin=0.05) # could use animation instead

