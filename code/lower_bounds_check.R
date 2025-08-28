# Test for feasible lower bounds for the probability and rate parameters

# Load libraries
library(tidyverse)

# Import data from Julia
bounds_data = read_csv("data/lower_bounds_test.csv") |> 
  # translate durations back to rates
  mutate(R = 1/r, G = 1/g)

# Coarse tests
bounds_test = filter(bounds_data, N_offspring >= 1) |> 
  filter(g == 489)

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

library(dplyr)
library(plotly)
library(mgcv)   # smoother, avoids akima scale errors

df <- bounds_data %>%
  mutate(x1 = p, x2 = r, x3 = g, y = N_offspring) |> 
  filter(p < 0.3)

# choose values of x3 to show
x3_vals <- pretty(range(df$x3), 20)

# build a combined dataframe of all slices
all_grids <- do.call(rbind, lapply(x3_vals, function(val) {
  df_sub <- df %>% filter(abs(x3 - val) < 51)
  if(nrow(df_sub) < 10) return(NULL)
  
  # fit smoother on subset
  gam_fit <- gam(y ~ s(x1, x2), data = df_sub)
  
  grid <- expand.grid(
    x1 = seq(min(df_sub$x1), max(df_sub$x1), length=80),
    x2 = seq(min(df_sub$x2), max(df_sub$x2), length=80)
  )
  grid$yhat <- predict(gam_fit, newdata=grid)
  grid$x3_val <- val   # tag the slice
  grid
}))

# build interactive plotly contour plot with slider
p <- plot_ly(
  data = all_grids,
  x = ~x1, y = ~x2, z = ~yhat,
  frame = ~x3_val,    # <-- key: makes slider
  type = "contour",
  contours = list(
    showlabels = TRUE,
    start = 1, end = 1, size = 1   # plot contour where y=1
  )
) %>%
  layout(
    title = "Contour where y = 1",
    xaxis = list(title = "x1"),
    yaxis = list(title = "x2")
  )

p  # last line so it renders
library(htmlwidgets)

# save to standalone HTML file
htmlwidgets::saveWidget(p, "contour_slider.html", selfcontained = TRUE)

  
