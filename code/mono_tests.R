# Load required library
library(mgcv)
library(readr)
library(trend)
library(mblm)
library(ggplot2)
library(dplyr)

# Read in data
data = read_csv("./data/julia_outputs.csv")
# data = read_csv("./data/julia_outputs.csv")

# Generate example data
set.seed(123)
subsamps = sample(1:length(data$sigma), 1000)
x <- data$sigma[subsamps]#seq(0, 100, length.out = 100)  # Predictor variable
y <- log10(data$R0[subsamps]) #x^2 - x + rnorm(100, sd = 20)  # Quadratic relationship with noise
df <- data.frame(x = x, y = y)


# Mann-Kendall Test for Monotonicity
mann_kendall_test <- function(x, y) {
  # Ensure data is sorted by x
  order_idx <- order(x)
  x <- x[order_idx]
  y <- y[order_idx]
  
  # Perform the Mann-Kendall test
  test_result <- mk.test(y)
  
  # Return results
  list(
    tau = test_result$statistic,
    p_value = test_result$p.value,
    trend = ifelse(test_result$statistic > 0, "Increasing", "Decreasing")
  )
}

result <- mann_kendall_test(x, y)
cat("Test Case:\n")
cat("Tau:", result$tau, "\n")
cat("P-value:", result$p_value, "\n")
cat("Trend:", result$trend, "\n\n")


# If fails the Mann-Kendall test, see if it is roughly constant (which is still "monotonic")
constant_relationship_test <- function(x, y, alpha = 0.05) {
  # Fit a LOESS model to smooth the relationship
  loess_fit <- loess(y ~ x, span = 0.75)
  y_pred <- predict(loess_fit)
  
  # Step 1: Compute deviations from the mean (constant function)
  y_mean <- mean(y)
  deviations <- y_pred - y_mean
  
  # Step 2: Perform an F-test for variance explained by the trend
  total_variance <- sum((y - y_mean)^2)
  explained_variance <- sum(deviations^2)
  residual_variance <- total_variance - explained_variance
  
  # Compute F-statistic
  df1 <- 1  # Degrees of freedom for the constant model
  df2 <- length(y) - 2  # Degrees of freedom for LOESS model
  f_statistic <- (explained_variance / df1) / (residual_variance / df2)
  
  # Compute p-value
  p_value <- pf(f_statistic, df1, df2, lower.tail = FALSE)
  
  # Conclusion
  result <- if (p_value < alpha) "Not Constant" else "Roughly Constant"
  
  # Return results
  list(
    explained_variance = explained_variance,
    residual_variance = residual_variance,
    f_statistic = f_statistic,
    p_value = p_value,
    result = result
  )
}

# Test Case 1: Constant function with noise
result <- constant_relationship_test(x, y)
cat("Test Case:\n")
cat("Explained Variance:", result$explained_variance, "\n")
cat("Residual Variance:", result$residual_variance, "\n")
cat("F-Statistic:", result$f_statistic, "\n")
cat("P-value:", result$p_value, "\n")
cat("Result:", result$result, "\n\n")

# Function to evaluate monotonicity and constancy
evaluate_relationships <- function(input_data, output_data) {
  results <- list()
  
  for (output_var in colnames(output_data)) {
    for (input_var in colnames(input_data)) {
      x <- input_data[[input_var]]
      y <- output_data[[output_var]]
      
      # Step 1: Mann-Kendall test for monotonicity
      mk_result <- mann_kendall_test(x, y)
      
      # If monotonic, record the result
      if (mk_result$p_value < 0.05) {
        results[[paste(input_var, "vs", output_var)]] <- list(
          monotonic = TRUE,
          trend = mk_result$trend,
          p_value = mk_result$p_value,
          constant = NA,
          constant_p_value = NA
        )
      } else {
        # Step 2: Constant relationship test
        cr_result <- constant_relationship_test(x, y)
        
        # Record the results for non-monotonic relationships
        results[[paste(input_var, "vs", output_var)]] <- list(
          monotonic = FALSE,
          trend = NA,
          p_value = mk_result$p_value,
          constant = cr_result$result == "Roughly Constant",
          constant_p_value = cr_result$p_value
        )
      }
    }
  }
  
  return(results)
}
subsamps = sample(1:length(data$sigma), 1000)
input_data = select(data, "sigma":"pG")[subsamps,]
output_data = select(data, "GCD":"R0")[subsamps,]

# Evaluate relationships
results <- evaluate_relationships(input_data, output_data)

# Print results
for (key in names(results)) {
  cat(key, "\n")
  cat("Monotonic:", results[[key]]$monotonic, "\n")
  if (!is.na(results[[key]]$trend)) cat("Trend:", results[[key]]$trend, "\n")
  cat("Monotonicity p-value:", results[[key]]$p_value, "\n")
  if (!is.na(results[[key]]$constant)) {
    cat("Roughly Constant:", results[[key]]$constant, "\n")
    cat("Constancy p-value:", results[[key]]$constant_p_value, "\n")
  }
  cat("\n")
}

# If all relationships are monotonic, calculate PRCCs