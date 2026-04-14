library(tidyverse)
library(sf)
library(terra)
library(raster)
library(zoo)


#### FUNCTIONS

## rescale function
scale_values <- function(x){(x-min(x))/(max(x)-min(x))}

## function to compute the first derivative
logistic_first_derivative <- function(x, beta_0, beta_1) {
  
  ## predicted probability
  y_hat <- 1 / (1 + exp(-(beta_0 + beta_1 * x)))
  
  ## first derivative
  dy_dx <- y_hat * (1 - y_hat) * beta_1
  
  return(dy_dx)
  
}


## function to smooth the temporal dimension of a SpatRaster
temporal_smooth <- function(x, window_size = 3) {
  
  if (sum(!is.na(x)) < 2) {
    return(rep(NA, length(x)))  # Return NA if insufficient data
  }
  return(rollapply(x, 
                   width = window_size, 
                   FUN = median, 
                   fill = "extend",  # Extend edges
                   align = "center"))
}

## function to fit logistic regression and find the midpoint (inflection)
fit_logistic_inflection <- function(segment_x, segment_y) {
  
  ## fit a logistic regression model
  logistic_model <- glm(segment_y ~ segment_x, family = quasibinomial)
  
  ## get coefficients
  coefficients <- coef(logistic_model)
  beta_0 <- coefficients[1]  # Intercept
  beta_1 <- coefficients[2]  # Slope
  
  # Create a sequence of x values
  x_vals <- seq(min(segment_x), max(segment_x), length.out = 100)
  
  # Calculate the first derivative for each x value
  first_derivative <- logistic_first_derivative(x_vals, beta_0, beta_1)
  
  max_derivative_index <- which.max((first_derivative))
  x_at_max_derivative <- x_vals[max_derivative_index]
  
  return(x_at_max_derivative)
}


## function to process segments and get the inflection points (2nd deriv)
process_segments_for_logistic_inflection <- function(data, condition_column) {
  
  ## empty object to hold first deriv
  irg_est <- c()
  
  ## split the data into groups based on the condition column
  data$group <- cumsum(c(TRUE, diff(data[,condition_column]) != 0))
  seg_list <- split(data, data$group)
  
  ## get only the condition column we want and require at least 4 rows of data
  seg_list <- Filter(function(chunk) {
    nrow(chunk) >= 4 && all(chunk[,condition_column])  ## at least 4 rows and all inc = TRUE
  }, seg_list)
  
  ##
  
  for(k in 1:length(seg_list)){
    
    current_segment <- seg_list[[k]]
    current_segment$rate_of_change <- c(NA, diff(current_segment$y) / diff(current_segment$x))
    segment_x <- current_segment[, 1]
    segment_y <- current_segment[, 2]
    
    max_derivative_index <- which.max(abs(current_segment$rate_of_change))
    tmp_est <- segment_x[max_derivative_index]
    irg_est <- c(irg_est, tmp_est)
    
  }
  
  return(irg_est)
}

## try catch the spline fitting
fit_spline_safe <- function(x, y, spar = 0.5) {
  # Check for missing or infinite values in the inputs
  if (any(is.na(x)) || any(is.na(y)) || any(is.infinite(x)) || any(is.infinite(y))) {
    warning("Missing or infinite values detected; skipping spline fitting.")
    return(NULL)
  }
  
  # Try to fit the spline
  spline_model <- try(smooth.spline(x, y, spar = spar, cv = TRUE), silent = TRUE)
  
  # Check if an error occurred
  if (inherits(spline_model, "try-error")) {
    warning("Error in smooth.spline fitting; returning NULL.")
    return(NULL)
  }
  
  # Return the fitted spline model
  return(spline_model)
}

## fit weighted spline
fit_weighted_spline <- function(x, y, spar = 0.5) {
  
  # Define weights (e.g., higher weights for regions with more change)
  weights <- 1 / (abs(diff(c(0, y))) + 0.01)  # Adjust weights based on y changes
  spline_model <- smooth.spline(x, y, spar = spar, w = weights)
  
  spline_model$y <- pmax(spline_model$y, 0)
  
  return(spline_model)
  
}
