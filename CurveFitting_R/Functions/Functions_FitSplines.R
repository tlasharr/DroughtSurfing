library(tidyverse)
library(sf)
library(terra)
library(raster)
library(zoo)
library(lubridate)
library(scales)


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

## temporal smoothing function
temporal_smooth <- function(x, window_size = 3) {
  if (sum(!is.na(x)) < 2) {
    return(rep(NA, length(x)))  # return NA if insufficient data
  }
  return(rollapply(x, 
                   width = window_size, 
                   FUN = function(y) {
                     if (sum(!is.na(y)) < 2) {
                       return(NA)
                     }
                     return(median(y, na.rm = TRUE))
                   }, 
                   fill = "extend",  # extend edges
                   align = "center"))
}

## apply the temporal smoothing function to each pixel in a time series
apply_temporal_smooth <- function(x, window_size = 3) {
  tryCatch(
    {
      return(temporal_smooth(x, window_size))
    },
    error = function(e) {
      message("An error occurred: ", e$message)
      return(rep(NA, length(x)))  # return NA if an error occurs
    }
  )
}


## function to fit logistic regression and find the midpoint (inflection)
fit_logistic_inflection <- function(segment_x, segment_y) {
  
  ## fit a logistic regression model
  logistic_model <- glm(segment_y ~ segment_x, family = quasibinomial)
  
  ## get coefficients
  coefficients <- coef(logistic_model)
  beta_0 <- coefficients[1]  # Intercept
  beta_1 <- coefficients[2]  # Slope
  
  ## create a sequence of x values
  x_vals <- seq(min(segment_x), max(segment_x), length.out = 100)
  
  ## calculate the first derivative for each x value
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
  
  ## get only the condition column we want and require at least 4 rows of data OR
  ## if there is a change of at least 0.1
  seg_list <- Filter(function(chunk) {
    if (all(chunk[, condition_column])) {
      if (nrow(chunk) >= 4) {
        return(TRUE)  # Keep segments with 3 or more rows
      } else {
      ## then check if the change between consecutive rows meets the threshold
        ndvi_change <- diff(chunk$y)
        return(any(abs(ndvi_change) >= 0.1))
      }
    }
    return(FALSE)
  }, seg_list)
  
  ##
  
  for(k in 1:length(seg_list)){
    
    current_segment <- seg_list[[k]]
    current_segment$rate_of_change <- c(NA, diff(current_segment$y) / diff(current_segment$jday))
    segment_x <- current_segment[, 1]
    segment_y <- current_segment[, 2]
    
    max_derivative_index <- which.max(abs(current_segment$rate_of_change))
    tmp_est <- current_segment[max_derivative_index,"jday"]
    irg_est <- c(irg_est, tmp_est)
    
  }
  
  return(irg_est)
}

## try catch the spline fitting
fit_spline_safe <- function(x, y, spar = 0.5) {
  
  ## check for missing or infinite values in the inputs
  if (any(is.na(x)) || any(is.na(y)) || any(is.infinite(x)) || any(is.infinite(y))) {
    warning("Missing or infinite values detected; skipping spline fitting.")
    return(NULL)
  }
  
  ## try to fit the spline
  spline_model <- try(smooth.spline(x, y, spar = spar, cv = TRUE), silent = TRUE)
  
  ## check if an error occurred
  if (inherits(spline_model, "try-error")) {
    warning("Error in smooth.spline fitting; returning NULL.")
    return(NULL)
  }
  
  ## return the fitted spline model
  return(spline_model)
}

## fit weighted spline
fit_weighted_spline <- function(x, y, spar = 0.5) {
  
  ## define weights (e.g., higher weights for regions with more change)
  weights <- 1 / (abs(diff(c(0, y))) + 0.01)  # adjust weights based on y changes
  spline_model <- smooth.spline(x, y, spar = spar, w = weights)
  
  spline_model$y <- pmax(spline_model$y, 0)
  
  return(spline_model)
  
}

## calculate green up periods
calc_greenup_metrics <- function(spl, data){
  
  data <- data[,c("date", "jday")]
  
  ## get the different segments of increasing or decreasing values
  segments <- data.frame(
    x = spl$x,
    y = spl$y,
    segment_id = c(NA, cumsum(c(1, diff(sign(diff(spl$y))) != 0)))  # Create segment ids based on changes in sign
  )
  
  
  ## get the total change for each segment
  segment_changes <- aggregate(y ~ segment_id, data = segments, function(x) {
    if (length(x) > 1) {
      return(x[length(x)] - x[1])
    } else {
      return(0)
    }
  })
  
  segment_changes$type <- ifelse(segment_changes$y < 0, "dec", "inc")
  
  ## identify significant segments
  threshold = 0.1  # threshold for detectable changes (10% change of the yearly NDVI values)
  significant_segments <- segment_changes$segment_id[abs(segment_changes$y) > threshold]

  if(length(significant_segments)==0){
    
    results_list <- NULL
    
    return(results_list)
    
    
  }else{
  
  ## get the points that are part of the significant segments
  significant_points <- segments[segments$segment_id %in% significant_segments,]
  significant_points <- merge(significant_points, segment_changes[,c("segment_id", "type")], by = "segment_id")

  ## create a dataframe
  result <- data.frame(x = spl$x, y = spl$y)
  result$inc <- FALSE
  result$dec <- FALSE
  
  if(all(significant_points$type=='inc')){
    result$dec <- FALSE
  }else{
    result[result$x %in% significant_points[significant_points$type=="dec",]$x,]$dec <- TRUE
  }
  
  if(all(significant_points$type=='dec')){
    result$inc <- FALSE
    
  }else{
    result[result$x %in% significant_points[significant_points$type=="inc",]$x,]$inc <- TRUE
    
  }
    

  ## get the dates and the julian date
  result$jday <- result$x
  
  ## bound them within 0 and 1
  result$y <- ifelse(result$y < 0, 0, result$y)
  result$y <- ifelse(result$y > 1, 1, result$y)
  
  inflection_points_increasing <- try(process_segments_for_logistic_inflection(result, "inc"), silent = TRUE)
  inflection_points_decreasing <- try(process_segments_for_logistic_inflection(result, "dec"), silent = TRUE)
  
  if(class(inflection_points_increasing)=='try-error'){
    inflection_points_increasing <- NA
  }
  if(class(inflection_points_decreasing)=='try-error'){
    inflection_points_decreasing <- NA
  }
  
  result <- result %>%
    mutate(inc_group = cumsum(!inc  | (lag(inc , default = FALSE) & !inc)))
  result$inc_group <- ifelse(result$inc==FALSE, NA, result$inc_group)
  
  inc_periods <- result %>%
    group_by(inc_group) %>%
    summarise(
      start_date = first(jday),
      end_date = last(jday),
      .groups = 'drop'  # This will drop the grouping so you don't need to ungroup later
    )
  
  ## remove NA groups
  inc_periods <- inc_periods[!is.na(inc_periods$inc_group),]
  
  ## now only get it for the inflection points we identified
  inc_periods <- inc_periods %>%
    rowwise() %>%
    filter(any(inflection_points_increasing >= start_date & inflection_points_increasing <= end_date))
  
  result <- result %>%
    mutate(dec_group = cumsum(!dec  | (lag(dec , default = FALSE) & !dec)))
  result$dec_group <- ifelse(result$dec==FALSE, NA, result$dec_group)
  
  dec_periods <- result %>%
    group_by(dec_group) %>%
    summarise(
      start_date = first(jday),
      end_date = last(jday),
    )
  
  dec_periods <- dec_periods[!is.na(dec_periods$dec_group),]
  
  ## now only get it for the inflection points we identified
  dec_periods <- dec_periods %>%
    rowwise() %>%
    filter(any(inflection_points_decreasing >= start_date & inflection_points_decreasing <= end_date))
  
  greenup_start <- if(nrow(inc_periods)==0){NA}else{inc_periods$start_date}
  greenup_end <- if(nrow(inc_periods)==0){NA}else{inc_periods$end_date}
  
  senesc_start <- if(nrow(dec_periods)==0){NA}else{dec_periods$start_date}
  senesc_end <- if(nrow(dec_periods)==0){NA}else{dec_periods$end_date}
  

  results_list <- list(inflection_points_increasing=inflection_points_increasing,
                       inflection_points_decreasing=inflection_points_decreasing,
                       greenup_start=greenup_start, greenup_end=greenup_end,
                       senesc_start=senesc_start, senesc_end=senesc_end)
  
  return(results_list)
  
  }
  
}
