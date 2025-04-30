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

## function to unscale the NDVI data
unscale_function <- function(x) {
  # x[1,] are the spline predictions, x[2] is min, and x[3] is max for each pixel
  if (any(is.na(x[1,]))) {
    return(rep(NA, length(x[1,])))
  } else {
    # Apply the unscale formula
    return(x[1] * (x[3] - x[2]) + x[2])
  }
}

## define a function to pad or trim elements to save rasters with the same number of layers
pad_or_trim <- function(x, length_out = 4) {
  x <- as.numeric(x)
  if (length(x) >= length_out) {
    return(x[1:length_out])
  } else {
    return(c(x, rep(NA, length_out - length(x))))
  }
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
                   fill = NA,  # extend edges
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

## function to process segments and get the inflection points
process_segments_for_logistic_inflection <- function(data,
                                                     condition_column,
                                                     ndvi_threshold = 0.05,
                                                     day_threshold = 14
                                                     ) {
  
  ## empty object 
  irg_est <- c()
  
  ## split the data into groups based on the condition column
  data$group <- cumsum(c(TRUE, diff(data[,condition_column]) != 0))
  seg_list <- split(data, data$group)
  
  ## get only the condition column we want and require at least 14 days of data OR
  ## if there is a change of at least 0.05
  seg_list <- Filter(function(chunk) {
    
    if (all(chunk[, condition_column])) {
      
      if (nrow(chunk) >= day_threshold) {
        
        return(TRUE)  # Keep segments with 14 or more days
      } else {
      
        ## then check if the change between consecutive rows meets the threshold
        ndvi_change <- diff(chunk$y)
        return(any(abs(ndvi_change) >= ndvi_threshold))
      }
    }
    return(FALSE)
  }, seg_list)
  
  ##
  
  for(k in 1:length(seg_list)){
    
    current_segment <- seg_list[[k]]
    current_segment$rate_of_change <- c(NA, diff(current_segment$y) / diff(current_segment$x))
    segment_x <- current_segment[, 1]
    segment_y <- current_segment[, 2]
    
    max_derivative_index <- which.max(abs(current_segment$rate_of_change))
    tmp_est <- current_segment[max_derivative_index,"x"]
    irg_est <- c(irg_est, tmp_est)
    
  }
  
  return(irg_est)
}

## fit weighted spline
fit_weighted_spline <- function(x, y, spar = 0.5) {
  
  ## check for missing or infinite values in the input
  if (any(is.na(x)) || any(is.na(y)) || any(is.infinite(x)) || any(is.infinite(y))) {
    warning("Missing or infinite values detected; skipping spline fitting.")
    return(NULL)
  }
  
  ## define weights based on changes in y, adjusting to avoid division by zero
  weights <- 1 / (abs(diff(c(0, y))) + 0.001)  # Adding a small number to avoid division by zero
  
  ## try to fit the spline with error handling
  spline_model <- try(smooth.spline(x, y, spar = spar, w = weights), silent = TRUE)
  
  ## check if an error occurred during the spline fitting
  if (inherits(spline_model, "try-error")) {
    warning("Error in smooth.spline fitting; returning NULL.")
    return(NULL)
  }
  
  ## ensure predictions don't go below 0 (if applicable)
  spline_model$y <- pmin(pmax(spline_model$y, 0), 1)  # Capping between 0 and 1
  
  return(spline_model)
}

## function to get splines for each raster cell
fit_spline_raster <- function(pixel_ts) {
  
  if (all(is.na(pixel_ts))) {
    
    ## return a vector of NAs with the length matching the number of output days * the number of summaries needed
    return(rep(NA, 365))  # Adjust based on what you decide to extract
  }
  
  ## assuming the time points (x-values) are 1:length(pixel_ts)
  time_points <- 1:length(pixel_ts)
  
  ## fit the spline
  spline_model <- fit_weighted_spline(time_points, pixel_ts, spar = 0.4)
  
  if (is.null(spline_model)) {
    return(rep(NA, 365))  # Adjust based on what you decide to extract
  }
  
  ## prediction for every day of the year
  predicted <- predict(spline_model, seq(from = 1, to = length(pixel_ts), length.out = 365))
  
  ## ensure no negative predictions, and adjust predictions if needed
  predicted$y <- pmin(pmax(predicted$y, 0), 1)  # Capping between 0 and 1
  
  return(predicted$y)  # Flatten the output to fit a vector
}


## calculate green up periods
calc_greenup_metrics <- function(pixel_values,
                                 ndvi_threshold = 0.05,
                                 day_threshold = 14){
  
  if (length(pixel_values) <= 1 || all(is.na(pixel_values))) {
    
    return(rep(NA, 32))  # Return NA if insufficient data
  }
  
  # Assuming the x values are just the indices of pixel_values
  x <- 1:length(pixel_values)
  y <- as.numeric(pixel_values)
  
  ## get the different segments of increasing or decreasing values
  segments <- data.frame(
    x = x,
    y = y,
    segment_id = c(NA, cumsum(c(1, diff(sign(diff(y))) != 0)))  # Create segment ids based on changes in sign
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
  sig_threshold = 0.05  # threshold for detectable changes (5% change of the yearly NDVI values)
  significant_segments <- segment_changes$segment_id[abs(segment_changes$y) > sig_threshold]
  
  ## get the points that are part of the significant segments
  significant_points <- segments[segments$segment_id %in% significant_segments,]
  significant_points <- merge(significant_points, segment_changes[,c("segment_id", "type")], by = "segment_id")
  
  ## create a dataframe
  result <- data.frame(x = x, y = y)
  result$inc <- FALSE
  result$dec <- FALSE
  
  if(nrow(significant_points[significant_points$type=="inc",])==0){
    
    result$inc <- FALSE
    
  }else{
    result[result$x %in% significant_points[significant_points$type=="inc",]$x,]$inc <- TRUE
  }
  
  if(nrow(significant_points[significant_points$type=="dec",])==0){
    
    result$dec <- FALSE
    
    
  }else{
    result[result$x %in% significant_points[significant_points$type=="dec",]$x,]$dec <- TRUE
  }
  
  ## get the dates and the julian date
  result$jday <- result$x
  
  ## bound them within 0 and 1
  result$y <- ifelse(result$y < 0, 0, result$y)
  result$y <- ifelse(result$y > 1, 1, result$y)
  ndvi_thresh <- ndvi_threshold
  day_thresh <- day_threshold
  
  inflection_points_increasing <- try(process_segments_for_logistic_inflection(result, "inc", ndvi_threshold = ndvi_thresh, day_threshold = day_thresh), silent = TRUE)
  inflection_points_decreasing <- try(process_segments_for_logistic_inflection(result, "dec", ndvi_threshold = ndvi_thresh, day_threshold = day_thresh), silent = TRUE)
  
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
      start_date = first(x),
      end_date = last(x),
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
      start_date = first(x),
      end_date = last(x),
    )
  
  dec_periods <- dec_periods[!is.na(dec_periods$dec_group),]
  
  ## now only get it for the inflection points we identified
  dec_periods <- dec_periods %>%
    rowwise() %>%
    filter(any(inflection_points_decreasing >= start_date & inflection_points_decreasing <= end_date))
  
  greenup_start <- if(nrow(inc_periods)==0){NA}else{inc_periods$start_date}
  greenup_end <- if(nrow(inc_periods)==0){NA}else{inc_periods$end_date}
  greenup_change <- if(nrow(inc_periods)==0){NA}else{result[result$x %in% inc_periods$end_date, "y"] - result[result$x %in% inc_periods$start_date, "y"]}
    
  senesc_start <- if(nrow(dec_periods)==0){NA}else{dec_periods$start_date}
  senesc_end <- if(nrow(dec_periods)==0){NA}else{dec_periods$end_date}
  senesc_change <- if(nrow(dec_periods)==0){NA}else{result[result$x %in% dec_periods$start_date, "y"] - result[result$x %in% dec_periods$end_date, "y"]}

  results_list <- list(greenup_peak=inflection_points_increasing,
                       senesc_peak=inflection_points_decreasing,
                       greenup_start=greenup_start,
                       greenup_end=greenup_end,
                       senesc_start=senesc_start,
                       senesc_end=senesc_end,
                       greenup_change = greenup_change,
                       senesc_change = senesc_change)

  
  ## apply to each element of your list and combine into one vector
  results_vector <- unlist(lapply(results_list, pad_or_trim, length_out = 4))
  results_vector <- as.numeric(results_vector)
  return(results_vector)
  
}




