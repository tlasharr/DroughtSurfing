#### OLD FUNCTIONS FROM PREVIOUS WORKFLOW,
#### CAN DELETE EVENTUALLY BUT WANT TO SAVE THEM FOR NOW


## function to compute the first derivative
logistic_first_derivative <- function(x, beta_0, beta_1) {
  
  ## predicted probability
  y_hat <- 1 / (1 + exp(-(beta_0 + beta_1 * x)))
  
  ## first derivative
  dy_dx <- y_hat * (1 - y_hat) * beta_1
  
  return(dy_dx)
  
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


## calculate green up periods
calc_greenup_metrics <- function(pixel_values){
  
  if (length(pixel_values) <= 1 || all(is.na(pixel_values))) {
    
    return(rep(NA, 6))  # Return NA if insufficient data
  }
  
  # Assuming the x values are just the indices of pixel_values
  x <- 1:length(pixel_values)
  y <- pixel_values
  
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
  threshold = 0.05  # threshold for detectable changes (10% change of the yearly NDVI values)
  significant_segments <- segment_changes$segment_id[abs(segment_changes$y) > threshold]
  
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
  
  senesc_start <- if(nrow(dec_periods)==0){NA}else{dec_periods$start_date}
  senesc_end <- if(nrow(dec_periods)==0){NA}else{dec_periods$end_date}
  
  results_list <- list(inflection_points_increasing=inflection_points_increasing,
                       inflection_points_decreasing=inflection_points_decreasing,
                       greenup_start=greenup_start, greenup_end=greenup_end,
                       senesc_start=senesc_start, senesc_end=senesc_end)
  
  
  ## save the names of the metrics we want layers for
  metric_names <- c("inflection_points_increasing", "inflection_points_decreasing", 
                    "greenup_start", "greenup_end", "senesc_start", "senesc_end")
  
  out <- NULL
  
  for (metric_name in metric_names) {
    
    ## check if there are more than green up or senescence period. If there are 
    ## more than 1, we are going to store them all together in the layer and will have to seperate them out
    ## at a later time. So for example, if you had green up on julian dates 150 and 322, it would be
    ## stored as 150322. Additionally, this pads the values so they all have 0s in front, but 
    ## because the values are stored as a number, the leading 0s for the very first one will not
    ## keep the 0s. So when we read them back out we will have to figure out if it's divisible by 3 and 
    ## add appropriate 0s there.
    
    if(length(unlist(results_list[metric_name])) > 1){
      
      tmpvals <- unlist(results_list[metric_name])
      tmpvals <- sprintf("%03d", tmpvals)
      tmpvals <- as.numeric(paste0(tmpvals, collapse = ""))
      
      
    }else{
      
      tmpvals <- as.numeric(unlist(results_list[metric_name]))
      
    }
    
    out <- c(out, tmpvals)
    out <-as.numeric(out)
  }
  
  return(out)
  
}
