library(tidyverse)
library(sf)
library(terra)
library(raster)
library(zoo)
library(lubridate)
library(scales)


#### functions ####

## rescale function -- accepts explicit min/max so buffer periods
## don't distort the scaling of the focal year
scale_values <- function(x, x_min = NULL, x_max = NULL) {
  if (is.null(x_min)) x_min <- min(x, na.rm = TRUE)
  if (is.null(x_max)) x_max <- max(x, na.rm = TRUE)
  if (is.infinite(x_min) || is.infinite(x_max) || x_max == x_min) return(rep(NA, length(x)))
  (x - x_min) / (x_max - x_min)
}

## unscale function
unscale_function <- function(x, x_min, x_max) {
  if (any(is.na(x))) return(rep(NA, length(x)))
  x * (x_max - x_min) + x_min
}

## pad or trim elements to a fixed length for saving rasters with the same number of layers
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
    return(rep(NA, length(x)))
  }
  return(rollapply(x,
                   width = window_size,
                   FUN = function(y) {
                     if (sum(!is.na(y)) < 2) return(NA)
                     return(median(y, na.rm = TRUE))
                   },
                   fill = NA,
                   align = "center"))
}

## apply temporal smoothing to each pixel
apply_temporal_smooth <- function(x, window_size = 3) {
  tryCatch(
    { return(temporal_smooth(x, window_size)) },
    error = function(e) {
      message("an error occurred: ", e$message)
      return(rep(NA, length(x)))
    }
  )
}

## remove short runs of the floor (fill) value that sit between real
## observations well above the floor; long runs (winter baseline) are kept
drop_floor_dips <- function(y, max_run = 2, frac = 0.25, tol = 1e-6) {
  ok <- which(!is.na(y))
  if (length(ok) < 5) return(y)
  v   <- y[ok]
  fl  <- min(v)
  amp <- max(v) - fl
  is_floor <- abs(v - fl) < tol
  r      <- rle(is_floor)
  ends   <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  for (k in which(r$values & r$lengths <= max_run)) {
    s <- starts[k]; e <- ends[k]
    if (s == 1 || e == length(v)) next
    if (v[s - 1] > fl + frac * amp && v[e + 1] > fl + frac * amp) {
      y[ok[s:e]] <- NA
    }
  }
  y
}

## process segments and get inflection points
process_segments_for_logistic_inflection <- function(data,
                                                     condition_column,
                                                     ndvi_threshold = 0.05,
                                                     day_threshold  = 14) {
  irg_est <- c()
  
  data$group <- cumsum(c(TRUE, diff(data[, condition_column]) != 0))
  seg_list   <- split(data, data$group)
  
  seg_list <- Filter(function(chunk) {
    if (all(chunk[, condition_column])) {
      if (nrow(chunk) >= day_threshold) return(TRUE)
      ndvi_change <- diff(chunk$y)
      return(any(abs(ndvi_change) >= ndvi_threshold))
    }
    return(FALSE)
  }, seg_list)
  
  for (k in seq_along(seg_list)) {
    current_segment <- seg_list[[k]]
    current_segment$rate_of_change <- c(NA, diff(current_segment$y) /
                                          diff(current_segment$x))
    max_derivative_index <- which.max(abs(current_segment$rate_of_change))
    tmp_est  <- current_segment[max_derivative_index, "x"]
    irg_est  <- c(irg_est, tmp_est)
  }
  
  return(irg_est)
}

## fit weighted spline
fit_weighted_spline <- function(x, y, spar = 0.5) {
  if (any(is.na(x)) || any(is.na(y)) ||
      any(is.infinite(x)) || any(is.infinite(y))) {
    warning("missing or infinite values detected; skipping spline fitting.")
    return(NULL)
  }

  spline_model <- try(smooth.spline(x, y, spar = spar), silent = TRUE)
  
  if (inherits(spline_model, "try-error")) {
    warning("error in smooth.spline fitting; returning NULL.")
    return(NULL)
  }
  spline_model$y <- pmin(pmax(spline_model$y, 0), 1)
  return(spline_model)
}

## fit spline and return 365-day predictions for a single pixel
fit_spline_raster <- function(pixel_ts) {
  
  fit_weighted_spline <- function(x, y, spar = 0.5) {
    if (any(is.na(x)) || any(is.na(y)) ||
        any(is.infinite(x)) || any(is.infinite(y))) return(NULL)

    spline_model <- try(smooth.spline(x, y, spar = spar), silent = TRUE)
    
    if (inherits(spline_model, "try-error")) return(NULL)
    spline_model$y <- pmin(pmax(spline_model$y, 0), 1)
    return(spline_model)
  }
  
  if (all(is.na(pixel_ts))) return(rep(NA, 365))
  
  time_points  <- 1:length(pixel_ts)
  spline_model <- fit_weighted_spline(time_points, as.numeric(pixel_ts), spar = 0.4)
  
  if (is.null(spline_model)) return(rep(NA, 365))
  
  predicted    <- predict(spline_model,
                          seq(from = 1, to = length(pixel_ts), length.out = 365))
  predicted$y  <- pmin(pmax(predicted$y, 0), 1)
  
  return(predicted$y)
}

## compute rss from spline fit
rss_spline_func <- function(pixel_ts) {
  
  if (all(is.na(pixel_ts))) return(NA)
  
  fit_weighted_spline <- function(x, y, spar = 0.5) {
    if (any(is.na(x)) || any(is.na(y)) ||
        any(is.infinite(x)) || any(is.infinite(y))) return(NULL)

    spline_model <- try(smooth.spline(x, y, spar = spar), silent = TRUE)
    
    if (inherits(spline_model, "try-error")) return(NULL)
    spline_model$y <- pmin(pmax(spline_model$y, 0), 1)
    return(spline_model)
  }
  
  time_points  <- 1:length(pixel_ts)
  spline_model <- fit_weighted_spline(time_points, pixel_ts, spar = 0.4)
  
  if (is.null(spline_model)) return(NA)
  
  rss_spline <- try(sum(resid(spline_model)^2))
  return(rss_spline)
}


#### calc_greenup_metrics ####

## output structure (36 elements):
##   greenup_peak    [1:4]   -- doy of peak rate of increase for each green-up event
##   senesc_peak     [5:8]   -- doy of peak rate of decrease for each senescence event
##   greenup_start   [9:12]  -- doy of start of increasing segment
##   greenup_end     [13:16] -- doy of end of increasing segment (where increase stops)
##   green_period_end[17:20] -- doy of first senescence start AFTER each green-up peak
##                              (shared across green-up events with no senescence between them)
##   senesc_start    [21:24] -- doy of start of decreasing segment
##   senesc_end      [25:28] -- doy of end of decreasing segment
##   greenup_change  [29:32] -- total ndvi change during green-up segment
##   senesc_change   [33:36] -- total ndvi change during senescence segment

calc_greenup_metrics <- function(pixel_values,
                                 ndvi_threshold = 0.05,
                                 day_threshold  = 14) {
  
  require(tidyverse)
  require(terra)
  
  ## inner copies of helper functions needed for parallel workers
  process_segments_for_logistic_inflection <- function(data,
                                                       condition_column,
                                                       ndvi_threshold = 0.05,
                                                       day_threshold  = 14) {
    irg_est  <- c()
    data$group <- cumsum(c(TRUE, diff(data[, condition_column]) != 0))
    seg_list   <- split(data, data$group)
    seg_list   <- Filter(function(chunk) {
      if (all(chunk[, condition_column])) {
        if (nrow(chunk) >= day_threshold) return(TRUE)
        ndvi_change <- diff(chunk$y)
        return(any(abs(ndvi_change) >= ndvi_threshold))
      }
      return(FALSE)
    }, seg_list)
    for (k in seq_along(seg_list)) {
      current_segment <- seg_list[[k]]
      current_segment$rate_of_change <- c(NA, diff(current_segment$y) /
                                            diff(current_segment$x))
      max_derivative_index <- which.max(abs(current_segment$rate_of_change))
      tmp_est  <- current_segment[max_derivative_index, "x"]
      irg_est  <- c(irg_est, tmp_est)
    }
    return(irg_est)
  }
  
  pad_or_trim <- function(x, length_out = 4) {
    x <- as.numeric(x)
    if (length(x) >= length_out) return(x[1:length_out])
    c(x, rep(NA, length_out - length(x)))
  }
  
  if (length(pixel_values) <= 1 || all(is.na(pixel_values))) {
    return(rep(NA, 36))
  }
  
  x <- 1:length(pixel_values)
  y <- as.numeric(pixel_values)
  
  ## identify increasing and decreasing segments
  segments <- data.frame(
    x          = x,
    y          = y,
    segment_id = c(NA, cumsum(c(1, diff(sign(diff(y))) != 0)))
  )
  
  delta_y <- tapply(segments$y, segments$segment_id, function(vals) {
    if (length(vals) > 1) vals[length(vals)] - vals[1] else 0
  })
  segment_changes       <- data.frame(segment_id = as.numeric(names(delta_y)), y = delta_y)
  segment_changes$type  <- ifelse(segment_changes$y < 0, "dec", "inc")
  
  sig_threshold         <- 0.05
  significant_segments  <- segment_changes$segment_id[abs(segment_changes$y) > sig_threshold]
  significant_points    <- segments[segments$segment_id %in% significant_segments, ]
  significant_points    <- merge(significant_points,
                                 segment_changes[, c("segment_id", "type")],
                                 by = "segment_id")
  
  result      <- data.frame(x = x, y = y)
  result$inc  <- FALSE
  result$dec  <- FALSE
  
  result$inc[result$x %in% significant_points$x[significant_points$type == "inc"]] <- TRUE
  result$dec[result$x %in% significant_points$x[significant_points$type == "dec"]] <- TRUE
  result$jday <- result$x
  result$y    <- pmax(0, pmin(1, result$y))
  
  ndvi_thresh <- ndvi_threshold
  day_thresh  <- day_threshold
  
  inflection_points_increasing <- try(
    process_segments_for_logistic_inflection(result, "inc",
                                             ndvi_threshold = ndvi_thresh,
                                             day_threshold  = day_thresh),
    silent = TRUE
  )
  inflection_points_decreasing <- try(
    process_segments_for_logistic_inflection(result, "dec",
                                             ndvi_threshold = ndvi_thresh,
                                             day_threshold  = day_thresh),
    silent = TRUE
  )
  
  if (inherits(inflection_points_increasing, "try-error")) inflection_points_increasing <- NA
  if (inherits(inflection_points_decreasing, "try-error")) inflection_points_decreasing <- NA
  
  ## green-up periods
  result <- result %>%
    mutate(inc_group = cumsum(!inc | (lag(inc, default = FALSE) & !inc)))
  result$inc_group <- ifelse(result$inc == FALSE, NA, result$inc_group)
  
  inc_periods <- result %>%
    group_by(inc_group) %>%
    summarise(start_date = first(x), end_date = last(x), .groups = "drop") %>%
    filter(!is.na(inc_group)) %>%
    filter(apply(across(c(start_date, end_date)), 1, function(row)
      any(inflection_points_increasing >= row[1] &
            inflection_points_increasing <= row[2])))
  
  ## senescence periods
  result <- result %>%
    mutate(dec_group = cumsum(!dec | (lag(dec, default = FALSE) & !dec)))
  result$dec_group <- ifelse(result$dec == FALSE, NA, result$dec_group)
  
  dec_periods <- result %>%
    group_by(dec_group) %>%
    summarise(start_date = first(x), end_date = last(x), .groups = "drop") %>%
    filter(!is.na(dec_group)) %>%
    filter(apply(across(c(start_date, end_date)), 1, function(row)
      any(inflection_points_decreasing >= row[1] &
            inflection_points_decreasing <= row[2])))
  
  ## basic green-up metrics
  greenup_start  <- if (nrow(inc_periods) == 0) NA else inc_periods$start_date
  greenup_end    <- if (nrow(inc_periods) == 0) NA else inc_periods$end_date
  greenup_change <- if (nrow(inc_periods) == 0) NA else {
    result$y[match(inc_periods$end_date,   result$x)] -
      result$y[match(inc_periods$start_date, result$x)]
  }
  
  senesc_start  <- if (nrow(dec_periods) == 0) NA else dec_periods$start_date
  senesc_end    <- if (nrow(dec_periods) == 0) NA else dec_periods$end_date
  senesc_change <- if (nrow(dec_periods) == 0) NA else {
    result$y[match(dec_periods$start_date, result$x)] -
      result$y[match(dec_periods$end_date,   result$x)]
  }
  
  ## green_period_end: for each green-up peak, find the start of the first
  ## senescence event that occurs after that peak.
  ## if two green-up peaks share the same next senescence (no brown-down
  ## between them), both get that senescence start as their green_period_end.
  ## if there is a senescence event between two peaks, the first peak uses
  ## that intermediate senescence start and the second peak gets the next one.
  ## if no senescence event follows a green-up peak, green_period_end is NA.
  green_period_end <- if (length(inflection_points_increasing) == 0 ||
                          all(is.na(inflection_points_increasing)) ||
                          nrow(dec_periods) == 0) {
    rep(NA, length(inflection_points_increasing))
  } else {
    map_dbl(inflection_points_increasing, function(gu_peak) {
      if (is.na(gu_peak)) return(NA_real_)
      ## find senescence events that start after this green-up peak
      candidate_senesc <- senesc_start[senesc_start > gu_peak]
      if (length(candidate_senesc) == 0 || all(is.na(candidate_senesc))) {
        return(NA_real_)
      }
      ## return the earliest one
      min(candidate_senesc, na.rm = TRUE)
    })
  }
  
  ## pack everything into the output vector
  results_list <- list(
    greenup_peak     = inflection_points_increasing,
    senesc_peak      = inflection_points_decreasing,
    greenup_start    = greenup_start,
    greenup_end      = greenup_end,
    green_period_end = green_period_end,
    senesc_start     = senesc_start,
    senesc_end       = senesc_end,
    greenup_change   = greenup_change,
    senesc_change    = senesc_change
  )
  
  results_vector <- unlist(lapply(results_list, pad_or_trim, length_out = 4))
  results_vector <- as.numeric(results_vector)
  return(results_vector)   ## 36 elements
}


#### fit_spline_buffered ####

## number of 8-day composites to use as buffer on each side of the focal year
## 6 images * 8 days = 48 days of context; events flagged within 21 days of
## year boundary are set to NA as potentially artifactual
n_buffer_images   <- 6
buffer_flag_days  <- 7

## fits the spline using a buffered window (year Y-1 tail + year Y + year Y+1 head)
## to avoid edge artifacts near Jan 1 / Dec 31.
## scale factors are computed from year Y only so the buffer doesn't distort scaling.
## returns a 37-element vector: 36 green-up metrics + rss_spline

fit_spline_buffered <- function(pixel_ts_prev,
                                pixel_ts_curr,
                                pixel_ts_next,
                                n_buffer       = n_buffer_images,
                                flag_days      = buffer_flag_days,
                                spar           = 0.4,
                                ndvi_threshold = 0.05,
                                day_threshold  = 14) {
  
  ## inner copy of fit_weighted_spline for parallel workers
  fit_weighted_spline_inner <- function(x, y, spar = 0.5) {
    if (any(is.na(x)) || any(is.na(y)) ||
        any(is.infinite(x)) || any(is.infinite(y))) return(NULL)
 
    spline_model <- try(smooth.spline(x, y, spar = spar), silent = TRUE)
    
     if (inherits(spline_model, "try-error")) return(NULL)
    spline_model$y <- pmin(pmax(spline_model$y, 0), 1)
    return(spline_model)
  }
  
  if (all(is.na(pixel_ts_curr))) return(rep(NA, 37))
  
  n_curr <- length(pixel_ts_curr)
  
  ## extract buffer from neighboring years, pad with NA if unavailable
  prev_buf <- if (!is.null(pixel_ts_prev) && length(pixel_ts_prev) >= n_buffer) {
    tail(pixel_ts_prev, n_buffer)
  } else {
    rep(NA, n_buffer)
  }
  
  next_buf <- if (!is.null(pixel_ts_next) && length(pixel_ts_next) >= n_buffer) {
    head(pixel_ts_next, n_buffer)
  } else {
    rep(NA, n_buffer)
  }
  
  ts_buffered <- c(prev_buf, pixel_ts_curr, next_buf)
  ts_buffered <- drop_floor_dips(ts_buffered)
  
  ## compute scale factors from (cleaned) year Y only -- not the buffer
  curr_idx <- (n_buffer + 1):(n_buffer + n_curr)
  curr_min <- min(ts_buffered[curr_idx], na.rm = TRUE)
  curr_max <- max(ts_buffered[curr_idx], na.rm = TRUE)
  
  if (is.infinite(curr_min) || is.infinite(curr_max) ||
      curr_max == curr_min) return(rep(NA, 37))
  
  ts_scaled <- scale_values(ts_buffered, x_min = curr_min, x_max = curr_max)
  ts_scaled  <- pmin(pmax(ts_scaled, 0), 1)
  
  ## time points in julian days
  ## each of the n_curr images covers ~8 days; use midpoints
  days_per_image <- 365 / n_curr
  tp_curr <- seq(days_per_image / 2, 365 - days_per_image / 2, length.out = n_curr)
  tp_prev <- tp_curr[1]      - rev(seq_len(n_buffer)) * days_per_image
  tp_next <- tp_curr[n_curr] +      seq_len(n_buffer)  * days_per_image
  
  time_points <- c(tp_prev, tp_curr, tp_next)
  
  ## fit spline across full buffered window
  ## drop missing observations, keeping time/value pairs aligned
  keep <- !is.na(ts_scaled)
  if (sum(keep) < 5) return(rep(NA, 37))
  
  ## fit spline across full buffered window
  spline_model <- fit_weighted_spline_inner(time_points[keep], ts_scaled[keep], spar = spar)
  
  if (is.null(spline_model)) return(rep(NA, 37))
  
  ## predict on daily grid for year Y only (DOY 1-365)
  ## this is what gets passed to calc_greenup_metrics so all output
  ## peak/start/end values are in DOY 1-365
  daily_grid <- 1:365
  pred_y     <- pmin(pmax(predict(spline_model, daily_grid)$y, 0), 1)
  
  ## extract green-up metrics
  metrics <- calc_greenup_metrics(pred_y,
                                  ndvi_threshold = ndvi_threshold,
                                  day_threshold  = day_threshold)
  
  ## output structure (36 elements):
  ##   greenup_peak     [1:4]
  ##   senesc_peak      [5:8]
  ##   greenup_start    [9:12]
  ##   greenup_end      [13:16]
  ##   green_period_end [17:20]
  ##   senesc_start     [21:24]
  ##   senesc_end       [25:28]
  ##   greenup_change   [29:32]
  ##   senesc_change    [33:36]
  
  ## flag events with peaks within flag_days of year boundary as potentially
  ## artifactual -- set peak and its associated start/end/change to NA
  ## green-up peaks: indices 1:4; senescence peaks: indices 5:8
  for (peak_group_start in c(1L, 5L)) {
    for (j in 1:4) {
      idx <- peak_group_start + j - 1L
      val <- metrics[idx]
      if (!is.na(val) && (val < flag_days || val > (365 - flag_days))) {
        ## zero out this event: peak, start, end, green_period_end (if greenup),
        ## and change
        if (peak_group_start == 1L) {
          ## green-up event
          metrics[idx]      <- NA   ## greenup_peak
          metrics[8L  + j]  <- NA   ## greenup_start
          metrics[12L + j]  <- NA   ## greenup_end
          metrics[16L + j]  <- NA   ## green_period_end
          metrics[28L + j]  <- NA   ## greenup_change
        } else {
          ## senescence event
          metrics[idx]      <- NA   ## senesc_peak
          metrics[20L + j]  <- NA   ## senesc_start
          metrics[24L + j]  <- NA   ## senesc_end
          metrics[32L + j]  <- NA   ## senesc_change
        }
      }
    }
  }
  
  ## compute rss over year Y time points only for consistency with existing pipeline
  rss <- tryCatch({
    yr_pred <- pmin(pmax(predict(spline_model, tp_curr)$y, 0), 1)
    yr_obs  <- ts_scaled[(n_buffer + 1):(n_buffer + n_curr)]
    sum((yr_obs - yr_pred)^2, na.rm = TRUE)
  }, error = function(e) NA_real_)
  
  c(metrics, rss)   ## 37 elements total
}