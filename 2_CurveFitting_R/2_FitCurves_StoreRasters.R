library(tidyverse)
library(sf)
library(terra)
library(raster)
library(zoo)
library(lubridate)
library(scales)

folder_of_interest <- "Utah"

## source functions
source("C:/Users/Tayler/Documents/GitHub/DroughtSurfing/2_CurveFitting_R/Functions/Functions_FitSplines_Raster.R")
source("C:/Users/Tayler/Documents/GitHub/DroughtSurfing/2_CurveFitting_R/Functions/Functions_FitDLC.R")


#### paths and parameters ####

ndvi_base_path  <- paste0("G:/My Drive/NDVI_Merged/", folder_of_interest)
out_path        <- paste0("G:/Shared drives/wyo-coop-lasharr/PostDoc/CurveFitting/Wyoming_Updated/")

## number of 8-day composites to use as buffer on each side of the focal year
## 6 images * ~8 days = ~48 days of context on each side
n_buffer_images  <- 6
buffer_flag_days <- 7   ## events within 7 doy of year boundary get flagged as NA

## band names for the output raster -- 36 metric bands + rss_spline + dlc bands
## note: green_period_end is new (between greenup_end and senesc_start)
metric_band_names <- c(
  "greenup_peak_1",     "greenup_peak_2",     "greenup_peak_3",     "greenup_peak_4",
  "senescence_peak_1",  "senescence_peak_2",  "senescence_peak_3",  "senescence_peak_4",
  "greenup_start_1",    "greenup_start_2",    "greenup_start_3",    "greenup_start_4",
  "greenup_end_1",      "greenup_end_2",      "greenup_end_3",      "greenup_end_4",
  "green_period_end_1", "green_period_end_2", "green_period_end_3", "green_period_end_4",
  "senesc_start_1",     "senesc_start_2",     "senesc_start_3",     "senesc_start_4",
  "senesc_end_1",       "senesc_end_2",       "senesc_end_3",       "senesc_end_4",
  "greenup_change_1",   "greenup_change_2",   "greenup_change_3",   "greenup_change_4",
  "senesc_change_1",    "senesc_change_2",    "senesc_change_3",    "senesc_change_4",
  "rss_spline"
)

dlc_band_names <- c(
  "xmidS", "scalS", "xmidA", "scalA",
  "xmidS_SE", "scalS_SE", "xmidA_SE", "scalA_SE",
  "maxNDVI", "maxNDVIdate",
  "maxIRG", "maxIRGdate",
  "sumIRG", "springStart", "springEnd", "rss_dlc"
)

all_band_names <- c(metric_band_names, dlc_band_names)


#### load all available rasters ####

list.rast <- list.files(ndvi_base_path)
## adjust index range as needed for the years you want to process
list.rast <- list.rast[17:23]

## extract year strings for easy lookup
yr_strings <- substr(list.rast, nchar(list.rast) - 7, nchar(list.rast) - 4)

print(paste0("processing years: ", paste(yr_strings, collapse = ", ")))


#### main processing loop ####

## we process each year but load the neighboring years' rasters to get
## the buffer time series avoids edge artifacts near jan 1 / dec 31

NDVI_splines <- lapply(seq_along(list.rast), function(j) {
  
  yr <- yr_strings[j]
  print(paste0("starting year ", yr, " (", j, " of ", length(list.rast), ")"))
  
  ## ---- load focal year raster ----
  ndvi_curr <- rast(file.path(ndvi_base_path, list.rast[[j]]))
  dates_curr <- names(ndvi_curr)
  n_images_curr <- nlyr(ndvi_curr)
  
  ## ---- load previous year raster for buffer ----
  ndvi_prev <- if (j > 1) {
    rast(file.path(ndvi_base_path, list.rast[[j - 1]]))
  } else {
    NULL
  }
  
  ## ---- load next year raster for buffer ----
  ndvi_next <- if (j < length(list.rast)) {
    rast(file.path(ndvi_base_path, list.rast[[j + 1]]))
  } else {
    NULL
  }
  
  
  #### filter for valid pixels (require > 4 non-na observations) ####
  
  modisdata_curr <- app(ndvi_curr, function(x) {
    if (sum(!is.na(x)) <= 4) return(rep(NA, length(x)))
    else return(x)
  }, cores = parallel::detectCores() - 2)
  
  
  #### scale year Y using its own min/max ####
  ## scale_values now accepts explicit x_min and x_max, but here we pass
  ## the full pixel vector and let it compute min/max from year Y only
  ## (the buffer scaling happens inside fit_spline_buffered)
  
  scaled_curr <- app(modisdata_curr, scale_values,
                     cores = parallel::detectCores() - 2)
  message(paste0(yr, ": finished scaling"))
  
  
  #### fit DLC (unchanged -- operates on scaled year Y only) ####
  
  dlcparams <- app(scaled_curr, fun = fundlc,
                   cores = parallel::detectCores() - 2)
  message(paste0(yr, ": finished fitting DLC"))
  
  
  #### extract buffer layers from neighboring years ####
  ## we extract the last n_buffer_images from year Y-1 and the first
  ## n_buffer_images from year Y+1 as plain matrices so we can pass them
  ## per-pixel into fit_spline_buffered
  
  ## previous year buffer: last n_buffer_images layers
  prev_mat <- if (!is.null(ndvi_prev)) {
    n_prev <- nlyr(ndvi_prev)
    buf_idx <- max(1, n_prev - n_buffer_images + 1):n_prev
    values(ndvi_prev[[buf_idx]])   ## n_pixels x n_buffer matrix
  } else {
    matrix(NA, nrow = ncell(modisdata_curr), ncol = n_buffer_images)
  }
  
  ## next year buffer: first n_buffer_images layers
  next_mat <- if (!is.null(ndvi_next)) {
    buf_idx <- 1:n_buffer_images
    values(ndvi_next[[buf_idx]])   ## n_pixels x n_buffer matrix
  } else {
    matrix(NA, nrow = ncell(modisdata_curr), ncol = n_buffer_images)
  }
  
  ## current year values as matrix
  curr_mat <- values(modisdata_curr)   ## n_pixels x n_images matrix
  
  message(paste0(yr, ": buffer matrices ready"))
  
  
  #### fit buffered spline per pixel ####
  ## terra::app doesn't easily support multi-raster inputs, so we use
  ## a manual pixel loop here. for large rasters you may want to
  ## parallelise this with parallel::mclapply or future.apply
  
  n_pixels   <- nrow(curr_mat)
  n_out      <- 37   ## 36 metrics + rss_spline
  
  spline_metrics_mat <- matrix(NA, nrow = n_pixels, ncol = n_out)
  
  ## use parallel processing via parallel::mclapply
  n_cores <- parallel::detectCores() - 2
  
  results_list <- parallel::mclapply(seq_len(n_pixels), function(i) {
    fit_spline_buffered(
      pixel_ts_prev  = as.numeric(prev_mat[i, ]),
      pixel_ts_curr  = as.numeric(curr_mat[i, ]),
      pixel_ts_next  = as.numeric(next_mat[i, ]),
      n_buffer       = n_buffer_images,
      flag_days      = buffer_flag_days,
      spar           = 0.4,
      ndvi_threshold = 0.05,
      day_threshold  = 14
    )
  }, mc.cores = n_cores)
  
  ## pack results into matrix
  for (i in seq_len(n_pixels)) {
    res <- results_list[[i]]
    if (!is.null(res) && length(res) == n_out) {
      spline_metrics_mat[i, ] <- res
    }
  }
  
  message(paste0(yr, ": finished fitting buffered splines"))
  
  
  #### convert metrics matrix back to raster ####
  
  ## use the current year raster as a spatial template
  template_r <- modisdata_curr[[1]]
  
  ## build a raster stack from the metrics matrix (one layer per metric)
  greenup_metrics_rast <- rast(lapply(seq_len(n_out), function(k) {
    r        <- template_r
    values(r) <- spline_metrics_mat[, k]
    r
  }))
  
  ## separate rss_spline (last column) from the metric bands
  rss_spline_rast      <- greenup_metrics_rast[[n_out]]
  greenup_metrics_rast <- greenup_metrics_rast[[1:(n_out - 1)]]   ## 36 metric bands
  
  
  #### min/max and unscaled predictions ####
  ## unscaling uses the year Y min/max, consistent with how we scaled
  
  min_vals    <- app(ndvi_curr, fun = min, na.rm = TRUE,
                     cores = parallel::detectCores() - 2)
  max_vals    <- app(ndvi_curr, fun = max, na.rm = TRUE,
                     cores = parallel::detectCores() - 2)
  
  ## refit spline predictions for unscaling (same as before -- 365-day output)
  spline_out  <- app(scaled_curr, fit_spline_raster,
                     cores = parallel::detectCores() - 2)
  message(paste0(yr, ": finished fitting spline for unscaled predictions"))
  
  unscaled_pred <- (spline_out * (max_vals - min_vals)) + min_vals
  message(paste0(yr, ": finished unscaling predictions"))
  
  
  #### combine all output bands ####
  
  greenup_metrics_all <- c(greenup_metrics_rast,
                           rss_spline_rast,
                           dlcparams)
  
  ## assign band names
  names(greenup_metrics_all) <- all_band_names
  
  ## save outputs
  saveRDS(greenup_metrics_all,
          file.path(out_path, paste0(folder_of_interest, "_SplineMetrics_Raster_BufferedTest_", yr, ".rds")))
  writeRaster(greenup_metrics_all,
              file.path(out_path, paste0(folder_of_interest, "_SplineMetrics_Raster_BufferedTest_", yr, ".tif")),
              overwrite = TRUE)
  
  message(paste0("FINISHED WITH RASTER ", list.rast[[j]]))
  
  ## clean up large objects before next iteration
  rm(ndvi_curr, ndvi_prev, ndvi_next, modisdata_curr, scaled_curr,
     prev_mat, curr_mat, next_mat, spline_metrics_mat, results_list,
     greenup_metrics_rast, spline_out, unscaled_pred)
  gc()
  
  return(greenup_metrics_all)
})
