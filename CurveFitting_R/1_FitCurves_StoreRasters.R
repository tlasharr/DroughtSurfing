library(tidyverse)
library(sf)
library(terra)
library(raster)
library(zoo)
library(lubridate)
library(scales)

<<<<<<< Updated upstream
folder_of_interest <- "Wyoming"
=======
folder_of_interest <- "Yellowstone"
>>>>>>> Stashed changes

####### Source functions

source("C:/Users/tlasharr/Documents/GitHub/DroughtSurfing/CurveFitting_R/Functions/Functions_FitSplines_Raster.R")
source("C:/Users/tlasharr/Documents/GitHub/DroughtSurfing/CurveFitting_R/Functions/Functions_FitDLC.R")

########## read in the data 

list.rast <- list.files(paste0("H:/My Drive/NDVI_Merged/", folder_of_interest))

NDVI_splines <- lapply(seq_along(list.rast), function(j) {
  
  ## get the raster we care about
  ndvi.rast <- rast(paste0("H:/My Drive/NDVI_Merged/", folder_of_interest, "/", list.rast[[j]]))
  
  #### get the dates
  dates <- names(ndvi.rast)
  ydates <- yday(as.Date(substr(dates, 1, 10), "%Y_%m_%d"))
  
  #### filter for valid pixels
  Sys.time()
  modisdata2 <- app(ndvi.rast, function(x) {
    if (sum(!is.na(x)) <= 4) return(rep(NA, length(x)))
    else return(x)
  },
  cores = parallel::detectCores() - 4)
  
  ## get julian dates
  jul2 <- yday(as.Date(substr(names(ndvi.rast), 1, 10), "%Y_%m_%d"))
  
  ## scale the raster
  Sys.time()
  scaled_r <- app(modisdata2, scale_values,
                  cores = parallel::detectCores() - 4)
  
  ## fit the DLC
  Sys.time()
  dlcparams <- app(scaled_r, fun = fundlc,
                   cores = parallel::detectCores() - 4)  

  #### fit spline and predict 365 days
  Sys.time()
  spline_out <- app(scaled_r, fit_spline_raster,
                    cores = parallel::detectCores() - 4)

  #### min/max values
  min_vals <- app(ndvi.rast, fun = min, na.rm = TRUE,
                  cores = parallel::detectCores() - 4)
  max_vals <- app(ndvi.rast, fun = max, na.rm = TRUE,
                  cores = parallel::detectCores() - 4)
  
  #### unscale predictions
  Sys.time()
  unscaled_pred <- (spline_out * (max_vals - min_vals)) + min_vals
  
  ## calculate green-up metrics
  # Sys.time()
  # greenup_metrics_list <- app(unscaled_pred, function(pixel_values) {
  #   calc_greenup_metrics(pixel_values, ndvi_threshold = 0.05, day_threshold = 14)
  # }, cores = parallel::detectCores() - 4)
  # 
  # Sys.time()
  
  greenup_metrics_list <- app(unscaled_pred,
                              fun = calc_greenup_metrics,
                              cores = parallel::detectCores() - 4)
  Sys.time()
  
  greenup_metrics_all <- c(greenup_metrics_list, dlcparams)
  
  ## assign names to the output
  names(greenup_metrics_all)  <- c("greenup_peak_1", "greenup_peak_2",
                                    "greenup_peak_3", "greenup_peak_4",
                                    "senescence_peak_1", "senescence_peak_2",
                                    "senescence_peak_3", "senescence_peak_4",
                                    "greenup_start_1", "greenup_start_2",
                                    "greenup_start_3", "greenup_start_4",
                                    "greenup_end_1", "greenup_end_2",
                                    "greenup_end_3", "greenup_end_4",
                                    "senesc_start_1", "senesc_start_2",
                                    "senesc_start_3", "senesc_start_4",
                                    "senesc_end_1", "senesc_end_2",
                                    "senesc_end_3", "senesc_end_4",
                                    "greenup_change_1", "greenup_change_2",
                                    "greenup_change_3", "greenup_change_4",
                                    "senesc_change_1", "senesc_change_2",
                                    "senesc_change_3", "senesc_change_4",
                                    "xmidS", "scalS", "xmidA", "scalA",
                                    "xmidS_SE", "scalS_SE",
                                    "xmidA_SE", "scalA_SE",
                                    "maxNDVI", "maxNDVIdate",
                                    "maxIRG", "maxIRGdate",
                                    "sumIRG", "springStart", "springEnd")
  
  ## save output
  yr <- substr(list.rast[[j]], nchar(list.rast[[j]])-7,  nchar(list.rast[[j]])-4)
  saveRDS(greenup_metrics_list, paste0("G:/Shared drives/wyo-coop-lasharr/PostDoc/CurveFitting/", folder_of_interest, "_SplineMetrics_Raster_", yr, ".rds"))
  
  print(paste0("FINISHED WITH RASTER ", list.rast[[j]]))
  
  return(greenup_metrics_list)
})

