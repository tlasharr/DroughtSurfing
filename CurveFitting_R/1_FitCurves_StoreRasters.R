library(tidyverse)
library(sf)
library(terra)
library(raster)
library(zoo)
library(lubridate)
library(scales)

folder_of_interest <- "Sublette"

####### Source functions

source("C:/Users/tlasharr/Documents/GitHub/DroughtSurfing/CurveFitting_R/Functions/Functions_FitSplines_Raster.R")
source("C:/Users/tlasharr/Documents/GitHub/DroughtSurfing/CurveFitting_R/Functions/Functions_FitDLC.R")

########## read in the data 

list.rast <- list.files(paste0("H:/My Drive/NDVI_Merged/", folder_of_interest))


########## create an empty list to store the output and loop through each year

NDVI_splines <- list()

for(j in 1:length(list.rast)){
  
  ## get the raster we care about
  ndvi.rast <- rast(paste0("H:/My Drive/NDVI_Merged/", folder_of_interest, "/", list.rast[[j]]))
  
  #### apply temporal smoothing across layers
  # smoothed_irg_rast <- app(ndvi.rast, function(x) apply_temporal_smooth(x, window_size = 3))

  #### get the dates
  dates <- list(names(ndvi.rast))
  ydates <- yday(as.Date(substr(dates[[1]], 1, 10), "%Y_%m_%d"))
 
  #### scale the raster
  scaled_r <- app(ndvi.rast, scale_values)
  
  #### fit a spline and predict out the values to 365 days a year
  spline_out <- app(scaled_r, fit_spline_raster)
  
  #### calculate min and max for each pixel
  min_vals <- app(ndvi.rast, fun = min, na.rm = TRUE)
  max_vals <- app(ndvi.rast, fun = max, na.rm = TRUE)
  
  #### get the unscaled values
  unscaled_pred <- (spline_out * (max_vals - min_vals)) + min_vals
  
  ## calculate green up metrics
  greenup_metrics_list <- app(unscaled_pred, calc_greenup_metrics)
  
  ## started at 12:31
  
  ## save the names of the metrics as the layer names
  names(greenup_metrics_list)  <- c("inflection_points_increasing",
                                    "inflection_points_decreasing", 
                                    "greenup_start", "greenup_end",
                                    "senesc_start", "senesc_end")
  
  greenup_metrics_list
  
  yr <- substr(list.rast[[j]], nchar(list.rast[[j]])-7,  nchar(list.rast[[j]])-4)
  
  saveRDS(greenup_metrics_list, paste0("G:/Shared drives/wyo-coop-lasharr/PostDoc/CurveFitting/", folder_of_interest, "_SplineMetrics_Raster_", yr, ".rds"))

  print(paste0("FINISHED WITH RASTER ", list.rast[[j]]))
  
}


