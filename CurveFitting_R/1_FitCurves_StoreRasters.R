library(tidyverse)
library(sf)
library(terra)
library(raster)
library(zoo)
library(lubridate)
library(scales)

folder_of_interest <- "Sublette"

####### Source functions

source("G:/Shared drives/wyo-coop-lasharr/PostDoc/Exploration_R/Functions/Functions_FitSplines.R")
source("G:/Shared drives/wyo-coop-lasharr/PostDoc/Exploration_R/Functions/Functions_FitDLC.R")

########## read in the data 

list.rast <- list.files(paste0("H:/My Drive/NDVI_Merged/", folder_of_interest))


########## create an empty list to store the output and loop through each year

NDVI_splines <- list()

for(j in 1:length(list.rast)){
  
  ## get the raster we care about
  ndvi.rast <- rast(paste0("H:/My Drive/NDVI_Merged/", folder_of_interest, "/", list.rast[[j]]))
  
  #### apply temporal smoothing across layers
  # smoothed_irg_rast <- app(ndvi.rast, function(x) apply_temporal_smooth(x, window_size = 3))

  smoothed_irg_rast <- ndvi.rast
  
  #### get the dates
  dates <- list(names(ndvi.rast))
  ydates <- yday(as.Date(substr(dates[[1]], 1, 10), "%Y_%m_%d"))
 
  #### get an empty rasters to store things in
  unscaled_rast <-rast(ndvi.rast, nlyrs = 365)
  scaled_rast <-rast(ndvi.rast, nlyrs = 365)
  values(unscaled_rast) <- NA
  values(scaled_rast) <- NA
  
  #### scale the raster
  scaled_r <- app(smoothed_irg_rast, scale_values)
  
  #### fit a spline and predict out the values to 365 days a year
  spline_out <- app(scaled_r, fit_spline_raster)
  
  #### calculate min and max for each pixel
  min_vals <- app(smoothed_irg_rast, fun = min, na.rm = TRUE)
  max_vals <- app(smoothed_irg_rast, fun = max, na.rm = TRUE)
  
  #### get the unscaled values
  unscaled_pred <- (spline_out * (max_vals - min_vals)) + min_vals
  
  ## calculate green up metrics
  greenup_metrics_list <- app(unscaled_pred, calc_greenup_metrics)
  
  ## save the names of the metrics we want layers for
  metric_names <- c("inflection_points_increasing", "inflection_points_decreasing", 
                    "greenup_start", "greenup_end", "senesc_start", "senesc_end")
  
  
  ## create a new raster with a layer for each metric
  metrics_rast <- rast(scaled_r, nlyr=length(metric_names))
  names(metrics_rast) <- metric_names
  
  ## fill the raster with data from greenup_metrics_list
  
  for (i in 1:length(greenup_metrics_list)) {
    
    metrics_item <- greenup_metrics_list[[i]]
    
    # Check if the entire item is NA (assuming it's a list with a single NA)
    if (length(metrics_item) == 1 && is.na(unlist(metrics_item))) {
      
      values(metrics_rast)[i, ] <- NA  # Assign NA to all layers for this cell
    
      } else {
        
      ## loop through and fill each layer with the respective metric
        
      for (metric_name in metric_names) {
        
        ## check if there are more than green up or senescence period. If there are 
        ## more than 1, we are going to store them all together in the layer and will have to seperate them out
        ## at a later time. So for example, if you had green up on julian dates 150 and 322, it would be
        ## stored as 150322. Additionally, this pads the values so they all have 0s in front, but 
        ## because the values are stored as a number, the leading 0s for the very first one will not
        ## keep the 0s. So when we read them back out we will have to figure out if it's divisible by 3 and 
        ## add appropriate 0s there.
        
        if(length(unlist(greenup_metrics_list[[i]][metric_name])) > 1){
          
          tmpvals <- unlist(greenup_metrics_list[[i]][metric_name])
          tmpvals <- sprintf("%03d", tmpvals)
          tmpvals <- paste0(tmpvals, collapse = "")
          
          values(metrics_rast)[i, ][[which(names(metrics_rast) == metric_name)]] <- as.numeric(tmpvals)
          
        }else{
          
          values(metrics_rast)[i, ][[which(names(metrics_rast) == metric_name)]] <- unlist(greenup_metrics_list[[i]][metric_name])
          
        }
        
      }
    }
  }
  
  print(paste0("FINISHED WITH RASTER ", list.rast[[j]]))
  
  
  yr <- substr(list.rast[[j]], nchar(list.rast[[j]])-7,  nchar(list.rast[[j]])-4)
  
  saveRDS(metrics_rast, paste0("G:/Shared drives/wyo-coop-lasharr/PostDoc/CurveFitting/", folder_of_interest, "_SplineMetrics_Raster_", yr, ".rds"))

}


