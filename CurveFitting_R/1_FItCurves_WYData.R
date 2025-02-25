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
  
  ## apply temporal smoothing across layers
  smoothed_irg_rast <- app(ndvi.rast, function(x) apply_temporal_smooth(x, window_size = 3))

  spat.ext <- ext(smoothed_irg_rast)
  ext.poly <- as.polygons(spat.ext, crs = crs(smoothed_irg_rast))
  
  random.sample <- spatSample(ext.poly, 500, method= "random")
  random.sample.df <- as.data.frame(random.sample, geom = "XY")
  
  random.sample.df <- data.frame(longitude = random.sample.df$x,
                                 latitude = random.sample.df$y,
                                 date = rep(as.POSIXct("2018-01-01 00:00:01"), nrow(random.sample.df)))
  
  random.sample_obj <- st_as_sf(random.sample.df, coords = c("longitude", "latitude"), crs = crs(smoothed_irg_rast))
  random.sample_obj <- st_transform(random.sample_obj, crs = 4326)
  random.sample_obj$point <- 1:500 + 500*(j-1)
  
  plot(smoothed_irg_rast[[2]])
  plot(st_geometry(st_transform(random.sample_obj, crs=crs(smoothed_irg_rast, proj=TRUE))), add = T)
  
  ## smoothed NDVI
  ndvi.values <- terra::extract(smoothed_irg_rast, 
                                st_transform(random.sample_obj,
                                             crs=crs(smoothed_irg_rast, proj=TRUE)),
                                na.rm = T)
  
  dates <- list(names(ndvi.rast))
  dates <- data.frame(date = dates[[1]],
                      ID = 1:length(dates[[1]]))
  
  ndvi.values <- reshape(ndvi.values,
                         direction = "long", 
                         varying = list(names(ndvi.values)[2:ncol(ndvi.values)]),
                         v.names = "NDVI",
                         timevar = "ID")
  
  ndvi.values$point <- ndvi.values$id + 500*(j-1)
  ndvi.values$id <- NULL
  
  ndvi.values <- ndvi.values[order(ndvi.values$point, ndvi.values$ID),]
  ndvi.values <- merge(ndvi.values, dates, by = "ID")
  ndvi.values$date <- substr(ndvi.values$date, 1, 10)
  
  ndvi.values <- ndvi.values[order(ndvi.values$point, ndvi.values$ID),]
  ndvi.values$date <- as.Date(ndvi.values$date, "%Y_%m_%d")
  ndvi.values$date <- as.Date(substr(names(ndvi.rast), 1, 10), "%Y_%m_%d")
  head(ndvi.values)
  
  ## raw NDVI
  ndvi.values_raw <- terra::extract(ndvi.rast, 
                                    st_transform(random.sample_obj,
                                                 crs=crs(ndvi.rast, proj=TRUE)),
                                    na.rm = T)
  
  ndvi.values_raw <- reshape(ndvi.values_raw,
                             direction = "long", 
                             varying = list(names(ndvi.values_raw)[2:ncol(ndvi.values_raw)]),
                             v.names = "NDVI_raw",
                             timevar = "ID")
  
  ndvi.values_raw$point <- ndvi.values_raw$id + 500*(j-1)
  ndvi.values_raw <- merge(ndvi.values_raw, dates, by = "ID")
  ndvi.values_raw$date <- substr(ndvi.values_raw$date, 1, 10)
  ndvi.values_raw$date <- as.Date(ndvi.values_raw$date, "%Y_%m_%d")
  ndvi.values_raw$ID <- NULL
  ndvi.values_raw$id <- NULL
  ndvi.values_raw <- ndvi.values_raw[order(ndvi.values_raw$point, ndvi.values_raw$date),]
  
  ## merge sf object and NDVI
  random.sample_obj <- merge(random.sample_obj, ndvi.values, by = "point")
  
  ## merge raw and smoothed
  ndvi.values <- merge(ndvi.values, ndvi.values_raw, by = c("point", "date"))
  ndvi.values <- ndvi.values[order(ndvi.values$point, ndvi.values$date),]

  ## test with a bunch of random points
  
  for(i in unique(ndvi.values$point)){
    
    ## get a point and rescale the values
    tmp <- ndvi.values[ndvi.values$point == i,]
    tmp$jday <- yday(tmp$date)
    
    if(all(is.na(tmp$NDVI))){
      
      print(paste0("Skipped point ", i, " pixel values null"))
    } else{
      
      tmp$rescale <- scale_values(tmp$NDVI)
      
      ## lists to store results
      second_derivatives_increasing <- list()
      second_derivatives_decreasing <- list()
      
      ## fit a spline to the time series
      spl <- fit_spline_safe(tmp$ID, tmp$rescale, spar = 0.4)
      
      if(is.null(spl)){
        
        png(paste0("G:/Shared drives/wyo-coop-lasharr/PostDoc/CurveFitting/", folder_of_interest, "/", folder_of_interest, "_RandomPoint_", i, '.png'), 1600, 1000)
        
        par(mfrow = c(1,2))
        plot(tmp$ID, tmp$NDVI,
             type = "n", ylim = c(0,1), ylab = "NDVI", xlab = "Date",
             main = paste(random.sample_obj[random.sample_obj$point==i,]$geometry[1]), axes = F)
        axis(2, seq(0, 1, 0.2))
        axis(1, seq(1, 50, 10),
             labels = tmp$date[seq(1, 50, 10)])
        box(bty = "l")
        points(tmp$ID, tmp$NDVI_raw, col = alpha("blue", 0.5), pch = 16)
        lines(tmp$ID, tmp$NDVI_raw, col = alpha("blue", 0.5))
        points(tmp$ID, tmp$NDVI, col = alpha("gray50", 0.5), pch = 16)
        lines(tmp$ID, tmp$NDVI, col = alpha("gray50", 0.5))
        points(tmp$ID, tmp$rescale, col = "black", pch = 16)
        lines(tmp$ID, tmp$rescale, col = "black")
        
        legend('topleft', c("Processed", "Smoothed"),
               col = c("blue", 'gray'), pch = 16, bty = "n")
        
        plot(ndvi.rast[[22]], main = i, col = terrain.colors(100))
        plot(st_geometry(random.sample_obj[random.sample_obj$point==i,]), add = T,
             pch = 16, col = "red")
        
        dev.off()

        print(paste0("Point ", i, " could not fit a spline."))
        
        
      }else{
        
        yr <- year(tmp[1,]$date)
        
        # days <- if(yr %in% c(2024, 2020, 2016, 2012, 2008, 2004, 2000)){1:366}else{1:365}
        days <- 1:365
        
        ## fit the spline with julian date and scaled NDVI values
        spl <- fit_weighted_spline(tmp$jday, tmp$rescale, spar = 0.4)
        
        ## get the greenup metrics
        greenup_metrics <- calc_greenup_metrics(spl, tmp)
        
        if(is.null(greenup_metrics)){
          
          print(paste0("No significant changes in spline for point ", i))
          
        }else{
        
        ## get the daily values and prevent predictions from going below 0
        predictions <- predict(spl, x = days)
        predictions$y <- pmax(predictions$y, 0)  ## ensure predictions don't go below 0
        
        ## then unscale it so it's closer to the raw data
        min_val <- min(tmp$NDVI)
        max_val <- max(tmp$NDVI)   
        
        unscaled_predictions <- predictions$y * (max_val - min_val) + min_val

        dlc_estimates <- dlc_coefs(tmp)
        
        ## then plot it
        png(paste0("G:/Shared drives/wyo-coop-lasharr/PostDoc/CurveFitting/", folder_of_interest, "/", folder_of_interest, "_RandomPoint_", i, '.png'), 1600, 1000)
        
        par(mfrow = c(1,2))
        plot(tmp$jday, tmp$NDVI,
             type = "n", ylim = c(0,1), ylab = "NDVI", xlab = "Date",
             main = paste(random.sample_obj[random.sample_obj$point==i,]$geometry[1]), axes = F)
        axis(2, seq(0, 1, 0.2))
        axis(1, tmp$jday[seq(1, nrow(tmp), 5)],
             labels = tmp$jday[seq(1, nrow(tmp), 5)])
        box(bty = "l")

        points(tmp$jday, tmp$NDVI_raw, col = alpha("blue", 0.5), pch = 16)
        lines(tmp$jday, tmp$NDVI_raw, col = alpha("blue", 0.5))
        points(tmp$jday, tmp$NDVI, col = alpha("gray50", 0.5), pch = 16)
        lines(tmp$jday, tmp$NDVI, col = alpha("gray50", 0.5))
        points(tmp$jday, tmp$rescale, col = "black", pch = 16)
        lines(tmp$jday, tmp$rescale, col = "black")
        lines(spl, col = "red")
        abline(v = greenup_metrics$inflection_points_increasing, col = "purple", pch = 19, lwd = 2)
        
        lines(unscaled_predictions, col = "blue")
        
        legend('topleft', c("Processed", "Smoothed", "Scaled", "Spline", "Green up Periods"),
               col = c("blue", 'gray', 'black', 'red', 'purple'), pch = 16, bty = "n")
        
        plot(ndvi.rast[[22]], main = i, col = terrain.colors(100))
        plot(st_geometry(random.sample_obj[random.sample_obj$point==i,]), add = T,
             pch = 16, col = "red")
        
        dev.off()
        
        ## save the output
        
        ## get the location to save for later
        tmp.point <- random.sample_obj[random.sample_obj$point==i,]
        tmp.point <- tmp.point[1, "geometry"]
        
        ## create a dataframe of the daily estimate predictions
        tmp_spline_output <- data.frame(jday = 1:365,
                                        raw_NDVI = unscaled_predictions,
                                        scaled_NDVI = predictions$y)
        
        ## get all of the greenup metrics
        max_greenup <- greenup_metrics$inflection_points_increasing
        greenup_start <- greenup_metrics$greenup_start
        greenup_end <- greenup_metrics$greenup_end
        
        ## get all of the senescence metrics
        max_sen <- greenup_metrics$inflection_points_decreasing
        sen_start <- greenup_metrics$senesc_start
        sen_end <- greenup_metrics$senesc_end
        
        ## get the year
        yr <- year(tmp$date)[1]
        
        ## get the unscaled predictions to see the actual change in NDVI with green up and senescence
        greenup_change <- unscaled_predictions[greenup_end] - unscaled_predictions[greenup_start]
        
        sen_change <- unscaled_predictions[sen_start] - unscaled_predictions[sen_end]
        
        
        tmp_spline_output <- list(tmp_spline_output, tmp, tmp.point,
                                  max_greenup, greenup_start, greenup_end, greenup_change,
                                  max_sen, sen_start, sen_end, sen_change,
                                  yr, dlc_estimates)
        
        NDVI_splines[[i]] <- tmp_spline_output

        rm(tmp, second_derivatives_decreasing, second_derivatives_increasing, result, 
           inflection_points_increasing, inflection_points_decreasing, spl, 
           max_greenup, greenup_start, greenup_end, greenup_change,
           max_sen, sen_start, sen_end, sen_change, yr, dlc_estimates)
        
        }
        
      }
      
    }
    
    print(i)
    
  }
  
  print(paste0("FINISHED WITH RASTER ", list.rast[[j]]))
  
}

saveRDS(NDVI_splines, paste0("G:/Shared drives/wyo-coop-lasharr/PostDoc/CurveFitting/", folder_of_interest, "_NDVIExtracted_MultiYear.rds"))

