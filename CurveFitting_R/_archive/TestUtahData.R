library(tidyverse)
library(sf)
library(terra)
library(raster)
library(zoo)

## read in the data
irg.rast1 <- rast("H:/My Drive/NDVI_Exports/Merged/UT_NDVI_2021.tif")
irg.rast2 <- rast("H:/My Drive/NDVI_Exports/Merged/UT_NDVI_2022.tif")
irg.rast3 <- rast("H:/My Drive/NDVI_Exports/Merged/UT_NDVI_2023.tif")

irg.rast1 <- irg.rast1[[43:46]]
irg.rast3 <- irg.rast3[[1:4]]

irg.rast <- c(irg.rast1, irg.rast2, irg.rast3)


# Apply temporal smoothing across layers
smoothed_irg_rast <- app(irg.rast, function(x) temporal_smooth(x, window_size = 3))

spat.ext <- ext(smoothed_irg_rast)
ext.poly <- as.polygons(spat.ext, crs = crs(smoothed_irg_rast))

random.sample <- spatSample(ext.poly, 500, method= "random")
random.sample.df <- as.data.frame(random.sample, geom = "XY")

random.sample.df <- data.frame(longitude = random.sample.df$x,
                               latitude = random.sample.df$y,
                               date = rep(as.POSIXct("2018-01-01 00:00:01"), nrow(random.sample.df)))

random.sample_obj <- st_as_sf(random.sample.df, coords = c("longitude", "latitude"), crs = crs(smoothed_irg_rast))
random.sample_obj <- st_transform(random.sample_obj, crs = 4326)
random.sample_obj$point <- 1:500

plot(smoothed_irg_rast[[2]])
plot(st_geometry(st_transform(random.sample_obj, crs=crs(smoothed_irg_rast, proj=TRUE))), add = T)

## smoothed NDVI
ndvi.values <- terra::extract(smoothed_irg_rast, 
                              st_transform(random.sample_obj,
                                           crs=crs(smoothed_irg_rast, proj=TRUE)),
                              na.rm = T)

dates <- list(names(ndvi.values)[2:ncol(ndvi.values)])
dates <- data.frame(date = dates[[1]],
                    ID = 1:length(dates[[1]]))

ndvi.values <- reshape(ndvi.values,
                       direction = "long", 
                       varying = list(names(ndvi.values)[2:ncol(ndvi.values)]),
                       v.names = "NDVI",
                       timevar = "ID")
ndvi.values$point <- ndvi.values$id
ndvi.values$id <- NULL

ndvi.values <- ndvi.values[order(ndvi.values$point, ndvi.values$ID),]
ndvi.values <- merge(ndvi.values, dates, by = "ID")
ndvi.values$date <- substr(ndvi.values$date, 1, 10)

ndvi.values <- ndvi.values[order(ndvi.values$point, ndvi.values$ID),]
ndvi.values$date <- as.Date(ndvi.values$date, "%Y_%m_%d")
head(ndvi.values)

## raw NDVI
ndvi.values_raw <- terra::extract(irg.rast, 
                                  st_transform(random.sample_obj,
                                               crs=crs(irg.rast, proj=TRUE)),
                                  na.rm = T)

dates <- list(names(ndvi.values_raw)[2:ncol(ndvi.values_raw)])
dates <- data.frame(date = dates[[1]],
                    ID = 1:length(dates[[1]]))

ndvi.values_raw <- reshape(ndvi.values_raw,
                           direction = "long", 
                           varying = list(names(ndvi.values_raw)[2:ncol(ndvi.values_raw)]),
                           v.names = "NDVI_raw",
                           timevar = "ID")

ndvi.values_raw$point <- ndvi.values_raw$id
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

setwd("G:/Shared drives/wyo-coop-lasharr/PostDoc/CurveFitting/Paunsagunt")

for(i in unique(ndvi.values$point)[1:50]){
  
  ## get a point and rescale the values
  tmp <- ndvi.values[ndvi.values$point == i,]
  
  if(all(is.na(tmp$NDVI))){
    
    print(paste0("Skipped point ", i, " pixel values null"))
  }
  else{
    
    tmp$rescale <- scale_values(tmp$NDVI)
    
    ## lists to store results
    second_derivatives_increasing <- list()
    second_derivatives_decreasing <- list()
    
    ## fit a spline to the time series
    spl <- fit_spline_safe(tmp$ID, tmp$rescale, spar = 0.5)
    
    if(is.null(spl)){
      
      png(paste0("Paun_RandomPoint_", i, '.png'), 1600, 1000)
      
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
      
      plot(irg.rast[[22]], main = i)
      plot(st_geometry(random.sample_obj[random.sample_obj$point==i,]), add = T,
           pch = 16, col = "red")
      
      dev.off()
      
      
    }else{
      
      spl <- fit_weighted_spline(tmp$ID, tmp$rescale)
      
      pred <- predict(spl, x = tmp$ID, deriv = 1)
      
      ## figure out if it's increasing, decreasing, or relatively flat
      inc <- diff(spl$y) > 0
      inc <- c(NA, inc)
      
      dec <- diff(spl$y) < 0
      dec <- c(NA, dec)
      
      flat <- abs(diff(spl$y)) < 0.005 
      flat <- c(NA, flat)
      
      ## create a dataframe
      result <- data.frame(x = spl$x, y = spl$y, inc = inc, dec = dec, flat = flat)
      result[1,"inc"] <- result[2, "inc"] ## set the increasing in point 1 as what it is in point 2
      result[1,"dec"] <- result[2, "dec"] ## same as above
      result[1,"flat"] <- result[2, "flat"] ## same as above
      
      ## remove all the 'flat' pieces
      flat <- result[result$flat==TRUE,]
      result <- result[!result$flat==TRUE,]
      result$y <- ifelse(result$y < 0, 0, result$y)
      result$y <- ifelse(result$y > 1, 1, result$y)
      
      ## Find inflection points for increasing and decreasing segments
      inflection_points_increasing <- process_segments_for_logistic_inflection(result, "inc")
      inflection_points_decreasing <- process_segments_for_logistic_inflection(result, "dec")
      
      ## then plot
      png(paste0("Paun_RandomPoint_", i, '.png'), 1600, 1000)
      
      par(mfrow = c(1,2))
      plot(tmp$ID, tmp$NDVI,
           type = "n", ylim = c(0,1), ylab = "NDVI", xlab = "Date",
           main = paste(random.sample_obj[random.sample_obj$point==i,]$geometry[1]), axes = F)
      axis(2, seq(0, 1, 0.2))
      axis(1, seq(1, 50, 10),
           labels = tmp$date[seq(1, 50, 10)])
      box(bty = "l")
      
      abline(v = 5, lwd = 2, lty = 2)
      abline(v = 51, lwd = 2, lty = 2)
      
      points(tmp$ID, tmp$NDVI_raw, col = alpha("blue", 0.5), pch = 16)
      lines(tmp$ID, tmp$NDVI_raw, col = alpha("blue", 0.5))
      points(tmp$ID, tmp$NDVI, col = alpha("gray50", 0.5), pch = 16)
      lines(tmp$ID, tmp$NDVI, col = alpha("gray50", 0.5))
      points(tmp$ID, tmp$rescale, col = "black", pch = 16)
      lines(tmp$ID, tmp$rescale, col = "black")
      lines(spl, col = "red")
      
      abline(v = inflection_points_increasing, col = "purple", pch = 19, lwd = 2)
      #abline(v = inflection_points_decreasing, col = "purple", pch = 19, lwd = 2)
      
      points(result[result$inc==TRUE,]$x, result[result$inc==TRUE,]$y, col = "purple", pch = 16)
      points(result[result$dec==TRUE,]$x, result[result$dec==TRUE,]$y, col = "green", pch = 16)
      points(flat$x, flat$y, col = "red", pch = 16)
      
      
      legend('topleft', c("Processed", "Smoothed", "Scaled", "No Change", "Increasing", "Decreasing"),
             col = c("blue", 'gray', 'black', 'red', 'purple', 'green'), pch = 16, bty = "n")
      
      plot(irg.rast[[22]], main = i)
      plot(st_geometry(random.sample_obj[random.sample_obj$point==i,]), add = T,
           pch = 16, col = "red")
      
      dev.off()
      
      rm(tmp, second_derivatives_decreasing, second_derivatives_increasing, result, 
         inflection_points_increasing, inflection_points_decreasing, spl)
      
    }
    
  }
  
}

