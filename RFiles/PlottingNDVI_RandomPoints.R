library(terra)
library(sf)
library(scales)

setwd("H:/My Drive/DroughtSurfing")
subfolders <- c("WY", "UT", "AZ")

## keep values
values.list <- list()

for(j in 2:length(subfolders)){
  
  ## read in raw and filtered
  bb <- st_read(paste0(subfolders[j], "/bounding_box.shp"))
  
  raw <- rast(paste0(subfolders[j], "/modis_ndvi.tif"))
  filtered <- rast(paste0(subfolders[j], "/modis_ndvi_filtered.tif"))
  
  ## read in smooth
  smooth_files <- list.files(path = paste0(getwd(), "/", subfolders[j]), pattern = "smoothed_ndvi", full.names = TRUE)
  smoothed <- rast(smooth_files)
  
  ## read in corrected1
  cor_files <- list.files(path = paste0(getwd(), "/", subfolders[j]), pattern = "corrected_modis_ndvi", full.names = TRUE)
  corrected <- rast(cor_files)
  
  nlayers <- nlyr(raw)
  
  bbox <- st_bbox(bb)  # Get the bounding box coordinates
  sample_points <- st_as_sf(data.frame(x = runif(20, bbox["xmin"], bbox["xmax"]),
                                       y = runif(20, bbox["ymin"], bbox["ymax"])),
                            coords = c("x", "y"), crs = st_crs(bb))
  
  sample_points_vect <- vect(sample_points)
  
  # Initialize lists to store the values at sample locations
  values_raw <- vector("list", nlayers)
  values_smoothed <- vector("list", nlayers)
  values_corrected <- vector("list", nlayers)
  dates <- names(raw)  # Extract layer names which contain the dates
  
  # Extract values at sample locations for each layer
  for (i in 1:nlayers) {
    values_raw[[i]] <- extract(raw[[i]], sample_points_vect, na.rm = TRUE)
    values_smoothed[[i]] <- extract(smoothed[[i]], sample_points_vect, na.rm = TRUE)
    values_corrected[[i]] <- extract(corrected[[i]], sample_points_vect, na.rm = TRUE)
    
    cat("Processed layer", dates[i], "\n")  # Print the processed date
  }
  
  # Combine extracted values into a single dataframe for each raster and date
  values_df <- data.frame(
    date = rep(as.Date(dates, format = "%Y_%m_%d_NDVI"), each = 20),  # Repeat dates for each point
    location_id = rep(1:20, nlayers),  # Identifier for each sampled location
    x = rep(st_coordinates(sample_points)[, 1], nlayers),  # X coordinates of the sample points
    y = rep(st_coordinates(sample_points)[, 2], nlayers),  # Y coordinates of the sample points
    raw_ndvi = unlist(lapply(values_raw, function(x) x[, 2])),  # NDVI values from raw
    smoothed_ndvi = unlist(lapply(values_smoothed, function(x) x[, 2])),  # NDVI values from smoothed
    corrected_ndvi = unlist(lapply(values_corrected, function(x) x[, 2]))  # NDVI values from corrected
  )
  
  values.list[[j]] <- values_df 
  
}

combined_values_df <- do.call(rbind, values.list)


## arizona drought years: 2019, 2020
titles <- c("Sublette", "Paunsaugunt", "Rainbow Valley")

for(i in 1:3){
  
  means_df <- values.list[[i]]
  
  cor_files <- list.files(path = paste0(getwd(), "/", subfolders[i]), pattern = "corrected_modis_ndvi", full.names = TRUE)
  corrected <- rast(cor_files)
  
  points <- st_as_sf(means_df, coords = c("x", "y"), crs = st_crs(bb))  # Use the CRS of your original bounding box
  
  # Plot the results
  for(j in 1:20){
    
  title <- paste0(titles[i], " point ", j)
  par(mfrow=c(1,2), mai=c(2,2,1.1,1.1), mgp=c(6.5,2.5,0))
  plot(means_df$date, means_df$smoothed_ndvi, 
       type = "n", xlab = "Date",
       ylab = "NDVI",
       xlim = c(as.Date("2018-01-01", "%Y-%m-%d"), as.Date("2021-01-01", "%Y-%m-%d")),
       main = title, lty = 2, lwd = 2,
       cex.axis = 2, cex.main = 2, cex.lab = 2, col = alpha("white", 0), ylim = c(0,1))
  
  # polygon(c(as.Date("2020-01-01"), as.Date("2020-01-01"), as.Date("2020-12-31"), as.Date("2020-12-31")),
  #          c(-1, 1, 1, -1),
  #          col = alpha("red", 0.2), border = NA)
  # 
  #  polygon(c(as.Date("2017-01-01"), as.Date("2017-01-01"), as.Date("2017-12-31"), as.Date("2017-12-31")),
  #          c(-1, 1, 1, -1),
  #          col = alpha("blue", 0.2), border = NA)
  
  abline(v = as.Date("2018-01-01", "%Y-%m-%d"), col = "gray")
  abline(v = as.Date("2019-01-01", "%Y-%m-%d"), col = "gray")
  abline(v = as.Date("2020-01-01", "%Y-%m-%d"), col = "gray")
  abline(v = as.Date("2021-01-01", "%Y-%m-%d"), col = "gray")
  abline(v = as.Date("2022-01-01", "%Y-%m-%d"), col = "gray")
  abline(v = as.Date("2023-01-01", "%Y-%m-%d"), col = "gray")
  
    
    lines(means_df[means_df$location_id==j,]$date, means_df[means_df$location_id==j,]$raw_ndvi, type = "l", col = alpha("black", 0.3), lty = 2, lwd = 2)
    lines(means_df[means_df$location_id==j,]$date, means_df[means_df$location_id==j,]$smoothed_ndvi, type = "l", col = alpha("blue", 0.3), lty = 2, lwd = 2)
    lines(means_df[means_df$location_id==j,]$date, means_df[means_df$location_id==j,]$corrected_ndvi, type = "l", col = alpha("red", 0.8), lty = 1, lwd = 2)
   
  plot(corrected[[25]])
  plot(st_geometry(points[j,]), pch = 16, col = "black", cex = 2, add = TRUE)
  
  }

}


# saveRDS(means.list, "G:/Shared drives/wyo-coop-lasharr/PostDoc/Exploration_R/meanNDVIlist.rds")
