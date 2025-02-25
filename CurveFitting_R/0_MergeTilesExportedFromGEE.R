library(terra)
library(raster)
library(sf)

folder_of_interest <- "Paun"

## list of NDVI rasters
ndvi.list <- list.files(paste0("H:/My Drive/NDVI_Exports/", folder_of_interest))

## get all the years
uni.years <- unique(sapply(ndvi.list, function(x) { unlist(strsplit(x, "_")) [3]}))

## combine them together
setwd(paste0("H:/My Drive/NDVI_Exports/", folder_of_interest))

for(i in 1:length(uni.years)){
  
  tmp.year <- grep(uni.years[[i]], ndvi.list, value = TRUE)

  rasters <- lapply(tmp.year, rast)
  
  s <- sprc(rasters)
  mosaic_raster <- merge(s)
  plot(mosaic_raster[[10]], main = tmp.year[[1]])
  
  writeRaster(mosaic_raster, paste0('H:/My Drive/NDVI_Merged/', folder_of_interest, "/", folder_of_interest, "_", uni.years[[i]], ".tif"), overwrite = TRUE)
  
}
