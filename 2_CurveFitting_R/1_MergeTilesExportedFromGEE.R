library(terra)
library(raster)
library(sf)


## list of NDVI rasters
ndvi.list <- list.files(paste0("G:/My Drive/NDVI_Exports/"))

## get all the years
uni.years <- unique(sapply(ndvi.list, function(x) unlist(strsplit(x, "_"))[4]))
uni.years 


## combine them together
setwd(paste0("G:/My Drive/NDVI_Exports/"))

ser <- c(2,3,4,5,6,7)

for(i in ser){ #1:length(uni.years)
  #i=1
  tmp.year <- grep(uni.years[[i]], ndvi.list, value = TRUE)

  rasters <- lapply(tmp.year, rast)
  
  s <- sprc(rasters)
  mosaic_raster <- merge(s)
  
  #test only
  #plot(mosaic_raster[[10]], main = tmp.year[[1]])
  
  varnames(mosaic_raster) <- paste0("FULL_WESTERN_US_", uni.years[[i]])
  
    writeRaster(mosaic_raster, paste0("G:/Shared drives/wyo-coop-WestWideMuleDeerProject/GIS/NDVI/FULL_WESTERN_US/FULL_WESTERN_US_", uni.years[[i]], ".tif"), overwrite = TRUE)
  
}
