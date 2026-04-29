library(terra)
library(raster)
library(sf)
library(parallel)
library(foreach)
library(doParallel)

## list of NDVI rasters
ndvi.list <- list.files(paste0("C:/Users/lwilde2/Desktop/tmp NDVI"))

## get all the years
uni.years <- unique(sapply(ndvi.list, function(x) unlist(strsplit(x, "_"))[4]))[-1]
uni.years 



# Sequential version below ####

## combine them together

setwd("C:/Users/lwilde2/Desktop/tmp NDVI")


for(i in 1:length(uni.years)){
  #i=2
  tmp.year <- grep(uni.years[[i]], ndvi.list, value = TRUE)

  if(length(tmp.year) < 9){
    # Handle case where there are less than 9 files
    message("Fewer than 9 tiles found for year ", uni.years[[i]])
    rm(tmp.year)
    next

  } else{


  rasters <- lapply(tmp.year, rast)


  message("Started production:", uni.years[[i]])
  s <- sprc(rasters)
  mosaic_raster <- merge(s)

  ## for one-time test only!
  #plot(mosaic_raster[[10]], main = tmp.year[[1]])

  varnames(mosaic_raster) <- paste0("FULL_WESTERN_US_", uni.years[[i]])

    writeRaster(mosaic_raster, paste0("G:/Shared drives/wyo-coop-WestWideMuleDeerProject/GIS/NDVI/FULL_WESTERN_US/FULL_WESTERN_US_", uni.years[[i]], ".tif"), overwrite = TRUE)

    rm(mosaic_raster, s, rasters)
    
    #file.remove(tmp.year); message("Source files removed for year ", uni.years[[i]])
    rm(tmp.year)


  }
  
}


# Parallel version ####

ndvi.list <- list.files("C:/Users/lwilde2/Desktop/tmp NDVI")
uni.years <- unique(sapply(ndvi.list, function(x) unlist(strsplit(x, "_"))[4]))[1]

local_tmp_base <- "C:/Users/lwilde2/AppData/Local/Temp/terra_ndvi"
local_out_base <- "C:/Users/lwilde2/Desktop/ndvi_merged"
dir.create(local_tmp_base, recursive = TRUE, showWarnings = FALSE)
dir.create(local_out_base, recursive = TRUE, showWarnings = FALSE)

# Conservative core count — these are large rasters
n.cores <- min(4, length(uni.years), detectCores() - 2)
cl <- makeCluster(n.cores, type = "PSOCK")
registerDoParallel(cl)

results <- foreach(
  i = 1:length(uni.years),
  .packages = c("terra", "sf"),
  .export = c("ndvi.list", "uni.years", "local_tmp_base", "local_out_base"),
  .errorhandling = "pass"  # returns the error object instead of crashing
) %dopar% {
  
  yr <- uni.years[[i]]
  yr_tmp <- file.path(local_tmp_base, yr)
  dir.create(yr_tmp, showWarnings = FALSE, recursive = TRUE)
  
  # Cleanup runs no matter what within the worker
  on.exit({
    terra::tmpFiles(remove = TRUE)
    unlink(yr_tmp, recursive = TRUE)
    gc()
  }, add = TRUE, after = TRUE)
  
  tryCatch({
    
    terraOptions(memfrac = 0.6 / n.cores, tempdir = yr_tmp)
    
    setwd("C:/Users/lwilde2/Desktop/tmp NDVI")
    tmp.year <- grep(yr, ndvi.list, value = TRUE)
    
    if(length(tmp.year) < 9){
      return(list(year = yr, status = "skipped", msg = "fewer than 9 tiles"))
    }
    
    local_out <- file.path(local_out_base, paste0("FULL_WESTERN_US_", yr, ".tif"))
    
    rasters <- lapply(tmp.year, rast)
    s <- sprc(rasters)
    mosaic_raster <- merge(s, filename = local_out, overwrite = TRUE)
    varnames(mosaic_raster) <- paste0("FULL_WESTERN_US_", yr)
    
    rm(mosaic_raster, s, rasters)
    return(list(year = yr, status = "success", local_out = local_out))
    
  }, error = function(e){
    return(list(year = yr, status = "error", msg = e$message))
  })
}

stopCluster(cl)

# Report parallel results
for(r in results){
  message(r$year, " — ", r$status, if(!is.null(r$msg)) paste0(": ", r$msg))
}


