library(MerkleLabGIS)
library(terra)
library(sf)
library(dplyr)
library(zoo)
library(lubridate)

folder_of_interest <- "Sublette"

## list of random points and their predicted NDVI curves
NDVI_splines_list <- readRDS(paste0("G:/Shared drives/wyo-coop-lasharr/PostDoc/CurveFitting/", folder_of_interest, "_NDVIExtracted_MultiYear.rds"))

## precip stuff
precip.files <- list.files(paste0("H:/My Drive/GEE_Precipitation_Exports/", folder_of_interest))
setwd(paste0("H:/My Drive/GEE_Precipitation_Exports/", folder_of_interest))

precip <- rast(precip.files[1:365])
precip2 <- rast(precip.files[365:729])
precip3 <- rast(precip.files[729:1093])
precip4 <- rast(precip.files[1094:1458])
precip5 <- rast(precip.files[1457:1821])

precip[[1]]
names(precip[[1:4]])

prcp2018 <- precip[[names(precip)=="prcp"]][[1:364]]
prcp2019 <- precip2[[names(precip2)=="prcp"]][[1:364]]
prcp2020 <- precip3[[names(precip3)=="prcp"]][[1:364]]
prcp2021 <- precip4[[names(precip4)=="prcp"]][[1:364]]
prcp2022 <- precip5[[names(precip5)=="prcp"]][[1:364]]

swe2018 <- precip[[names(precip)=="swe"]][[1:364]]
swe2019 <- precip2[[names(precip2)=="swe"]][[1:364]]
swe2020 <- precip3[[names(precip3)=="swe"]][[1:364]]
swe2021 <- precip4[[names(precip4)=="swe"]][[1:364]]
swe2022 <- precip5[[names(precip5)=="swe"]][[1:364]]

tmin2018 <- precip[[names(precip)=="tmin"]][[1:364]]
tmin2019 <- precip2[[names(precip2)=="tmin"]][[1:364]]
tmin2020 <- precip3[[names(precip3)=="tmin"]][[1:364]]
tmin2021 <- precip4[[names(precip4)=="tmin"]][[1:364]]
tmin2022 <- precip5[[names(precip5)=="tmin"]][[1:364]]

# elev <- rast("/vsicurl/https://pathfinder.arcc.uwyo.edu/merklelab/DEM/DEM_30m.tif")
# aspect <- rast("/vsicurl/https://pathfinder.arcc.uwyo.edu/merklelab/DEM/Aspect_30m.tif")
# slope <- rast("/vsicurl/https://pathfinder.arcc.uwyo.edu/merklelab/DEM/Slope_30m.tif")

## store the empty df
final.df <- list()

for(i in 1:length(NDVI_splines_list)){
  
  ## get the list and the data of the curves
  tmp.list <- NDVI_splines_list[[i]]
  
  ## get the df with daily estimates, df with the actual NDVI values, and the location
  loc <- tmp.list[[3]]
  yr <- year(tmp.list[[2]]$date)[1]
  
  if(is.null(loc)){
    
    print(paste0("Iteration ", i, " is empty."))
    
  }else{
    
    ## get the correct precip layers
    prcp <-if(yr == 2018){prcp2018}else if(yr == 2019){prcp2019}else if(yr == 2020){prcp2020}else if(yr == 2021){prcp2021}else{prcp2022}
    tmin <-if(yr == 2018){tmin2018}else if(yr == 2019){tmin2019}else if(yr == 2020){tmin2020}else if(yr == 2021){tmin2021}else{tmin2022}
    swe <-if(yr == 2018){swe2018}else if(yr == 2019){swe2019}else if(yr == 2020){swe2020}else if(yr == 2021){swe2021}else{swe2022}
    
    ## precip
    tmp.prcp <- t(terra::extract(prcp, st_transform(loc, crs = crs(prcp, proj=TRUE)),
                                 df=TRUE, na.rm=TRUE))
    
    if(all(is.na(tmp.prcp[-1,1]))){
    
      print(paste0("Iteration ", i, " has no precipitation data."))
      
    }else{
    
    tmp.prcp <- tmp.prcp[2:length(tmp.prcp)]
      
    ## snow water equivalent
    tmp.swe <- t(terra::extract(swe, st_transform(loc, crs = crs(swe, proj=TRUE)),
                                df=TRUE, na.rm=TRUE))
    tmp.swe <- tmp.swe[2:length(tmp.swe)]
    
    ## temperature
    tmp.tmin <- t(terra::extract(tmin, st_transform(loc, crs = crs(tmin, proj=TRUE)),
                                 df=TRUE, na.rm=TRUE))
    
    ## get the cumulative precip
    tmp.df <- data.frame(prcp = tmp.prcp, tmin = tmp.tmin[2:365,])
    tmp.df[365,] <- 0
    tmp.df$prcpwarm <- ifelse(tmp.df$tmin > 0, tmp.df$prcp, 0)
    tmp.df$cumprecp <- cumsum(tmp.df$prcpwarm)
    tmp.df <- tmp.df %>% mutate(is_decreasing = c(NA, diff(prcpwarm) < 0))
    
    cumprcp <- (tmp.df$cumprecp - min(c(tmp.df$cumprecp, tmp.swe), na.rm = TRUE)) /
      (max(c(tmp.df$cumprecp, tmp.swe), na.rm = TRUE) - min(c(tmp.df$cumprecp, tmp.swe), na.rm = TRUE))
    
    swe.plot <- (tmp.swe - min(c(tmp.df$cumprecp, tmp.swe), na.rm = TRUE)) /
      (max(c(tmp.df$cumprecp, tmp.swe), na.rm = TRUE) - min(c(tmp.df$cumprecp, tmp.swe), na.rm = TRUE))
    
    
    # elev.extract <- terra::extract(elev, st_transform(loc, crs = crs(elev, proj=TRUE)))[,2]
    # aspect.extract <- terra::extract(aspect, st_transform(loc, crs = crs(aspect, proj=TRUE)))[,2]
    # slope.extract <- terra::extract(slope, st_transform(loc, crs = crs(slope, proj=TRUE)))[,2]
    
    
    ### identify precip events ####
    
    ## precip start date
    precip_spl <- smooth.spline(tmp.df$cumprecp, spar = 0.4)
    
    inc <- c(NA, diff(precip_spl$y) > 0)
    dec <- c(NA, diff(precip_spl$y) < 0)
    flat <- c(NA, abs(diff(precip_spl$y)) < 0.005)
    
    ## create a dataframe
    precip_spl <- data.frame(x = precip_spl$x, y = precip_spl$y, inc = inc, dec = dec, flat = flat)
    precip_spl$y <- ifelse(precip_spl$y < 0, 0, precip_spl$y)
    
    precip_spl[1,"inc"] <- precip_spl[2, "inc"] ## set the increasing in point 1 as what it is in point 2
    precip_spl[1,"dec"] <- precip_spl[2, "dec"] ## same as above
    precip_spl[1,"flat"] <- precip_spl[2, "flat"] ## same as above
    
    precip_spl$inc_val <- c(0, abs(diff(precip_spl$y)))
    
    precip_spl <- precip_spl %>%
      mutate(inc_group = cumsum(!inc | (lag(inc, default = FALSE) & !inc)))
    
    precip_spl$inc_group <- ifelse(precip_spl$inc==FALSE, NA, precip_spl$inc_group)
    
    total_change <- max(precip_spl$y)
    
    absolute_threshold = 8  # change in mm considered significant
    percentage_threshold = 1  # 1% of total annual change
    
    ## get the minimum mm increased based oin the percentage value
    percentage_value = (percentage_threshold / 100) * total_change
    
    uni <- unique(precip_spl$inc_group)
    uni <- uni[!is.na(uni)]
    
    tmp.precip.events <- NULL
    
    for(k in uni){
      
      ## get each increasing group
      tmp.precip <- precip_spl[precip_spl$inc_group==k,]
      tmp.precip <- tmp.precip[!is.na(tmp.precip$y),]
      
      tmp.precip$rolling_cumulative_inc = rollapply(tmp.precip$inc_val, width = 5, FUN = sum, partial = TRUE, align = "right", fill = NA)
      
      tmp.precip$significant  <- ifelse(tmp.precip$rolling_cumulative_inc > absolute_threshold | tmp.precip$rolling_cumulative_inc > percentage_value, TRUE, FALSE)
      
      tmp.precip <- tmp.precip %>%
        mutate(inc_group_tmp = cumsum(!significant  | (lag(significant , default = FALSE) & !inc)))
      tmp.precip$inc_group_tmp <- ifelse(tmp.precip$significant==FALSE, NA, tmp.precip$inc_group_tmp)
      tmp.precip <- tmp.precip[!is.na(tmp.precip$inc_group_tmp),]
      
      periods <- tmp.precip %>%
        group_by(inc_group_tmp) %>%
        summarise(
          start_date = first(x),
          end_date = last(x),
          total_precip =  last(y) - first(y)
        )
      
      periods <- data.frame(periods)
      
      tmp.precip.events <- rbind(tmp.precip.events, periods)
      
    }
    
    ### identify snow melt events ####
    
    ## put it in a dataframe and smooth
    snow <- data.frame(raw = swe.plot)
    snow$smooth <- rollapply(swe.plot, width = 5, FUN = mean, partial = TRUE, align = "right", fill = NA)
    snow$smooth <- ifelse(snow$smooth < 0.0001, 0, snow$smooth)
    
    ## add a decreasing column
    snow <- snow %>% mutate(is_decreasing = c(NA, diff(smooth) < 0))
    
    ## add a column to see if it reached 0
    snow <- snow %>% mutate(is_zero_after_decrease = smooth == 0 & lag(smooth) > 0 & lag(is_decreasing))
    snow[1,]$is_zero_after_decrease <- FALSE ## add FALSE for the first row
    snow$jday <- 1:nrow(snow) ## get the julian dates
      
    ## get the end of snow melt
    end_melt_indices <- snow[snow$is_zero_after_decrease==TRUE,"jday"]
    
    if(length(end_melt_indices)>=1){
      
      ## get the start melt indices
      start_melt_indices <- sapply(end_melt_indices, function(end_idx) {
        ## Start searching from the end index backwards to find the first 'not decreasing' row
        rev_idx <- rev(1:end_idx)  # Reverse index to search backwards
        for (idx in rev_idx) {
          if (!snow$is_decreasing[idx] | idx == 1) {  # When it finds not decreasing or reaches the start
            return(idx + 1)  # Return the start of the decreasing period
          }
        }
      })
      
      
      tmp.snowmelt.events <- data.frame(start_date = start_melt_indices,
                                        end_date = end_melt_indices)
      
      tmp.snowmelt.events$snowmelt_amount <- snow[start_melt_indices, "smooth"] - snow[end_melt_indices, "smooth"] 
      
    }

    ## make sure that they are listed and named appropriately
    ## also deals with pixels that had no snow during a given year
    

    ## combine them together
    tmp_moisture_info <- list(precip = tmp.precip.events, snowmelt = tmp.snowmelt.events,
                              cumulative_precip = cumprcp,
                              daily_precip = tmp.df$prcpwarm,
                              swe = c(tmp.swe, tmp.swe[length(tmp.swe)]))
    
    ## add them to the larger list
    final.df[[i]] <- tmp_moisture_info
    
    print(i)
    
    }
  }
  
}


saveRDS(final.df, paste0("G:/Shared drives/wyo-coop-lasharr/PostDoc/CurveFitting/", folder_of_interest, "_multiyear_dataframe.rds"))



