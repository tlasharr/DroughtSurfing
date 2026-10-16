library(terra)
library(sf)
library(tidyverse)
library(lubridate)
library(furrr)


########## set up paths and parameters

base_path_rasters <- "G:/Shared drives/wyo-coop-lasharr/PostDoc/CurveFitting"
base_path_deer    <- "G:/Shared drives/wyo-coop-WestWideMuleDeerProject/HerdData"

rss_threshold <- 2   # rss threshold above which spline fits are considered poor
years_of_interest <- 2018:2022


########## load rasters into a named list ########## 

## define a helper to load rasters for a given region -- returns a named list by year
load_greenup_rasters <- function(base_path, region, years) {
  
  rast_list <- setNames(
    lapply(years, function(yr) {
      path <- file.path(base_path, region, paste0(region, "_SplineMetrics_Raster_", yr, ".rds"))
      if (file.exists(path)) rast(path) else NULL # skip years that haven't been processed yet
    }),
    as.character(years)
  )
  
  Filter(Negate(is.null), rast_list) # drop any nulls from missing files
}

utah_rasters <- load_greenup_rasters(base_path_rasters, region = "Utah",    years = years_of_interest)
wy_rasters   <- load_greenup_rasters(base_path_rasters, region = "Wyoming",  years = years_of_interest)


########## define the bands to extract by name

## using names instead of indices so this doesn't break if band order changes
bands_greenup <- c("greenup_peak_1", "greenup_peak_2",
                   "greenup_peak_3", "greenup_peak_4",
                   "rss_spline", "maxIRGdate")


########## load gps and migration data

## helper to load a population's gps and migration data from standardized file paths
load_population <- function(base_path, herd_id) {
  list(
    deer = readRDS(file.path(base_path, herd_id, "FinalCleanData",
                             paste0(herd_id, "_locs_formatted.rds"))),
    mig  = readRDS(file.path(base_path, herd_id, "FinalCleanData",
                             paste0(herd_id, "_migtime_formatted.rds")))
  )
}

## define which populations to run and which raster stack they map to
population_meta <- tribble(
  ~herd_id,                   ~rast_region,
  "UT_Paunsaugunt",           "Utah",
  "UT_LaSal",                 "Utah",
  "UT_BookCliffs",            "Utah",
  "UT_Cache",                 "Utah",
  "UT_MorganSouthRich",       "Utah",
  "UT_EagleMountain",         "Utah",
  "UT_Beaver",                "Utah",
  "UT_WasatchEast",           "Utah",
  "UT_SanJuan",               "Utah",
  "UT_PineValley",            "Utah"
  # "WY_RedDesert",             "Wyoming"
)

## load all populations into a named list
populations <- setNames(
  lapply(population_meta$herd_id, load_population, base_path = base_path_deer),
  population_meta$herd_id
)


########## core extraction function

## for a single population, extracts green-up metrics at each gps location,
## joins migration timing, and computes surfing metrics
extract_greenup_for_deer <- function(deer_sf,
                                     mig_df,
                                     rast_list,
                                     bands,
                                     rss_threshold = 2) {
  
  ## add year, julian day, and migration season to the gps data
  deer_sf <- deer_sf %>%
    mutate(year = year(dt.utc),
           jday = yday(dt.utc)) %>%
    left_join(
      mig_df %>%
        ungroup() %>%
        select(aid.yr, spring.start, spring.end, fall.start, fall.end),
      by = "aid.yr"
    ) %>%
    mutate(season = case_when(
      dt.utc >= spring.start & dt.utc <= spring.end ~ "spring",
      dt.utc >= fall.start   & dt.utc <= fall.end   ~ "fall",
      TRUE ~ NA_character_
    ))
  
  ## loop over years and extract raster values at gps locations
  map_dfr(names(rast_list), function(yr) {
    
    r   <- rast_list[[yr]]
    sub <- deer_sf %>% filter(year == as.integer(yr))
    
    if (nrow(sub) == 0) return(NULL) # skip years with no data for this population
    
    ## reproject points to match raster crs -- do this once per year, not per band
    sub_proj <- st_transform(sub, crs = crs(r))
    
    ## extract all bands in a single call rather than one call per band
    extracted <- terra::extract(r[[bands]], sub_proj)
    
    bind_cols(sub, extracted[, -1]) # drop the id column terra adds
  }) %>%
    
    ## compute surfing metrics
    mutate(
      
      ## calculate difference between each green-up peak and the deer's current julian day
      delta_peak_1 = greenup_peak_1 - jday,
      delta_peak_2 = greenup_peak_2 - jday,
      delta_peak_3 = greenup_peak_3 - jday,
      delta_peak_4 = greenup_peak_4 - jday,
      
      ## find which peak is closest in time to when the deer is there
      closest_peak_idx = apply(
        cbind(abs(delta_peak_1), abs(delta_peak_2), abs(delta_peak_3), abs(delta_peak_4)),
        1,
        function(x) if (all(is.na(x))) NA_integer_ else which.min(x)
      ),
      
      ## pull out the actual date of the closest peak and how far the deer is from it
      closest_peak_day = case_when(
        closest_peak_idx == 1 ~ greenup_peak_1,
        closest_peak_idx == 2 ~ greenup_peak_2,
        closest_peak_idx == 3 ~ greenup_peak_3,
        closest_peak_idx == 4 ~ greenup_peak_4
      ),
      days_from_peak = closest_peak_day - jday,
      
      ## set surfing metrics to na where spline fit is poor
      closest_peak_day = if_else(rss_spline > rss_threshold, NA_real_, closest_peak_day),
      days_from_peak   = if_else(rss_spline > rss_threshold, NA_real_, days_from_peak)
    )
}


########## run extraction across all populations

## pull the raster stack for each population based on its region
raster_lookup <- list(
  Utah    = utah_rasters
  # Wyoming = wy_rasters
)

## iterate sequentially over herds -- parallelism already happens inside app() during curve fitting,
## and the rasters are too large to ship across parallel workers here
all_deer_greenup <- map_dfr(
  population_meta$herd_id,
  function(herd_id) {
    
    region    <- population_meta$rast_region[population_meta$herd_id == herd_id]
    rast_list <- raster_lookup[[region]]
    
    message("extracting green-up for: ", herd_id)
    
    extract_greenup_for_deer(
      deer_sf       = populations[[herd_id]]$deer,
      mig_df        = populations[[herd_id]]$mig,
      rast_list     = rast_list,
      bands         = bands_greenup,
      rss_threshold = rss_threshold
    ) %>%
      mutate(herd = herd_id)
  }
)


########## quality checks and summaries

## check how many locations have valid green-up extractions per herd
summaries <- all_deer_greenup %>%
  group_by(herd, year, season) %>%
  summarise(
    n_total      = n(),
    n_valid      = sum(!is.na(closest_peak_day)),
    pct_valid    = round(n_valid / n_total * 100, 1),
    mean_days_from_peak = round(mean(days_from_peak, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  print(n = Inf)

summaries

########## save output

saveRDS(all_deer_greenup,
        file.path(base_path_rasters, "GreenupSurfing_AllHerds_AllYears.rds"))

message("done -- output saved to: ",
        file.path(base_path_rasters, "GreenupSurfing_AllHerds_AllYears.rds"))