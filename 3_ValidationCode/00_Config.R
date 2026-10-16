#### libraries ####

library(tidyverse)
library(terra)
library(sf)
library(lubridate)
library(zoo)
library(patchwork)
library(spdep)
library(ggtext)

## functions for fitting splines and the dlc - needed for scripts 03 and 04
source("C:/Users/Tayler/Documents/GitHub/DroughtSurfing/2_CurveFitting_R/Functions/Functions_FitSplines_Raster.R")
source("C:/Users/Tayler/Documents/GitHub/DroughtSurfing/2_CurveFitting_R/Functions/Functions_FitDLC.R")


#### paths ####

## merged spline + dlc metrics rasters (one per year, full western US extent)
base_path_rasters <- "G:/Shared drives/wyo-coop-WestWideMuleDeerProject/GIS/Merged_SplineSurfing"

## raw ndvi rasters, used for re-fitting splines/dlc in scripts 03 and 04
ndvi_base_path <- "G:/My Drive/NDVI_Merged/Utah"

## elevation, used in script 04 to stratify agreement by elevation quartile
dem_path <- "G:/Shared drives/wyo-coop-WestWideMuleDeerProject/GIS/Elevation/Elevation_EPSG4269_Res30m_UnitsMeters.tif"

## where all validation outputs (rds, png) get saved
base_path_saveproducts <- "G:/Shared drives/wyo-coop-WestWideMuleDeerProject/SplineSurfingValidation/Data"
base_path_savefigures <- "G:/Shared drives/wyo-coop-WestWideMuleDeerProject/SplineSurfingValidation/Figures"


#### study area ####
## CAN REPLACE THIS BOUNDING BOX WITH THE HERD POLYGONS EVENTUALLY

## bounding box spanning a north-south transect from northern Montana
## to southern Arizona, capturing MT, WY, UT, and AZ
## kept narrow (east-west) to stay a representative latitudinal slice
## rather than the full width of each state
md_poly <- st_bbox(
  c(xmin = -113.5,   # just west of wasatch front/I-15
    xmax = -109.0,   # utah's eastern border
    ymin =   36.9,   # utah-arizona border
    ymax =   42.1),  # utah-idaho border
  crs = st_crs(4326)
) %>%
  st_as_sfc() %>%
  st_as_sf()

md_vect_latlon <- vect(md_poly)   # reproject to match each raster as needed


#### shared parameters ####

## years with merged spline/dlc rasters available
years_of_interest <- 2003:2025

## rss cutoff to retain a pixel for analysis -- same threshold for both methods
rss_threshold <- 2

## julian day window defining "spring" green-up
spring_start_jday <- 60
spring_end_jday   <- 180

## baseline spline parameters (chosen after running 03_Validation_Sensitivity)
spar_baseline      <- 0.4
threshold_baseline <- 0.05

## minimum number of valid years a pixel needs for temporal consistency checks
min_valid_years <- 5


#### plotting ####

method_colors <- c("Spline" = "#2ca25f", "DLC" = "#756bb1")


#### get md_vect in raster crs from the first available raster ####
## this gets reused at the top of each validation script so the cropped
## extent is identical across all of them

ref_rast_path <- NULL
for (yr in years_of_interest) {
  p <- file.path(base_path_rasters, paste0("FULL_SplineMetrics_Raster_", yr, ".tif"))
  if (file.exists(p)) { ref_rast_path <- p; break }
}
if (is.null(ref_rast_path)) stop("no raster files found in base_path_rasters")

ref_rast <- rast(ref_rast_path)[[1]]
md_vect  <- project(md_vect_latlon, crs(ref_rast))
rm(ref_rast); gc()
