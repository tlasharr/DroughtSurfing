library(sf)
library(terra)

## read in the home ranges and migration routes
setwd("G:/Shared drives/wyo-coop-WestWideMuleDeerProject/HumanFootprintProject/02_Outputs/Objects/HomeRanges")

summer <- readRDS("SummerKDES_95cont_250m.rds")
winter <- readRDS("WinterKDES_95cont_250m.rds")
mig <- readRDS("MigrationLocations_AllHerds_CoarseDataDropped.rds")

proj_aea <- paste0(
  "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 ",
  "+x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs ",
  "+ellps=GRS80 +towgs84=0,0,0"
)

## convert to AEA
summer <- st_transform(summer, proj_aea)
winter <- st_transform(winter, proj_aea)
mig <- st_transform(mig, proj_aea)

## add a 1000 m buffer to the migration route
mig <- st_buffer(mig, 1000)
# mig$level <- "1000m"

## get the herds to make the shapefile
herds_of_interest <- c(
  "KaibabNorth", "SanFranciscoPeaks", "18A", "RainbowValley",
  "Modoc", "Mendocino", "KernRiver", "CasaDiablo",
  "TetonRiver", "BoiseRiver", "NorthForkSalmonRiver", "AlbionMountains",
  "RockyMtnFront", "Gardiner", "CarbonCounty", "GlendiveSidney",
  "RioChama", "JemezSprings", "SanCristobal", "CrowMesa",
  "CanyonDeChelly", "RamPasture", "SouthCarrizo",        # McCrackenMesa missing
  "RubyMountains", "Washoe", "SpringMountains", "Area17Toiyabe",
  "Steens", "DevilsGarden", "Elkhorn", "KlamathBasin",
  "BookCliffs", "Paunsaugunt", "Cache", "PineValley",
  "Methow", "Klickitat", "Chelan",                       # BlueMountains missing
  "Jackson", "RedDesert", "UpperPowderRiver", "SouthWindRiver"
)

summer <- summer[,c("state", "herd", "aid.yr")]
winter <- winter[,c("state", "herd", "aid.yr")]
mig_tmp <- mig[,c("state", "herd", "aid.yr")]

summer <- summer[summer$herd %in% herds_of_interest,]
winter <- winter[winter$herd %in% herds_of_interest,]
mig_tmp <- mig_tmp[mig_tmp$herd %in% herds_of_interest,]

rm(mig)
gc()

## combine into a single sf object
all <- st_union(summer, winter)
rm(summer, winter)
gc()

all <- st_union(all, mig_tmp)
rm(mig_tmp)
gc()

## check if it works
# plot(st_geometry(all))

saveRDS(all, "G:/Shared drives/wyo-coop-WestWideMuleDeerProject/SplineSurfingValidation/Data/muledeerarea.rds")

all <- readRDS("G:/Shared drives/wyo-coop-WestWideMuleDeerProject/SplineSurfingValidation/Data/muledeerarea.rds")
