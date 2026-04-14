library(tidyverse)
library(ggplot2)
library(scales)

## folder of interest
folder_of_interest <- "Sublette"

## read in the data
Moisture_list <- readRDS(paste0("G:/Shared drives/wyo-coop-lasharr/PostDoc/CurveFitting/", folder_of_interest, "_multiyear_dataframe.rds"))
NDVI_splines_list <- readRDS(paste0("G:/Shared drives/wyo-coop-lasharr/PostDoc/CurveFitting/", folder_of_interest, "_NDVIExtracted_MultiYear.rds"))

## look at the two lists:
Moisture_list[[1]]
## [1] rainfall events including start date, end date, and total precip that fell during the event
## [2] snowmelt events, including start date, end date, and total amount of snow melt
## [3] daily cumulative precipitation
## [4] daily precipitation
## [5] daily swe

NDVI_splines_list[[1]] 
## [1] predicted daily NDVI values
## [2] raw, smoothed and rescaled NDVI values (8 day sequence)
## [3] location
## [4] max green up dates
## [5] start of green up period
## [6] end of green up period
## [7] change in NDVI during green up period
## [8] max senescence data
## [9] start of senescence period
## [10] end of senescence period
## [11] change in NDVI during senescence period
## [12] year
## [13] DLC parameters

analyses_df <- NULL

for(i in 1:length(NDVI_splines_list)){

  tmp.df <- data.frame(greenup_start = as.numeric(NDVI_splines_list[[i]][[5]]),
                       greenup_end = as.numeric(NDVI_splines_list[[i]][[6]]),
                       NDVI_change = as.numeric(NDVI_splines_list[[i]][[7]]))
  
  tmp.df <- tmp.df[complete.cases(tmp.df$greenup_start),]
  
  tmp.rain <- Moisture_list[[i]][[1]]
  tmp.snowmelt <- Moisture_list[[i]][[2]]
  
  if(nrow(tmp.df)==0){
    
    print(paste0("No NDVI data for point ", i))
    
    }else if (is.null(tmp.rain)){
    
      print(paste0("No Precip data for point ", i))
    
    }else{
  
    tmp.df$greenup_date <- as.numeric(NDVI_splines_list[[i]][[4]])

    ## now calculate the closest moisture event that is occuring
    closest.moisture.event <- data.frame(moisture_date = c(tmp.rain$start_date,
                                                           tmp.snowmelt$end_date),
                                         moisture_amt = c(tmp.rain$total_precip,
                                                          tmp.snowmelt$snowmelt_amount),
                                         type = c(rep("rain", nrow(tmp.rain)), rep("snow", nrow(tmp.snowmelt))))
    
    # Prepare expanded data frame
    expanded_df <- merge(tmp.df, closest.moisture.event, by = NULL)
    expanded_df$days_from_start <- expanded_df$moisture_date - (expanded_df$greenup_start - 7)
    expanded_df$days_from_end <- expanded_df$moisture_date - expanded_df$greenup_end
    
    # Filter events
    valid_rain_events <- subset(expanded_df, type == "rain" & days_from_start >= 0 & days_from_end <= 0)
    valid_snow_events <- subset(expanded_df, type == "snow")
    valid_events <- rbind(valid_rain_events, valid_snow_events)
    valid_events$abs_diff <- abs(valid_events$greenup_date - valid_events$moisture_date)
    
    # Find closest events
    closest_events <- valid_events[valid_events$abs_diff == ave(valid_events$abs_diff, valid_events$greenup_date, FUN=min),]
    closest_events <- closest_events[!duplicated(closest_events$greenup_date),]
    
    # Merge back
    closest_events <- closest_events[, c("greenup_date", "moisture_date", "moisture_amt", "type")]
    final_df <- merge(tmp.df, closest_events, by = "greenup_date", all.x = TRUE)
    
    long <- as.numeric(st_coordinates(NDVI_splines_list[[i]][[3]])[,1])
    lat <- as.numeric(st_coordinates(NDVI_splines_list[[i]][[3]])[,2])
    
    final_df$long <- long
    final_df$lat <- lat
    
    final_df$point_id <- i
    final_df$year <- NDVI_splines_list[[i]][[12]]
    
    dlc_predicted_values <- NDVI_splines_list[[i]][[13]]$predicted_values
    dlc_irg_values <- NDVI_splines_list[[i]][[13]]$irg
    dlc_sd <- NDVI_splines_list[[i]][[13]]$dlcparams
    
    
    if(is.null(dlc_irg_values)){
      
      final_df$dlc_greenup <- NA
      final_df$dlc_sd <- NA
      
    }else{
      
      final_df$dlc_greenup <- which.max(dlc_irg_values)
      final_df$dlc_sd <- dlc_sd[5]
      
      tmp <- NDVI_splines_list[[i]][[2]]
      tmp$jday <- yday(tmp$date)
      
      
      png(paste0("G:/Shared drives/wyo-coop-lasharr/PostDoc/CurveFitting/Tmp_Fits/", folder_of_interest, "_", i, ".png"), 800, 800)
      
      par(mfrow=c(1,1), mai=c(2,2,1.1,1.1), mgp=c(6.5,2.5,0))
      
      plot(tmp$jday, tmp$NDVI , type = "l", col = "gray50",
           ylim = c(0,1), ylab = "NDVI", xlab = "Julian Date", cex.axis = 1.5, cex.lab = 1.5, lwd = 2, main = i)
      points(tmp$jday, tmp$NDVI , pch = 16, col = "gray50")
      lines(NDVI_splines_list[[i]][[1]]$jday, NDVI_splines_list[[i]][[1]]$raw_NDVI, type = "l", lwd = 2)
      lines(NDVI_splines_list[[i]][[1]]$jday, NDVI_splines_list[[i]][[13]]$predicted_values, type = "l", col = "red", lwd = 2)
      abline(v = final_df$greenup_date, col = "cyan3", lwd = 2)
      abline(v = final_df$dlc_greenup, col = "green3", lwd = 2)
      legend("topright", c("Raw NDVI", "Spline", "DLC", "Spline Green-up", "DLC Green-up"),
              col = c("gray50", "black", "red", "cyan3", "green3"), bty= "n", cex = 1.25,
              pch = 16)
      legend("bottom", paste0("Total Change: ",  round(final_df$NDVI_change, digits = 3)), bty= "n", cex = 1.25)
      
      dev.off()
    }
    
    final_df$pop <- folder_of_interest
    analyses_df <- rbind(final_df, analyses_df)
    
  }

}


head(analyses_df)

analyses_df$sd_col <- ifelse(analyses_df$dlc_sd < 10, "blue", "gray70")

table(analyses_df$dlc_sd < 10, analyses_df$pop)

plot(analyses_df$dlc_greenup, analyses_df$greenup_date, pch = 16, 
     col = alpha(analyses_df$sd_col, 0.6),
     xlab = "Green-up (DLC)",
     ylab = "Green-up (Spline)")



plot(analyses_df[analyses_df$dlc_sd >= 10,]$dlc_greenup,
     analyses_df[analyses_df$dlc_sd >= 10,]$greenup_date, pch = 16, 
     col = alpha(analyses_df[analyses_df$dlc_sd >= 10, ]$sd_col, 0.6),
     xlab = "Green-up (DLC)",
     ylab = "Green-up (Spline)")

points(analyses_df[analyses_df$dlc_sd < 10,]$dlc_greenup,
       analyses_df[analyses_df$dlc_sd < 10,]$greenup_date, pch = 16, 
       col = alpha(analyses_df[analyses_df$dlc_sd < 10, ]$sd_col, 0.2))


tmp <- analyses_df %>% 
   group_by(point_id) %>%
    filter(greenup_start ==first(greenup_start ))


## make greenup numeric
analyses_df$greenup_date <- as.numeric(analyses_df$greenup_date)

## add in season and color for plots
analyses_df$season <- NA
analyses_df$season <- ifelse(analyses_df$greenup_date >= 90 & analyses_df$greenup_date <= 227, 'spring', NA)
analyses_df$season <- ifelse(analyses_df$greenup_date > 227, 'fall', analyses_df$season)
analyses_df$season <- dplyr::if_else(is.na(analyses_df$season), 'earlyspring', analyses_df$season)


table(analyses_df$season)

analyses_df$col <- ifelse(analyses_df$season=="earlyspring", "lightgreen", 
                          ifelse(analyses_df$season=="spring", "green4", "darkorange"))

analyses_df$col <- 'gray40'

head(analyses_df)
nrow(analyses_df)

## model some things
mod1 <- lm(greenup_date  ~ moisture_date, data = analyses_df)
summary(mod1)

######### PLOTTING #################

#### snow melt ####
par(mfrow=c(1,1), mai=c(2,2,1.1,1.1), mgp=c(6.5,2.5,0))

analyses_df_tmp <- analyses_df[analyses_df$NDVI_change > 0.025,]
analyses_df_tmp <- analyses_df_tmp[complete.cases(analyses_df_tmp$NDVI_change),]

analyses_df_tmp2 <- analyses_df[analyses_df$dlc_sd < 10,]
analyses_df_tmp2 <- analyses_df_tmp2[complete.cases(analyses_df_tmp2$dlc_sd),]
# analyses_df_tmp2$col <- ifelse(analyses_df_tmp2$pop=="Sublette", "blue2", "green3")


mod1 <- lm(greenup_date  ~ moisture_date, data = analyses_df_tmp2)
summary(mod1)

mod2 <- lm(dlc_greenup ~ moisture_date, data = analyses_df_tmp2)
summary(mod2)

plot(analyses_df$moisture_date, analyses_df$greenup_date, xlim = c(0, 375), ylim = c(0, 375),
     col = alpha(analyses_df$col, 0.5), pch = 16,
     xlab = "Moisture Event (Julian Day)",
     ylab = 'Date of Peak Green-Up (Julian Day)', cex = 1.25, cex.axis = 2, cex.lab = 2,
     bty = "l")
box(bty = "o", lwd = 2)


plot(analyses_df_tmp$moisture_date, analyses_df_tmp$greenup_date, xlim = c(0, 375), ylim = c(0, 375),
     col = alpha(analyses_df_tmp$col, 0.5), pch = 16,
     xlab = "Moisture Event (Julian Day)",
     ylab = 'Date of Peak Green-Up (Julian Day)', cex = 1.25, cex.axis = 2, cex.lab = 2,
     bty = "l")
box(bty = "o", lwd = 2)
abline(a=0, b=1, lwd = 4, lty = 2)

plot(analyses_df_tmp$moisture_date, analyses_df_tmp$dlc_greenup, xlim = c(0, 375), ylim = c(0, 375),
     col = alpha(analyses_df_tmp$col, 0.5), pch = 16,
     xlab = "Moisture Event (Julian Day)",
     ylab = 'Date of Peak Green-Up (Julian Day)', cex = 1.25, cex.axis = 2, cex.lab = 2,
     bty = "l")
box(bty = "o", lwd = 2)
abline(a=0, b=1, lwd = 4, lty = 2)


## predicted dataframe
newdat <- data.frame(moisture_date = seq(min(analyses_df_tmp$moisture_date, na.rm=T),
                                         max(analyses_df_tmp$moisture_date, na.rm=T),
                                         len=1000))
preddat <- data.frame(predict(mod1, 
                              newdata = newdat, 
                              se.fit = TRUE, interval="confidence"))

preddat$moisture_date <- seq(min(analyses_df_tmp$moisture_date, na.rm=T),
                             max(analyses_df_tmp$moisture_date, na.rm=T),
                             len=1000)

preddat$lci <- preddat$fit.fit - (preddat$se.fit*1.96)
preddat$uci <- preddat$fit.fit + (preddat$se.fit*1.96)

clip(-10, 400, -10, 400)
abline(a=0, b=1, lwd = 4, lty = 2)



## clip to the season
clip(min(analyses_df_tmp$moisture_date, na.rm = T)-10,
     max(analyses_df_tmp$moisture_date, na.rm = T)+10,
     min(analyses_df_tmp$greenup_date, na.rm = T)-10,
     max(analyses_df_tmp$greenup_date, na.rm = T)+10)

## plot the preddat
polygon(c(rev(preddat$moisture_date), preddat$moisture_date),
        c(rev(preddat$uci), preddat$lci), col = alpha('green4', 0.3), border = NA)

lines(preddat$moisture_date, preddat$fit.fit, lwd=3, col = "green4") ## this is the line


clip(-10, 400, -10, 400)
legend('topleft', c('Early Spring (before March 15)',
                     'Spring/Summer (March 15 - July 15',
                     'Autumn (after July 15)'),
       col = c('lightgreen', 'green4', 'darkorange'),
       pch = 16, bty = 'n', cex = 1.25)

