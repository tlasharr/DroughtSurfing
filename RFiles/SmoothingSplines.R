library(segmented)
library(strucchange)

ndvi_data <- readRDS("G:/Shared drives/wyo-coop-lasharr/PostDoc/Exploration_R/meanNDVIlist.rds")

par(mfrow=c(1,1))

## WY
wy_data <- ndvi_data[[1]]
wy_data$Date_num <- as.numeric(as.Date(wy_data$date))

breakpoints_model <- breakpoints(mean.corrected ~ Date_num, data = wy_data)

spline_fit <- smooth.spline(wy_data$date, wy_data$mean.corrected, spar = 0.25)
plot(wy_data$date, wy_data$mean.corrected, type = "l", main = "Breakpoints in NDVI Time Series")
lines(spline_fit, col = "red", lwd = 2)  # lwd = line width

## UT
ut_data <- ndvi_data[[2]]
ut_data$Date_num <- as.numeric(as.Date(ut_data$date))


spline_fit <- smooth.spline(ut_data$date, ut_data$mean.corrected, spar = 0.25)
plot(ut_data$date, ut_data$mean.corrected, type = "l", main = "Breakpoints in NDVI Time Series")
lines(spline_fit, col = "red", lwd = 2)  # lwd = line width


## AZ
az_data <- ndvi_data[[3]]
az_data$Date_num <- as.numeric(as.Date(az_data$date))

spline_fit <- smooth.spline(az_data$date, az_data$mean.corrected, spar = 0.25)
plot(az_data$date, az_data$mean.corrected, type = "l", main = "Breakpoints in NDVI Time Series")
lines(spline_fit, col = "red", lwd = 2)  # lwd = line width

# Example of logistic growth
logistic_model <- nls(mean.corrected ~ Asym / (1 + exp((xmid - Date_num) / scal)), 
                      data = az_data, 
                      start = list(Asym = 0.3, xmid = mean(az_data$Date_num), scal = 10))

# Add predicted logistic curve to the plot
plot(az_data$Date_num, az_data$mean.corrected, type = "l", main = "Breakpoints in NDVI Time Series")
lines(az_data$Date_num, predict(logistic_model), col = "blue", lwd = 2)

