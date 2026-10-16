## for pixels where both spline and dlc pass rss < 2, compare the closest
## spline green-up peak to the dlc maxIRGdate using linear regression
## this evaluates whether the two methods agree on spring green-up timing
## in pixels where both methods fit well

## run 00_Config.R first -- uses base_path_rasters, base_path_saveproducts,
## dem_path, md_vect, years_of_interest, rss_threshold, spring_start_jday,
## spring_end_jday

bands_needed <- c("greenup_peak_1", "greenup_peak_2",
                  "greenup_peak_3", "greenup_peak_4",
                  "rss_spline", "maxIRGdate", "rss_dlc")


#### step 1: extract spring green-up dates for all pixels x years ####

## only retain pixels where both methods pass rss < 2

print("extracting green-up dates across years...")

comparison_data <- map_dfr(years_of_interest, function(yr) {

  path <- file.path(base_path_rasters,
                    paste0("FULL_SplineMetrics_Raster_", yr, ".tif"))
  if (!file.exists(path)) {
    print(paste0(yr, " -- file not found, skipping"))
    return(NULL)
  }

  r      <- rast(path)[[bands_needed]]
  r_crop <- crop(r, md_vect)
  r_mask <- mask(r_crop, md_vect)
  rm(r, r_crop); gc()

  v <- values(r_mask)
  rm(r_mask); gc()

  ## column indices
  i_p1   <- 1L; i_p2 <- 2L; i_p3 <- 3L; i_p4 <- 4L
  i_rss  <- 5L; i_irg <- 6L; i_rdlc <- 7L

  ## spline: find closest peak to dlc maxIRGdate within spring window
  ## we want the spline peak that is most comparable to the dlc spring
  ## estimate, so we find the spring peak (within days 60-180) closest
  ## to maxIRGdate
  dlc_date <- v[, i_irg]
  rss_sp   <- v[, i_rss]
  rss_dlc  <- v[, i_rdlc]

  ## for each pixel, find the spline peak closest to maxIRGdate
  ## among peaks within the spring window
  p_mat <- v[, c(i_p1, i_p2, i_p3, i_p4), drop = FALSE]
  p_mat[p_mat < spring_start_jday | p_mat > spring_end_jday] <- NA

  closest_spline_peak <- map_dbl(seq_len(nrow(p_mat)), function(i) {
    peaks <- p_mat[i, ]
    peaks <- peaks[!is.na(peaks)]
    if (length(peaks) == 0) return(NA_real_)
    dlc   <- dlc_date[i]
    if (is.na(dlc)) return(NA_real_)
    ## return the peak closest to dlc date
    peaks[which.min(abs(peaks - dlc))]
  })

  tibble(
    pixel_id            = seq_len(nrow(v)),
    year                = yr,
    spline_spring_peak  = closest_spline_peak,
    dlc_spring_peak     = dlc_date,
    rss_spline          = rss_sp,
    rss_dlc             = rss_dlc
  ) %>%
    ## keep only pixels where both methods pass rss threshold
    ## and both have valid spring green-up dates
    filter(
      !is.na(rss_spline),
      !is.na(rss_dlc),
      rss_spline < rss_threshold,
      rss_dlc    < rss_threshold,
      !is.na(spline_spring_peak),
      !is.na(dlc_spring_peak),
      dlc_spring_peak >= spring_start_jday,
      dlc_spring_peak <= spring_end_jday
    )
})

print(paste0("total pixel-year records passing both rss filters: ",
             nrow(comparison_data)))
print(paste0("unique years: ", n_distinct(comparison_data$year)))

saveRDS(comparison_data,
        file.path(base_path_saveproducts,
                  "Validation4_GreenupDateComparison_Raw.rds"))


#### step 2: linear regression -- spline peak ~ dlc peak ####

print("running linear regression...")

## overall model
lm_overall <- lm(spline_spring_peak ~ dlc_spring_peak,
                 data = comparison_data)

lm_summary <- summary(lm_overall)

print("overall regression: spline peak ~ dlc peak")
print(paste0("n         = ", nrow(comparison_data)))
print(paste0("intercept = ", round(coef(lm_overall)[1], 2)))
print(paste0("slope     = ", round(coef(lm_overall)[2], 3)))
print(paste0("r2        = ", round(lm_summary$r.squared, 3)))
print(paste0("rmse      = ",
             round(sqrt(mean(residuals(lm_overall)^2)), 2), " days"))

## mean absolute difference
mad_overall <- mean(abs(comparison_data$spline_spring_peak -
                          comparison_data$dlc_spring_peak),
                    na.rm = TRUE)
print(paste0("mad       = ", round(mad_overall, 2), " days"))

## by year
lm_by_year <- comparison_data %>%
  group_by(year) %>%
  summarise(
    n            = n(),
    intercept    = round(coef(lm(spline_spring_peak ~ dlc_spring_peak))[1], 2),
    slope        = round(coef(lm(spline_spring_peak ~ dlc_spring_peak))[2], 3),
    r_squared    = round(summary(lm(spline_spring_peak ~ dlc_spring_peak))$r.squared, 3),
    rmse         = round(sqrt(mean(lm(spline_spring_peak ~ dlc_spring_peak)$residuals^2)), 2),
    mad          = round(mean(abs(spline_spring_peak - dlc_spring_peak)), 2),
    .groups = "drop"
  )

print("by year:")
print(lm_by_year, n = Inf)

saveRDS(lm_by_year,
        file.path(base_path_saveproducts,
                  "Validation4_GreenupDateComparison_ByYear.rds"))


#### step 3: plots ####

## overall scatter with regression line, downsample for plotting speed
plot_sample <- comparison_data %>%
  sample_n(min(100000, n()))

p_scatter <- ggplot(plot_sample,
                    aes(x = dlc_spring_peak,
                        y = spline_spring_peak)) +
  geom_point(alpha = 0.05, size = 0.4, color = "gray40") +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "gray20", linewidth = 0.8) +
  geom_smooth(method = "lm", se = TRUE,
              color = "#2ca25f", linewidth = 1) +
  annotate("text",
           x = spring_start_jday + 5,
           y = spring_end_jday   - 5,
           hjust = 0, vjust = 1, size = 3.5,
           label = paste0("n = ", formatC(nrow(comparison_data),
                                          format = "d", big.mark = ","),
                          "\nSlope = ", round(coef(lm_overall)[2], 3),
                          "\nR² = ",    round(lm_summary$r.squared, 3),
                          "\nRMSE = ",  round(sqrt(mean(residuals(lm_overall)^2)), 1),
                          " days\nMAD = ", round(mad_overall, 1), " days")) +
  labs(
    title    = "Spring green-up date: spline vs. DLC",
    subtitle = paste0("Pixels where both methods pass RSS < ", rss_threshold,
                      " | Dashed = 1:1 line | Green = regression fit"),
    x = "DLC spring peak (maxIRGdate; Julian day)",
    y = "Spline spring peak (Julian day)"
  ) +
  theme_bw(base_size = 12)

p_scatter

ggsave(file.path(base_path_savefigures,
                 "Validation4_GreenupDate_Scatter.png"),
       p_scatter, width = 7, height = 7, dpi = 200)

## r2 and rmse by year
p_by_year <- lm_by_year %>%
  dplyr::select(year, r_squared, rmse, mad) %>%
  pivot_longer(cols = c(r_squared, rmse, mad),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
                         "r_squared" = "R²",
                         "rmse"      = "RMSE (days)",
                         "mad"       = "MAD (days)")) %>%
  ggplot(aes(x = year, y = value, color = metric, group = metric)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c("R²"         = "#2ca25f",
                                "RMSE (days)" = "#d95f02",
                                "MAD (days)"  = "#756bb1")) +
  labs(
    title    = "Agreement between spline and DLC spring green-up dates by year",
    subtitle = paste0("Pixels passing RSS < ", rss_threshold, " in both methods"),
    x = "Year", y = "Value"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

p_by_year
ggsave(file.path(base_path_savefigures,
                 "Validation4_GreenupDate_ByYear.png"),
       p_by_year, width = 7, height = 9, dpi = 200)

## residual distribution
p_resid <- ggplot(comparison_data,
                  aes(x = spline_spring_peak - dlc_spring_peak)) +
  geom_histogram(binwidth = 2, fill = "#2ca25f",
                 color = "white", alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "gray20", linewidth = 0.8) +
  labs(
    title    = "Residuals: spline peak - DLC peak",
    subtitle = paste0("Positive = spline estimates later green-up than DLC\n",
                      "Mean = ", round(mean(comparison_data$spline_spring_peak -
                                              comparison_data$dlc_spring_peak,
                                            na.rm = TRUE), 1), " days"),
    x = "Difference (Julian days)",
    y = "Count"
  ) +
  theme_bw(base_size = 12)

p_resid
ggsave(file.path(base_path_savefigures,
                 "Validation4_GreenupDate_ResidualDist.png"),
       p_resid, width = 7, height = 5, dpi = 200)

saveRDS(list(lm_overall = lm_overall,
             lm_summary = lm_summary),
        file.path(base_path_saveproducts,
                  "Validation4_GreenupDate_LMResults.rds"))


#### step 4: stratify agreement by elevation and max ndvi productivity ####

print("extracting elevation and max ndvi for stratification...")

## load reference raster to get target extent, crs, resolution
ref_r    <- rast(file.path(base_path_rasters,
                           paste0("FULL_SplineMetrics_Raster_",
                                  years_of_interest[1], ".tif")))[[1]]
ref_crop <- mask(crop(ref_r, md_vect), md_vect)

## elevation -- reproject and resample to match spline raster
dem_raw  <- rast(dem_path)
dem_proj <- project(dem_raw, ref_crop, method = "bilinear")
dem_crop <- mask(crop(dem_proj, md_vect), md_vect)
elev_vals <- values(dem_crop)[, 1]
rm(dem_raw, dem_proj); gc()

print(paste0("elevation range: ", round(min(elev_vals, na.rm = TRUE)),
             " - ", round(max(elev_vals, na.rm = TRUE)), " m"))

## max ndvi -- average across years for a stable productivity index
print("computing mean max ndvi across years...")

maxndvi_stack <- map(years_of_interest, function(yr) {
  path <- file.path(base_path_rasters,
                    paste0("FULL_SplineMetrics_Raster_", yr, ".tif"))
  if (!file.exists(path)) return(NULL)
  r      <- rast(path)[["maxNDVI"]]
  r_crop <- mask(crop(r, md_vect), md_vect)
  values(r_crop)[, 1]
}) %>%
  Filter(Negate(is.null), .) %>%
  do.call(cbind, .)

mean_maxndvi <- rowMeans(maxndvi_stack, na.rm = TRUE)
mean_maxndvi[is.nan(mean_maxndvi)] <- NA
rm(maxndvi_stack); gc()

print(paste0("max ndvi range: ", round(min(mean_maxndvi, na.rm = TRUE), 3),
             " - ", round(max(mean_maxndvi, na.rm = TRUE), 3)))

## pack into a data frame keyed by pixel_id (row index of cropped raster)
pixel_strata <- tibble(
  pixel_id  = seq_len(ncell(ref_crop)),
  elevation = elev_vals,
  max_ndvi  = mean_maxndvi
) %>%
  filter(!is.na(elevation) | !is.na(max_ndvi))

rm(ref_r, ref_crop, dem_crop, elev_vals, mean_maxndvi); gc()


#### step 4b: join strata to comparison_data and assign quartiles ####

comparison_strat <- comparison_data %>%
  left_join(pixel_strata, by = "pixel_id") %>%
  filter(!is.na(elevation), !is.na(max_ndvi)) %>%
  mutate(
    elev_quartile   = ntile(elevation, 4),
    ndvi_quartile   = ntile(max_ndvi,  4),
    elev_label      = case_when(
      elev_quartile == 1 ~ "Low elev (Q1)",
      elev_quartile == 2 ~ "Mid-low elev (Q2)",
      elev_quartile == 3 ~ "Mid-high elev (Q3)",
      elev_quartile == 4 ~ "High elev (Q4)"
    ),
    ndvi_label      = case_when(
      ndvi_quartile == 1 ~ "Low productivity (Q1)",
      ndvi_quartile == 2 ~ "Mid-low productivity (Q2)",
      ndvi_quartile == 3 ~ "Mid-high productivity (Q3)",
      ndvi_quartile == 4 ~ "High productivity (Q4)"
    ),
    elev_label = factor(elev_label,
                        levels = c("Low elev (Q1)", "Mid-low elev (Q2)",
                                   "Mid-high elev (Q3)", "High elev (Q4)")),
    ndvi_label = factor(ndvi_label,
                        levels = c("Low productivity (Q1)",
                                   "Mid-low productivity (Q2)",
                                   "Mid-high productivity (Q3)",
                                   "High productivity (Q4)"))
  )

print(paste0("records after joining strata: ", nrow(comparison_strat)))


#### step 4c: regression statistics by stratum ####

compute_lm_stats <- function(df) {
  if (nrow(df) < 10) return(NULL)
  fit <- lm(spline_spring_peak ~ dlc_spring_peak, data = df)
  sm  <- summary(fit)
  tibble(
    n         = nrow(df),
    intercept = round(coef(fit)[1], 2),
    slope     = round(coef(fit)[2], 3),
    r_squared = round(sm$r.squared, 3),
    rmse      = round(sqrt(mean(residuals(fit)^2)), 2),
    mad       = round(mean(abs(df$spline_spring_peak -
                                 df$dlc_spring_peak)), 2),
    mean_diff = round(mean(df$spline_spring_peak -
                             df$dlc_spring_peak), 2)
  )
}

## by elevation quartile
lm_by_elev <- comparison_strat %>%
  group_by(elev_quartile, elev_label) %>%
  group_modify(~ compute_lm_stats(.x)) %>%
  ungroup()

print("agreement by elevation quartile:")
print(lm_by_elev)

## by max ndvi quartile
lm_by_ndvi <- comparison_strat %>%
  group_by(ndvi_quartile, ndvi_label) %>%
  group_modify(~ compute_lm_stats(.x)) %>%
  ungroup()

print("agreement by max ndvi quartile:")
print(lm_by_ndvi)

saveRDS(lm_by_elev,
        file.path(base_path_saveproducts,
                  "Validation4_GreenupDate_ByElevation.rds"))
saveRDS(lm_by_ndvi,
        file.path(base_path_saveproducts,
                  "Validation4_GreenupDate_ByNDVI.rds"))


#### step 4d: scatter plots by elevation and ndvi quartile ####

## downsample for plotting
plot_strat <- comparison_strat %>%
  group_by(elev_quartile) %>%
  sample_n(min(25000, n())) %>%
  ungroup()

p_elev_scatter <- ggplot(plot_strat,
                         aes(x = dlc_spring_peak,
                             y = spline_spring_peak)) +
  geom_point(alpha = 0.05, size = 0.4, color = "gray40") +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "gray20", linewidth = 0.7) +
  geom_smooth(method = "lm", se = TRUE,
              color = "#2ca25f", linewidth = 1) +
  geom_text(data = lm_by_elev,
            aes(x = 65, y = 175,
                label = paste0("R² = ", r_squared,
                               "\nRMSE = ", rmse, " d",
                               "\nMAD = ", mad, " d")),
            hjust = 0, vjust = 1, size = 3, inherit.aes = FALSE) +
  facet_wrap(~ elev_label, ncol = 2) +
  labs(
    title    = "Spring green-up agreement by elevation quartile",
    subtitle = "Spline peak vs. DLC maxIRGdate | Both methods RSS < 2",
    x = "DLC spring peak (Julian day)",
    y = "Spline spring peak (Julian day)"
  ) +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(size = 9))

ggsave(file.path(base_path_savefigures,
                 "Validation4_GreenupDate_ByElevation_Scatter.png"),
       p_elev_scatter, width = 9, height = 8, dpi = 200)

## scatter by max ndvi quartile
plot_strat_ndvi <- comparison_strat %>%
  group_by(ndvi_quartile) %>%
  sample_n(min(25000, n())) %>%
  ungroup()

p_ndvi_scatter <- ggplot(plot_strat_ndvi,
                         aes(x = dlc_spring_peak,
                             y = spline_spring_peak)) +
  geom_point(alpha = 0.05, size = 0.4, color = "gray40") +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "gray20", linewidth = 0.7) +
  geom_smooth(method = "lm", se = TRUE,
              color = "#2ca25f", linewidth = 1) +
  geom_text(data = lm_by_ndvi,
            aes(x = 65, y = 175,
                label = paste0("R² = ", r_squared,
                               "\nRMSE = ", rmse, " d",
                               "\nMAD = ", mad, " d")),
            hjust = 0, vjust = 1, size = 3, inherit.aes = FALSE) +
  facet_wrap(~ ndvi_label, ncol = 2) +
  labs(
    title    = "Spring green-up agreement by productivity quartile (max NDVI)",
    subtitle = "Spline peak vs. DLC maxIRGdate | Both methods RSS < 2",
    x = "DLC spring peak (Julian day)",
    y = "Spline spring peak (Julian day)"
  ) +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(size = 9))

ggsave(file.path(base_path_savefigures,
                 "Validation4_GreenupDate_ByNDVI_Scatter.png"),
       p_ndvi_scatter, width = 9, height = 8, dpi = 200)

## summary bar chart: r2 and mad by quartile for both strata
summary_long <- bind_rows(
  lm_by_elev %>%
    dplyr::select(label = elev_label, r_squared, mad, rmse) %>%
    mutate(stratum = "Elevation"),
  lm_by_ndvi %>%
    dplyr::select(label = ndvi_label, r_squared, mad, rmse) %>%
    mutate(stratum = "Max NDVI (productivity)")
) %>%
  pivot_longer(cols = c(r_squared, mad, rmse),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
                         "r_squared" = "R²",
                         "mad"       = "MAD (days)",
                         "rmse"      = "RMSE (days)"))

p_summary_bars <- ggplot(summary_long,
                         aes(x = label, y = value, fill = stratum)) +
  geom_col(position = "dodge", alpha = 0.85) +
  facet_grid(metric ~ stratum, scales = "free") +
  scale_fill_manual(values = c("Elevation"              = "#2c7bb6",
                               "Max NDVI (productivity)" = "#2ca25f")) +
  labs(
    title    = "Agreement between spline and DLC by elevation and productivity",
    subtitle = "Both methods RSS < 2",
    x = NULL, y = "Value"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x  = element_text(angle = 30, hjust = 1, size = 8),
    legend.position = "none"
  )

ggsave(file.path(base_path_savefigures,
                 "Validation4_GreenupDate_StrataSummary.png"),
       p_summary_bars, width = 10, height = 8, dpi = 200)

print("step 4 complete")
print("validation 4 (green-up date comparison) complete")
