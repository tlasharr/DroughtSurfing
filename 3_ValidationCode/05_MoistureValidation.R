# =============================================================================
# 05_Validation_MoistureGreenup.R
#
# evaluate how moisture events (snowmelt and precipitation) relate to
# green-up events identified by the spline and double logistic methods.
# tests whether the spline method identifies green-up events that are
# more tightly associated with moisture events than the dlc.
# =============================================================================

## run 00_Config.R first

library(tidyverse)
library(terra)
library(sf)
library(lubridate)
library(zoo)
library(patchwork)

daymet_path           <- "G:/My drive/DAYMET_merged"
base_path_savefigures <- "G:/Shared drives/wyo-coop-WestWideMuleDeerProject/SplineSurfingValidation/Figures"
n_sample              <- 3000

#### step 1: helper functions to identify moisture events ####

## identify snowmelt events from swe time series
find_snowmelt_events <- function(swe_vec, jdays, min_melt_mm = 10) {
  
  if (all(is.na(swe_vec))) return(NULL)
  
  ## smooth swe with rolling mean to reduce daily variability
  swe_smooth <- rollmean(swe_vec, k = 5, fill = NA, align = "right")
  
  events     <- list()
  in_melt    <- FALSE
  melt_start <- NA
  
  for (i in seq_along(swe_smooth)) {
    if (is.na(swe_smooth[i])) next
    if (!in_melt && swe_smooth[i] > 0) {
      in_melt    <- TRUE
      melt_start <- i
    }
    if (in_melt && swe_smooth[i] <= 0) {
      total_melt <- swe_vec[melt_start] - swe_vec[i]
      if (!is.na(total_melt) && total_melt >= min_melt_mm) {
        events[[length(events) + 1]] <- tibble(
          event_type     = "snowmelt",
          start_jday     = jdays[melt_start],
          end_jday       = jdays[i],
          total_moisture = total_melt
        )
      }
      in_melt    <- FALSE
      melt_start <- NA
    }
  }
  bind_rows(events)
}

## identify rainfall events from cumulative precip time series
find_rainfall_events <- function(prcp_vec, tmin_vec, jdays,
                                 min_precip_pct = 0.01,
                                 min_precip_mm  = 10) {
  
  if (all(is.na(prcp_vec))) return(NULL)
  
  ## only count precip on days above freezing as rain
  rain_vec <- ifelse(!is.na(tmin_vec) & tmin_vec > 0, prcp_vec, 0)
  rain_vec[is.na(rain_vec)] <- 0
  
  cum_rain   <- cumsum(rain_vec)
  total_rain <- max(cum_rain, na.rm = TRUE)
  if (total_rain == 0) return(NULL)
  
  valid  <- !is.na(cum_rain)
  sp_fit <- tryCatch(
    smooth.spline(jdays[valid], cum_rain[valid], spar = 0.4),
    error = function(e) NULL
  )
  if (is.null(sp_fit)) return(NULL)
  
  pred_y  <- predict(sp_fit, jdays)$y
  delta_y <- diff(pred_y)
  
  rle_dir    <- rle(sign(delta_y))
  seg_ends   <- cumsum(rle_dir$lengths)
  seg_starts <- c(1, head(seg_ends, -1) + 1)
  
  events <- list()
  for (i in seq_along(rle_dir$values)) {
    if (is.na(rle_dir$values[i]) || rle_dir$values[i] != 1) next
    idx          <- seg_starts[i]:seg_ends[i]
    total_precip <- pred_y[max(idx)] - pred_y[idx[1]]
    if (total_precip < max(min_precip_mm, total_rain * min_precip_pct)) next
    events[[length(events) + 1]] <- tibble(
      event_type     = "rain",
      start_jday     = jdays[idx[1]],
      end_jday       = jdays[max(idx)],
      total_moisture = total_precip
    )
  }
  bind_rows(events)
}

## match a green-up event to the nearest moisture event
## rain: must start before greenup_start - 7 and end before greenup_end
## snow: use end date of melt with no temporal constraint
## returns lag (abs diff to peak), moisture date, event type, and match flag
get_moisture_match <- function(greenup_peak, greenup_start, greenup_end,
                               moisture_events) {
  
  empty <- tibble(lag           = NA_real_,
                  moisture_jday = NA_real_,
                  event_type    = NA_character_,
                  matched       = FALSE)
  
  if (is.null(moisture_events) || nrow(moisture_events) == 0) return(empty)
  if (is.na(greenup_peak) || is.na(greenup_start) || is.na(greenup_end)) return(empty)
  
  ## rain must start at or before greenup_start - 7 and end before greenup_end
  rain_events <- moisture_events %>%
    filter(event_type == "rain",
           start_jday <= (greenup_start - 7),
           end_jday   <= greenup_end) %>%
    mutate(ref_jday = start_jday)
  
  ## snowmelt: use end date, no temporal constraint
  snow_events <- moisture_events %>%
    filter(event_type == "snowmelt") %>%
    mutate(ref_jday = end_jday)
  
  candidates <- bind_rows(rain_events, snow_events) %>%
    mutate(abs_diff = abs(greenup_peak - ref_jday)) %>%
    arrange(abs_diff)
  
  if (nrow(candidates) == 0) return(empty)
  
  best <- candidates %>% slice(1)
  tibble(lag           = best$abs_diff,
         moisture_jday = best$ref_jday,
         event_type    = best$event_type,
         matched       = TRUE)
}


#### step 2: sample pixels by coordinates ####

print("loading data and sampling pixels...")

ref_r    <- rast(ref_rast_path)[[1]]
ref_c    <- mask(crop(ref_r, md_vect), md_vect)
rss_vals <- values(ref_c)[, 1]
eligible <- which(!is.na(rss_vals))

set.seed(42)
sampled_cells <- sample(eligible, min(n_sample, length(eligible)))
sampled_xy    <- xyFromCell(ref_c, sampled_cells)
rm(ref_r, ref_c, rss_vals); gc()

print(paste0("sampled ", nrow(sampled_xy), " pixels"))


#### step 3: loop over years -- extract moisture events and green-up dates ####

print("extracting moisture events and green-up dates...")

moisture_greenup <- map_dfr(years_of_interest, function(yr) {
  
  daymet_file  <- file.path(daymet_path, paste0("DAYMET_", yr, ".tif"))
  metrics_file <- file.path(base_path_rasters,
                            paste0("FULL_SplineMetrics_Raster_", yr, ".tif"))
  
  if (!file.exists(daymet_file) || !file.exists(metrics_file)) {
    print(paste0(yr, " -- skipping, file not found"))
    return(NULL)
  }
  
  ## load and crop daymet
  r_day <- rast(daymet_file)
  r_dc  <- mask(crop(r_day, md_vect), md_vect)
  rm(r_day); gc()
  
  ## parse band names -- format is 20180101_prcp_2018_01_01
  bnames    <- names(r_dc)
  variables <- sub("^\\d{8}_([a-z]+)_.*$", "\\1", bnames)
  dates_str <- sub("^.*_(\\d{4}_\\d{2}_\\d{2})$", "\\1", bnames)
  jdays_all <- yday(as.Date(dates_str, "%Y_%m_%d"))
  
  vals   <- terra::extract(r_dc, sampled_xy)
  has_id <- names(vals)[1] == "ID"
  rm(r_dc); gc()
  
  ## load spline metrics -- all four peaks with start/end windows
  bands_gu <- c("greenup_peak_1",  "greenup_peak_2",
                "greenup_peak_3",  "greenup_peak_4",
                "greenup_start_1", "greenup_start_2",
                "greenup_start_3", "greenup_start_4",
                "greenup_end_1",   "greenup_end_2",
                "greenup_end_3",   "greenup_end_4",
                "rss_spline",      "maxIRGdate",
                "springStart",     "springEnd",
                "rss_dlc")
  
  r_met   <- rast(metrics_file)[[bands_gu]]
  r_metc  <- mask(crop(r_met, md_vect), md_vect)
  gu_vals <- terra::extract(r_metc, sampled_xy)
  rm(r_met, r_metc); gc()
  
  print(paste0(yr, " -- extracted, processing pixels..."))
  
  map_dfr(seq_len(nrow(sampled_xy)), function(i) {
    
    v    <- if (has_id) as.numeric(vals[i, -1]) else as.numeric(vals[i, ])
    prcp <- v[variables == "prcp"]
    swe  <- v[variables == "swe"]
    tmin <- v[variables == "tmin"]
    jd   <- jdays_all[variables == "prcp"]
    
    if (length(prcp) < 10 || all(is.na(prcp))) return(NULL)
    
    ## identify moisture events
    snow_events <- find_snowmelt_events(swe, jd)
    rain_events <- find_rainfall_events(prcp, tmin, jd)
    all_events  <- bind_rows(snow_events, rain_events)
    
    ## green-up metrics -- column order matches bands_gu above
    g          <- as.numeric(gu_vals[i, ])
    sp_peak_1  <- g[1];  sp_peak_2  <- g[2];  sp_peak_3  <- g[3];  sp_peak_4  <- g[4]
    gu_start_1 <- g[5];  gu_start_2 <- g[6];  gu_start_3 <- g[7];  gu_start_4 <- g[8]
    gu_end_1   <- g[9];  gu_end_2   <- g[10]; gu_end_3   <- g[11]; gu_end_4   <- g[12]
    rss_sp     <- g[13]
    irg_dat    <- g[14]
    dlc_start  <- g[15]; dlc_end    <- g[16]
    rss_dlc    <- g[17]
    
    ## spline: use first spring peak with its start/end window
    sp_peaks       <- c(sp_peak_1, sp_peak_2, sp_peak_3, sp_peak_4)
    sp_peaks_valid <- sp_peaks[!is.na(sp_peaks) 
                               # &
                               #   sp_peaks >= spring_start_jday &
                               #   sp_peaks <= spring_end_jday
                               ]
    
    spline_match <- if (length(sp_peaks_valid) > 0 && !is.na(rss_sp)) {
      get_moisture_match(min(sp_peaks_valid), gu_start_1, gu_end_1, all_events)
    } else tibble(lag = NA_real_, moisture_jday = NA_real_,
                  event_type = NA_character_, matched = FALSE)
    
    ## dlc: use springStart and springEnd as the green-up window
    dlc_match <- if (!is.na(irg_dat) && !is.na(rss_dlc) 
                     # &&
                     # irg_dat >= spring_start_jday &&
                     # irg_dat <= spring_end_jday
                     ) {
      get_moisture_match(irg_dat, dlc_start, dlc_end, all_events)
    } else tibble(lag = NA_real_, moisture_jday = NA_real_,
                  event_type = NA_character_, matched = FALSE)
    
    tibble(
      pixel_id             = i,
      year                 = yr,
      ## all four spline peaks with start/end windows
      spline_peak_1        = as.numeric(sp_peak_1),
      spline_peak_2        = as.numeric(sp_peak_2),
      spline_peak_3        = as.numeric(sp_peak_3),
      spline_peak_4        = as.numeric(sp_peak_4),
      spline_start_1       = as.numeric(gu_start_1),
      spline_start_2       = as.numeric(gu_start_2),
      spline_start_3       = as.numeric(gu_start_3),
      spline_start_4       = as.numeric(gu_start_4),
      spline_end_1         = as.numeric(gu_end_1),
      spline_end_2         = as.numeric(gu_end_2),
      spline_end_3         = as.numeric(gu_end_3),
      spline_end_4         = as.numeric(gu_end_4),
      ## convenience column: earliest spring peak
      spline_peak_jday     = if (length(sp_peaks_valid) > 0) min(sp_peaks_valid) else NA_real_,
      spline_start_jday    = as.numeric(gu_start_1),
      spline_end_jday      = as.numeric(gu_end_1),
      ## dlc
      dlc_peak_jday        = as.numeric(irg_dat),
      dlc_start_jday       = as.numeric(dlc_start),
      dlc_end_jday         = as.numeric(dlc_end),
      ## moisture matching results
      spline_lag_days      = spline_match$lag,
      dlc_lag_days         = dlc_match$lag,
      spline_moisture_jday = spline_match$moisture_jday,
      dlc_moisture_jday    = dlc_match$moisture_jday,
      spline_event_type    = spline_match$event_type,
      dlc_event_type       = dlc_match$event_type,
      spline_matched       = spline_match$matched,
      dlc_matched          = dlc_match$matched,
      rss_spline           = as.numeric(rss_sp),
      rss_dlc              = as.numeric(rss_dlc),
      n_rain_events        = if (!is.null(rain_events)) nrow(rain_events) else 0L,
      n_snow_events        = if (!is.null(snow_events)) nrow(snow_events) else 0L
    )
  })
})

print(paste0("done. rows: ", nrow(moisture_greenup)))
head(moisture_greenup)

saveRDS(moisture_greenup,
        file.path(base_path_saveproducts, "Validation5_MoistureGreenup.rds"))


#### step 4: compare lag distributions between methods ####

lag_long <- moisture_greenup %>%
  filter(!is.na(spline_lag_days) | !is.na(dlc_lag_days)) %>%
  pivot_longer(cols = c(spline_lag_days, dlc_lag_days),
               names_to = "method", values_to = "lag_days") %>%
  mutate(method = case_match(method,
                             "spline_lag_days" ~ "Spline",
                             "dlc_lag_days"    ~ "DLC")) %>%
  filter(!is.na(lag_days))

lag_summary <- lag_long %>%
  group_by(method) %>%
  summarise(
    n              = n(),
    mean_lag       = round(mean(lag_days,   na.rm = TRUE), 1),
    median_lag     = round(median(lag_days, na.rm = TRUE), 1),
    sd_lag         = round(sd(lag_days,     na.rm = TRUE), 1),
    pct_within_7d  = round(mean(lag_days <= 7,  na.rm = TRUE) * 100, 1),
    pct_within_14d = round(mean(lag_days <= 14, na.rm = TRUE) * 100, 1),
    pct_matched    = round(mean(!is.na(lag_days)) * 100, 1),
    .groups = "drop"
  )

print("lag to nearest moisture event by method:")
print(lag_summary)

## paired wilcoxon on pixels where both methods have a lag
paired <- moisture_greenup %>%
  filter(!is.na(spline_lag_days), !is.na(dlc_lag_days))

if (nrow(paired) > 10) {
  wt <- wilcox.test(paired$spline_lag_days, paired$dlc_lag_days, paired = TRUE)
  print(paste0("wilcoxon: W = ", wt$statistic, "  p = ", round(wt$p.value, 4)))
  print(paste0("median difference (spline - dlc): ",
               round(median(paired$spline_lag_days - paired$dlc_lag_days), 1),
               " days"))
}

saveRDS(lag_summary,
        file.path(base_path_saveproducts, "Validation5_MoistureLag_Summary.rds"))


#### step 4b: how many green-up events occur without a moisture match? ####

match_rate <- moisture_greenup %>%
  mutate(
    has_spline_greenup = !is.na(spline_peak_jday) & rss_spline < rss_threshold,
    has_dlc_greenup    = !is.na(dlc_peak_jday)    & rss_dlc    < rss_threshold
  ) %>%
  summarise(
    n_spline_greenup   = sum(has_spline_greenup,  na.rm = TRUE),
    n_spline_matched   = sum(has_spline_greenup & spline_matched, na.rm = TRUE),
    pct_spline_matched = round(mean(spline_matched[has_spline_greenup],
                                    na.rm = TRUE) * 100, 1),
    n_dlc_greenup      = sum(has_dlc_greenup,     na.rm = TRUE),
    n_dlc_matched      = sum(has_dlc_greenup & dlc_matched, na.rm = TRUE),
    pct_dlc_matched    = round(mean(dlc_matched[has_dlc_greenup],
                                    na.rm = TRUE) * 100, 1)
  )

print("proportion of green-up events with a detectable moisture driver:")
print(match_rate)

## event type breakdown for matched spline green-ups
event_type_summary <- moisture_greenup %>%
  filter(!is.na(spline_peak_jday), rss_spline < rss_threshold, spline_matched) %>%
  count(spline_event_type) %>%
  mutate(pct = round(n / sum(n) * 100, 1))

print("spline green-up event type breakdown:")
print(event_type_summary)

## lag distribution summary
lag_dist <- moisture_greenup %>%
  filter(!is.na(spline_peak_jday), rss_spline < rss_threshold, spline_matched) %>%
  summarise(
    median_lag     = round(median(spline_lag_days, na.rm = TRUE), 1),
    pct_lag_lt_14d = round(mean(spline_lag_days <= 14, na.rm = TRUE) * 100, 1),
    pct_lag_lt_30d = round(mean(spline_lag_days <= 30, na.rm = TRUE) * 100, 1),
    pct_lag_gt_30d = round(mean(spline_lag_days >  30, na.rm = TRUE) * 100, 1)
  )

print("spline green-up lag distribution:")
print(lag_dist)


#### step 5: plots ####

## lag distribution comparison
p_lag <- ggplot(lag_long, aes(x = lag_days, fill = method, color = method)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
  scale_fill_manual(values  = method_colors) +
  scale_color_manual(values = method_colors) +
  labs(
    title    = "Lag from moisture event to spring green-up peak",
    subtitle = paste0("Spline: n = ", lag_summary$n[lag_summary$method == "Spline"],
                      "  |  DLC: n = ", lag_summary$n[lag_summary$method == "DLC"]),
    x     = "Lag (days; absolute difference between moisture event and green-up peak)",
    y     = "Density",
    fill  = "Method",
    color = "Method"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(base_path_savefigures, "Validation5_MoistureLag_Distribution.png"),
       p_lag, width = 8, height = 5, dpi = 200)

## proportion matched by method
p_match_rate <- lag_summary %>%
  dplyr::select(method, pct_matched, pct_within_7d, pct_within_14d) %>%
  pivot_longer(cols = c(pct_matched, pct_within_7d, pct_within_14d),
               names_to = "metric", values_to = "pct") %>%
  mutate(metric = case_match(metric,
                             "pct_matched"    ~ "% with moisture match",
                             "pct_within_7d"  ~ "% within 7 days",
                             "pct_within_14d" ~ "% within 14 days")) %>%
  ggplot(aes(x = metric, y = pct, fill = method)) +
  geom_col(position = "dodge", alpha = 0.85) +
  geom_text(aes(label = paste0(pct, "%")),
            position = position_dodge(width = 0.9),
            vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = method_colors) +
  labs(
    title = "Green-up association with moisture events: spline vs. DLC",
    x = NULL, y = "Percent of green-up events (%)",
    fill = "Method"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(base_path_savefigures, "Validation5_MoistureMatch_Rate.png"),
       p_match_rate, width = 7, height = 5, dpi = 200)

## lag by year
p_lag_year <- moisture_greenup %>%
  filter(!is.na(spline_lag_days) | !is.na(dlc_lag_days)) %>%
  pivot_longer(cols = c(spline_lag_days, dlc_lag_days),
               names_to = "method", values_to = "lag_days") %>%
  mutate(method = case_match(method,
                             "spline_lag_days" ~ "Spline",
                             "dlc_lag_days"    ~ "DLC")) %>%
  filter(!is.na(lag_days)) %>%
  group_by(year, method) %>%
  summarise(median_lag = median(lag_days, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = year, y = median_lag, color = method, group = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_color_manual(values = method_colors) +
  labs(
    title    = "Median lag to nearest moisture event by year",
    subtitle = "Lower = green-up more tightly coupled to moisture",
    x = "Year", y = "Median lag (days)",
    color = "Method"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(base_path_savefigures, "Validation5_MoistureLag_ByYear.png"),
       p_lag_year, width = 8, height = 5, dpi = 200)

## scatter: green-up date vs moisture event date
scatter_long <- bind_rows(
  moisture_greenup %>%
    filter(!is.na(dlc_peak_jday), !is.na(dlc_moisture_jday)) %>%
    transmute(greenup_jday  = dlc_peak_jday,
              moisture_jday = dlc_moisture_jday,
              method        = "DLC (A)"),
  moisture_greenup %>%
    filter(!is.na(spline_peak_jday), !is.na(spline_moisture_jday)) %>%
    transmute(greenup_jday  = spline_peak_jday,
              moisture_jday = spline_moisture_jday,
              method        = "Spline (B)")
) %>%
  mutate(method = factor(method, levels = c("DLC (A)", "Spline (B)")))

lm_stats <- scatter_long %>%
  group_by(method) %>%
  summarise(
    r2    = round(cor(moisture_jday, greenup_jday, use = "complete.obs")^2, 2),
    slope = round(coef(lm(greenup_jday ~ moisture_jday))[2], 2),
    .groups = "drop"
  )

p_scatter_moisture <- ggplot(scatter_long,
                             aes(x = moisture_jday, y = greenup_jday)) +
  geom_point(alpha = 0.15, size = 0.6, color = "gray30") +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_smooth(method = "lm", se = FALSE,
              color = "#2ca25f", linewidth = 1) +
  geom_text(data = lm_stats,
            aes(x = 25, y = 330,
                label = paste0("R² = ", r2, "\nβ = ", slope)),
            hjust = 0, vjust = 1, size = 4, inherit.aes = FALSE) +
  facet_wrap(~ method, ncol = 2) +
  labs(
    title    = "Green-up date vs. nearest preceding moisture event",
    subtitle = "Dashed = 1:1 line",
    x = "Moisture event (Julian day)",
    y = "Green-up date (Julian day)"
  ) +
  theme_bw(base_size = 12) +
  theme(strip.text = element_text(size = 11))

ggsave(file.path(base_path_savefigures,
                 "Validation5_GreenupVsMoisture_Scatter.png"),
       p_scatter_moisture, width = 10, height = 6, dpi = 200)


#### step 6: example pixel time series plots ####

## pick pixels with rss_spline < 2 that have a moisture match
## sample across different lag lengths for variety
set.seed(42)

## count how many peaks each pixel has
moisture_greenup_counts <- moisture_greenup %>%
  filter(rss_spline < rss_threshold,
         !is.na(spline_peak_jday),
         spline_matched) %>%
  mutate(
    n_greenup_events = rowSums(!is.na(across(c(spline_peak_1, spline_peak_2,
                                               spline_peak_3, spline_peak_4))))
  )

## 3 single-event pixels
single_event <- moisture_greenup_counts %>%
  filter(n_greenup_events == 1) %>%
  slice_sample(n = 6)

## 3 multi-event pixels
multi_event <- moisture_greenup_counts %>%
  filter(n_greenup_events > 1) %>%
  slice_sample(n = 6)

example_pixels <- bind_rows(single_event, multi_event)
# 
# example_pixels <- moisture_greenup %>%
#   filter(rss_spline < rss_threshold,
#          !is.na(spline_peak_jday),
#          spline_matched) %>%
#   mutate(lag_bin = cut(spline_lag_days,
#                        breaks = c(0, 14, 30, 60, Inf),
#                        labels = c("0-14d", "15-30d", "31-60d", ">60d"))) %>%
#   group_by(lag_bin, year) %>%
#   slice_sample(n = 1) %>%
#   ungroup() %>%
#   slice_sample(n = 12)


example_plots <- map(seq_len(nrow(example_pixels)), function(j) {
  
  px  <- example_pixels[j, ]
  yr  <- px$year
  pid <- px$pixel_id
  
  daymet_file <- file.path(daymet_path, paste0("DAYMET_", yr, ".tif"))
  if (!file.exists(daymet_file)) return(NULL)
  
  r_day <- rast(daymet_file)
  r_dc  <- mask(crop(r_day, md_vect), md_vect)
  rm(r_day); gc()
  
  bnames    <- names(r_dc)
  variables <- sub("^\\d{8}_([a-z]+)_.*$", "\\1", bnames)
  dates_str <- sub("^.*_(\\d{4}_\\d{2}_\\d{2})$", "\\1", bnames)
  jdays_all <- yday(as.Date(dates_str, "%Y_%m_%d"))
  
  ## extract this pixel only
  xy_px <- sampled_xy[pid, , drop = FALSE]
  v_raw <- terra::extract(r_dc, xy_px)
  has_id_ex <- names(v_raw)[1] == "ID"
  v_raw <- if (has_id_ex) as.numeric(v_raw[1, -1]) else as.numeric(v_raw[1, ])
  rm(r_dc); gc()
  
  prcp <- v_raw[variables == "prcp"]
  swe  <- v_raw[variables == "swe"]
  jd   <- jdays_all[variables == "prcp"]
  
  df <- tibble(jday = jd, prcp = prcp, swe = swe)
  
  ## build green-up events data frame from all four peaks + windows
  greenup_events_df <- bind_rows(
    tibble(peak  = px$spline_peak_1,
           start = px$spline_start_1,
           end   = px$spline_end_1),
    tibble(peak  = px$spline_peak_2,
           start = px$spline_start_2,
           end   = px$spline_end_2),
    tibble(peak  = px$spline_peak_3,
           start = px$spline_start_3,
           end   = px$spline_end_3),
    tibble(peak  = px$spline_peak_4,
           start = px$spline_start_4,
           end   = px$spline_end_4)
  ) %>%
    filter(!is.na(peak), !is.na(start), !is.na(end))
  
  ggplot(df) +
    ## shaded green-up windows
    geom_rect(data = greenup_events_df,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              fill = "#2ca25f", alpha = 0.15, inherit.aes = FALSE) +
    ## swe as dashed blue line
    geom_line(aes(x = jday, y = swe),
              color = "#2171b5", linetype = "dashed", linewidth = 0.8,
              na.rm = TRUE) +
    ## daily precip as black bars
    geom_segment(aes(x = jday, xend = jday, y = 0, yend = prcp),
                 color = "gray20", linewidth = 0.4, na.rm = TRUE) +
    ## green-up peak vertical lines
    geom_vline(data = greenup_events_df,
               aes(xintercept = peak),
               color = "#2ca25f", linewidth = 1.2) +
    labs(
      title    = paste0("Pixel ", pid, " | ", yr),
      subtitle = paste0("Spline RSS = ", round(px$rss_spline, 2),
                        " | lag = ", round(px$spline_lag_days, 0), " days",
                        " | ", px$spline_event_type),
      x = "Julian day",
      y = "Moisture (mm)"
    ) +
    scale_x_continuous(limits = c(0, 365)) +
    theme_bw(base_size = 10) +
    theme(plot.title    = element_text(size = 9),
          plot.subtitle = element_text(size = 7))
})

example_plots <- Filter(Negate(is.null), example_plots)

if (length(example_plots) > 0) {
  p_examples <- wrap_plots(example_plots, ncol = 4)
  ggsave(file.path(base_path_savefigures,
                   "Validation5_ExamplePixels_MoistureGreenup.png"),
         p_examples, width = 12, height = 12, dpi = 200)
  print("saved example pixel plots")
}

print("moisture-greenup analysis complete")