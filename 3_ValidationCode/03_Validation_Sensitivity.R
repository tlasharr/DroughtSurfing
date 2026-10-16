## test sensitivity of spline green-up estimates to the two key parameters:
##   1. smoothing parameter (spar) -- current value 0.4
##   2. ndvi change threshold      -- current value 0.05
## runs on a sample of pixels from the utah latitudinal bounding box

## run 00_Config.R first -- uses base_path_rasters, base_path_saveproducts,
## ndvi_base_path, md_vect, years_of_interest, rss_threshold,
## spring_start_jday, spring_end_jday, spar_baseline, threshold_baseline,
## method_colors -- also needs fundlc() sourced in config

spar_values      <- c(0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6)
threshold_values <- c(0.025, 0.05, 0.075, 0.10, 0.125, 0.15)
n_sample_pixels  <- 2000

set.seed(42)


#### step 1: sample representative pixels ####

print("sampling pixels...")

r_crop   <- crop(rast(ref_rast_path)[[1]], md_vect)
r_mask   <- mask(r_crop, md_vect)
rss_vals <- values(r_mask)[, 1]
eligible <- which(!is.na(rss_vals))
sampled_idx <- sample(eligible, min(n_sample_pixels, length(eligible)))
rm(r_crop, r_mask, rss_vals); gc()

print(paste0("sampled ", length(sampled_idx), " pixels"))


#### step 1b: refit dlc from raw ndvi for sampled pixels, per year ####

## compute dlc rss directly so there's no pixel mismatch with the spline fits

print("refitting dlc for sampled pixels...")

dlc_rss_by_year <- map_dfr(years_of_interest, function(yr) {

  ndvi_path <- file.path(ndvi_base_path, paste0("Utah_", yr, ".tif"))
  if (!file.exists(ndvi_path)) {
    print(paste0("skipping ", yr, " -- ndvi raster not found"))
    return(NULL)
  }

  ndvi_rast  <- rast(ndvi_path)
  ndvi_crop  <- crop(ndvi_rast, md_vect)
  ndvi_mask  <- mask(ndvi_crop, md_vect)
  rm(ndvi_rast, ndvi_crop); gc()

  layer_dates <- names(ndvi_mask)

  all_vals       <- values(ndvi_mask)
  n_valid        <- rowSums(!is.na(all_vals))
  eligible       <- which(n_valid >= 5)

  if (length(eligible) == 0) {
    rm(ndvi_mask, all_vals); gc()
    return(NULL)
  }

  set.seed(42)
  sampled_idx_yr <- sample(eligible, min(n_sample_pixels, length(eligible)))
  pixel_vals     <- all_vals[sampled_idx_yr, , drop = FALSE]
  rm(ndvi_mask, all_vals); gc()

  scale_values <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) return(rep(NA, length(x)))
    (x - rng[1]) / (rng[2] - rng[1])
  }

  map_dfr(seq_len(nrow(pixel_vals)), function(i) {

    ts_raw <- as.numeric(pixel_vals[i, ])
    if (sum(!is.na(ts_raw)) < 5) return(NULL)

    ts_scaled  <- scale_values(ts_raw)
    if (all(is.na(ts_scaled))) return(NULL)

    ## name the vector with dates so fundlc can use them
    ndvi_named <- setNames(ts_scaled, layer_dates)

    ## fit dlc and compute rss
    dlc_params <- tryCatch(fundlc(ndvi_named), error = function(e) rep(NA, 16))

    dlc_rss <- if (!any(is.na(dlc_params[1:4]))) {
      ## reconstruct the fitted dlc curve at the observed time points
      ## fundlc returns xmidS, scalS, xmidA, scalA as first 4 elements
      jdays_obs <- yday(as.Date(substr(layer_dates, 1, 10), "%Y_%m_%d"))
      valid_idx  <- which(!is.na(ts_scaled))
      t_obs      <- jdays_obs[valid_idx]
      y_obs      <- ts_scaled[valid_idx]

      y_fit <- (1 / (1 + exp((dlc_params[1] - t_obs) / dlc_params[2]))) -
        (1 / (1 + exp((dlc_params[3] - t_obs) / dlc_params[4])))

      sum((y_obs - y_fit)^2, na.rm = TRUE)
    } else {
      NA_real_
    }

    tibble(
      pixel_id = sampled_idx_yr[i],
      year     = yr,
      rss_dlc  = dlc_rss
    )
  })
})

print(paste0("dlc rss computed. rows: ", nrow(dlc_rss_by_year)))
print(paste0("pixels with valid dlc rss: ",
             sum(!is.na(dlc_rss_by_year$rss_dlc))))
print(paste0("mean dlc rss: ",
             round(mean(dlc_rss_by_year$rss_dlc, na.rm = TRUE), 3)))


#### step 2: sensitivity function -- refit spline across parameter grid ####

run_sensitivity_for_pixel <- function(scaled_ts, jdays, spar_vals, threshold_vals,
                                      spring_start, spring_end) {
  valid_idx   <- which(!is.na(scaled_ts))
  if (length(valid_idx) < 5) return(NULL)

  ## use actual julian days as time points, not sequential indices
  time_points <- jdays[valid_idx]
  scaled_ts   <- scaled_ts[valid_idx]

  map_dfr(spar_vals, function(sp) {
    map_dfr(threshold_vals, function(thresh) {

      weights      <- 1 / (abs(diff(c(0, scaled_ts))) + 0.001)
      spline_model <- tryCatch(
        smooth.spline(time_points, scaled_ts, spar = sp, w = weights),
        error = function(e) NULL
      )
      if (is.null(spline_model)) {
        return(tibble(spar = sp, threshold = thresh,
                      spring_peak = NA_real_, n_greenup = NA_integer_, rss = NA_real_))
      }

      ## predict on a daily grid so peak detection is at day resolution
      daily_grid <- seq(min(time_points), max(time_points))
      pred_y     <- pmin(pmax(predict(spline_model, daily_grid)$y, 0), 1)
      rss        <- sum((scaled_ts - predict(spline_model, time_points)$y)^2,
                        na.rm = TRUE)

      delta_y <- diff(pred_y)
      if (all(is.na(delta_y)) || all(delta_y == 0)) {
        return(tibble(spar = sp, threshold = thresh,
                      spring_peak = NA_real_, n_greenup = 0L, rss = rss))
      }

      rle_dir    <- rle(sign(delta_y))
      seg_ends   <- cumsum(rle_dir$lengths)
      seg_starts <- c(1, head(seg_ends, -1) + 1)

      n_greenup    <- 0L
      spring_peaks <- numeric(0)

      for (i in seq_along(rle_dir$values)) {
        if (is.na(rle_dir$values[i]) || rle_dir$values[i] != 1) next

        idx          <- seg_starts[i]:seg_ends[i]
        end_idx      <- min(max(idx) + 1, length(pred_y))
        total_change <- pred_y[end_idx] - pred_y[idx[1]]
        if (is.na(total_change) || total_change < thresh) next

        n_greenup <- n_greenup + 1L
        if (length(idx) < 2) next

        roc      <- diff(pred_y[idx])
        if (length(roc) == 0 || all(is.na(roc))) next

        ## peak_day is now an actual julian day
        peak_day <- daily_grid[idx[which.max(roc)]]

        if (!is.na(peak_day) && peak_day >= spring_start && peak_day <= spring_end)
          spring_peaks <- c(spring_peaks, peak_day)
      }

      tibble(spar = sp, threshold = thresh,
             spring_peak = if (length(spring_peaks) == 0) NA_real_ else min(spring_peaks),
             n_greenup = n_greenup, rss = rss)
    })
  })
}


#### step 3: run across sampled pixels x years ####

print("running sensitivity analysis...")

sensitivity_results <- map_dfr(years_of_interest, function(yr) {

  ndvi_path <- file.path(ndvi_base_path, paste0("Utah_", yr, ".tif"))
  if (!file.exists(ndvi_path)) {
    print(paste0("skipping ", yr, " -- ndvi raster not found"))
    return(NULL)
  }

  ndvi_rast  <- rast(ndvi_path)
  ndvi_crop  <- crop(ndvi_rast, md_vect)
  ndvi_mask  <- mask(ndvi_crop, md_vect)
  rm(ndvi_rast, ndvi_crop); gc()

  ## parse julian days from layer names
  layer_dates <- names(ndvi_mask)
  jdays <- yday(as.Date(substr(layer_dates, 1, 10), "%Y_%m_%d"))
  print(paste0(yr, " jday range: ", min(jdays, na.rm = TRUE),
               " - ", max(jdays, na.rm = TRUE)))

  ## sample pixels from this raster
  all_vals       <- values(ndvi_mask)
  n_valid        <- rowSums(!is.na(all_vals))
  eligible       <- which(n_valid >= 5)

  if (length(eligible) == 0) {
    print(paste0("skipping ", yr, " -- no eligible pixels"))
    rm(ndvi_mask, all_vals); gc()
    return(NULL)
  }

  set.seed(42)
  sampled_idx_yr <- sample(eligible, min(n_sample_pixels, length(eligible)))
  pixel_vals     <- all_vals[sampled_idx_yr, , drop = FALSE]
  rm(ndvi_mask, all_vals); gc()

  scale_values <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) return(rep(NA, length(x)))
    (x - rng[1]) / (rng[2] - rng[1])
  }

  map_dfr(seq_len(nrow(pixel_vals)), function(i) {
    ts_raw    <- as.numeric(pixel_vals[i, ])
    if (sum(!is.na(ts_raw)) < 5) return(NULL)
    ts_scaled <- scale_values(ts_raw)
    if (all(is.na(ts_scaled))) return(NULL)

    result          <- run_sensitivity_for_pixel(
      ts_scaled, jdays,
      spar_values, threshold_values,
      spring_start_jday, spring_end_jday
    )
    if (is.null(result)) return(NULL)
    result$pixel_id <- sampled_idx_yr[i]
    result$year     <- yr
    result
  })
})

saveRDS(sensitivity_results,
        file.path(base_path_saveproducts, "Validation3_SensitivityResults_Raw.rds"))
print(paste0("done. rows: ", nrow(sensitivity_results)))


#### step 4: summarise sensitivity ####

baseline <- sensitivity_results %>%
  filter(spar == spar_baseline, threshold == threshold_baseline) %>%
  dplyr::select(pixel_id, year, spring_peak_baseline = spring_peak,
         n_greenup_baseline = n_greenup)

sensitivity_delta <- sensitivity_results %>%
  left_join(baseline, by = c("pixel_id", "year")) %>%
  mutate(delta_peak    = spring_peak - spring_peak_baseline,
         delta_greenup = n_greenup   - n_greenup_baseline)

param_sensitivity_summary <- sensitivity_delta %>%
  filter(!is.na(delta_peak)) %>%
  group_by(spar, threshold) %>%
  summarise(n_pixels       = n(),
            mean_delta     = round(mean(delta_peak,       na.rm = TRUE), 2),
            mad_delta      = round(mean(abs(delta_peak),  na.rm = TRUE), 2),
            sd_delta       = round(sd(delta_peak,         na.rm = TRUE), 2),
            pct_within_7d  = round(mean(abs(delta_peak) <= 7,  na.rm = TRUE) * 100, 1),
            pct_within_14d = round(mean(abs(delta_peak) <= 14, na.rm = TRUE) * 100, 1),
            .groups = "drop")

print("sensitivity summary:")
print(param_sensitivity_summary, n = Inf)


#### step 4b: compare spline rss across parameter grid to dlc rss ####

## summarise dlc rss across all pixels and years
dlc_rss_summary <- dlc_rss_by_year %>%
  summarise(
    n_pixels        = sum(!is.na(rss_dlc)),
    mean_rss_dlc    = round(mean(rss_dlc,   na.rm = TRUE), 3),
    median_rss_dlc  = round(median(rss_dlc, na.rm = TRUE), 3),
    pct_fit_dlc     = round(mean(rss_dlc < rss_threshold, na.rm = TRUE) * 100, 1)
  )

print("dlc rss summary:")
print(dlc_rss_summary)

## for each parameter combination, what proportion of pixels have
## spline rss < dlc rss (i.e., spline fits better than dlc)?
sensitivity_with_dlc <- sensitivity_results %>%
  left_join(dlc_rss_by_year, by = c("pixel_id", "year"))

rss_comparison <- sensitivity_with_dlc %>%
  filter(!is.na(rss_dlc), !is.na(rss)) %>%
  group_by(spar, threshold) %>%
  summarise(
    n_pixels            = n(),
    mean_rss_spline     = round(mean(rss,     na.rm = TRUE), 3),
    mean_rss_dlc        = round(mean(rss_dlc, na.rm = TRUE), 3),
    pct_spline_fits     = round(mean(rss     < rss_threshold, na.rm = TRUE) * 100, 1),
    pct_dlc_fits        = round(mean(rss_dlc < rss_threshold, na.rm = TRUE) * 100, 1),
    pct_spline_beats_dlc = round(mean(rss < rss_dlc, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  )

print("rss comparison -- spline vs dlc across parameter grid:")
print(rss_comparison %>%
        filter(threshold == threshold_baseline) %>%
        dplyr::select(spar, threshold, mean_rss_spline, mean_rss_dlc,
               pct_spline_fits, pct_dlc_fits, pct_spline_beats_dlc),
      n = Inf)

## heatmap: proportion of pixels where spline rss < dlc rss
p_rss_comparison <- rss_comparison %>%
  ggplot(aes(x = factor(spar), y = factor(threshold),
             fill = pct_spline_beats_dlc)) +
  geom_tile(color = "white") +
  geom_text(aes(label = paste0(pct_spline_beats_dlc, "%")), size = 3.5) +
  scale_fill_gradient(low = "#f7fbff", high = "#2ca25f") +
  labs(
    title    = "Proportion of pixels where spline RSS < DLC RSS",
    subtitle = paste0("RSS threshold = ", rss_threshold,
                      " | Green = spline fits better than DLC"),
    x    = "Smoothing parameter (spar)",
    y    = "NDVI change threshold",
    fill = "% spline\nbeats DLC"
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(base_path_savefigures,
                 "Validation3_RSS_SplineVsDLC_Heatmap.png"),
       p_rss_comparison, width = 10, height = 5, dpi = 200)

## line plot: mean rss by spar (threshold fixed at baseline)
p_rss_by_spar <- rss_comparison %>%
  filter(threshold == threshold_baseline) %>%
  dplyr::select(spar, mean_rss_spline, mean_rss_dlc) %>%
  pivot_longer(cols = c(mean_rss_spline, mean_rss_dlc),
               names_to = "method", values_to = "mean_rss") %>%
  mutate(method = recode(method,
                         "mean_rss_spline" = "Spline",
                         "mean_rss_dlc"    = "DLC")) %>%
  ggplot(aes(x = spar, y = mean_rss, color = method, group = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = rss_threshold,
             linetype = "dashed", color = "red") +
  geom_vline(xintercept = spar_baseline,
             linetype = "dashed", color = "gray40") +
  scale_color_manual(values = method_colors) +
  labs(
    title    = "Mean RSS by smoothing parameter: spline vs. DLC",
    subtitle = paste0("Threshold fixed at ", threshold_baseline,
                      " | Red dashed = RSS threshold (",
                      rss_threshold, ")"),
    x     = "spar",
    y     = "Mean RSS",
    color = "Method"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(base_path_savefigures,
                 "Validation3_RSS_MeanBySpline_vs_DLC.png"),
       p_rss_by_spar, width = 7, height = 5, dpi = 200)

saveRDS(rss_comparison,
        file.path(base_path_saveproducts, "Validation3_RSS_Comparison.rds"))
saveRDS(dlc_rss_summary,
        file.path(base_path_saveproducts, "Validation3_DLC_RSS_Summary.rds"))


#### step 5: sensitivity plots ####

p_heatmap <- param_sensitivity_summary %>%
  ggplot(aes(x = factor(spar), y = factor(threshold), fill = mad_delta)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(mad_delta, 1)), size = 3.5) +
  scale_fill_gradient(low = "#f7fbff", high = "#08306b") +
  labs(title    = "Sensitivity of spring green-up peak date to parameter choices",
       subtitle = paste0("MAD from baseline (spar = ", spar_baseline,
                         ", threshold = ", threshold_baseline, ")"),
       x = "Smoothing parameter (spar)", y = "NDVI change threshold",
       fill = "MAD (days)") +
  theme_bw(base_size = 12)

ggsave(file.path(base_path_savefigures, "Validation3_Sensitivity_Heatmap.png"),
       p_heatmap, width = 8, height = 5, dpi = 200)

p_spar <- param_sensitivity_summary %>%
  filter(threshold == threshold_baseline) %>%
  ggplot(aes(x = spar, y = mad_delta)) +
  geom_line(linewidth = 1, color = "#2ca25f") +
  geom_point(size = 3, color = "#2ca25f") +
  geom_vline(xintercept = spar_baseline, linetype = "dashed", color = "gray40") +
  labs(title    = "Sensitivity to smoothing parameter (spar)",
       subtitle = paste0("Threshold fixed at ", threshold_baseline),
       x = "spar", y = "MAD from baseline (days)") +
  theme_bw(base_size = 12)

p_thresh <- param_sensitivity_summary %>%
  filter(spar == spar_baseline) %>%
  ggplot(aes(x = threshold, y = mad_delta)) +
  geom_line(linewidth = 1, color = "#756bb1") +
  geom_point(size = 3, color = "#756bb1") +
  geom_vline(xintercept = threshold_baseline, linetype = "dashed", color = "gray40") +
  labs(title    = "Sensitivity to NDVI change threshold",
       subtitle = paste0("spar fixed at ", spar_baseline),
       x = "NDVI change threshold", y = "MAD from baseline (days)") +
  theme_bw(base_size = 12)

ggsave(file.path(base_path_savefigures, "Validation3_Sensitivity_Spar.png"),
       p_spar, width = 6, height = 5, dpi = 200)
ggsave(file.path(base_path_savefigures, "Validation3_Sensitivity_Threshold.png"),
       p_thresh, width = 6, height = 5, dpi = 200)

print("sensitivity analysis complete")
saveRDS(param_sensitivity_summary,
        file.path(base_path_saveproducts, "Validation3_SensitivitySummary.rds"))


#### step 6: visual validation -- spline and dlc fits with green-up events ####

## plot a sample of pixels at the baseline parameter settings, with the
## fitted curves and identified green-up events overlaid

print("building green-up visualization...")

viz_yr       <- 2020
n_viz_pixels <- 16

ndvi_path_viz <- file.path(ndvi_base_path, paste0("Utah_", viz_yr, ".tif"))
ndvi_rast_viz <- rast(ndvi_path_viz)
ndvi_crop_viz <- crop(ndvi_rast_viz, md_vect)
ndvi_mask_viz <- mask(ndvi_crop_viz, md_vect)
rm(ndvi_rast_viz, ndvi_crop_viz); gc()

layer_dates_viz <- names(ndvi_mask_viz)
jdays_viz       <- yday(as.Date(substr(layer_dates_viz, 1, 10), "%Y_%m_%d"))

all_vals_viz <- values(ndvi_mask_viz)
n_valid_viz  <- rowSums(!is.na(all_vals_viz))
eligible_viz <- which(n_valid_viz >= 5)

set.seed(99)
sampled_viz    <- sample(eligible_viz, min(n_viz_pixels, length(eligible_viz)))
pixel_vals_viz <- all_vals_viz[sampled_viz, , drop = FALSE]
rm(ndvi_mask_viz, all_vals_viz); gc()

scale_values <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(NA, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

identify_greenup_events <- function(pred_y, daily_grid,
                                    thresh, spring_start, spring_end) {
  delta_y <- diff(pred_y)
  if (all(is.na(delta_y)) || all(delta_y == 0)) return(NULL)

  rle_dir    <- rle(sign(delta_y))
  seg_ends   <- cumsum(rle_dir$lengths)
  seg_starts <- c(1, head(seg_ends, -1) + 1)

  map_dfr(seq_along(rle_dir$values), function(i) {
    if (is.na(rle_dir$values[i]) || rle_dir$values[i] != 1) return(NULL)
    idx          <- seg_starts[i]:seg_ends[i]
    end_idx      <- min(max(idx) + 1, length(pred_y))
    total_change <- pred_y[end_idx] - pred_y[idx[1]]
    if (is.na(total_change) || total_change < thresh) return(NULL)
    if (length(idx) < 2) return(NULL)
    roc      <- diff(pred_y[idx])
    if (length(roc) == 0 || all(is.na(roc))) return(NULL)
    peak_day <- daily_grid[idx[which.max(roc)]]
    tibble(
      start_jday = daily_grid[idx[1]],
      end_jday   = daily_grid[min(end_idx, length(daily_grid))],
      peak_jday  = peak_day,
      is_spring  = !is.na(peak_day) &
        peak_day >= spring_start &
        peak_day <= spring_end
    )
  })
}

## run the loop and store as a list -- each pixel gets raw data, fitted
## lines, and identified green-up events for both methods
viz_results <- map(seq_len(nrow(pixel_vals_viz)), function(i) {

  ts_raw <- as.numeric(pixel_vals_viz[i, ])
  if (sum(!is.na(ts_raw)) < 5) return(NULL)

  ts_scaled <- scale_values(ts_raw)
  if (all(is.na(ts_scaled))) return(NULL)

  valid_idx   <- which(!is.na(ts_scaled))
  time_points <- jdays_viz[valid_idx]
  ts_fit      <- ts_scaled[valid_idx]

  ## raw data frame
  raw_df <- tibble(
    pixel_id    = sampled_viz[i],
    jday        = time_points,
    ndvi_scaled = ts_fit
  )

  #### spline ####
  weights      <- 1 / (abs(diff(c(0, ts_fit))) + 0.001)
  spline_model <- tryCatch(
    smooth.spline(time_points, ts_fit,
                  spar = spar_baseline, w = weights),
    error = function(e) NULL
  )

  spline_line   <- NULL
  spline_events <- NULL

  if (!is.null(spline_model)) {
    daily_grid <- seq(min(time_points), max(time_points))
    pred_y     <- pmin(pmax(predict(spline_model, daily_grid)$y, 0), 1)
    rss_sp     <- sum((ts_fit - predict(spline_model, time_points)$y)^2,
                      na.rm = TRUE)

    spline_line <- tibble(
      pixel_id    = sampled_viz[i],
      jday        = daily_grid,
      ndvi_scaled = pred_y,
      method      = "Spline",
      rss         = rss_sp
    )

    events <- identify_greenup_events(
      pred_y, daily_grid,
      threshold_baseline,
      spring_start_jday, spring_end_jday
    )
    if (!is.null(events) && nrow(events) > 0) {
      spline_events <- events %>%
        mutate(pixel_id = sampled_viz[i], method = "Spline")
    }
  }

  #### dlc ####
  ndvi_named <- setNames(ts_fit, layer_dates_viz[valid_idx])
  dlc_params <- tryCatch(fundlc(ndvi_named), error = function(e) rep(NA, 16))

  dlc_line   <- NULL
  dlc_events <- NULL

  if (!any(is.na(dlc_params[1:4]))) {
    daily_grid_dlc <- seq(min(time_points), max(time_points))
    dlc_pred <- (1 / (1 + exp((dlc_params[1] - daily_grid_dlc) /
                                dlc_params[2]))) -
      (1 / (1 + exp((dlc_params[3] - daily_grid_dlc) /
                      dlc_params[4])))
    rss_dlc  <- sum(
      (ts_fit -
         ((1 / (1 + exp((dlc_params[1] - time_points) / dlc_params[2]))) -
            (1 / (1 + exp((dlc_params[3] - time_points) / dlc_params[4]))))) ^ 2,
      na.rm = TRUE
    )

    dlc_line <- tibble(
      pixel_id    = sampled_viz[i],
      jday        = daily_grid_dlc,
      ndvi_scaled = dlc_pred,
      method      = "DLC",
      rss         = rss_dlc
    )

    ## dlc green-up window: +/- 30 days around peak irg
    irg      <- diff(dlc_pred)
    peak_idx <- which.max(irg)
    if (length(peak_idx) > 0) {
      start_idx <- max(1, peak_idx - 30)
      end_idx   <- min(length(daily_grid_dlc), peak_idx + 30)
      peak_day  <- daily_grid_dlc[peak_idx]

      dlc_events <- tibble(
        pixel_id   = sampled_viz[i],
        method     = "DLC",
        start_jday = daily_grid_dlc[start_idx],
        end_jday   = daily_grid_dlc[end_idx],
        peak_jday  = peak_day,
        is_spring  = peak_day >= spring_start_jday &
          peak_day <= spring_end_jday
      )
    }
  }

  list(
    raw    = raw_df,
    lines  = bind_rows(spline_line, dlc_line),
    events = bind_rows(spline_events, dlc_events)
  )
})

## remove null results
viz_results <- Filter(Negate(is.null), viz_results)

## unpack each component separately
raw_data_viz <- map_dfr(viz_results, ~ .x$raw)
fitted_viz   <- map_dfr(viz_results, ~ .x$lines)
events_viz   <- map_dfr(viz_results, ~ .x$events)

## build panel labels with rss for both methods, colored green if < threshold
## and red if >= threshold using html formatting via ggtext
rss_labels <- fitted_viz %>%
  group_by(pixel_id, method) %>%
  summarise(rss = round(unique(rss), 2), .groups = "drop") %>%
  pivot_wider(names_from   = method,
              values_from  = rss,
              names_prefix = "rss_") %>%
  mutate(
    spline_color = if_else(rss_Spline < rss_threshold, "#2ca25f", "#d7191c"),
    dlc_color    = if_else(rss_DLC    < rss_threshold, "#2ca25f", "#d7191c"),
    label = paste0(
      "Pixel ", pixel_id, "<br>",
      "<span style='color:", spline_color, "'>Spline RSS = ", rss_Spline, "</span>",
      " | ",
      "<span style='color:", dlc_color,    "'>DLC RSS = ",    rss_DLC,    "</span>"
    )
  )

raw_data_viz <- raw_data_viz %>%
  left_join(rss_labels %>% dplyr::select(pixel_id, label), by = "pixel_id")

fitted_viz <- fitted_viz %>%
  left_join(rss_labels %>% dplyr::select(pixel_id, label), by = "pixel_id") %>%
  mutate(method = factor(method, levels = c("Spline", "DLC")))

events_viz <- events_viz %>%
  left_join(rss_labels %>% dplyr::select(pixel_id, label), by = "pixel_id") %>%
  mutate(method = factor(method, levels = c("Spline", "DLC"))) %>%
  filter(
    !is.na(label),
    !is.na(start_jday),
    !is.na(end_jday),
    !is.na(peak_jday),
    start_jday >= 1,
    end_jday   <= 365,
    end_jday   > start_jday
  )

p_viz <- ggplot() +
  geom_rect(
    data = events_viz,
    aes(xmin = start_jday, xmax = end_jday,
        ymin = 0, ymax = 1, fill = method),
    alpha = 0.15, inherit.aes = FALSE
  ) +
  geom_point(
    data = raw_data_viz,
    aes(x = jday, y = ndvi_scaled),
    color = "gray55", size = 1.2, alpha = 0.7
  ) +
  geom_line(
    data = fitted_viz,
    aes(x = jday, y = ndvi_scaled, color = method),
    linewidth = 0.9
  ) +
  geom_vline(
    data = events_viz,
    aes(xintercept = peak_jday, color = method),
    linetype = "dashed", linewidth = 0.7
  ) +
  facet_wrap(~ label, ncol = 4) +
  scale_color_manual(values = method_colors, name = "Method") +
  scale_fill_manual(values  = method_colors, name = "Method") +
  ylim(0, 1) +
  labs(
    title    = paste0("Spline vs. DLC: fitted curves and green-up events (",
                      viz_yr, ")"),
    subtitle = paste0("spar = ", spar_baseline,
                      " | threshold = ", threshold_baseline,
                      "\nShaded = green-up window; dashed = peak date"),
    x = "Julian day",
    y = "Scaled NDVI"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.text      = element_markdown(size = 7),   ## needs ggtext for html color
    panel.spacing   = unit(0.4, "cm"),
    axis.text       = element_text(size = 7)
  )

p_viz

ggsave(
  file.path(base_path_savefigures,
            "Validation3_GreenupViz_SplineVsDLC.png"),
  p_viz, width = 14, height = 14, dpi = 200
)
print("saved green-up visualization")
