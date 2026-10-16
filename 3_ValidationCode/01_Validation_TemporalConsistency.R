## for pixels where both methods fit well across most years, compare the
## year-to-year stability (sd) of estimated spring green-up dates
##   - for each pixel, compute the sd of spring peak date across years,
##     separately for spline and dlc
##   - only include years where the pixel passed rss < threshold for that method
##   - compare sd distributions between methods overall
##   - also compute the cv of annual peak date to account for pixels that green
##     up very early vs. very late (early pixels naturally have lower variance)
##   - flag pixels where spline sd >> dlc sd as potential instability cases
##     and inspect them to make sure they're biologically sensible rather than
##     fitting noise

## run 00_Config.R first -- uses base_path_rasters, base_path_saveproducts,
## md_vect, years_of_interest, rss_threshold, spring_start_jday,
## spring_end_jday, min_valid_years


#### step 1: extract spring peak date and rss for every pixel x year ####

## memory strategy:
##   - crop and mask the raster to md_vect immediately on load so we only
##     ever hold the clipped extent in memory
##   - extract only the 7 bands we need (not the full ~48-band stack)
##   - accumulate running sufficient statistics (n, sum, sum-of-squares) per
##     pixel rather than storing the full year x pixel table -- memory use is
##     O(n_pixels) regardless of how many years we process, not O(n_pixels x n_years)
##   - remove each raster from memory with rm() + gc() after processing

print("extracting spring peak dates across all years...")

## get stable pixel ids from the first available raster
## pixel ids are assigned from the cropped extent of the first raster found,
## and every subsequent year is cropped to the same extent so ids stay aligned
ref_rast_path <- NULL
for (yr in years_of_interest) {
  p <- file.path(base_path_rasters, paste0("FULL_SplineMetrics_Raster_", yr, ".tif"))
  if (file.exists(p)) { ref_rast_path <- p; break }
}
if (is.null(ref_rast_path)) stop("no raster files found in base_path_rasters")

ref_rast    <- rast(ref_rast_path)[[1]]
ref_cropped <- crop(ref_rast, md_vect)
ref_masked  <- mask(ref_cropped, md_vect)
n_pixels    <- ncell(ref_masked)
rm(ref_rast, ref_cropped); gc()

print(paste0("pixels within mule deer range: ", n_pixels))

## initialise running accumulators (one slot per pixel)
## spline spring peak: n, sum, sum_sq (to compute mean and sd at the end)
## dlc spring peak: n, sum, sum_sq
acc <- list(
  sp_n    = integer(n_pixels),
  sp_sum  = numeric(n_pixels),
  sp_ss   = numeric(n_pixels),   ## sum of squares
  dlc_n   = integer(n_pixels),
  dlc_sum = numeric(n_pixels),
  dlc_ss  = numeric(n_pixels)
)

bands_needed <- c("greenup_peak_1", "greenup_peak_2",
                  "greenup_peak_3", "greenup_peak_4",
                  "rss_spline", "maxIRGdate", "rss_dlc")

## loop over years, update accumulators, never store the full table
for (yr in years_of_interest) {

  path <- file.path(base_path_rasters, paste0("FULL_SplineMetrics_Raster_", yr, ".tif"))
  if (!file.exists(path)) {
    print(paste0(yr, " -- file not found, skipping"))
    next
  }

  ## load only the bands we need, crop and mask immediately
  r       <- rast(path)[[bands_needed]]
  r_crop  <- crop(r, md_vect)
  r_mask  <- mask(r_crop, md_vect)
  rm(r, r_crop); gc()

  ## pull values as a plain matrix -- faster than as.data.frame for arithmetic
  v <- values(r_mask)   ## n_pixels x 7 matrix
  rm(r_mask); gc()

  ## column indices (bands_needed order defined above)
  i_p1  <- 1L; i_p2 <- 2L; i_p3 <- 3L; i_p4 <- 4L
  i_rss <- 5L; i_irg <- 6L; i_rdlc <- 7L

  ## spline spring peak: earliest of the 4 peaks within the spring window,
  ## set to na where rss >= threshold
  p_mat <- v[, c(i_p1, i_p2, i_p3, i_p4), drop = FALSE]

  ## mask peaks outside spring window
  p_mat[p_mat < spring_start_jday | p_mat > spring_end_jday] <- NA

  ## row-wise minimum (earliest spring peak per pixel)
  ## using apply here instead of pmin since pmin propagates na if any column is na
  sp_peak <- apply(p_mat, 1, function(x) {
    v <- x[!is.na(x)]
    if (length(v) == 0) NA_real_ else min(v)
  })

  ## zero out poor fits
  sp_peak[!is.na(v[, i_rss]) & v[, i_rss] >= rss_threshold] <- NA

  ## dlc spring peak: use maxIRGdate if within spring window and rss < threshold
  dlc_peak <- v[, i_irg]
  dlc_peak[is.na(dlc_peak) |
             dlc_peak < spring_start_jday |
             dlc_peak > spring_end_jday  |
             (!is.na(v[, i_rdlc]) & v[, i_rdlc] >= rss_threshold)] <- NA

  ## update running accumulators for valid (non-na) pixels
  sp_valid  <- !is.na(sp_peak)
  dlc_valid <- !is.na(dlc_peak)

  acc$sp_n[sp_valid]    <- acc$sp_n[sp_valid]    + 1L
  acc$sp_sum[sp_valid]  <- acc$sp_sum[sp_valid]  + sp_peak[sp_valid]
  acc$sp_ss[sp_valid]   <- acc$sp_ss[sp_valid]   + sp_peak[sp_valid]^2

  acc$dlc_n[dlc_valid]    <- acc$dlc_n[dlc_valid]    + 1L
  acc$dlc_sum[dlc_valid]  <- acc$dlc_sum[dlc_valid]  + dlc_peak[dlc_valid]
  acc$dlc_ss[dlc_valid]   <- acc$dlc_ss[dlc_valid]   + dlc_peak[dlc_valid]^2

  rm(v, p_mat, sp_peak, dlc_peak); gc()
  print(paste0(yr, " done"))
}

## compute mean and sd from accumulators (welford-style from sum/ss)
## sd = sqrt( (sum_sq - n * mean^2) / (n - 1) )  -- sample sd
compute_sd <- function(n, s, ss) {
  mean_val <- ifelse(n > 0, s / n, NA_real_)
  var_val  <- ifelse(n > 1, (ss - n * mean_val^2) / (n - 1L), NA_real_)
  list(mean = mean_val, sd = sqrt(pmax(var_val, 0)))
}

sp_stats  <- compute_sd(acc$sp_n,  acc$sp_sum,  acc$sp_ss)
dlc_stats <- compute_sd(acc$dlc_n, acc$dlc_sum, acc$dlc_ss)

## pack into a tidy per-pixel data frame -- this is the only large object kept
spring_peaks_raw <- tibble(
  pixel_id       = seq_len(n_pixels),
  n_valid_spline = acc$sp_n,
  mean_spline    = sp_stats$mean,
  sd_spline      = sp_stats$sd,
  n_valid_dlc    = acc$dlc_n,
  mean_dlc       = dlc_stats$mean,
  sd_dlc         = dlc_stats$sd
)

rm(acc, sp_stats, dlc_stats); gc()

saveRDS(spring_peaks_raw,
        file.path(base_path_saveproducts, "Validation1_SpringPeaks_Raw.rds"))
print(paste0("done. pixels with data: ", sum(spring_peaks_raw$n_valid_spline > 0)))


#### step 2: finalise per-pixel stability metrics ####

print("finalising per-pixel temporal stability metrics...")

## spring_peaks_raw already has mean and sd per pixel from the accumulators
## just add cv, sd_diff, and apply the min_valid_years filter here
pixel_stability <- spring_peaks_raw %>%
  mutate(
    cv_spline = (sd_spline / mean_spline) * 100,
    cv_dlc    = (sd_dlc    / mean_dlc)    * 100,
    sd_diff   = sd_spline - sd_dlc        ## positive = spline more variable
  ) %>%
  filter(n_valid_spline >= min_valid_years,
         n_valid_dlc    >= min_valid_years)

print(paste0("pixels with >= ", min_valid_years,
             " valid years in both methods: ", nrow(pixel_stability)))


#### step 3: overall summary ####

stability_summary <- pixel_stability %>%
  summarise(
    n_pixels               = n(),
    mean_sd_spline         = round(mean(sd_spline,  na.rm = TRUE), 2),
    mean_sd_dlc            = round(mean(sd_dlc,     na.rm = TRUE), 2),
    median_sd_spline       = round(median(sd_spline, na.rm = TRUE), 2),
    median_sd_dlc          = round(median(sd_dlc,    na.rm = TRUE), 2),
    mean_cv_spline         = round(mean(cv_spline,  na.rm = TRUE), 2),
    mean_cv_dlc            = round(mean(cv_dlc,     na.rm = TRUE), 2),
    pct_spline_more_stable = round(mean(sd_diff < 0, na.rm = TRUE) * 100, 1)
  )

print("stability summary:")
print(stability_summary)


#### step 4: statistical test -- are sd distributions different? ####

wilcox_results <- wilcox.test(pixel_stability$sd_spline,
                              pixel_stability$sd_dlc,
                              paired = TRUE)

print("paired wilcoxon test: spline sd vs dlc sd per pixel")
print(paste0("W = ", wilcox_results$statistic,
             "  p = ", round(wilcox_results$p.value, 6)))
print(paste0("median difference (spline - dlc): ",
             round(median(pixel_stability$sd_diff, na.rm = TRUE), 2), " days"))


#### step 5: flag unstable spline pixels ####

sd_absolute_threshold <- 14
sd_relative_threshold <- 10

unstable_pixels <- pixel_stability %>%
  filter(sd_spline > sd_absolute_threshold,
         sd_diff   > sd_relative_threshold) %>%
  arrange(desc(sd_diff))

print(paste0("pixels flagged as potentially unstable (spline >> dlc): ",
             nrow(unstable_pixels)))


#### step 6: distribution plots ####

stability_long <- pixel_stability %>%
  dplyr::select(pixel_id, sd_spline, sd_dlc, cv_spline, cv_dlc) %>%
  pivot_longer(
    cols      = c(sd_spline, sd_dlc, cv_spline, cv_dlc),
    names_to  = c("metric", "method"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  mutate(
    method = recode(method, "spline" = "Spline", "dlc" = "DLC"),
    metric = recode(metric,
                    "sd" = "SD of spring peak (days)",
                    "cv" = "CV of spring peak (%)")
  )

p_dist <- ggplot(stability_long %>% filter(!is.na(value)),
                 aes(x = value, fill = method, color = method)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  facet_wrap(~ metric, scales = "free") +
  scale_fill_manual(values  = method_colors) +
  scale_color_manual(values = method_colors) +
  labs(
    title    = "Year-to-year stability of spring green-up estimates",
    subtitle = paste0("Pixels with >= ", min_valid_years,
                      " valid years in both methods (RSS < ", rss_threshold, ")"),
    x = "Value", y = "Density", fill = "Method", color = "Method"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(base_path_savefigures, "Validation1_StabilityDistributions.png"),
       p_dist, width = 8, height = 6, dpi = 200)

p_scatter <- ggplot(pixel_stability, aes(x = sd_dlc, y = sd_spline)) +
  geom_point(alpha = 0.2, size = 0.6, color = "gray40") +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "gray20") +
  labs(
    title    = "Per-pixel SD: Spline vs. DLC",
    subtitle = "Points above dashed line = spline more variable than DLC",
    x = "DLC spring peak SD (days)",
    y = "Spline spring peak SD (days)"
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(base_path_savefigures, "Validation1_SDScatter_SplineVsDLC.png"),
       p_scatter, width = 7, height = 6, dpi = 200)

print("validation 1 complete")


#### save outputs ####

saveRDS(pixel_stability,
        file.path(base_path_saveproducts, "Validation1_PixelStability.rds"))
saveRDS(stability_summary,
        file.path(base_path_saveproducts, "Validation1_StabilitySummary.rds"))
saveRDS(unstable_pixels,
        file.path(base_path_saveproducts, "Validation1_UnstablePixels.rds"))

gc()
