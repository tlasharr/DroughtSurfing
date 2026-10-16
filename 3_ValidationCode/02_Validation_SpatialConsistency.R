## for pixels where both methods fit well across most years, compare the
## spatial consistency of estimated spring green-up dates -- does fit
## quality and green-up timing cluster spatially in ecologically
## interpretable ways, and do both methods show similar spatial structure

## run 00_Config.R first -- uses base_path_rasters, base_path_saveproducts,
## md_vect, years_of_interest, rss_threshold, spring_start_jday,
## spring_end_jday, method_colors

n_neighbors <- 8

#### step 1: build per-pixel fit-proportion raster ####

print("building fit-proportion raster...")

template   <- mask(crop(rast(ref_rast_path)[[1]], md_vect), md_vect)
n_pixels   <- ncell(template)

print(paste0("pixels in bounding box: ", n_pixels))

fit_count_spline <- integer(n_pixels)
fit_count_dlc    <- integer(n_pixels)
year_count       <- integer(n_pixels)

bands_needed <- c("rss_spline", "maxIRGdate", "rss_dlc")

for (yr in years_of_interest) {
  path <- file.path(base_path_rasters, paste0("FULL_SplineMetrics_Raster_", yr, ".tif"))
  if (!file.exists(path)) { print(paste0(yr, " skipped")); next }

  r      <- rast(path)[[bands_needed]]
  r_crop <- crop(r, md_vect)
  r_mask <- mask(r_crop, md_vect)
  rm(r, r_crop); gc()

  v       <- values(r_mask)
  rm(r_mask); gc()

  rss_sp  <- v[, 1]
  irg     <- v[, 2]
  rss_dlc <- v[, 3]

  ## a pixel "attempted" if it has any non-na value, "fit" if rss < threshold
  attempted        <- !is.na(rss_sp) | !is.na(rss_dlc)
  year_count       <- year_count + as.integer(attempted)
  fit_count_spline <- fit_count_spline +
    as.integer(!is.na(rss_sp) & rss_sp < rss_threshold)
  fit_count_dlc    <- fit_count_dlc +
    as.integer(!is.na(rss_dlc) & rss_dlc < rss_threshold &
                 !is.na(irg) & irg >= spring_start_jday & irg <= spring_end_jday)

  rm(v); gc()
  print(paste0(yr, " done"))
}

prop_spline <- fit_count_spline / pmax(year_count, 1L)
prop_dlc    <- fit_count_dlc    / pmax(year_count, 1L)
prop_spline[year_count == 0] <- NA
prop_dlc[year_count    == 0] <- NA

fit_rast              <- c(template, template, template)
values(fit_rast[[1]]) <- prop_spline
values(fit_rast[[2]]) <- prop_dlc
values(fit_rast[[3]]) <- as.numeric(year_count)
names(fit_rast)       <- c("prop_fit_spline", "prop_fit_dlc", "n_years")

writeRaster(fit_rast,
            file.path(base_path_saveproducts, "Validation2_FitProportion_SplineVsDLC.tif"),
            overwrite = TRUE)

## positive = spline fits more years, negative = dlc fits more years
diff_rast        <- fit_rast[["prop_fit_spline"]] - fit_rast[["prop_fit_dlc"]]
names(diff_rast) <- "fit_prop_diff"
writeRaster(diff_rast,
            file.path(base_path_saveproducts, "Validation2_FitProportionDiff_SplineMinusDLC.tif"),
            overwrite = TRUE)

rm(fit_count_spline, fit_count_dlc, year_count, prop_spline, prop_dlc); gc()
print("fit-proportion rasters saved")


#### step 2: moran's i on fit proportion ####

print("computing moran's i on fit proportion...")

compute_morans_I <- function(r, layer_name, n_neighbors = 8) {

  r_single <- r[[layer_name]]

  vals   <- values(r_single)[, 1]
  non_na <- which(!is.na(vals))

  if (length(non_na) < 100) { print("too few non-na pixels"); return(NULL) }
  if (length(non_na) > 50000) { set.seed(42); non_na <- sample(non_na, 50000) }

  ## xyFromCell uses cell numbers directly -- no alignment issues
  xy <- xyFromCell(r_single, non_na)
  z  <- vals[non_na]

  if (nrow(xy) != length(z)) {
    print(paste0("dimension mismatch: ", nrow(xy), " coords vs ", length(z), " values"))
    return(NULL)
  }

  nb <- knn2nb(knearneigh(xy, k = n_neighbors))
  lw <- nb2listw(nb, style = "W")
  mi <- moran.test(z, lw)

  tibble(
    layer      = layer_name,
    moran_I    = mi$estimate["Moran I statistic"],
    expected_I = mi$estimate["Expectation"],
    z_score    = mi$statistic,
    p_value    = mi$p.value,
    n_cells    = length(z)
  )
}

morans_results <- bind_rows(
  compute_morans_I(fit_rast, "prop_fit_spline", n_neighbors) %>% mutate(method = "Spline"),
  compute_morans_I(fit_rast, "prop_fit_dlc",    n_neighbors) %>% mutate(method = "DLC")
)

print("moran's i results:")
print(morans_results)
saveRDS(morans_results,
        file.path(base_path_saveproducts, "Validation2_MoransI_Results.rds"))


#### step 3: moran's i on spring green-up date ####

## rather than asking whether fit quality clusters spatially (step 2),
## here we ask whether the spring green-up date itself is spatially
## autocorrelated, and whether both methods produce similar spatial structure
##
## for each year we build a mean spring peak date surface (averaging across
## pixels that were fit in that year), then compute moran's i on that surface
## for both methods. we summarise across years and compare
##
## high and similar moran's i for both methods = both capturing the same
## landscape-scale phenology gradient (e.g., low elevation greens up first,
## high elevation later). a big difference in i suggests the methods are
## producing different spatial structure in their phenology estimates

print("computing moran's i on spring green-up date...")

bands_date <- c("greenup_peak_1", "greenup_peak_2",
                "greenup_peak_3", "greenup_peak_4",
                "rss_spline", "maxIRGdate", "rss_dlc")

morans_date_results <- map_dfr(years_of_interest, function(yr) {

  path <- file.path(base_path_rasters,
                    paste0("FULL_SplineMetrics_Raster_", yr, ".tif"))
  if (!file.exists(path)) { print(paste0(yr, " skipped")); return(NULL) }

  r      <- rast(path)[[bands_date]]
  r_crop <- crop(r, md_vect)
  r_mask <- mask(r_crop, md_vect)
  rm(r, r_crop); gc()

  v <- values(r_mask)
  rm(r_mask); gc()

  ## spline: earliest spring peak within window, rss filter applied
  p_mat <- v[, 1:4, drop = FALSE]
  p_mat[p_mat < spring_start_jday | p_mat > spring_end_jday] <- NA
  sp_peak <- apply(p_mat, 1, function(x) {
    x <- x[!is.na(x)]; if (length(x) == 0) NA_real_ else min(x)
  })
  sp_peak[!is.na(v[, 5]) & v[, 5] >= rss_threshold] <- NA

  ## dlc: maxIRGdate within spring window, rss filter applied
  dlc_peak <- v[, 6]
  dlc_peak[is.na(dlc_peak) |
             dlc_peak < spring_start_jday |
             dlc_peak > spring_end_jday   |
             (!is.na(v[, 7]) & v[, 7] >= rss_threshold)] <- NA

  rm(v, p_mat); gc()

  ## re-load a single band as coordinate template for xyFromCell
  r_tmp   <- rast(path)[[1]]
  r_tmp_c <- mask(crop(r_tmp, md_vect), md_vect)
  rm(r_tmp); gc()

  compute_morans_date <- function(peak_vals, template_rast, method_name, yr) {
    non_na <- which(!is.na(peak_vals))
    if (length(non_na) < 100) {
      print(paste0(yr, " ", method_name, " -- too few valid pixels"))
      return(NULL)
    }
    if (length(non_na) > 50000) {
      set.seed(42); non_na <- sample(non_na, 50000)
    }
    xy <- xyFromCell(template_rast, non_na)
    z  <- peak_vals[non_na]
    if (nrow(xy) != length(z)) return(NULL)

    nb <- knn2nb(knearneigh(xy, k = n_neighbors))
    lw <- nb2listw(nb, style = "W")
    mi <- moran.test(z, lw)

    tibble(
      year       = yr,
      method     = method_name,
      moran_I    = mi$estimate["Moran I statistic"],
      expected_I = mi$estimate["Expectation"],
      z_score    = mi$statistic,
      p_value    = mi$p.value,
      n_cells    = length(z)
    )
  }

  bind_rows(
    compute_morans_date(sp_peak,  r_tmp_c, "Spline", yr),
    compute_morans_date(dlc_peak, r_tmp_c, "DLC",    yr)
  )
})

print("moran's i on spring green-up date, by year:")
print(morans_date_results, n = Inf)

## summary across years
morans_date_summary <- morans_date_results %>%
  group_by(method) %>%
  summarise(
    n_years        = n(),
    mean_moran_I   = round(mean(moran_I,  na.rm = TRUE), 3),
    median_moran_I = round(median(moran_I, na.rm = TRUE), 3),
    sd_moran_I     = round(sd(moran_I,    na.rm = TRUE), 3),
    min_moran_I    = round(min(moran_I,   na.rm = TRUE), 3),
    max_moran_I    = round(max(moran_I,   na.rm = TRUE), 3),
    .groups = "drop"
  )

print("summary of moran's i on spring green-up date:")
print(morans_date_summary)

## plot moran's i by year for both methods
p_morans_date <- morans_date_results %>%
  filter(!is.na(moran_I)) %>%
  ggplot(aes(x = year, y = moran_I, color = method, group = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_color_manual(values = method_colors) +
  labs(
    title    = "Spatial autocorrelation of spring green-up date",
    subtitle = "Moran's I by year -- higher = stronger spatial structure in phenology",
    x        = "Year",
    y        = "Moran's I",
    color    = "Method"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

p_morans_date

ggsave(file.path(base_path_savefigures, "Validation2_MoransI_GreenupDate_ByYear.png"),
       p_morans_date, width = 8, height = 5, dpi = 200)

saveRDS(morans_date_results,
        file.path(base_path_saveproducts, "Validation2_MoransI_GreenupDate.rds"))
saveRDS(morans_date_summary,
        file.path(base_path_saveproducts, "Validation2_MoransI_GreenupDate_Summary.rds"))


#### step 4: spatial maps of spring green-up date for example years ####

print("building spatial maps of spring green-up date...")

## years to map -- adjust based on the moran's i by year plot
plot_years <- c(2016, 2019, 2021, 2023)

for (yr in plot_years) {

  path <- file.path(base_path_rasters,
                    paste0("FULL_SplineMetrics_Raster_", yr, ".tif"))
  if (!file.exists(path)) { print(paste0(yr, " skipped")); next }

  r      <- rast(path)[[bands_date]]
  r_crop <- crop(r, md_vect)
  r_mask <- mask(r_crop, md_vect)
  rm(r, r_crop); gc()

  v <- values(r_mask)

  ## spline spring peak
  p_mat <- v[, 1:4, drop = FALSE]
  p_mat[p_mat < spring_start_jday | p_mat > spring_end_jday] <- NA
  sp_peak <- apply(p_mat, 1, function(x) {
    x <- x[!is.na(x)]; if (length(x) == 0) NA_real_ else min(x)
  })
  sp_peak[!is.na(v[, 5]) & v[, 5] >= rss_threshold] <- NA

  ## dlc spring peak
  dlc_peak <- v[, 6]
  dlc_peak[is.na(dlc_peak) |
             dlc_peak < spring_start_jday |
             dlc_peak > spring_end_jday   |
             (!is.na(v[, 7]) & v[, 7] >= rss_threshold)] <- NA

  rm(v, p_mat); gc()

  ## build rasters for plotting
  sp_rast  <- dlc_rast <- r_mask[[1]]
  values(sp_rast)  <- sp_peak
  values(dlc_rast) <- dlc_peak
  rm(r_mask); gc()

  ## convert to data frames for ggplot
  sp_df <- as.data.frame(sp_rast,   xy = TRUE) %>%
    rename(jday = 3) %>% mutate(method = "Spline")
  dlc_df <- as.data.frame(dlc_rast, xy = TRUE) %>%
    rename(jday = 3) %>% mutate(method = "DLC")

  plot_df <- bind_rows(sp_df, dlc_df) %>%
    filter(!is.na(jday)) %>%
    mutate(method = factor(method, levels = c("Spline", "DLC")))

  rm(sp_df, dlc_df, sp_rast, dlc_rast); gc()

  ## shared color scale: julian days 60-180 (spring window)
  p_map <- ggplot(plot_df, aes(x = x, y = y, fill = jday)) +
    geom_raster() +
    facet_wrap(~ method, ncol = 2) +
    scale_fill_gradientn(
      colors = c("#ffffcc", "#78c679", "#238443", "#004529"),
      limits = c(spring_start_jday, spring_end_jday),
      na.value = "white",
      name = "Spring green-up\n(Julian day)"
    ) +
    coord_equal() +
    labs(
      title    = paste0("Spring green-up date: ", yr),
      subtitle = "Earlier = darker green; white = no valid estimate (RSS >= 2)",
      x = "Longitude", y = "Latitude"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position  = "right",
      strip.text       = element_text(size = 12, face = "bold"),
      panel.spacing    = unit(0.5, "cm"),
      axis.text        = element_text(size = 8)
    )

  fname <- paste0("Validation2_SpringGreenupMap_", yr, ".png")
  ggsave(file.path(base_path_savefigures, fname),
         p_map, width = 10, height = 7, dpi = 200)
  plot(p_map)
  print(paste0("saved ", fname))

  rm(plot_df); gc()
}


#### step 5: local moran's i (lisa) maps for example years ####

## decomposes global moran's i into per-pixel cluster membership:
##   high-high = pixel and neighbors both have late green-up
##   low-low   = pixel and neighbors both have early green-up
##   high-low  = late pixel surrounded by early neighbors (spatial outlier)
##   low-high  = early pixel surrounded by late neighbors (spatial outlier)
##   ns        = not significant (p > 0.05)

print("building lisa cluster maps...")

lisa_plot_years <- c(2016, 2019, 2021, 2023)
## lisa_plot_years <- c(2015)

for (yr in lisa_plot_years) {

  path <- file.path(base_path_rasters,
                    paste0("FULL_SplineMetrics_Raster_", yr, ".tif"))
  if (!file.exists(path)) { print(paste0(yr, " skipped")); next }

  r      <- rast(path)[[bands_date]]
  r_crop <- crop(r, md_vect)
  r_mask <- mask(r_crop, md_vect)
  rm(r, r_crop); gc()

  v <- values(r_mask)

  ## spline spring peak
  p_mat <- v[, 1:4, drop = FALSE]
  p_mat[p_mat < spring_start_jday | p_mat > spring_end_jday] <- NA
  sp_peak <- apply(p_mat, 1, function(x) {
    x <- x[!is.na(x)]; if (length(x) == 0) NA_real_ else min(x)
  })
  sp_peak[!is.na(v[, 5]) & v[, 5] >= rss_threshold] <- NA

  ## dlc spring peak
  dlc_peak <- v[, 6]
  dlc_peak[is.na(dlc_peak) |
             dlc_peak < spring_start_jday |
             dlc_peak > spring_end_jday   |
             (!is.na(v[, 7]) & v[, 7] >= rss_threshold)] <- NA

  rm(v, p_mat); gc()

  ## template raster for coordinates
  r_tmp   <- rast(path)[[1]]
  r_tmp_c <- mask(crop(r_tmp, md_vect), md_vect)
  rm(r_tmp); gc()

  compute_lisa <- function(peak_vals, template_rast, method_name, yr,
                           n_neighbors = 8, subsample = 150000) {

    non_na <- which(!is.na(peak_vals))
    if (length(non_na) < 100) return(NULL)

    if (length(non_na) > subsample) {
      set.seed(42)
      non_na <- sample(non_na, subsample)
    }

    xy <- xyFromCell(template_rast, non_na)
    z  <- peak_vals[non_na]
    if (nrow(xy) != length(z)) return(NULL)

    nb <- knn2nb(knearneigh(xy, k = n_neighbors))
    lw <- nb2listw(nb, style = "W")
    local_mi <- localmoran(z, lw)

    z_scaled   <- scale(z)[, 1]
    lag_scaled <- lag.listw(lw, z_scaled)

    cluster <- case_when(
      local_mi[, "Pr(z != E(Ii))"] > 0.05  ~ "NS",
      z_scaled > 0 & lag_scaled > 0         ~ "High-High",
      z_scaled < 0 & lag_scaled < 0         ~ "Low-Low",
      z_scaled > 0 & lag_scaled < 0         ~ "High-Low",
      z_scaled < 0 & lag_scaled > 0         ~ "Low-High",
      TRUE                                   ~ "NS"
    )

    tibble(
      x       = xy[, 1],
      y       = xy[, 2],
      lisa_I  = local_mi[, "Ii"],
      p_value = local_mi[, "Pr(z != E(Ii))"],
      cluster = cluster,
      method  = method_name,
      year    = yr
    )
  }

  sp_lisa  <- compute_lisa(sp_peak,  r_tmp_c, "Spline", yr)
  dlc_lisa <- compute_lisa(dlc_peak, r_tmp_c, "DLC",    yr)

  rm(r_tmp_c, sp_peak, dlc_peak); gc()

  if (is.null(sp_lisa) || is.null(dlc_lisa)) {
    print(paste0(yr, " -- insufficient data for lisa"))
    next
  }

  lisa_df <- bind_rows(sp_lisa, dlc_lisa) %>%
    mutate(
      cluster = factor(cluster,
                       levels = c("High-High", "Low-Low",
                                  "High-Low",  "Low-High", "NS")),
      method  = factor(method, levels = c("Spline", "DLC"))
    )

  ## lisa cluster colours: standard convention
  lisa_colors <- c(
    "High-High" = "#d7191c",
    "Low-Low"   = "#2c7bb6",
    "High-Low"  = "#fdae61",
    "Low-High"  = "#abd9e9"
  )

  ## resolution of the background gray layer -- match modis 250m pixel size
  ## in degrees (approximately 0.002 degrees at utah latitudes)
  res_deg <- 0.002

  ## full background grid from all non-na pixels (both methods combined)
  ## so gray ns fills the whole valid extent
  all_xy <- bind_rows(
    tibble(x = sp_lisa$x,  y = sp_lisa$y),
    tibble(x = dlc_lisa$x, y = dlc_lisa$y)
  ) %>% distinct()

  background_df <- crossing(
    method = factor(c("Spline", "DLC"), levels = c("Spline", "DLC"))
  ) %>%
    mutate(data = list(all_xy)) %>%
    tidyr::unnest(data)

  ## significant pixels only
  sig_df <- lisa_df %>%
    filter(cluster != "NS") %>%
    mutate(cluster = factor(cluster,
                            levels = c("High-High", "Low-Low",
                                       "High-Low",  "Low-High")))

  p_lisa <- ggplot() +
    ## gray background = ns pixels (all valid pixels)
    geom_point(data = background_df,
               aes(x = x, y = y),
               color = "gray88", size = 0.05, shape = 15) +
    ## significant lisa clusters on top
    geom_point(data = sig_df,
               aes(x = x, y = y, color = cluster),
               size = 0.15, shape = 15) +
    facet_wrap(~ method, ncol = 2) +
    scale_color_manual(values = lisa_colors,
                       name   = "LISA cluster",
                       guide  = guide_legend(override.aes = list(size = 4,
                                                                 shape = 15))) +
    coord_equal() +
    labs(
      title    = paste0("LISA clusters of spring green-up date: ", yr),
      subtitle = "Red = late green-up clusters; Blue = early green-up clusters;\nOrange/Lt.Blue = spatial outliers; Gray = not significant",
      x = "Longitude", y = "Latitude"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "right",
      strip.text      = element_text(size = 12, face = "bold"),
      panel.spacing   = unit(0.8, "cm"),
      axis.text       = element_text(size = 8),
      panel.grid      = element_blank()
    )

  fname <- paste0("Validation2_LISA_SpringGreenup_", yr, ".png")
  ggsave(file.path(base_path_savefigures, fname),
         p_lisa, width = 12, height = 7, dpi = 250)
  print(paste0("saved ", fname))
  p_lisa
  rm(lisa_df); gc()
}

print("validation 2 complete")
