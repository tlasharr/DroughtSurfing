## functions for extracting and classifying green-up metrics from spline raster outputs
## and matching them to gps locations
##
## functions are organized into three groups:
##   1. data loading helpers
##   2. raster extraction
##   3. pulse inventory and surfing metric classification


########## data loading helpers

## loads rasters for a given region into a named list by year
## skips years where the file doesn't exist yet
load_greenup_rasters <- function(base_path, region, years) {
  
  rast_list <- setNames(
    lapply(years, function(yr) {
      path <- file.path(base_path, region, paste0(region, "_SplineMetrics_Raster_", yr, ".tif"))
      if (file.exists(path)) rast(path) else NULL # skip years that haven't been processed yet
    }),
    as.character(years)
  )
  
  Filter(Negate(is.null), rast_list) # drop any nulls from missing files
}

## loads a population's gps and migration data from standardized file paths
load_population <- function(base_path, herd_id) {
  list(
    deer = readRDS(file.path(base_path, herd_id, "FinalCleanData",
                             paste0(herd_id, "_locs_formatted.rds"))),
    mig  = readRDS(file.path(base_path, herd_id, "FinalCleanData",
                             paste0(herd_id, "_migtime_formatted.rds")))
  )
}


########## raster extraction

## for a single population, extracts green-up metrics at each gps location,
## joins migration timing, assigns season, and then runs pulse classification
extract_greenup_for_deer <- function(deer_sf, mig_df, rast_list, bands, rss_threshold = 2) {
  
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
    
    ## classify pulse status and compute all surfing metrics for each fix
    add_pulse_metrics(rss_threshold = rss_threshold)
}


########## pulse inventory

## counts the number of non-na green-up pulses for a single fix
count_pulses <- function(peak_1, peak_2, peak_3, peak_4) {
  sum(!is.na(c(peak_1, peak_2, peak_3, peak_4)))
}

## assembles a small dataframe of valid pulses for a single fix
## returns one row per valid pulse with start, peak, end, duration, and intensity
## returns an empty dataframe if no pulses exist at that pixel
inventory_pulses <- function(peak_1,   peak_2,   peak_3,   peak_4,
                             start_1,  start_2,  start_3,  start_4,
                             end_1,    end_2,    end_3,    end_4,
                             change_1, change_2, change_3, change_4) {
  
  ## bundle into vectors for clean iteration
  peaks   <- c(peak_1,   peak_2,   peak_3,   peak_4)
  starts  <- c(start_1,  start_2,  start_3,  start_4)
  ends    <- c(end_1,    end_2,    end_3,    end_4)
  changes <- c(change_1, change_2, change_3, change_4)
  
  ## only keep pulses where the peak is non-na -- peak is the anchor
  valid <- !is.na(peaks)
  
  if (!any(valid)) return(data.frame()) # early exit if no pulses at all
  
  data.frame(
    pulse_idx       = which(valid),
    peak_day        = peaks[valid],
    start_day       = starts[valid],
    end_day         = ends[valid],
    pulse_intensity = changes[valid],
    pulse_duration  = ends[valid] - starts[valid] # width of the green-up window in days
  )
}


########## window membership

## checks whether a given julian day falls inside any green-up pulse window
## returns the pulse index of the matching window, or na if none
## if the day falls in two overlapping windows, takes the more intense pulse
check_window_membership <- function(jday, pulse_inventory) {
  
  if (nrow(pulse_inventory) == 0) return(NA_integer_) # no pulses, nothing to check
  
  ## find pulses where jday falls between start and end
  in_window <- with(pulse_inventory,
                    !is.na(start_day) & !is.na(end_day) &
                      jday >= start_day & jday <= end_day)
  
  if (!any(in_window)) return(NA_integer_) # not inside any window
  
  ## if in multiple overlapping windows, take the one with the largest greenup_change
  matching <- pulse_inventory[in_window, ]
  matching$pulse_idx[which.max(matching$pulse_intensity)]
}


########## proximity to nearest window boundary

## when a deer is not inside any pulse window, find the nearest window edge
## returns a list with:
##   days_to_next_start  -- days until the next pulse window opens (na if none upcoming)
##   days_since_last_end -- days since the last pulse window closed (na if none past)
##   proximity_status    -- categorical label for the deer's phenological situation
nearest_window_boundary <- function(jday, pulse_inventory) {
  
  ## default return for no-pulse pixels
  empty <- list(days_to_next_start  = NA_real_,
                days_since_last_end = NA_real_,
                proximity_status    = "no_pulse")
  
  if (nrow(pulse_inventory) == 0) return(empty)
  
  ## upcoming pulse starts -- windows that haven't opened yet
  upcoming_starts <- pulse_inventory$start_day[!is.na(pulse_inventory$start_day) &
                                                 pulse_inventory$start_day > jday]
  
  ## past pulse ends -- windows that have already closed
  past_ends <- pulse_inventory$end_day[!is.na(pulse_inventory$end_day) &
                                         pulse_inventory$end_day < jday]
  
  days_to_next    <- if (length(upcoming_starts) > 0) min(upcoming_starts) - jday else NA_real_
  days_since_last <- if (length(past_ends) > 0) jday - max(past_ends) else NA_real_
  
  ## assign a status label based on what's available
  status <- case_when(
    !is.na(days_to_next) & is.na(days_since_last)  ~ "pre_pulse",       # before any pulse has occurred
    is.na(days_to_next)  & !is.na(days_since_last) ~ "post_pulse",      # after all pulses have ended
    !is.na(days_to_next) & !is.na(days_since_last) ~ "between_pulses",  # between two pulses
    TRUE                                             ~ "no_pulse"
  )
  
  list(days_to_next_start  = days_to_next,
       days_since_last_end = days_since_last,
       proximity_status    = status)
}


########## main classification function

## takes a single fix's extracted raster values and julian day
## runs all three stages (inventory, window membership, boundary proximity)
## returns a named list of pulse metrics for that fix
## designed to be called inside add_pulse_metrics() row by row
classify_pulse_status <- function(jday,
                                  peak_1,   peak_2,   peak_3,   peak_4,
                                  start_1,  start_2,  start_3,  start_4,
                                  end_1,    end_2,    end_3,    end_4,
                                  change_1, change_2, change_3, change_4,
                                  rss_spline, rss_threshold = 2) {
  
  ## return all na if spline fit is too poor to trust
  if (!is.na(rss_spline) && rss_spline > rss_threshold) {
    return(list(
      n_pulses            = NA_integer_,
      pulse_status        = "poor_fit",
      active_pulse_idx    = NA_integer_,
      days_into_window    = NA_real_,
      window_duration     = NA_real_,
      pulse_intensity     = NA_real_,
      days_from_peak      = NA_real_,
      days_to_next_start  = NA_real_,
      days_since_last_end = NA_real_,
      proximity_status    = NA_character_
    ))
  }
  
  ## stage 1: build pulse inventory for this pixel
  inv      <- inventory_pulses(peak_1,   peak_2,   peak_3,   peak_4,
                               start_1,  start_2,  start_3,  start_4,
                               end_1,    end_2,    end_3,    end_4,
                               change_1, change_2, change_3, change_4)
  n_pulses <- nrow(inv) # 0 if no valid pulses detected at this pixel
  
  ## early exit for no-pulse pixels
  if (n_pulses == 0) {
    return(list(
      n_pulses            = 0L,
      pulse_status        = "no_pulse",
      active_pulse_idx    = NA_integer_,
      days_into_window    = NA_real_,
      window_duration     = NA_real_,
      pulse_intensity     = NA_real_,
      days_from_peak      = NA_real_,
      days_to_next_start  = NA_real_,
      days_since_last_end = NA_real_,
      proximity_status    = "no_pulse"
    ))
  }
  
  ## stage 2: check if deer is inside a green-up window right now
  active_idx <- check_window_membership(jday, inv)
  
  if (!is.na(active_idx)) {
    
    ## deer is inside a pulse window -- fill in all in-window metrics
    active_pulse <- inv[inv$pulse_idx == active_idx, ]
    
    return(list(
      n_pulses            = n_pulses,
      pulse_status        = "in_greenup",
      active_pulse_idx    = active_idx,
      days_into_window    = jday - active_pulse$start_day,  # 0 = first day of window
      window_duration     = active_pulse$pulse_duration,
      pulse_intensity     = active_pulse$pulse_intensity,
      days_from_peak      = jday - active_pulse$peak_day,   # signed: negative = before peak
      days_to_next_start  = NA_real_,
      days_since_last_end = NA_real_,
      proximity_status    = NA_character_
    ))
  }
  
  ## stage 3: deer is not in any window -- find nearest boundary
  boundary <- nearest_window_boundary(jday, inv)
  
  ## find the nearest peak outside any window, keeping the sign
  ## negative = deer is before the peak, positive = deer is after
  peak_deltas      <- jday - inv$peak_day
  nearest_peak_row <- inv[which.min(abs(peak_deltas)), ]
  
  return(list(
    n_pulses            = n_pulses,
    pulse_status        = boundary$proximity_status,
    active_pulse_idx    = NA_integer_,
    days_into_window    = NA_real_,
    window_duration     = NA_real_,
    pulse_intensity     = NA_real_,
    days_from_peak      = peak_deltas[which.min(abs(peak_deltas))], # nearest peak even if outside window
    days_to_next_start  = boundary$days_to_next_start,
    days_since_last_end = boundary$days_since_last_end,
    proximity_status    = boundary$proximity_status
  ))
}


########## wrapper: apply classify_pulse_status across a full dataframe

## runs classify_pulse_status row-wise on a dataframe that already has raster
## values extracted. expects all greenup_peak/start/end/change columns to be present.
## returns the input dataframe with pulse metric columns appended.
## intentionally uses a plain for loop -- the conditional logic across 0-4 pulses
## is complex enough that clarity matters more than speed here
add_pulse_metrics <- function(df, rss_threshold = 2) {
  
  message("classifying pulse status for ", nrow(df), " fixes...")
  
  results <- vector("list", nrow(df))
  
  for (i in seq_len(nrow(df))) {
    results[[i]] <- classify_pulse_status(
      jday      = df$jday[i],
      peak_1    = df$greenup_peak_1[i],   peak_2   = df$greenup_peak_2[i],
      peak_3    = df$greenup_peak_3[i],   peak_4   = df$greenup_peak_4[i],
      start_1   = df$greenup_start_1[i],  start_2  = df$greenup_start_2[i],
      start_3   = df$greenup_start_3[i],  start_4  = df$greenup_start_4[i],
      end_1     = df$greenup_end_1[i],    end_2    = df$greenup_end_2[i],
      end_3     = df$greenup_end_3[i],    end_4    = df$greenup_end_4[i],
      change_1  = df$greenup_change_1[i], change_2 = df$greenup_change_2[i],
      change_3  = df$greenup_change_3[i], change_4 = df$greenup_change_4[i],
      rss_spline    = df$rss_spline[i],
      rss_threshold = rss_threshold
    )
  }
  
  bind_cols(df, bind_rows(results))
}