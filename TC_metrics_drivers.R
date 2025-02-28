#------------------------------------------------------------------------------
# Compute TC metrics 
#------------------------------------------------------------------------------
compute_tc_metrics <- function(data) {
  # Compute the number of unique event (SID) per year for intensity = "cat35"
  cat35_count <- data %>%
    filter(intensity == "cat35") %>%
    group_by(modified_year, BASIN) %>%
    summarise(Cat35_count = n_distinct(SID), .groups = 'drop')
  
  # Compute the number of unique event (SID) per year for intensity = "cat15" and "cat35"
  cat15_count <- data %>%
    filter(intensity %in% c("cat15", "cat35")) %>%
    group_by(modified_year, BASIN) %>%
    summarise(Cat15_count = n_distinct(SID), .groups = 'drop')
  
  # Compute the annual average ACE (summed over each storm first)
  Check_Ace <- data %>%
    arrange(date) %>%  # Ensure the data is ordered by date
    mutate(hour = format(date, "%H")) %>%
    filter(hour %in% c("00", "06", "12", "18")) %>%  # Keep only 6-hour intervals
    filter(USA_WIND >= 34) %>%  # Filter to include only storms with wind speed >= 34 knots
    mutate(ACE = (USA_WIND^2) * 10^(-4))
  
  ACE_per_storm_year <- Check_Ace %>%
    group_by(SID, modified_year, BASIN) %>%
    summarise(total_ACE = sum(ACE, na.rm = TRUE), .groups = 'drop')
  
  ACE_count <- ACE_per_storm_year %>%
    group_by(modified_year, BASIN) %>%
    summarise(avg_ACE = mean(total_ACE, na.rm = TRUE), .groups = 'drop')
  
  # Compute the Latitude of maximum wind for each storm
  lat_maxwind_per_storm <- data %>%
    group_by(SID, modified_year, BASIN) %>%
    filter(!all(is.na(USA_WIND))) %>%  # Remove groups where all USA_WIND values are NA
    filter(USA_WIND == max(USA_WIND, na.rm = TRUE)) %>%  # Get the rows with max USA_WIND per storm
    slice(1) %>%  # In case there are ties, select the first row
    summarise(lat_maxwind = LAT, .groups = 'drop')
  
  lat_maxwind_per_year <- lat_maxwind_per_storm %>%
    group_by(modified_year, BASIN) %>%
    summarise(avg_lat_maxwind = mean(lat_maxwind, na.rm = TRUE), .groups = 'drop')
  
  # Join all metrics together
  TC_metrics <- ACE_count %>%
    left_join(lat_maxwind_per_year, by = c("modified_year", "BASIN")) %>%
    left_join(cat15_count, by = c("modified_year", "BASIN")) %>%
    left_join(cat35_count, by = c("modified_year", "BASIN"))
  
  # Add the Hemisphere information
  TC_metrics <- TC_metrics %>%
    mutate(HEMISPHERE = if_else(BASIN %in% c("SI", "SP"), "SH", "NH"))
  
  # Rename the columns as required
  names(TC_metrics) <- c("YR", "BASIN", "avg_ACE", "avg_lat_maxwind", "Cat15_count", "Cat35_count", "HEMISPHERE")
  
  return(TC_metrics)
}

#------------------------------------------------------------------------------
# Compute TC drivers by fetching and processing ENSO and AMO data
#------------------------------------------------------------------------------

extract_tc_drivers <- function(enso_url, amo_url) {
  
  # Fetch and process ENSO data
  enso <- tryCatch({
    enso_data <- read_table(enso_url, col_names = TRUE) %>%
      as_tibble()  # Ensure data is in tibble format
    
    # Rename columns
    colnames(enso_data) <- c("YR", "MON", "NINO1+2", "ANOM1+2", "NINO3", "ANOM3", "NINO4", "ANOM4", "NINO3.4", "ANOM3.4")
    
    # Filter and select specific columns
    enso_data <- enso_data %>%
      dplyr::filter(YR <= 2022) %>%  # Filter years of interest
      dplyr::select(YR, MON, ANOM3.4)  # Select relevant columns
    
    message("ENSO data successfully processed.")

    # Explicitly return the processed ENSO data
    enso_data
    
  }, error = function(e) {
    message("Error fetching ENSO data: ", e)
    return(NULL)
  })
  
  if (is.null(enso)) {
    message("ENSO data could not be fetched or processed.")
    return(NULL)  # Exit if ENSO data can't be fetched
  }
  
  # Fetch and process AMO data using your specific method
  amo <- tryCatch({
    # Read the data from the AMO URL
    amo_data <- read_table(amo_url, skip = 127, col_names = FALSE)
    
    # Convert to a data frame
    amo_data <- as.data.frame(amo_data)
    
    # Remove the last 5 rows of the data
    amo_data <- amo_data[1:(nrow(amo_data) - 5), ]
    
    # Rename columns for year and months
    names(amo_data) <- c('YR', 1:12)
    
    # Convert character columns to numeric if needed
    amo_data %<>% mutate_if(is.character, as.numeric)
    
    # Reshape the data: pivot longer to have months as a column and values as "AMO"
    amo_data <- amo_data %>%
      pivot_longer(
        cols = -YR,          # All columns except 'YR'
        names_to = "MON",    # New column name for months
        values_to = "AMO"    # New column name for values
      )
    
    # Ensure that month values are numeric
    amo_data %<>% mutate_if(is.character, as.numeric)
    
    # Print a success message and display a sample of the data
    print("AMO data successfully processed.")
    
    # Explicitly return the processed AMO data
    amo_data
    
  }, error = function(e) {
    # Handle errors by printing a message and returning NULL
    message("Error fetching AMO data: ", e)
    NULL
  })
  
  # Check if AMO data is successfully processed
  if (is.null(amo)) {
    message("AMO data could not be fetched or processed.")
  } else {
    print("AMO data ready for merging.")
  }
  
  # Merge ENSO and AMO datasets by year (YR) and month (MON)
  merged_data <- enso %>% 
    full_join(amo, by = c("YR", "MON"))
  
  # Create a Date column in the merged dataset
  merged_data$Date <- as.Date(paste(merged_data$YR, merged_data$MON, "01", sep = "-"), "%Y-%m-%d")
  
  return(merged_data)
}

#------------------------------------------------------------------------------
# Count tracks per grid cells
#------------------------------------------------------------------------------

compute_TC_counts <- function(grid_size, basin, ALL_df, nh_seasonal, sh_seasonal) {
  # Determine the hemisphere based on the basin
  hemisphere <- ifelse(basin %in% c("WP", "EP", "NA", "NI"), "NH", "SH")
  
  # Select the relevant data from ALL_df based on the basin
  track_location_ <- ALL_df %>%
    dplyr::select(SID, BASIN, LAT, LON, USA_WIND, modified_year, month) %>%
    dplyr::filter(BASIN == basin)
  
  # Choose the appropriate seasonal data based on the hemisphere
  seasonal_data <- if (hemisphere == "NH") {
    nh_seasonal %>%
      filter(HEMISPHERE == 'NH') %>%
      mutate(NINO34 = ifelse(seasonal_nino34_index >= 0.5, "Nino",
                             ifelse(seasonal_nino34_index <= -0.5, "Nina", 
                                    "Neutral")))
  } else {
    sh_seasonal %>%
      filter(HEMISPHERE == 'SH') %>%
      mutate(NINO34 = ifelse(seasonal_nino34_index >= 0.5, "Nino",
                             ifelse(seasonal_nino34_index <= -0.5, "Nina", 
                                    "Neutral")))
  }
  
  # Join track location data with the seasonal data
  df_tracks_enso <- track_location_ %>%
    full_join(seasonal_data, by = c("modified_year" = "YR"))
  
  # Filter the data to exclude "Neutral" NINO34 conditions
  df_tracks_filtered <- df_tracks_enso %>%
    filter(NINO34 != "Neutral")
  
  # Create a grid with the specified resolution and count tracks
  grid_counts <- df_tracks_filtered %>%
    mutate(
      lon_bin = floor(LON / grid_size) * grid_size,
      lat_bin = floor(LAT / grid_size) * grid_size
    ) %>%
    group_by(lon_bin, lat_bin, NINO34, modified_year) %>%
    dplyr::summarize(track_count = n_distinct(SID), .groups = 'drop') %>%
    group_by(lon_bin, lat_bin, NINO34) %>%
    dplyr::summarize(avg_track_count = mean(track_count, na.rm = TRUE), .groups = 'drop') %>%
    ungroup()
  
  # Prepare data for the difference plot (Nino - Nina)
  grid_counts_diff <- grid_counts %>%
    tidyr::spread(NINO34, avg_track_count) %>%
    mutate(nino_nina_diff = Nino - Nina)
  
  # Return the resulting data
  return(list(grid_counts = grid_counts, grid_counts_diff = grid_counts_diff))
}

#------------------------------------------------------------------------------
# Compute seasonal indices
#------------------------------------------------------------------------------
compute_seasonal_indices <- function(tc_drivers, nh_season_months, sh_season_months) {
  # Filter and aggregate for the Northern Hemisphere
  nh_seasonal <- tc_drivers %>%
    filter(MON %in% nh_season_months) %>%
    group_by(YR) %>%
    dplyr::summarize(seasonal_nino34_index = mean(ANOM3.4, na.rm = TRUE),
                     seasonal_amo_index = mean(AMO, na.rm = TRUE)) %>%
    mutate(HEMISPHERE = "NH")
  
  # Renaming columns for consistency
  names(nh_seasonal) <- c("YR", "seasonal_nino34_index", "seasonal_amo_index", "HEMISPHERE")
  
  # Filter and aggregate for the Southern Hemisphere
  sh_seasonal <- tc_drivers %>%
    filter(MON %in% sh_season_months) %>%
    mutate(season_year = if_else(MON %in% c(10, 11, 12), YR + 1, YR)) %>%  # Adjust season year for months Oct-Dec
    group_by(season_year) %>%
    dplyr::summarize(seasonal_nino34_index = mean(ANOM3.4, na.rm = TRUE),
                     seasonal_amo_index = mean(AMO, na.rm = TRUE)) %>%
    mutate(HEMISPHERE = "SH")
  
  # Renaming columns for consistency
  names(sh_seasonal) <- c("YR", "seasonal_nino34_index", "seasonal_amo_index", "HEMISPHERE")
  
  # Return the results as a list
  return(list(nh_seasonal = nh_seasonal, sh_seasonal = sh_seasonal))
}

#------------------------------------------------------------------------------
# Compute correlation matrix
#------------------------------------------------------------------------------
library(dplyr)
library(tidyr)

process_tc_metrics <- function(TC_metrics, seasonal_indices, cutoff_year = 2022) {
  
  # Step 1: Merge TC metrics with seasonal indices
  df_index_metrics <- TC_metrics %>%
    full_join(seasonal_indices, by = c("YR", "HEMISPHERE"))
  
  # Step 2: Handle missing values and filter based on the year
  df_index_metrics <- df_index_metrics %>%
    mutate(Cat15_count = ifelse(is.na(Cat15_count), 0, Cat15_count),
           Cat35_count = ifelse(is.na(Cat35_count), 0, Cat35_count),
           avg_ACE = ifelse(is.na(avg_ACE), 0, avg_ACE)) %>%
    filter(YR <= cutoff_year) %>%
    as.data.frame()
  
  # Step 3: Drop "seasonal_" prefix from column names
  df_index_metrics <- df_index_metrics %>%
    rename_with(~ sub("^seasonal_", "", .x))
  
  # Step 4: Pivot the data for different metrics
  
  # Pivot the Cat15_count column
  cat15_wide <- df_index_metrics %>%
    dplyr::select(YR, BASIN, Cat15_count) %>%
    pivot_wider(names_from = BASIN, values_from = Cat15_count, names_prefix = "Cat15_count_")
  
  # Pivot the Cat35_count column
  cat35_wide <- df_index_metrics %>%
    dplyr::select(YR, BASIN, Cat35_count) %>%
    pivot_wider(names_from = BASIN, values_from = Cat35_count, names_prefix = "Cat35_count_")
  
  # Pivot the avg_ACE column
  ACE_wide <- df_index_metrics %>%
    dplyr::select(YR, BASIN, avg_ACE) %>%
    pivot_wider(names_from = BASIN, values_from = avg_ACE, names_prefix = "avg_ACE_")
  
  # Pivot the nino34_index column
  nino34_wide <- df_index_metrics %>%
    dplyr::select(YR, HEMISPHERE, nino34_index) %>%
    distinct(YR, HEMISPHERE, .keep_all = TRUE) %>%
    pivot_wider(names_from = HEMISPHERE, values_from = nino34_index, names_prefix = "nino34_index_")
  
  # Pivot the amo_index column
  amo_wide <- df_index_metrics %>%
    dplyr::select(YR, HEMISPHERE, amo_index) %>%
    distinct(YR, HEMISPHERE, .keep_all = TRUE) %>%
    pivot_wider(names_from = HEMISPHERE, values_from = amo_index, names_prefix = "amo_index_")
  
  # Step 5: Merge the wide tables for frequency and ACE correlation data
  
  # Create the data frame for frequency correlations
  df_freq <- cat15_wide %>%
    full_join(nino34_wide, by = "YR") %>%
    full_join(amo_wide, by = "YR")
  
  df_for_freq_corr <- df_freq %>% dplyr::select(-YR)
  
  # Create the data frame for ACE correlations
  df_ace <- ACE_wide %>%
    full_join(nino34_wide, by = "YR") %>%
    full_join(amo_wide, by = "YR")
  
  df_for_ace_corr <- df_ace %>% dplyr::select(-YR)
  
  # Return the two data frames as a list
  return(list(
    df_for_freq_corr = df_for_freq_corr,
    df_for_ace_corr = df_for_ace_corr
  ))
}

# Function to compute correlations
cors <- function(df) {
  M <- Hmisc::rcorr(as.matrix(df), type = "pearson")
  Mdf <- map(M, ~data.frame(.x))
  return(Mdf)
}

# Function to format correlations for plotting
formatted_cors <- function(df){
  cors(df) %>%
    map(~rownames_to_column(.x, var="measure1")) %>%
    map(~pivot_longer(.x, -measure1, names_to = "measure2", values_to = "value")) %>%
    bind_rows(.id = "id") %>%
    pivot_wider(names_from = id, values_from = value) %>%
    rename(p = P) %>%
    mutate(sig_p = ifelse(p < .05, TRUE, FALSE),
           p_if_sig = ifelse(sig_p, p, NA),
           r_if_sig = ifelse(sig_p, r, NA)) 
}
