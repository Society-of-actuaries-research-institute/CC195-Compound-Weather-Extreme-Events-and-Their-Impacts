###############################################################################
# Tropical cyclone compound event project
#
# Description:
# This R script processes and analyzes the IBTrACS shapefile data, which contains
# global tropical cyclone (TC) track information. The script uses the 'sf' package 
# for geospatial data handling, and provides insights on cyclone patterns 
# such as frequency, intensity, and geographic distribution.
#
# Key Features:
# - Data Source: IBTrACS Version 04 shapefile (global tropical cyclones)
#     In the data/TropicalCyclone folder
# - Data visualization functions : TC_figures.R
# - Data manipulation functions : TC_metrics_drivers.R
#
# Usage:
# 1. Download the libraries
# 2. Run the script in R RStudio
#
# Date: [10/08/2024]
#
###############################################################################
#------------------------------------------------------------------------------
#--- Step 0 : LOAD LIBRARIES
#------------------------------------------------------------------------------
library(data.table)
library(dplyr)
library(fields)
library(ggplot2)
library(magrittr)
library(patchwork)
library(raster)
library(readr)
library(rworldmap)
library(sf)
library(sp)
library(statmod)
library(terra)
library(tidyr)
library(viridis)


source("TC_metrics_drivers.R")
source('TC_figures.R')

#------------------------------------------------------------------------------
#--- Step 1 : Read tropical cyclone data Ibtracs for all BASINS

# You can download the data in a shapefile format : https://www.ncei.noaa.gov/products/international-best-track-archive
# Kenneth R. Knapp, Scott Applequist, Howard J. Diamond, James P. Kossin, Michael Kruk, and Carl Schreck (2010). NCDC International Best Track Archive for Climate Stewardship (IBTrACS) Project, Version 3. [indicate subset used]. NOAA National Centers for Environmental Information. 
# Knapp, K. R., M. C. Kruk, D. H. Levinson, H. J. Diamond, and C. J. Neumann, 2010: The International Best Track Archive for Climate Stewardship (IBTrACS): Unifying tropical cyclone best track data. Bulletin of the American Meteor. Society, 91, 363-376.
#------------------------------------------------------------------------------
getwd()
ALL   <- st_read("../data/TropicalCyclone/IBTrACS.ALL.list.v04r01.lines.shp")

#------------------------------------------------------------------------------
#--- Step 2 : Format the spatial data frame
#------------------------------------------------------------------------------

# Set date format
ALL$date <- as.POSIXct(ALL$ISO_TIME, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
# Filter years of interest
ALL_filtered <- ALL %>% dplyr::filter(year >=1982 & year <= 2023)
# BASINS = EP NA NI SI SP WP. Remove the "SA" South Atlantic basin that is not relevant 
ALL_filtered %<>% dplyr::filter(BASIN != "SA")

# Change format to data.frame
ALL_df <- st_drop_geometry(ALL_filtered)

# Classification of TC Categories : Cat1+, Cat3+, ALL
ALL_df %<>%   mutate(intensity = case_when(
       USA_WIND >= 96 ~ "cat35",
       USA_WIND >= 64 ~ "cat15",
       TRUE ~ "ALL"))  

# Add a south hemisphere TC season October to September the following year
ALL_df %<>% mutate(modified_year = if_else(BASIN %in% c("SI", "SP") & month %in% c(10, 11, 12), year + 1, year))

#------------------------------------------------------------------------------
#--- Step 3 : Compute TC metrics (frequency, ACE, Latitude of max wind)
#------------------------------------------------------------------------------

TC_metrics <- compute_tc_metrics(ALL_df)

#------------------------------------------------------------------------------
#--- Step 4 : Compute TC drivers ENSO & AMO over north and south hemisphere
#------------------------------------------------------------------------------

enso_url <- "https://www.cpc.ncep.noaa.gov/data/indices/sstoi.indices"
amo_url <- "https://psl.noaa.gov/data/correlation/amon.us.long.data"

tc_drivers <- extract_tc_drivers(enso_url, amo_url)

# Plot time series
timeseries_nino34(tc_drivers,Date)


# Define the Northern Hemisphere and Southern Hemisphere cyclone season months
nh_season_months <- 6:11                       # Jun-Nov for the Northern Hemisphere
sh_season_months <- c(10, 11, 12, 1, 2, 3, 4)  # Oct-Apr for the Southern Hemisphere

# Call the function with the defined months and your data
result <- compute_seasonal_indices(tc_drivers, nh_season_months, sh_season_months)

# Access the results
nh_seasonal <- result$nh_seasonal
sh_seasonal <- result$sh_seasonal

seasonal_indices <- rbind(nh_seasonal,sh_seasonal)

#--------------------------------------------------------------------------------------------
# Compute the number of tracks per grid_size (= 2°) in the chosen basin (NI,SI,NA,WP,EP,SP)
#--------------------------------------------------------------------------------------------

result <- compute_TC_counts(grid_size = 2,  # chose grid resolution
                              basin = "NI", # chose basin
                              ALL_df = ALL_df, 
                              nh_seasonal = nh_seasonal, 
                              sh_seasonal = sh_seasonal)

grid_counts <- result$grid_counts
grid_counts_diff <- result$grid_counts_diff


#-- Display the results
#--- Define the domain for your figure (add+ or remove- a °Lat: y and °Lon: x)
world <- map_data("world")
lat_min <- min(grid_counts$lat_bin) #-+ y
lat_max <- max(grid_counts$lat_bin) #-+ y
lon_min <- min(grid_counts$lon_bin) #-+ x
lon_max <- max(grid_counts$lon_bin) #-+ x
world_expanded <- world %>%  filter(long >= lon_min & long <= lon_max & lat >= lat_min & lat <= lat_max)
ggplot() + geom_polygon(data = world_expanded, aes(x = long, y = lat, group = group), fill = NA, color = "black")

#--- Display thre results
plot_nino_nina_difference(grid_counts = grid_counts,
                          grid_counts_diff = grid_counts_diff,
                          world_expanded = world_expanded,
                          grid_size = 2,
                          lon_min = -180, lon_max = 180,
                          lat_min = -90, lat_max = 90)

#--------------------------------------------------------------------------------------------
# Merge TC Metrics and Seasonal Indices and Plot the time series of TC metrics and drivers
#--------------------------------------------------------------------------------------------

df_index_metrics <- TC_metrics %>%
  full_join(seasonal_indices, by = c("YR", "HEMISPHERE"))
df_index_metrics %<>% dplyr::filter(YR <= 2022)

# Replace NA values in Cat15_count and Cat35_count with 0
df_index_metrics %<>%
  mutate(Cat15_count = ifelse(is.na(Cat15_count), 0, Cat15_count),
         Cat35_count = ifelse(is.na(Cat35_count), 0, Cat35_count))
df_index_metrics <- as.data.frame(df_index_metrics)

plot_enso_freq(df_index_metrics,YR)
plot_enso_ACE(df_index_metrics,YR)

primary_y_limits <- c(-1, 1) # Set the y-axis limits
plot_AMO_freq(df_index_metrics,YR)
plot_AMO_ACE(df_index_metrics,YR)

#--------------------------------------------------------------------------------------------
# Multivariate XY plot index and normalized (remove mean) 
#--------------------------------------------------------------------------------------------

# Compute correlations for each basin
correlations <- df_index_metrics %>%
  group_by(BASIN) %>%
  dplyr::summarize(
    cor_enso_Cat15 = cor(seasonal_nino34_index, Cat15_count, method = "pearson"),
    cor_enso_Cat35 = cor(seasonal_nino34_index, Cat35_count, method = "pearson"),
    cor_enso_ACE   = cor(seasonal_nino34_index, avg_ACE, method = "pearson"),
    cor_enso_maxlat= cor(seasonal_nino34_index, avg_lat_maxwind, method = "pearson"),    
    
    cor_amo_Cat15  = cor(seasonal_amo_index, Cat15_count, method = "pearson"),
    cor_amo_Cat35  = cor(seasonal_amo_index, Cat35_count, method = "pearson"),
    cor_amo_ACE    = cor(seasonal_amo_index, avg_ACE, method = "pearson"),
    cor_amo_maxlat= cor(seasonal_amo_index, avg_lat_maxwind, method = "pearson"), 
    cor_ACE_Cat15  = cor(avg_ACE, Cat15_count, method = "pearson")
  )
df_index_metrics <- df_index_metrics %>%
  left_join(correlations, by = "BASIN")


# Compute the anomalies of the Cat15_count
df_index_metrics %<>%
  group_by(BASIN) %>%
  mutate(Cat15_mean = mean(Cat15_count, na.rm = TRUE),
         Cat15_anomaly = Cat15_count - Cat15_mean,
         Cat35_mean = mean(Cat35_count, na.rm = TRUE),
         Cat35_anomaly = Cat35_count - Cat35_mean,
         ACE_mean = mean(avg_ACE, na.rm = TRUE),
         ACE_anomaly = avg_ACE - ACE_mean,
         Maxlat_mean = mean(avg_lat_maxwind, na.rm = TRUE),
         Maxlat_anomaly = avg_lat_maxwind - Maxlat_mean) %>%
  ungroup()

# Compare frequency and intensity metrics
plot_freq_int(df_index_metrics)

# Plot correlation between ENSO and anomalies of Cat15_count and Cat35_count
cor_ENSO_freq(df_index_metrics)

# Plot correlation between ENSO and anomalies of ACE
cor_ENSO_ACE(df_index_metrics) 

# Plot correlation between ENSO and anomalies of Latitude of Maximum Winds 
cor_ENSO_maxlat(df_index_metrics)

# Plot correlation between AMO and anomalies of Cat15_count and Cat35_count
cor_AMO_freq(df_index_metrics)

# Plot correlation between AMO and anomalies of ACE
cor_AMO_ACE(df_index_metrics)


#--------------------------------------------------------------------------------------------
# CORRELATION matrix
#--------------------------------------------------------------------------------------------

# Load necessary libraries
library(knitr)
library(gt)
library(tidyverse, warn.conflict=F)
library(Hmisc)
library(purrr)
library(gridExtra)

# Format data into matrix
result <- process_tc_metrics(TC_metrics, seasonal_indices, cutoff_year = 2022)

# Access the frequency correlation data frame
df_for_freq_corr <- result$df_for_freq_corr

# Access the ACE correlation data frame
df_for_ace_corr <- result$df_for_ace_corr

# Subset for Southern Hemisphere (SH) and BASINS SI, SP
df_for_freq_corr_sh <- df_for_freq_corr %>%
  dplyr::select(starts_with("Cat15_count_SI"), starts_with("Cat15_count_SP"), starts_with("nino34_index_SH"), starts_with("amo_index_SH"))

df_for_ace_corr_sh <- df_for_ace_corr %>%
  dplyr::select(starts_with("avg_ACE_SI"), starts_with("avg_ACE_SP"), starts_with("nino34_index_SH"), starts_with("amo_index_SH"))

# Subset for Northern Hemisphere (NH) and BASINS EP, WP, NA, NI
df_for_freq_corr_nh <- df_for_freq_corr %>%
  dplyr::select(starts_with("Cat15_count_EP"), starts_with("Cat15_count_WP"), starts_with("Cat15_count_NA"), starts_with("Cat15_count_NI"), starts_with("nino34_index_NH"), starts_with("amo_index_NH"))

df_for_ace_corr_nh <- df_for_ace_corr %>%
  dplyr::select(starts_with("avg_ACE_EP"), starts_with("avg_ACE_WP"), starts_with("avg_ACE_NA"), starts_with("avg_ACE_NI"), starts_with("nino34_index_NH"), starts_with("amo_index_NH"))


# Calculate correlations for SH metrics
correlations_freq_sh <- formatted_cors(df_for_freq_corr_sh)
correlations_ace_sh <- formatted_cors(df_for_ace_corr_sh)

# Calculate correlations for NH metrics
correlations_freq_nh <- formatted_cors(df_for_freq_corr_nh)
correlations_ace_nh <- formatted_cors(df_for_ace_corr_nh)

create_correlation_plot_grid(correlations_freq_sh, correlations_ace_sh, correlations_freq_nh, correlations_ace_nh)



