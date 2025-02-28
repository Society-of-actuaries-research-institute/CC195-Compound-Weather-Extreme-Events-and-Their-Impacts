
#------------------------------------------------------------------------------
# Plot NINO3.4 time series
#------------------------------------------------------------------------------

timeseries_nino34 <- function(data, X) {ggplot(data, aes(x = {{X}})) +
  # Bar plot for ANOM3.4 anomalies
  geom_bar(aes(y = ANOM3.4 , fill = ANOM3.4  > 0), stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"), guide = FALSE) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = -0.5, linetype = "dashed") + 
  scale_y_continuous(
    name = "NINO3.4 Index") + 
  labs(title = "NINO3.4 Anomalies",
       x = "Time",
       color = "Series") +
  theme_minimal()
}

#------------------------------------------------------------------------------
# Plot Time series of ENSO and TC frequency and ACE
#------------------------------------------------------------------------------

plot_enso_freq <- function(data, X) {ggplot(df_index_metrics, aes(x = YR)) +
  # Bar plot for NINO3.4 anomalies
  geom_bar(aes(y = seasonal_nino34_index, fill = seasonal_nino34_index > 0), stat = "identity") +
  scale_fill_manual(values = c("TRUE" = rgb(255, 85, 113, maxColorValue = 255),  # Pink
                               "FALSE" = rgb(2, 77, 124, maxColorValue = 255)),  # Blue
                    guide = FALSE) +
  
  # Horizontal dashed lines
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = -0.5, linetype = "dashed") +
  
  # Line plot for Cat15 count (dark green, dotted line)
  geom_line(aes(y = Cat15_count / 10), color = "darkgreen", size = 1, linetype = "dotted") +
  
  # Line plot for Cat35 count (med-dark gray)
  geom_line(aes(y = Cat35_count / 10), color = rgb(127, 127, 127, maxColorValue = 255), size = 1) +
  
  # Customize y-axes
  scale_y_continuous(
    name = "NINO3.4 seasonal mean Index",
    sec.axis = sec_axis(~ . * 10, name = "Tropical cyclone frequency")  # Adjust the scale factor to match the line plots
  ) +
  
  # Labels and theme
  labs(
    title = "NINO3.4 Anomalies and Tropical Cyclone Count per Season",
    x = "Time",
    color = "Series"
  ) +
  theme_minimal() +
  facet_wrap(~ BASIN, scales = "free_y")
}


plot_enso_ACE <- function(data, X) {ggplot(df_index_metrics, aes(x = YR)) +
  geom_bar(aes(y = seasonal_nino34_index, fill = seasonal_nino34_index > 0), stat = "identity") +
  scale_fill_manual(values = c("TRUE" = rgb(255, 85, 113, maxColorValue = 255), 
                               "FALSE" = rgb(2, 77, 124, maxColorValue = 255)), guide = FALSE) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = -0.5, linetype = "dashed") +
  geom_line(aes(y = avg_ACE / 10), color = rgb(200, 50, 75, maxColorValue = 255), size = 1) +  # Darker pink line
  scale_y_continuous(
    name = "NINO3.4 seasonal mean Index",
    sec.axis = sec_axis(~ . * 10, name = "Tropical cyclone mean ACE")  # Adjust the scale factor to match the line plots
  ) +
  labs(
    title = "NINO3.4 Anomalies and Tropical cyclone mean ACE per season",
    x = "Time",
    color = "Series"
  ) +
  theme_minimal() +
  facet_wrap(~ BASIN, scales = "free_y")
}


#------------------------------------------------------------------------------
# Plot Time series of AMO and TC frequency and ACE
#------------------------------------------------------------------------------

plot_AMO_freq <- function(data, X) {ggplot(df_index_metrics, aes(x = YR)) +
  # Bar plot for AMO anomalies
  geom_bar(aes(y = seasonal_amo_index, fill = seasonal_amo_index > 0), stat = "identity") +
  scale_fill_manual(values = c("TRUE" = rgb(255, 85, 113, maxColorValue = 255),  # Pink
                               "FALSE" = rgb(2, 77, 124, maxColorValue = 255)),  # Blue
                    guide = FALSE) +
  # Line for Cat15 count (now using dotted line)
  geom_line(aes(y = Cat15_count / 25), color = "darkgreen", size = 1, linetype = "dotted") +  # Dotted line
  # Line for Cat35 count (using med-dark gray)
  geom_line(aes(y = Cat35_count / 25), color = rgb(127, 127, 127, maxColorValue = 255), size = 1) +  
  scale_y_continuous(
    name = "AMO seasonal mean Index",
    limits = primary_y_limits,
    sec.axis = sec_axis(~ . * 25, name = "Tropical cyclone frequency")  # Adjust the scale factor to match the line plots
  ) +
  labs(
    title = "AMO Anomalies and Tropical cyclone count per season",
    x = "Time",
    color = "Series"
  ) +
  theme_minimal() +
  facet_wrap(~ BASIN, scales = "free_y")
}

plot_AMO_ACE <- function(data, X) {ggplot(df_index_metrics, aes(x = YR)) +
  # Bar plot for AMO anomalies
  geom_bar(aes(y = seasonal_amo_index, fill = seasonal_amo_index > 0), stat = "identity") +
  scale_fill_manual(values = c("TRUE" = rgb(255, 85, 113, maxColorValue = 255),  # Pink
                               "FALSE" = rgb(2, 77, 124, maxColorValue = 255)),  # Blue
                    guide = FALSE) +
  
  # Line plot for avg_ACE (darker pink)
  geom_line(aes(y = avg_ACE / 20), color = rgb(200, 50, 75, maxColorValue = 255), size = 1) +  
  
  # Dotted trend line for avg_ACE (darker pink)
  geom_smooth(aes(y = avg_ACE / 20), method = "lm", linetype = "dotted", color = rgb(200, 50, 75, maxColorValue = 255), se = FALSE, size = 0.7) +
  
  # Customize y-axes
  scale_y_continuous(
    name = "AMO seasonal mean Index",
    limits = primary_y_limits,
    sec.axis = sec_axis(~ . * 20, name = "Tropical cyclone mean ACE")  # Adjust the scale factor to match the line plots
  ) +
  
  # Labels and theme
  labs(
    title = "AMO Anomalies and Tropical Cyclone Mean ACE per Season",
    x = "Time",
    color = "Series"
  ) +
  theme_minimal() +
  facet_wrap(~ BASIN, scales = "free_y")
}

#------------------------------------------------------------------------------
# Plot scatter diagram of the intensity vs frequency metrics
#------------------------------------------------------------------------------
plot_freq_int<- function(data) {ggplot(data) +
    geom_point(aes(x = Cat15_anomaly, y = ACE_anomaly), color = "grey30") +  # Scatter plot points for ACE_anomaly vs. Cat15_anomaly
    geom_smooth(aes(x = Cat15_anomaly, y = ACE_anomaly), 
                method = "lm", color = "grey30", linetype = "dotted", size = 0.7, se = FALSE) +  # Add dark grey trend line
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Horizontal zero line
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Vertical zero line
    labs(
      title = "Relationship between Category 1-5 count and Mean ACE Anomalies",
      x = "Cat15 Anomalies",
      y = "Mean ACE Anomalies"
    ) +
    theme_minimal() +
    facet_wrap(~ BASIN, scales = "free") +  # Create a separate plot for each basin
    geom_text(
      aes(
        label = paste("", round(cor_ACE_Cat15, 2)),  # Update correlation calculation
        x = -Inf, y = Inf
      ),
      hjust = -0.1, vjust = 1.1,  # Position next to the title
      size = 4,
      color = "grey30",
      inherit.aes = FALSE
    )
}

#------------------------------------------------------------------------------
# Plot scatter diagram of NINO3.4 Index vs Tropical Cyclone Count Anomalies
#------------------------------------------------------------------------------
cor_ENSO_freq <- function(data) {ggplot(data) +
    # Scatter plot for Cat15_anomaly (pink)
    geom_point(aes(x = seasonal_nino34_index, y = Cat15_anomaly), color = rgb(255, 85, 113, maxColorValue = 255)) +  
    # Scatter plot for Cat35_anomaly (blue)
    geom_point(aes(x = seasonal_nino34_index, y = Cat35_anomaly), color = rgb(2, 77, 124, maxColorValue = 255)) +  
    
    # Trend line for Cat15_anomaly (pink, dotted)
    geom_smooth(aes(x = seasonal_nino34_index, y = Cat15_anomaly), 
                method = "lm", color = rgb(255, 85, 113, maxColorValue = 255), 
                linetype = "dotted", size = 0.7, se = FALSE) + 
    # Trend line for Cat35_anomaly (blue, dotted)
    geom_smooth(aes(x = seasonal_nino34_index, y = Cat35_anomaly), 
                method = "lm", color = rgb(2, 77, 124, maxColorValue = 255), 
                linetype = "dotted", size = 0.7, se = FALSE) + 
    
    # Horizontal and vertical lines
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  
    geom_vline(xintercept = 0.5, linetype = "dotted", color = "black") +  
    geom_vline(xintercept = -0.5, linetype = "dotted", color = "black") +  
    
    # Labels for the plot
    labs(
      title = "Relationship between NINO3.4 Index and Tropical Cyclone Count Anomalies",
      x = "Average NINO3.4 Anomalies",
      y = "Cyclone Count Anomalies"
    ) +
    
    # Minimal theme for the plot
    theme_minimal() +
    
    # Facet by BASIN
    facet_wrap(~ BASIN, scales = "free") +  
    
    # Add correlation values for Cat15_anomaly (pink)
    geom_text(
      aes(
        label = paste("", round(cor_enso_Cat15, 2)),
        x = -Inf, y = Inf
      ),
      hjust = -0.1, vjust = 1.1,  # Adjust position next to the title
      size = 4,
      color = rgb(255, 85, 113, maxColorValue = 255),
      inherit.aes = FALSE
    ) +
    
    # Add correlation values for Cat35_anomaly (blue)
    geom_text(
      aes(
        label = paste("", round(cor_enso_Cat35, 2)),
        x = -Inf, y = Inf
      ),
      hjust = -0.1, vjust = 2.3,  # Adjust position next to the title
      size = 4,
      color = rgb(2, 77, 124, maxColorValue = 255),
      inherit.aes = FALSE
    )
}

cor_ENSO_ACE <- function(data) {ggplot(data) +
  geom_point(aes(x = seasonal_nino34_index, y = ACE_anomaly), color = "grey30") +  # Scatter plot points for ACE_anomaly
  geom_smooth(aes(x = seasonal_nino34_index, y = ACE_anomaly), 
              method = "lm", color = "grey30", linetype = "dotted", size = 0.7, se = FALSE) +  # Add dark grey trend line
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Horizontal zero line
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Vertical zero line
  geom_vline(xintercept = 0.5, linetype = "dotted", color = "black") +  # Vertical ENSO line
  geom_vline(xintercept = -0.5, linetype = "dotted", color = "black") +  # Vertical ENSO line
  labs(
    title = "Relationship between NINO3.4 Index and Mean ACE Anomalies",
    x = "Average NINO3.4 Anomalies",
    y = "Mean ACE Anomalies"
  ) +
  theme_minimal() +
  facet_wrap(~ BASIN, scales = "free") +  # Create a separate plot for each basin
  geom_text(
    aes(
      label = paste("", round(cor_enso_ACE, 2)),
      x = -Inf, y = Inf
    ),
    hjust = -0.1, vjust = 1.1,  # Position next to the title
    size = 4,
    color = "grey30",
    inherit.aes = FALSE
  )
}

cor_ENSO_maxlat <- function(data) {ggplot(data) +
  geom_point(aes(x = seasonal_nino34_index, y = Maxlat_anomaly), color = "grey30") +  # Scatter plot points for Maxlat_anomaly
  geom_smooth(aes(x = seasonal_nino34_index, y = Maxlat_anomaly), 
              method = "lm", color = "grey30", linetype = "dotted", size = 0.7, se = FALSE) +  # Add dark grey trend line
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Horizontal zero line
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Vertical zero line
  geom_vline(xintercept = 0.5, linetype = "dotted", color = "black") +  # Vertical ENSO line
  geom_vline(xintercept = -0.5, linetype = "dotted", color = "black") +  # Vertical ENSO line
  labs(
    title = "Relationship between NINO3.4 Index and Latitude of Maximum Winds Anomalies",
    x = "Average NINO3.4 Anomalies",
    y = "Latitude Max Winds Anomalies"
  ) +
  theme_minimal() +
  facet_wrap(~ BASIN, scales = "free") +  # Create a separate plot for each basin
  geom_text(
    aes(
      label = paste("", round(cor_enso_maxlat, 2)),
      x = -Inf, y = Inf
    ),
    hjust = -0.1, vjust = 1.1,  # Position next to the title
    size = 4,
    color = "grey30",
    inherit.aes = FALSE
  )
}

cor_AMO_freq <- function(data) {ggplot(data) +
  # Scatter plot for Cat15_anomaly (pink)
  geom_point(aes(x = seasonal_amo_index, y = Cat15_anomaly), color = rgb(255, 85, 113, maxColorValue = 255)) +  
  # Scatter plot for Cat35_anomaly (blue)
  geom_point(aes(x = seasonal_amo_index, y = Cat35_anomaly), color = rgb(2, 77, 124, maxColorValue = 255)) +  
  
  # Trend line for Cat15_anomaly (pink, dotted)
  geom_smooth(aes(x = seasonal_amo_index, y = Cat15_anomaly), method = "lm", 
              color = rgb(255, 85, 113, maxColorValue = 255), linetype = "dotted", se = FALSE) +  
  # Trend line for Cat35_anomaly (blue, dotted)
  geom_smooth(aes(x = seasonal_amo_index, y = Cat35_anomaly), method = "lm", 
              color = rgb(2, 77, 124, maxColorValue = 255), linetype = "dotted", se = FALSE) +  
  
  # Horizontal and vertical lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  
  
  # Labels for the plot
  labs(
    title = "Relationship between AMO Index and Tropical Cyclone Count Anomalies",
    x = "Average AMO Anomalies",
    y = "Cyclone Count Anomalies"
  ) +
  
  # Minimal theme for the plot
  theme_minimal() +
  
  # Facet by BASIN
  facet_wrap(~ BASIN, scales = "free") +  
  
  # Add correlation values for Cat15_anomaly (pink)
  geom_text(
    aes(
      label = paste("", round(cor_amo_Cat15, 2)),
      x = -Inf, y = Inf
    ),
    hjust = -0.1, vjust = 1.1,  # Adjust position next to the title
    size = 4,
    color = rgb(255, 85, 113, maxColorValue = 255),
    inherit.aes = FALSE
  ) +
  
  # Add correlation values for Cat35_anomaly (blue)
  geom_text(
    aes(
      label = paste("", round(cor_amo_Cat35, 2)),
      x = -Inf, y = Inf
    ),
    hjust = -0.1, vjust = 2.3,  # Adjust position next to the title, slightly below the first text
    size = 4,
    color = rgb(2, 77, 124, maxColorValue = 255),
    inherit.aes = FALSE
  )
}

cor_AMO_ACE <- function(data) {ggplot(data) +
  geom_point(aes(x = seasonal_amo_index, y = ACE_anomaly), color = "black") +  # Scatter plot points for ACE_anomaly
  geom_smooth(aes(x = seasonal_amo_index, y = ACE_anomaly), method = "lm", color = "grey30", linetype = "dotted", se = FALSE) +  # Add a trend line for ACE_anomaly
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Add a horizontal zero line
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Add a vertical zero line
  labs(
    title = "Relationship between AMO Index and Tropical Cyclone Mean ACE Anomalies",
    x = "Average AMO Anomalies",
    y = "Mean ACE Anomalies"
  ) +
  theme_minimal() +
  facet_wrap(~ BASIN, scales = "free") +  # Create a separate plot for each basin
  # Add correlation values to the plot in grey30
  geom_text(
    aes(
      label = paste("", round(cor_amo_ACE, 2)),
      x = -Inf, y = Inf
    ),
    hjust = -0.1, vjust = 1.1,  # Adjust position next to the title
    size = 4,
    color = "grey30",
    inherit.aes = FALSE
  )
}

#------------------------------------------------------------------------------
# Plot number of TC tracks per grid cells during El Nino / La Nina / Nino-Nina
#------------------------------------------------------------------------------
plot_nino_nina_difference <- function(grid_counts, grid_counts_diff, world_expanded, grid_size, 
                                      lon_min, lon_max, lat_min, lat_max) {
  
  # Determine the maximum track count for color scale limits
  max_track_count <- max(grid_counts$avg_track_count, na.rm = TRUE)
  
  # Plot Nino
  p_nino <- ggplot() +
    geom_tile(data = grid_counts %>% filter(NINO34 == "Nino" & !is.na(avg_track_count)), 
              aes(x = lon_bin + grid_size / 2, y = lat_bin + grid_size / 2, fill = avg_track_count)) +
    geom_polygon(data = world_expanded, aes(x = long, y = lat, group = group), fill = NA, color = "black") +
    coord_fixed(ratio = 1.3, xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max)) +
    scale_fill_gradientn(
      name = "Avg Track Count", 
      colors = c(
        rgb(255, 245, 247, maxColorValue = 255), # Almost white (very light pink)
        rgb(255, 190, 200, maxColorValue = 255), # Light pink
        rgb(255, 85, 113, maxColorValue = 255),  # Custom pink
        rgb(204, 51, 85, maxColorValue = 255)    # Darker pink
      ), 
      na.value = "transparent", 
      limits = c(0, max_track_count)
    ) +
    labs(
      title = "Nino Track Counts (2° x 2° Grid)",
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal() +
    theme(legend.position = "right")
  
  # Plot Nina
  p_nina <- ggplot() +
    geom_tile(data = grid_counts %>% filter(NINO34 == "Nina" & !is.na(avg_track_count)), 
              aes(x = lon_bin + grid_size / 2, y = lat_bin + grid_size / 2, fill = avg_track_count)) +
    geom_polygon(data = world_expanded, aes(x = long, y = lat, group = group), fill = NA, color = "black") +
    coord_fixed(ratio = 1.3, xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max)) +
    scale_fill_gradientn(
      name = "Avg Track Count", 
      colors = c(
        rgb(255, 245, 247, maxColorValue = 255), # Almost white (very light pink)
        rgb(255, 190, 200, maxColorValue = 255), # Light pink
        rgb(255, 85, 113, maxColorValue = 255),  # Custom pink
        rgb(204, 51, 85, maxColorValue = 255)    # Darker pink
      ), 
      na.value = "transparent", 
      limits = c(0, max_track_count)
    ) +
    labs(
      title = "Nina Track Counts (2° x 2° Grid)",
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal() +
    theme(legend.position = "right")
  
  # Plot Nino - Nina difference
  p_diff <- ggplot() +
    geom_tile(data = grid_counts_diff %>% filter(!is.na(nino_nina_diff)), 
              aes(x = lon_bin + grid_size / 2, y = lat_bin + grid_size / 2, fill = nino_nina_diff)) +
    geom_polygon(data = world_expanded, aes(x = long, y = lat, group = group), fill = NA, color = "black") +
    coord_fixed(ratio = 1.3, xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max)) +
    scale_fill_gradient2(
      name = "Nino - Nina Track Count", 
      low = rgb(2, 77, 124, maxColorValue = 255),  # Custom blue
      high = rgb(255, 85, 113, maxColorValue = 255),  # Custom pink
      mid = "white", 
      midpoint = 0, 
      na.value = "transparent"
    ) +
    labs(
      title = "Difference in Track Counts (Nino - Nina)",
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal() +
    theme(legend.position = "right")
  
  # Combine the three plots into one row
  combined_plot <- p_nino + p_nina + p_diff + plot_layout(ncol = 3)
  
  # Display the combined plot
  print(combined_plot)
}


#------------------------------------------------------------------------------
# Plot correlation matrix
#------------------------------------------------------------------------------

create_correlation_plot <- function(data, title) {
  plot <- data %>%
    ggplot(aes(measure1, measure2, fill = r, label = round(r_if_sig, 2))) +
    geom_tile() +
    labs(
      x = NULL, 
      y = NULL, 
      fill = "Pearson's\nCorrelation", 
      title = title
    ) +
    scale_fill_gradient2(
      mid = "#FBFEF9", 
      low = "#024D7C", 
      high = "#E1415D", 
      limits = c(-1, 1)
    ) +
    geom_text() +
    theme_classic() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme(
      text = element_text(family = "Roboto"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  return(plot)
}

create_correlation_plot_grid <- function(correlations_freq_sh, correlations_ace_sh, correlations_freq_nh, correlations_ace_nh) {
  # Create individual plots using the create_correlation_plot function
  plot_freq_sh <- create_correlation_plot(correlations_freq_sh, "Correlations Frequency Metrics (SH)")
  plot_ace_sh <- create_correlation_plot(correlations_ace_sh, "Correlations ACE Metrics (SH)")
  plot_freq_nh <- create_correlation_plot(correlations_freq_nh, "Correlations Frequency Metrics (NH)")
  plot_ace_nh <- create_correlation_plot(correlations_ace_nh, "Correlations ACE Metrics (NH)")
  
  # Arrange the four plots in a 2x2 grid
  grid.arrange(
    plot_freq_nh, plot_ace_nh,
    plot_freq_sh, plot_ace_sh,
    ncol = 2, nrow = 2
  )
}
