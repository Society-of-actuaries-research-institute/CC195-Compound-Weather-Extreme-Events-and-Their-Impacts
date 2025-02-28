# Compound Weather Extreme Events and their Impacts

Analysis of the effect of compound weather event of 3 types :
- Temporal Compounding event - Type 1
- Multivariate Compounding event - Type 2
- Spatially Compounding event - Type 3

We then project the future possible climates and analyse the evolution of single event (drought and positive temperature anomaly) and compound events (co-occurence of drought and positive temperature anomaly)
## Temporal Compounding event - Type 1

We study the effect of temporaly compounding drought and positive temperature anomaly on the occurence of wildfire in California

### Code
1- src/create_dataset_wildfire.py
2- src/wildfire_analysis.py

## Multivariate Compounding event - Type 2

We study the effect of multivariate (co-occurence) compounding drought and positive temperature anomaly on winter wheat yield loss in California.

### Code
1- src/yield_analysis.py

## Spatially Compounding event - Type 3

We study the effect of climate drivers (ENSO & AMO) on the cyclone activity in different bassins.

### Code
Main code : src/TC_Project_final.R
Data visualization functions : TC_figures.R
Data manipulation functions : TC_metrics_drivers.R

### Data
data/TropicalCyclone/IBTrACS Version 04 shapefile


## Future Climate

We study the evolution of climate in California. We analyse how single events (drought and positive temperature anomaly) and  multivariate (co-occurence) compounding drought and positive temperature anomaly are likely to evolve in a warming climate.

### Code
1- src/analysis_future_climate.py
