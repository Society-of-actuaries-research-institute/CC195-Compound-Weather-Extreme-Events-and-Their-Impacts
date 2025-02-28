# %%
import os
import warnings
from datetime import datetime

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from loguru import logger
from numpy import std
from pandas import Series
from scipy.stats import fisk, gamma, genextreme, norm
from tqdm import tqdm

warnings.filterwarnings("ignore")

# %%
OUTPUT_PATH = "./data/wildfire"
PATH_SHP_CALIFORNIA = "./data/California/wildfires-california.shp"
PATH_ERA5 = "./data/ERA5/download_10_06.nc"

if not os.path.exists(OUTPUT_PATH):
    os.makedirs(OUTPUT_PATH)
    logger.info(f"Directory '{OUTPUT_PATH}' created.")
else:
    logger.info(f"Directory '{OUTPUT_PATH}' already exists.")


# %%
def calculate_si(series_ref, series_calc, dist_type, by_month=True):
    si_values = Series(index=series_calc.index, dtype=float)
    if dist_type == "gamma":
        if by_month == True:
            for month in range(1, 13):
                data_ref = (
                    series_ref[series_ref.index.month == month].dropna().sort_values()
                )
                data_calc = (
                    series_calc[series_calc.index.month == month].dropna().sort_values()
                )
                *pars, loc, scale = gamma.fit(data_ref, scale=std(data_ref))
                cdf = gamma.cdf(data_calc, pars, loc=loc, scale=scale)
                ppf = norm.ppf(cdf)
                si_values.loc[data_calc.index] = ppf
        else:
            data_ref = series_ref.dropna().sort_values()
            data_calc = series_calc.dropna().sort_values()
            *pars, loc, scale = gamma.fit(data_ref, scale=std(data_ref))
            cdf = gamma.cdf(data_calc, pars, loc=loc, scale=scale)
            ppf = norm.ppf(cdf)
            si_values.loc[data_calc.index] = ppf

    elif dist_type == "log-logistic":
        if by_month == True:
            for month in range(1, 13):
                data_ref = (
                    series_ref[series_ref.index.month == month].dropna().sort_values()
                )
                data_calc = (
                    series_calc[series_calc.index.month == month].dropna().sort_values()
                )
                *pars, loc, scale = fisk.fit(data_ref, scale=std(data_ref))
                cdf = fisk.cdf(data_calc, pars, loc=loc, scale=scale)
                ppf = norm.ppf(cdf)
                si_values.loc[data_calc.index] = ppf
        else:
            data_ref = series_ref.dropna().sort_values()
            data_calc = series_calc.dropna().sort_values()
            *pars, loc, scale = fisk.fit(data_ref, scale=std(data_ref))
            cdf = fisk.cdf(data_calc, pars, loc=loc, scale=scale)
            ppf = norm.ppf(cdf)
            si_values.loc[data_calc.index] = ppf

    elif dist_type == "gev":
        if by_month == True:
            for month in range(1, 13):
                data_ref = (
                    series_ref[series_ref.index.month == month].dropna().sort_values()
                )
                data_calc = (
                    series_calc[series_calc.index.month == month].dropna().sort_values()
                )
                *pars, loc, scale = genextreme.fit(data_ref, scale=std(data_ref))
                cdf = genextreme.cdf(data_calc, pars, loc=loc, scale=scale)
                ppf = norm.ppf(cdf)
                si_values.loc[data_calc.index] = ppf
        else:
            data_ref = series_ref.dropna().sort_values()
            data_calc = series_calc.dropna().sort_values()
            *pars, loc, scale = genextreme.fit(data_ref, scale=std(data_ref))
            cdf = genextreme.cdf(data_calc, pars, loc=loc, scale=scale)
            ppf = norm.ppf(cdf)
            si_values.loc[data_calc.index] = ppf

    else:
        raise ValueError('dist_type must be either "gamma", "log-logistic", or "gev"')

    return pd.Series(si_values, index=series_calc.index)


def calculate_SPI(dataframe, dist_type, scale, period_reference):

    """Computing the SPI using a distribution type for precipitation and a scale in months"""

    dataframe_copy = dataframe.copy().set_index("time")
    dataframe_copy["tp_rolling"] = dataframe_copy["tp"].rolling(scale).mean()
    dataframe_copy["si"] = calculate_si(
        series_ref=dataframe_copy.loc[
            period_reference[0] : period_reference[1], "tp_rolling"
        ],
        series_calc=dataframe_copy["tp_rolling"],
        dist_type=dist_type,
        by_month=True,
    )
    return dataframe_copy.reset_index()


def normalized_temperature_anomaly(dataframe):
    dataframe["heatwaves_normalized"] = (
        dataframe["t2m"] - dataframe["t2m"].mean()
    ) / dataframe["t2m"].std()
    return dataframe


# %%
def read_wildfire() -> pd.DataFrame:

    """Read wildfires database """
    gdf = gpd.read_file(PATH_SHP_CALIFORNIA)
    only_wildfires = pd.DataFrame(gdf[gdf["Incid_Type"] == "Wildfire"])
    only_wildfires["Ig_Date"] = pd.to_datetime(only_wildfires["Ig_Date"])
    return only_wildfires


def read_era5() -> pd.DataFrame:

    """Read era5 data"""
    dataset = xr.open_dataset(PATH_ERA5)
    df = dataset.to_dataframe().reset_index()
    df["time"] = pd.to_datetime(df["time"])
    df.loc[:, "wildfire"] = 0
    df.loc[:, "BurnBndAc"] = 0
    return df

def add_wildfires_to_climate_dataframe(
    dataframe_wildfire, dataframe_climate
) -> pd.DataFrame:


    """Add wildifres events to the climate data.
    If two wildfires happen the same month at the same location,we only write the last one in the dataframe."""

    for wildfire in tqdm(dataframe_wildfire.iterrows()):
        wildfire_data = wildfire[1]
        input_latitude = float(wildfire_data.BurnBndLat)
        input_longitude = float(wildfire_data.BurnBndLon)
        input_date = wildfire_data["Ig_Date"]

        # Calculate the absolute differences
        dataframe_climate["latitude_diff"] = (
            dataframe_climate["latitude"] - input_latitude
        ).abs()
        dataframe_climate["longitude_diff"] = (
            dataframe_climate["longitude"] - input_longitude
        ).abs()

        # Find the rows with the smallest latitude and longitude differences
        closest_latitude_rows = dataframe_climate.loc[
            dataframe_climate["latitude_diff"]
            == dataframe_climate["latitude_diff"].min()
        ]

        closest_longitude_rows = closest_latitude_rows.loc[
            closest_latitude_rows["longitude_diff"]
            == closest_latitude_rows["longitude_diff"].min()
        ]

        # Among those rows, find the one with the closest date
        closest_longitude_rows["time_diff"] = (
            closest_longitude_rows["time"] - input_date
        ).abs()
        closest_row_index = closest_longitude_rows["time_diff"].idxmin()

        dataframe_climate.loc[closest_row_index, "wildfire"] = 1
        dataframe_climate.loc[closest_row_index, "BurnBndAc"] = wildfire_data[
            "BurnBndAc"
        ]

    return dataframe_climate


# %%
def compute_spi_temperature_anomaly(
    dataframe, period_reference, scale=6
) -> pd.DataFrame:

    "Computing standardized precipitation index (SPI) and standardized temperature anomaly "
    df_final = dataframe.groupby(["longitude", "latitude"], as_index=True).apply(
        lambda x: calculate_SPI(
            x, dist_type="gamma", scale=scale, period_reference=period_reference
        )
    )

    df_final = df_final.rename(columns={"si": f"spi_{scale}_months"})

    dataframe["month"] = [date.month for date in dataframe.time]

    df_heatwave = dataframe.groupby(["latitude", "longitude", "month"]).apply(
        lambda x: normalized_temperature_anomaly(x)
    )

    df_final_w_hw = df_final.set_index(["latitude", "longitude", "time"]).join(
        df_heatwave.set_index(["latitude", "longitude", "time"]), lsuffix="_df_spi"
    )

    columns_to_keep = (
        "t2m",
        "tp",
        "wildfire",
        "BurnBndAc",
        "spi_6_months",
        "month",
        "heatwaves_normalized"
    )

    output_df = df_final_w_hw.loc[:, columns_to_keep]

    return output_df


# %%
if __name__ == "__main__":

    logger.info("Reading wildfire database.")
    df_wildfire = read_wildfire()

    logger.info("Reading era5 data.")
    df_climate = read_era5()

    logger.info("Fusion of the two dataframes.")
    df_climate_wildfire = add_wildfires_to_climate_dataframe(df_wildfire, df_climate)

    period_reference = (
        datetime(min(df_wildfire["Ig_Date"]).year, 1, 1),
        datetime(max(df_wildfire["Ig_Date"]).year, 12, 31),
    )

    logger.info("Computing SPI and normalized temperature anomaly.")
    df_completed = compute_spi_temperature_anomaly(
        df_climate_wildfire, period_reference, scale=6
    )

    logger.info("Saving output dataframe")
    df_completed.to_parquet(
        f"{OUTPUT_PATH}/climate_data_w_wildfires_normalized.parquet"
    )
