# %%
import os
import warnings
from datetime import datetime

import cartopy.crs as crs
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import pymannkendall as mk
import seaborn as sns
import xarray as xr
from loguru import logger
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy import std
from pandas import Series
from scipy.stats import chi2_contingency, fisk, gamma, genextreme, norm
from shapely.geometry import mapping

warnings.filterwarnings("ignore")

# %%
OUTPUT_PATH = "./output/crop_loss"
PATH_ERA5 = "./data/ERA5/download_10_06.nc"
PATH_YIELD = "./data/Streamlit-yield-data/United States_gad_2.csv"
SHP_PATH = "./data/California/california_boundaries_4326.shp"

spi_scale: int = 2
crop, month_start, month_end = ("WHEAT_WINTER", 4, 5)

if not os.path.exists(OUTPUT_PATH):
    os.makedirs(OUTPUT_PATH)
    logger.info(f"Directory '{OUTPUT_PATH}' created.")
else:
    logger.info(f"Directory '{OUTPUT_PATH}' already exists.")


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

    """Computing the Standardized Precipitation Index (SPI) using a distribution type for precipitation and a scale in months"""

    df = dataframe.copy()
    df["tp_rolling"] = df["tp"].rolling(scale).mean()
    df["si"] = calculate_si(
        series_ref=df.loc[period_reference[0] : period_reference[1], "tp_rolling"],
        series_calc=df["tp_rolling"],
        dist_type=dist_type,
        by_month=True,
    )
    return df["si"]


def compute_cond_prob(
    df: pd.DataFrame,
    yield_threshold: float = 0,
    spi_threshold: float = np.inf,
    heat_threshold: float = np.NINF,
):
    """
    Computing the conditional probability of having a yield below the mean depending on different context.

    P(Yield|D ∩ H)=P(LowYield ∩ D ∩ H) / P(D ∩ H)
    where :
    Yield = Yield < yield_threshold
    D = Drought -> SPI <= spi_threhold
    H = Heatwave -> Standardized T° anomaly => heat_threshold
    to compute P(LowYield|H) or P(LowYield|D) , put the other driver threshold to np.inf or np.NINF
    """

    numerateur = df[
        (df.loc[:, f"spi_{spi_scale}_months"] <= spi_threshold)
        & (df.loc[:, "standardized T° anomaly"] >= heat_threshold)
        & (df.loc[:, "Z-score"] < yield_threshold)
    ].count()["Z-score"]

    denominateur = df[
        (df.loc[:, f"spi_{spi_scale}_months"] <= spi_threshold)
        & (df.loc[:, "standardized T° anomaly"] >= heat_threshold)
    ].count()["Z-score"]

    return numerateur, denominateur, numerateur / denominateur


def anomaly(x):
    """Computing the standardized Temperature Anomaly"""
    x["standardized T° anomaly"] = (x.t2m - x.t2m.mean()) / x.t2m.std()
    return x


def zscore(x):
    """Computing the Z-score, i.e, the standardized yield and the difference to the mean."""
    x["Z-score"] = (x["avgYield"] - x["avgYield"].mean()) / x["avgYield"].std()
    x["difference_to_the_mean_%"] = (
        (x["avgYield"] - x["avgYield"].mean()) / x["avgYield"].mean() * 100
    )
    return x


# %%


def read_climate_data() -> xr.Dataset:
    """Reading era5 data"""
    dataset = xr.open_dataset(PATH_ERA5)

    dataset = dataset.rio.write_crs("4326")

    dataset = dataset.rename({"latitude": "y", "longitude": "x"})

    return dataset


def read_yield_data() -> pd.DataFrame:

    """Reading yield data at county level for California."""

    us_yield = pd.read_csv(PATH_YIELD, parse_dates=["year"])

    us_yield_wheat_durum = us_yield[us_yield["original_name"] == crop]

    us_yield_wheat_durum["year"] = [time.year for time in us_yield_wheat_durum["year"]]

    us_yield_wheat_durum["gad2"] = us_yield_wheat_durum["gad2"].str.upper()

    return us_yield_wheat_durum


# %%


def detect_untrended_yield_data(dataframe):
    """
    Using mann kendall test, we figure that from 1985 to 2022 there is no trend detected.
    Hence, one does not need to detrend the yield data.
    For WHEAT_WINTER : 2005 ->2022
    We want to work only on a time period without a trend i.e 2005-2022"""

    for period in range(1950, 2022, 5):

        period_condition = dataframe.year > period

        yield_wheat = (
            dataframe[period_condition].groupby("year").mean(numeric_only=True)
        )

        result = mk.original_test(yield_wheat["avgYield"])

        logger.info(f"trend={result.trend},trend bool={result.h},p value={result.p}")

    return None


# %%
def computing_z_score(dataframe, period_reference) -> pd.DataFrame:
    "Computing Z-score on the dataframe"
    short_data = dataframe[dataframe.year >= period_reference[0].year]
    short_data = short_data.groupby("gad2").apply(lambda x: zscore(x))

    return short_data


# %%
def display_lineplot_yield(dataframe):

    """Display lineplot of the timeseries of the average yield."""
    plot = sns.lineplot(data=dataframe, x="year", y="avgYield")

    plot.set_ylabel("Average yield [T/ha]")

    plot.set_title(f"Yield evolution from 2005 to 2022 for winter wheat ")

    labels = ["2005", "2007", "2009", "2011", "2013", "2015", "2017", "2019", "2022"]

    plot.set_xticks(ticks=[int(l) for l in labels], labels=labels)

    plot.figure.savefig(f"{OUTPUT_PATH}/timeseries_yield.png")

    plot.figure.clear()

    return None


# %%
def display_lineplot_z_score(dataframe):

    """Display lineplot of the timeseries of the scaled average yield."""
    plot = sns.lineplot(data=dataframe, x="year", y="Z-score")

    plot.set_ylabel("Average yield scaled [Z-score]")

    plot.set_title(f"Yield evolution from 2005 to 2022 for winter wheat ")

    labels = ["2005", "2007", "2009", "2011", "2013", "2015", "2017", "2019", "2022"]

    plot.set_xticks(ticks=[int(l) for l in labels], labels=labels)

    plot.axhline(y=0)

    plot.figure.savefig(f"{OUTPUT_PATH}/timeseries_yield_scaled.png")

    plot.figure.clear()

    return None


# %%
def read_shapefile_california() -> gpd.GeoDataFrame:

    '''Reading the California shapefile'''
    gdf = gpd.read_file(SHP_PATH)

    gdf["gad2"] = gdf["NAME"].str.upper()

    return gdf


def average_over_county(
    geodataframe, yield_data, climate_data, period_reference
) -> pd.DataFrame:

    """
    Computing the averageing of climate data over counties to fit with the yield regional data.
    We only use data from 2005-2022 as this is the period for wich we observe no trend from yield data.
    """

    all_county_dataframe = pd.DataFrame()

    dataset_time_not_trend = climate_data.sel(
        time=slice(
            datetime(year=2005, month=1, day=1), datetime(year=2022, month=12, day=31)
        )
    )

    geojson_geometry = [mapping(value) for value in geodataframe.geometry.values]
    geojson_name = [name for name in geodataframe["gad2"]]

    for county_geometry, county_name in zip(geojson_geometry, geojson_name):

        if county_name in yield_data["gad2"].unique():

            tmp = dataset_time_not_trend.rio.clip(
                [county_geometry], geodataframe.crs, drop=True
            )
            tmp_mean = tmp.mean(dim=["x", "y"])
            df_mean_over_county = tmp_mean.to_dataframe()
            df_mean_over_county["month"] = [
                date.month
                for date in df_mean_over_county.index.get_level_values("time")
            ]
            df_mean_over_county["spi_2_months"] = calculate_SPI(
                df_mean_over_county,
                scale=spi_scale,
                dist_type="gamma",
                period_reference=period_reference,
            )

            df_mean_over_county["gad2"] = county_name

            df_temperature = (
                df_mean_over_county[
                    (df_mean_over_county.month >= month_start)
                    & (df_mean_over_county.month <= month_end)
                ]
                .groupby(pd.Grouper(freq="Y"))
                .mean(numeric_only=True)
            )

            df_temperature["standardized T° anomaly"] = (
                df_temperature.t2m - df_temperature.t2m.mean()
            ) / df_temperature.t2m.std()

            df_spi = (
                df_mean_over_county.loc[
                    df_mean_over_county.month == month_end, "spi_2_months"
                ]
                .groupby(pd.Grouper(freq="Y"))
                .mean(numeric_only=True)
            )

            df_spi_temperature = df_temperature.drop("spi_2_months", axis=1).join(
                df_spi
            )

            df_spi_temperature.index = [time.year for time in df_spi_temperature.index]

            df_spi_temperature["gad2"] = county_name

            df_spi_temperature = df_spi_temperature.drop(["spatial_ref"], axis=1)

            df_spi_temperature = df_spi_temperature.set_index("gad2", append=True)

            df_spi_temperature.index.names = ("year", "gad2")

            us_yield_poi = yield_data[(yield_data["gad2"] == county_name)].set_index(
                ["year", "gad2"]
            )
            join = us_yield_poi.join(df_spi_temperature)

            all_county_dataframe = pd.concat([all_county_dataframe, join])

        else:
            print(f"{str.upper(county_name)} is not among the county shortlisted")

    return all_county_dataframe


# %%


def display_jointplot(dataframe):

    """Display jointplot of standardized temperature anomaly and spi-2-months."""
    g = sns.jointplot(
        data=dataframe,
        x="standardized T° anomaly",
        y=f"spi_{spi_scale}_months",
        kind="scatter",
        ylim=[-3, 3],
        xlim=[-3, 3],
    )

    g.ax_joint.set_xlabel("Standardized T° anomaly")

    g.ax_joint.set_ylabel(f"SPI-{spi_scale}-month")

    g.plot_marginals(sns.kdeplot)

    g.refline(x=0, y=0, joint=True)

    g.figure.suptitle(
        t=f"Jointplot of SPI-{spi_scale}-month and standardized T° anomaly\n for winter wheat yield ",
        y=1.05,
    )

    g.savefig("/".join([OUTPUT_PATH, "jointplot_all_data.png"]))

    g.figure.clear()

    return None


# %%
def display_jointplot_differend_std(dataframe):

    """Display jointplot of standardized temperature anomaly and spi-2-months for different yield loss threshold."""
    for std in range(0, -3, -1):
        g = sns.jointplot(
            data=dataframe[dataframe.loc[:, "Z-score"] < std],
            x="standardized T° anomaly",
            y=f"spi_{spi_scale}_months",
            kind="scatter",
            ylim=[-3, 3],
            xlim=[-3, 3],
        )
        g.ax_joint.set_xlabel("Standardized T° anomaly")
        g.ax_joint.set_ylabel(f"SPI-{spi_scale}-month")

        g.plot_marginals(sns.kdeplot)
        g.refline(x=0, y=0, joint=True)
        g.figure.suptitle(
            t=f"Jointplot of SPI-{spi_scale}-month and standardized T° anomaly\n for winter wheat {std} std below average yield ",
            y=1.05,
        )
        g.savefig("/".join([OUTPUT_PATH, f"joinplot_yield_below_{std}_std.png"]))
        g.figure.clear()

    return None


# %%
def compute_probability_context(dataframe) -> pd.DataFrame:

    """Computing probability of yield loss depending on climate context.
    D : Drough -> SPI-2-months <=-1
    H : Heat -> Standardized Temperature anomaly >=1
    DH : Compound event -> D ∩ H"""

    data_dict = {}
    data_dict[crop] = {}

    "P(LowYield|D ∩ H)"
    _, _, Probability_DH = compute_cond_prob(
        dataframe, yield_threshold=0, spi_threshold=-1, heat_threshold=1
    )

    """
    P(LowYield|D)
    """
    _, _, Probability_D = compute_cond_prob(
        dataframe, yield_threshold=0, spi_threshold=-1
    )

    """
    P(LowYield|H)
    """
    _, _, Probability_H = compute_cond_prob(
        dataframe, yield_threshold=0, heat_threshold=1
    )

    """Median/Mean Yield|D ∩ H"""
    # %%
    data_dict[crop]["LowYield|D∩H"] = Probability_DH * 100
    data_dict[crop]["LowYield|D"] = Probability_D * 100
    data_dict[crop]["LowYield|H"] = Probability_H * 100

    output_df = pd.DataFrame(data=data_dict)
    long_df = output_df.stack().reset_index()
    long_df.columns = ["Metric", "Crop", "Value"]

    return long_df


def display_probability(dataframe):
    """Display probability of low yield depending on climate context."""
    plot = sns.catplot(
        data=dataframe,
        x="Metric",
        y="Value",
        kind="bar",
        palette=["#01608e", "#fe5974", "#77c4d5"],
        order=["LowYield|D", "LowYield|H", "LowYield|D∩H"],
    )
    plot.set_ylabels("Probability [%]")
    plot.set_xlabels("Context")
    plot.set_titles("Probability of low yield depending on context for whinter heat")
    plot.savefig("/".join([OUTPUT_PATH, "histogram_probability.png"]))
    plot.figure.clear()

    return None


# %%
def create_context_dataframe(dataframe) -> pd.DataFrame:

    """Create dataframe tagging months according to their climate."""
    condition_d = np.where(dataframe[f"spi_{spi_scale}_months"] <= -1, True, False)
    dataframe.loc[condition_d, "context"] = "D"

    condition_h = np.where(dataframe["standardized T° anomaly"] >= 1, True, False)
    dataframe.loc[condition_h, "context"] = "H"

    condition_d_h = condition_d & condition_h

    dataframe.loc[condition_d_h, "context"] = "D∩H"

    return dataframe


# %%


def display_violinplot(dataframe):

    """Display violinplot of observed yield according to their climate context. """
    plot = sns.violinplot(
        data=dataframe,
        x="context",
        y="difference_to_the_mean_%",
        order=["D", "H", "D∩H"],
        palette=["#01608e", "#fe5974", "#77c4d5"],
        inner="quart",
        split=False,
        density_norm="width",
        cut=0,
    )

    plot.set_ylabel("Difference to the mean [%]")

    plot.figure.suptitle(
        t="Distribution of yield difference to the mean depending on context"
    )
    contexts = ["D", "H", "D∩H"]

    for context in contexts:

        median_value = dataframe.loc[
            dataframe["context"] == context, "difference_to_the_mean_%"
        ].median()

        plot.text(x=context, y=-90, s=f"median = {median_value:.0f}%")

    plot.figure.savefig(f"{OUTPUT_PATH}/violinplot_yield.png")  # ,bbox_inches="tight"

    plot.figure.clear()

    return None


# %%

if __name__ == "__main__":

    period_reference = (
        datetime(2005, 1, 1),
        datetime(2022, 12, 31),
    )

    climate_df = read_climate_data()

    yield_df = computing_z_score(read_yield_data(), period_reference)

    display_lineplot_yield(yield_df)

    display_lineplot_z_score(yield_df)

    california_shapefile = read_shapefile_california()

    completed_df = average_over_county(
        california_shapefile, yield_df, climate_df, period_reference
    )

    display_jointplot(completed_df)

    display_jointplot_differend_std(completed_df)

    probability_df = compute_probability_context(completed_df)

    display_probability(probability_df)

    create_context_dataframe(completed_df)

    display_violinplot(completed_df)
