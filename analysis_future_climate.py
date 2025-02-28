# %%
import os
from datetime import datetime
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import xarray as xr
from loguru import logger
from numpy import std
from pandas import Series, to_datetime
from scipy.stats import fisk, gamma, genextreme, norm

# %%
OUTPUT_PATH = "./output/future"
INPUT_PATH = "./data/Climate-models"

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


def calculate_SPI(
    dataframe,
    dist_type,
    scale,
):

    """Computing the SPI using a distribution type for precipitation and a scale in months"""
    df = dataframe.copy().droplevel([1, 2])
    historic_index = df.timeframe == "historic"
    df["tp_rolling"] = df["tp"].rolling(scale).mean()
    df["si"] = calculate_si(
        series_ref=df.loc[historic_index, "tp_rolling"],
        series_calc=df["tp_rolling"],
        dist_type=dist_type,
        by_month=True,
    )
    return df["si"]


def anomaly_temperature(dataframe):
    """Compute the standardized temperature anomaly using the historic timeframe"""
    index_historic = dataframe.timeframe == "historic"
    mean_reference = dataframe.loc[index_historic, "tas"].mean()
    std_reference = dataframe.loc[index_historic, "tas"].std()
    dataframe["standardized T° anomaly"] = (
        dataframe.loc[:, "tas"] - mean_reference
    ) / std_reference
    return dataframe


def count_ce_se(dataframe):
    """Count the number of month with Compound Events (CE) or Single Events (SE)"""

    heatwave_condition = dataframe["standardized T° anomaly"] >= 1

    drought_condition = dataframe["spi_6_months"] <= -1

    compound_event_condition = heatwave_condition & drought_condition

    count_ce = dataframe.loc[compound_event_condition, "tas"].count()

    count_heatwave = dataframe.loc[heatwave_condition, "tas"].count()

    count_drought = dataframe.loc[drought_condition, "tas"].count()

    count_total = len(dataframe)

    return pd.DataFrame(
        index=["D ∩ H", "D", "H", "number_month_total"],
        data=[count_ce, count_drought, count_heatwave, count_total],
    )

def create_dataframe() -> pd.DataFrame:

    """Reading the different netcdfs from different models and computing the different indices.

    Climate scenarios used were from the NEX-GDDP-CMIP6 dataset, prepared by the Climate Analytics Group and NASA Ames Research Center using the NASA Earth Exchange and distributed by the NASA Center for Climate Simulation (NCCS). We acknowledge the World Climate Research Programme, which, through its Working Group on Coupled Modeling, coordinated and promoted CMIP6. We thank the climate modeling groups for producing and making available their model output, the Earth System Grid Federation (ESGF) for archiving the data and providing access, and the multiple funding agencies who support CMIP6 and ESGF
 """

    list_netcdf = Path(INPUT_PATH).glob("*.nc")

    list_df = list()

    for netcdf in list_netcdf:
        ds = xr.open_dataset(str(netcdf))
        temp = ds.to_dataframe()
        temp["model"] = netcdf.stem.split("_")[0]
        list_df.append(temp)

    df = pd.concat(list_df)

    df = df.dropna()

    df["timeframe"] = "2050"

    index_historic = df.index.get_level_values("time") <= pd.to_datetime("2014")

    df.loc[index_historic, "timeframe"] = "historic"

    nbr_second_day = 3600 * 24

    df["tp"] = df["pr"] * nbr_second_day

    df["Tavg [C°]"] = df["tas"] - 273

    df["month"] = [date.month for date in df.index.get_level_values("time")]

    spi_6 = df.groupby(["lon", "lat", "model"]).apply(
        lambda x: calculate_SPI(x, dist_type="gamma", scale=6)
    )

    df_ni = df.reset_index().set_index(["lon", "lat", "model", "time"])

    df_ni = df_ni.join([spi_6])

    df_ni = df_ni.rename(columns={"si": f"spi_6_months"})

    dataframe = df_ni.groupby(["lon", "lat", "month", "model"], as_index=False).apply(
        lambda x: anomaly_temperature(x)
    )
    return dataframe


def kde_temperature_historic_2050(dataframe):

    """Display kde plot of standardized temperature anomaly of historic and 2050 data"""
    fig, axes = plt.subplots(2, 3, figsize=(20, 20), sharey=True,sharex=True)
    f_axes = axes.ravel()
    """kde plot T° historic vs 2050"""
    for ax, (model, model_data) in zip(f_axes, dataframe.groupby("model")):
        plot = sns.kdeplot(
            data=model_data,
            x="standardized T° anomaly",
            hue="timeframe",
            fill=True,
            alpha=0.45,
            common_norm=False,
            cut=0,
            ax=ax,
            legend=True if model == "ACCESS-CM2" else None,
            hue_order=["historic", "2050"],
            palette=["#01608e", "#fe5974"],
        )
        plot.set_title(f"{model}",fontdict={"fontsize":"15"})
        plot.set_xlabel("Standardized temperature anomaly",fontdict={"fontsize":"15"})

    plt.savefig(f"{OUTPUT_PATH}/kde_temperature_historic_2050.png",dpi=1000)
    return None


def lineplot_temperature_historic_2050(dataframe):

    """Display line plot of standardized temperature anomaly of historic and 2050 data"""
    fig, axes = plt.subplots(2, 3, figsize=(20, 20), sharex=True,sharey=True)
    f_axes = axes.ravel()
    for ax, (model, model_data) in zip(f_axes, dataframe.groupby("model")):

        plot = sns.lineplot(
            data=model_data,
            y="standardized T° anomaly",
            x="month",
            hue="timeframe",
            ax=ax,
            legend=True if model == "ACCESS-CM2" else None,
            hue_order=["historic", "2050"],
            palette=["#01608e", "#fe5974"],
        )
        plot.set_title(f"{model}",fontdict={"fontsize":"15"})
        plot.set_xlabel("Months",fontdict={"fontsize":"15"})
        plot.set_ylabel("Standardized temperature anomaly",fontdict={"fontsize":"15"})

    plt.savefig(f"{OUTPUT_PATH}/lineplot_temperature_historic_2050.png",dpi=1000)
    return None


def kde_spi_historic_2050(dataframe):
    """Display kde plot of standardized precipitation index of historic and 2050 data"""
    fig, axes = plt.subplots(2, 3, figsize=(20, 20), sharey=True,sharex=True)
    f_axes = axes.ravel()
    for ax, (model, model_data) in zip(f_axes, dataframe.groupby("model")):

        plot = sns.kdeplot(
            data=model_data,
            x="spi_6_months",
            hue="timeframe",
            fill=True,
            alpha=0.45,
            common_norm=False,
            cut=0,
            legend=True if model == "ACCESS-CM2" else None,
            hue_order=["historic", "2050"],
            ax=ax,
            palette=["#01608e", "#fe5974"],
        )
        plot.set_title(f"{model}",fontdict={"fontsize":"15"})
        plot.set_xlabel("spi_6_months",fontdict={"fontsize":"15"})


    plt.savefig(f"{OUTPUT_PATH}/kde_spi_historic_2050.png",dpi=1000)
    return None


def lineplot_spi_historic_2050(dataframe):

    """Display lineplot of standardized precipitation index of historic and 2050 data"""

    fig, axes = plt.subplots(2, 3, figsize=(20, 20), sharex=True,sharey=True)
    f_axes = axes.ravel()
    for ax, (model, model_data) in zip(f_axes, dataframe.groupby("model")):

        plot = sns.lineplot(
            data=model_data,
            y="spi_6_months",
            x="month",
            hue="timeframe",
            ax=ax,
            legend=True if model == "ACCESS-CM2" else None,
            hue_order=["historic", "2050"],
            palette=["#01608e", "#fe5974"],
        )
        plot.set_title(f"{model}",fontdict={"fontsize":"15"})
        plot.set_xlabel("Months",fontdict={"fontsize":"15"})
        plot.set_ylabel("spi-6-month",fontdict={"fontsize":"15"})

    plt.savefig(f"{OUTPUT_PATH}/lineplot_spi_historic_2050.png",dpi=1000)
    return None


def jointplot_historic_2050(dataframe):

    """Display jointplot of standardized precipitation index and standardized precipitation index of historic and 2050 data"""
    for model, model_data in dataframe.groupby("model"):

        plot = sns.jointplot(
            data=model_data,
            y="spi_6_months",
            x="standardized T° anomaly",
            xlim=[-6, 6],
            ylim=[-6, 6],
            hue="timeframe",
            kind="kde",
            fill=True,
            alpha=0.5,
            hue_order=["historic", "2050"],
            palette=["#01608e", "#fe5974"],
        )
        plot.figure.suptitle(f"{model}")
        plot.refline(x=0, y=0)
        plt.savefig(f"{OUTPUT_PATH}/{model}_jointplot.png",dpi=1000)
    return None


# %%
def count_month_dataframe(dataframe) -> pd.DataFrame:

    """Creating the dataframe and computing probability of single events (SE) and compound events (CE) per months."""
    timeseries_count = dataframe.groupby(["model", "month", "timeframe"]).apply(
        lambda x: count_ce_se(x)
    )
    timeseries_count.index.names = ["model", "month", "timeframe", "event"]

    number_month_total = timeseries_count.unstack("event").iloc[0, 3]

    timeseries_probability = timeseries_count / number_month_total * 100

    probability_display = timeseries_probability.swaplevel(3, 0).drop(
        ["number_month_total"], axis=0
    )

    return probability_display


def catplot_month_historic_future_probabilities(dataframe):

    """Display catplot of probability of single events or compound events per month for historic and 2050 data"""
    plot = sns.catplot(
        data=dataframe,
        x="month",
        y=0,
        hue="timeframe",
        kind="bar",
        col="event",
        estimator="median",
        hue_order=["historic", "2050"],
        col_order=["D", "H", "D ∩ H"],
        palette=["#01608e", "#fe5974"])


    plot.set_ylabels("Probability [%]",fontdict={"fontsize":"15"})
    plot.savefig(f"{OUTPUT_PATH}/catplot_month_historic_future_probabilities.png",dpi=1000)
    return None


# %%
def count_dataframe(dataframe) -> pd.DataFrame:

    """Creating the dataframe and computing probability of single events (SE) and compound events (CE) overall."""

    ce_count = dataframe.groupby(["model", "timeframe"]).apply(lambda x: count_ce_se(x))

    ce_count.index.names = ["model", "timeframe", "event"]

    number_month_total = ce_count.iloc[3, 0]

    ce_probability = ce_count / number_month_total * 100

    ce_probability_display = ce_probability.swaplevel(0, 2).drop("number_month_total")

    return ce_probability_display


def catplot_event_historic_future_probabilities(dataframe):

    """Display catplot of probability of single events or compound events over the whole year for historic and 2050 data"""

    plot = sns.catplot(
        data=dataframe,
        x="event",
        y=0,
        hue="timeframe",
        kind="bar",
        estimator="median",
        hue_order=["historic", "2050"],
        order=["D", "H", "D ∩ H"],
        palette=["#01608e", "#fe5974"],
    )

    plot.set_ylabels("Probability [%]")
    plot.figure.savefig(
        f"{OUTPUT_PATH}/catplot_event_historic_future_probabilities.png",dpi=1000
    )
    return None


# %%
if __name__ == "__main__":

    logger.info("Creating the dataframe")
    df = create_dataframe()

    logger.info("Plotting kde plot on temperature")
    kde_temperature_historic_2050(df)

    logger.info("Plotting lineplot on temperature")
    lineplot_temperature_historic_2050(df)

    logger.info("Plotting kde plot on spi-6-month")
    kde_spi_historic_2050(df)

    logger.info("Plotting lineplot on spi-6-month")
    lineplot_spi_historic_2050(df)

    logger.info("Plotting jointplot on each model")
    jointplot_historic_2050(df)

    df_count_per_month = count_month_dataframe(df)

    logger.info("Plotting catplot for each month")
    catplot_month_historic_future_probabilities(df_count_per_month)

    df_count = count_dataframe(df)

    logger.info("Plotting catplot whole year round")
    catplot_event_historic_future_probabilities(df_count)
