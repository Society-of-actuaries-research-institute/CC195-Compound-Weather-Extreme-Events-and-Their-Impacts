# %%
import warnings
from math import log

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

warnings.filterwarnings("ignore")
import os
from datetime import datetime

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import xarray as xr
from loguru import logger
from numpy import linspace, std
from pandas import Series
from scipy.stats import (boschloo_exact, chi2_contingency, fisk, gamma,
                         genextreme, norm)
from tqdm import tqdm

# %%

OUTPUT_PATH = "./output/wildfires"
INPUT_PATH = "./data/wildfire/climate_data_w_wildfires_normalized.parquet"
if not os.path.exists(OUTPUT_PATH):
    os.makedirs(OUTPUT_PATH)
    logger.info(f"Directory '{OUTPUT_PATH}' created.")
else:
    logger.info(f"Directory '{OUTPUT_PATH}' already exists.")

# %%


def create_dataframe() -> pd.DataFrame:

    """Reading climate data with wildfires added and indices already computed"""

    df = pd.read_parquet(INPUT_PATH)

    index_ninf = df[df["spi_6_months"] == np.NINF].index

    df.loc[index_ninf, "spi_6_months"] = -4

    df = df.rename(
        columns={
            "heatwaves_normalized": "T° anomaly standardized",
        }
    )

    shifted_anomaly = df.groupby(["latitude", "longitude"])[
        "T° anomaly standardized"
    ].shift(1)

    df_w_lag = df.copy()

    df_w_lag["T° anomaly standardized lagged"] = shifted_anomaly

    return df_w_lag

def create_susceptibility_dataframe(dataframe):

    """Keep only places where there as been wildfires in the past as a proxy for places that are susceptible to burn."""

    df_indexed = dataframe.reset_index().set_index(["longitude", "latitude"])

    coords_w_wf = df_indexed[df_indexed.wildfire == 1].index

    df_w_wf_susceptibility = df_indexed.loc[coords_w_wf]

    return df_w_wf_susceptibility


def display_jointplot_where_susceptibility(dataframe):

    """Display jointplot of standardized temperature anomaly and standardized precipitation index """
    g = sns.jointplot(
        data=dataframe,
        x="T° anomaly standardized lagged",
        y="spi_6_months",
        kind="hex",
        xlim=[-3, 3],
        ylim=[-3, 3],
    )
    g.plot_marginals(sns.kdeplot)
    g.refline(x=0, y=0, joint=True)
    g.set_axis_labels(xlabel="T° anomaly standardized", ylabel="spi_6_month")
    g.figure.suptitle(
        t="Joint distribution of SPI 6 and heatwave over the 1950-2023 ", y=1.05
    )
    plt.colorbar(g.ax_joint.collections[0], ax=g.ax_joint, orientation="horizontal")

    g.savefig(f"{OUTPUT_PATH}/jointplot_where_susceptibility.png")

    g.figure.clear()
    return None


# %%
def display_jointplot_when_wildfires(dataframe):

    """Display jointplot of standardized temperature anomaly and standardized precipitation index for month where wildfire happens."""

    condition = np.where(dataframe["wildfire"] == 1, True, False)

    g = sns.jointplot(
        data=dataframe.loc[condition],
        x="T° anomaly standardized lagged",
        y="spi_6_months",
        kind="hex",
        xlim=[-3, 3],
        ylim=[-3, 3],
    )

    g.plot_marginals(sns.kdeplot, fill=False)
    g.refline(x=0, y=0, joint=True)
    g.set_axis_labels(xlabel="T° anomaly standardized", ylabel="spi_6_month")
    g.figure.suptitle(
        t="Joint distribution of SPI-6-month and heatwave over the 1950-2023\nconsidering wildfire events",
        y=1.05,
    )
    plt.colorbar(g.ax_joint.collections[0], ax=g.ax_joint, orientation="horizontal")
    g.savefig(f"{OUTPUT_PATH}/jointplot_when_wildfires.png")

    g.figure.clear()
    return None


# %%
def create_context(dataframe) -> pd.DataFrame:

    """Create dataframe with months tagged according to the climate context (D : Drought , H : Positive temperature anomaly , D∩H : Compound event)"""
    condition_D = np.where(dataframe["spi_6_months"] <= -1, True, False)
    dataframe.loc[condition_D, "context"] = "D"

    condition_H_lagged = np.where(
        dataframe["T° anomaly standardized lagged"] >= 1, True, False
    )
    dataframe.loc[condition_H_lagged, "context"] = "H"

    condition_D_H = condition_D & condition_H_lagged
    dataframe.loc[condition_D_H, "context"] = "D∩H"

    return dataframe


# %%
def create_season(dataframe) -> pd.DataFrame:

    """Discriminate months between wildfire season and outside wildfire season """

    wildfire_season = (dataframe.month >= 6) & (dataframe.month <= 10)
    outside_wildfire_season = (dataframe.month < 6) | (dataframe.month > 10)

    dataframe.loc[wildfire_season, "wildfire_season"] = True
    dataframe.loc[outside_wildfire_season, "wildfire_season"] = False

    return dataframe


def compute_probability_event_from_season(dataframe) -> pd.DataFrame:

    """Compute probability of wildifre depending on context and wildfire."""

    occurence_context = dataframe.groupby(["wildfire_season", "context"])[
        "wildfire"
    ].count()

    occurence_wildfire = dataframe.groupby(["wildfire_season", "context"])[
        "wildfire"
    ].sum()

    both = pd.DataFrame(occurence_context).join(
        occurence_wildfire, rsuffix="_occurence", lsuffix="_context"
    )

    both["probability"] = both["wildfire_occurence"] / both["wildfire_context"] * 100

    return both


def display_probability_season(dataframe):

    """Display probability of single events and compound events with regards to season (widlfire season vs not wildfire season)."""
    plot = sns.barplot(
        data=dataframe,
        x=dataframe.index.get_level_values("wildfire_season"),
        y="probability",
        hue=dataframe.index.get_level_values("context"),
        hue_order=["D", "H", "D∩H"],
        palette=["#01608e", "#fe5974", "#77c4d5"],
    )
    plot.set_xlabel("Season")
    plot.set_xticklabels(["Outside of wildfire season", "Wildfire season"])
    plot.set_ylabel("Probability [%]")
    plot.figure.savefig(f"{OUTPUT_PATH}/probability_wildfire_per_season.png")
    plot.figure.clear()

    return None


# %%
def compute_probability_event_from_context(dataframe) -> pd.DataFrame:

    """Compute probability of wildfire according to climate context."""
    occurence_context = dataframe.groupby("context")["wildfire"].count()
    occurence_wildfire = dataframe.groupby("context")["wildfire"].sum()
    both = pd.DataFrame(occurence_context).join(
        occurence_wildfire, rsuffix="_occurence", lsuffix="_context"
    )

    both["probability"] = both["wildfire_occurence"] / both["wildfire_context"] * 100

    return both


def display_probability_context(dataframe):

    """Display probability of wildfire according to climate context """

    plot = sns.barplot(
        data=dataframe,
        x=dataframe.index.get_level_values("context"),
        y="probability",
        order=["D", "H", "D∩H"],
        palette=["#01608e", "#fe5974", "#77c4d5"],
    )

    plot.set_ylabel("Probability [%]")

    plot.figure.savefig(f"{OUTPUT_PATH}/probability_wildfire_global.png")

    plot.figure.clear()

    return None


# %%
def display_month_occurence_wildfire(dataframe):

    """Display the occurence of wildfire per month."""

    wildfire_events = dataframe.iloc[np.where(dataframe.wildfire == 1)]

    months_event = [
        date.month for date in wildfire_events.index.get_level_values("time")
    ]

    plot = sns.histplot(x=months_event)

    plot.set_xticklabels(
        ["January", "February", "April", "June", "August", "October", "December"]
    )

    plot.set_title("Occurence of wildfire")

    plot.figure.savefig(f"{OUTPUT_PATH}/month_occurence_wildfire.png")

    plot.figure.clear()
    return None


if __name__ == "__main__":

    logger.info("Creating DataFrame")

    df = create_dataframe()

    display_month_occurence_wildfire(df)

    df_susceptibility = create_susceptibility_dataframe(df)

    logger.info("Saving jointplots")

    display_jointplot_where_susceptibility(df_susceptibility)

    display_jointplot_when_wildfires(df_susceptibility)

    logger.info("Creating susceptibility dataframe")

    df_context = create_context(df_susceptibility)

    df_probability_context = compute_probability_event_from_context(df_context)

    logger.info("Saving probability graph per context")

    display_probability_context(df_probability_context)

    logger.info("Creating season dataframe")

    df_season = create_season(df_context)

    df_probability_season = compute_probability_event_from_season(df_season)

    logger.info("Saving probability graph per season")

    display_probability_season(df_probability_season)
