import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

from tqdm import tqdm

rcParams["axes.facecolor"] = "FFFFFF"
rcParams["savefig.facecolor"] = "FFFFFF"
rcParams["xtick.direction"] = "in"
rcParams["ytick.direction"] = "in"

rcParams.update({"figure.autolayout": True})

import numpy as np

import pandas as pd
from pprint import pprint

import h5py
import json

import seaborn as sns

colors = sns.color_palette("Set2")

import argparse

import rioxarray as rxr

from skimage.draw import line

from bng_latlon import WGS84toOSGB36 as latlon_to_os
from bng_latlon import OSGB36toWGS84 as os_to_latlon

from pysolar.solar import get_altitude, get_azimuth
from pysolar.radiation import get_radiation_direct

from datetime import datetime, timezone, timedelta

from limao.limao import Limao


def intensityProjection(fileNameDSM, fileNameDTM, latLon, size, horizontal, vertical):

    nx = horizontal
    ny = vertical

    intensities = np.zeros((nx, ny))

    locOS = latlon_to_os(*latLon)

    leftOS = list(reversed([(locOS[0] - i, locOS[1]) for i in range(nx // 2)]))
    rightOS = [(locOS[0] + i, locOS[1]) for i in range(nx // 2)]

    osCoords = leftOS + rightOS
    llCoords = [os_to_latlon(*osCoord) for osCoord in osCoords]

    altitudes = [i for i in range(ny)]

    for i in tqdm(range(nx)):

        loc = llCoords[i]

        for j in range(ny):

            alt = altitudes[j]

            limao = Limao(fileNameDSM, fileNameDTM, loc, size, surfHeight=alt)

            table = limao.yearlyIntensityTable(progress=False)
            isNorth, _, _ = limao.intensityOnElevation(table)
            table["isNorth"] = isNorth

            northIntensity = table[~table["isNorth"]]["intensity_passed"]

            meanIntensity = northIntensity.mean()

            intensities[i][j] = meanIntensity

    plt.imshow(
        intensities.T,
        origin="lower",
        cmap="plasma",
        extent=(0, nx, 0, ny),
        interpolation="bilinear",
    )

    plt.colorbar().set_label(label="Average yearly intensity $(W/m^2)$", size=14)
    plt.ylabel("z height $(m)$", fontsize=14)
    plt.xlabel("x extent $(m)$", fontsize=14)

    plt.savefig("projIntensities.pdf")
    plt.clf()

    import pickle

    pickle.dump(intensities, open("intensities.pkl", "wb"))


def dailyAvgIntensity(fileNameDSM, fileNameDTM, latLon, size):

    limao = Limao(fileNameDSM, fileNameDTM, latLon, size, surfHeight=2.0)

    table = limao.yearlyIntensityTable()

    # plt.plot(table[table["altitude"] > 0]["azimuth"], ".")
    # plt.plot(table[table["altitude"] < 0]["azimuth"], ".")
    # plt.savefig("az.pdf")
    # plt.clf()
    #
    # plt.plot(table["altitude"], ".")
    # plt.savefig("alt.pdf")
    # plt.clf()

    isNorth, _, _ = limao.intensityOnElevation(table)
    table["isNorth"] = isNorth

    # plt.plot(table[table["isNorth"]]["azimuth"], ".")
    # plt.savefig("north_az.pdf")
    # plt.clf()
    #
    # plt.plot(table[~table["isNorth"]]["intensity_passed"], alpha=0.5, label="South")
    # plt.plot(table[table["isNorth"]]["intensity_passed"], alpha=0.5, label="North")
    #
    # plt.xlabel("Hours from 1/1", fontsize=14)
    # plt.ylabel("Direct sunlight intensity $(W/m^2)$", fontsize=14)
    # plt.legend(loc=0, fontsize=14)
    # plt.savefig("test.pdf")
    # plt.clf()

    tableDayAvg = (
        table.groupby(["day", "isNorth"])
        .agg(
            {
                "intensity_passed": ["mean", "std", "min", "max"],
                "altitude": ["mean", "std", "min", "max"],
                "azimuth": ["mean", "std", "min", "max"],
            }
        )
        .reset_index()
    )

    tableWeekAvg = (
        table.groupby(["week", "isNorth"])
        .agg(
            {
                "intensity_passed": ["mean", "std", "min", "max"],
                "altitude": ["mean", "std", "min", "max"],
                "azimuth": ["mean", "std", "min", "max"],
            }
        )
        .reset_index()
    )

    plt.plot(
        tableDayAvg[~tableDayAvg["isNorth"]]["day"],
        tableDayAvg[~tableDayAvg["isNorth"]][("intensity_passed", "mean")],
        alpha=0.5,
        label="South",
    )
    plt.plot(
        tableDayAvg[tableDayAvg["isNorth"]]["day"],
        tableDayAvg[tableDayAvg["isNorth"]][("intensity_passed", "mean")],
        alpha=0.5,
        label="North",
    )

    # plt.plot(tableDayAvg['day'], tableDayAvg['altitude'])
    # plt.plot(tableDayAvg["day"], tableDayAvg["azimuth"])

    plt.xlabel("Day", fontsize=14)
    plt.ylabel("Direct sunlight intensity $(W/m^2)$", fontsize=14)
    plt.legend(loc=0, fontsize=14)
    plt.savefig("dayAvg.pdf")
    plt.clf()

    plt.plot(
        tableWeekAvg[~tableWeekAvg["isNorth"]]["week"],
        tableWeekAvg[~tableWeekAvg["isNorth"]][("intensity_passed", "mean")],
        alpha=0.5,
        label="South",
    )
    plt.plot(
        tableWeekAvg[tableWeekAvg["isNorth"]]["week"],
        tableWeekAvg[tableWeekAvg["isNorth"]][("intensity_passed", "mean")],
        alpha=0.5,
        label="North",
    )
    plt.fill_between(
        tableWeekAvg[tableWeekAvg["isNorth"]]["week"],
        tableWeekAvg[tableWeekAvg["isNorth"]][("intensity_passed", "min")],
        tableWeekAvg[tableWeekAvg["isNorth"]][("intensity_passed", "max")],
        alpha=0.05,
        color="orange",
    )
    plt.fill_between(
        tableWeekAvg[~tableWeekAvg["isNorth"]]["week"],
        tableWeekAvg[~tableWeekAvg["isNorth"]][("intensity_passed", "min")],
        tableWeekAvg[~tableWeekAvg["isNorth"]][("intensity_passed", "max")],
        alpha=0.05,
        color="blue",
    )

    # plt.plot(tableDayAvg['day'], tableDayAvg['altitude'])
    # plt.plot(tableDayAvg["day"], tableDayAvg["azimuth"])

    plt.xlabel("Week number", fontsize=14)
    plt.ylabel("Direct sunlight intensity $(W/m^2)$", fontsize=14)
    plt.legend(loc=0, fontsize=14)
    plt.savefig("weekAvg.pdf", dpi = 300)
    plt.savefig("weekAvg.png", dpi = 300)
    plt.clf()

def run():

    # I'd like an argument, please
    argParser = argparse.ArgumentParser()

    argParser.add_argument(
        "-s", type=str, dest="fileNameDSM", default="", help="DSM input file."
    )

    argParser.add_argument(
        "-t", type=str, dest="fileNameDTM", default="", help="DTM input file."
    )

    argParser.add_argument(
        "--lat", type=float, dest="lat", default=None, help="Latitude."
    )

    argParser.add_argument(
        "--lon", type=float, dest="lon", default=None, help="Longitude."
    )

    argParser.add_argument(
        "--size",
        type=int,
        dest="size",
        default=100,
        help="Region Size around location.",
    )

    argParser.add_argument(
        "-proj",
        action="store_true",
        dest="proj",
        default=False,
        help="Plot a spatial 2d projection map of intensity.",
    )

    argParser.add_argument(
        "--vertical",
        type=int,
        dest="vertical",
        default=10,
        help="Vertical extent of the projection.",
    )

    argParser.add_argument(
        "--horizontal",
        type=int,
        dest="horizontal",
        default=10,
        help="Horizontal extent of the projection.",
    )

    args = argParser.parse_args()

    if args.proj:

        intensityProjection(
            args.fileNameDSM, args.fileNameDTM, (args.lat, args.lon), args.size, args.horizontal, args.vertical
        )

    else:

        dailyAvgIntensity(
            args.fileNameDSM, args.fileNameDTM, (args.lat, args.lon), args.size
        )


if __name__ == "__main__":

    run()
