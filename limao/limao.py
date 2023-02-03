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

import seaborn as sns

colors = sns.color_palette("Set2")

import argparse

import rioxarray as rxr
import earthpy as et

from skimage.draw import line

from bng_latlon import WGS84toOSGB36 as latlon_to_os
from bng_latlon import OSGB36toWGS84 as os_to_latlon

from pysolar.solar import get_altitude, get_azimuth
from pysolar.radiation import get_radiation_direct

from datetime import datetime, timezone, timedelta


class Limao(object):
    def __init__(self, fileNameDSM, fileNameDTM, locLL, size, buildingAltitude=0):

        self.fileNameDSM = fileNameDSM
        self.fileNameDTM = fileNameDTM
        self.locLL = locLL
        self.size = size

        # Does nothing, yet
        self.buildingAltitude = buildingAltitude

        self.startDate = datetime(2023, 1, 1, 0, 0, tzinfo=timezone.utc)

        self.locIdxX = self.size // 2
        self.locIdxY = self.size // 2

        self.data = rxr.open_rasterio(self.fileNameDSM, masked=True)
        self.dataDTM = rxr.open_rasterio(self.fileNameDTM, masked=True)

        self.setDataParams(self.data)

    def setDataParams(self, data):

        dx_min, dy_min, dx_max, dy_max = data.rio.bounds()

        self.dataReal = dx_max - dx_min
        self.dataPx = data.data.squeeze().shape[0]

        # Each real unit is this many pixels
        self.scaleToPx = self.dataReal / self.dataPx

        self.scaleToReal = self.dataPx / self.dataReal

        self.realBounds = data.rio.bounds()
        self.pxBounds = data.data.squeeze().shape

    def loc2px(self, data, loc):

        # Find where the location is in local real coordinates
        # and convert to pixels

        pxX = int((loc[0] - self.realBounds[0]) / self.scaleToPx)
        pxY = int((loc[1] - self.realBounds[1]) / self.scaleToPx)

        # Indexing is opposite - do this in a better way

        return self.dataPx - pxY, pxX

    def selectArea(self, data, loc, size=250):

        # Return a size x size region with loc in the centre

        # Location in real space

        pxLoc = self.loc2px(data, loc)

        # Get pixel bounds

        minIdxX = pxLoc[0] - size // 2
        maxIdxX = pxLoc[0] + size // 2

        minIdxY = pxLoc[1] - size // 2
        maxIdxY = pxLoc[1] + size // 2

        arr = data.data.squeeze()

        return arr[minIdxX:maxIdxX, :][:, minIdxY:maxIdxY].copy()

    def getPixelBoundary(self, region, azimuth):

        # Get the boundary region pixel on the line from the sun to the location

        # Azimuth in [-180, 180] clockwise (North = 0) IN DEGREES

        if azimuth > 180:
            azimuth -= 360

        sign = np.sign(azimuth)
        azimuth = np.abs(azimuth)

        angle = np.mod(azimuth, 45.0)
        octant = np.floor(azimuth / 45.0)

        regionSize = region.shape[0]
        quadrantSize = regionSize // 2

        maxIdx = regionSize - 1

        # Fraction along axis (from centre)
        ratio = np.tan(np.radians(angle))

        if octant == 0:
            r = 0
            c = int(quadrantSize + sign * ratio * quadrantSize)
        elif octant == 1:
            c = maxIdx if sign == 1 else 0
            r = int(ratio * quadrantSize)
        elif octant == 2:
            c = maxIdx if sign == 1 else 0
            r = int(quadrantSize + ratio * quadrantSize)
        elif octant == 3:
            c = (
                int(maxIdx - ratio * quadrantSize)
                if sign == 1
                else int(ratio * quadrantSize)
            )
            r = maxIdx

        return r, c

    def distanceMatrix(self, size, realDistance=True):

        # Create distance matrix with size, filled with pixel-wise distances
        # from the centre (query location)

        xx, yy = np.meshgrid(range(size), range(size))

        # Distance from loc
        dist = np.sqrt((xx - size // 2) ** 2 + (yy - size // 2) ** 2)

        # pixel to real distance
        dist *= self.scaleToPx

        # A hack so that the building at the location doesn't occlude itself
        dist[dist < 10] = 25

        return dist

    def observedSizes(self, region, regionDTM):

        size = region.shape[0]

        dist = self.distanceMatrix(size)

        # Only consider relative heights, using the terrain map for the 'ground'
        # floor
        mags = (region - regionDTM[self.locIdxX][self.locIdxY]) / dist

        return mags

    def loadRegion(self, data, locLL):

        locOS = latlon_to_os(*locLL)

        region = self.selectArea(data, locOS)

        return region

    def getRadiation(self, date, alt):
        return get_radiation_direct(date, alt) if alt > 0 else 0.0

    def getAzmuth(self, date, locLL):
        return get_azimuth(*locLL, date)

    def getAltitude(self, date, locLL):
        return get_altitude(*locLL, date)

    def getThreshold(self, alt):
        # Threshold of height / dist at which it occludes
        # alt IN DEGREES
        return np.tan(np.radians(alt))

    def yearlyIntensity(self):

        region = self.loadRegion(self.data, self.locLL)
        regionDTM = self.loadRegion(self.dataDTM, self.locLL)

        mags = self.observedSizes(region, regionDTM)

        hoursInYear = 24 * 7 * 52

        altitudes = np.zeros(hoursInYear)
        azimuths = np.zeros(hoursInYear)
        intensities = np.zeros(hoursInYear)
        occluded = np.zeros(hoursInYear)
        dates = []

        date = self.startDate

        for i in tqdm(range(hoursInYear)):

            altitude = self.getAltitude(date, self.locLL)
            azimuth = self.getAzmuth(date, self.locLL)

            radiation = self.getRadiation(date, altitude)

            threshold = self.getThreshold(altitude)

            ir, ic = line(
                self.locIdxX, self.locIdxY, *self.getPixelBoundary(region, azimuth)
            )

            if np.all(mags[ir, ic] < threshold):
                occluded[i] = True
            else:
                occluded[i] = False

            altitudes[i] = altitude
            azimuths[i] = azimuth
            intensities[i] = radiation
            dates.append(date)

            date = date + timedelta(hours=1)

        return altitudes, azimuths, intensities, occluded, np.array(dates)

    def yearlyIntensityTable(self):

        altitudes, azimuths, intensities, occluded, dates = self.yearlyIntensity()

        intensityTable = pd.DataFrame(
            {
                "date": dates,
                "altitude": altitudes,
                "azimuth": azimuths,
                "intensity": intensities,
                "occluded": occluded,
            }
        )

        intensityTable["intensity_passed"] = np.where(
            intensityTable["occluded"], 0, intensityTable["intensity"]
        )

        return intensityTable

    def intensityOnElevation(self, intensityTable, elevationOrientation):

        # North, elevationOrientation = 0, going clockwise
        # How to specify the 'normal'?
        left_az = np.mod(90 + elevationOrientation, 360)
        right_az = np.mod(270 + elevationOrientation, 360)

        az_north = (intensityTable["azimuth"] < left_az) | (
            intensityTable["azimuth"] > right_az
        )

        return (
            intensityTable[az_north]["intensity_passed"],
            intensityTable[~az_north]["intensity_passed"],
        )


if __name__ == "__main__":

    # I'd like an argument, please
    argParser = argparse.ArgumentParser()

    argParser.add_argument(
        "-s", type=str, dest="fileNameDSM", default="", help="DSM input file."
    )
    argParser.add_argument(
        "-t", type=str, dest="fileNameDTM", default="", help="DTM input file."
    )

    argParser.add_argument(
        "-lat", type=float, dest="lat", default=None, help="Latitutde."
    )
    argParser.add_argument(
        "-lon", type=float, dest="lon", default=None, help="Longitude."
    )

    args = argParser.parse_args()

    limao = Limao(args.fileNameDSM, args.fileNameDTM, (args.lat, args.lon), 100)

    table = limao.yearlyIntensityTable()

    print(table)
