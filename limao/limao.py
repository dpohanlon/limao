from tqdm import tqdm

import numpy as np

import pandas as pd
from pprint import pprint

import h5py
import json

import functools

import rioxarray as rxr
import earthpy as et

from skimage.draw import line

from bng_latlon import WGS84toOSGB36 as latlon_to_os
from bng_latlon import OSGB36toWGS84 as os_to_latlon

from pysolar.solar import get_altitude, get_azimuth
from pysolar.radiation import get_radiation_direct

from datetime import datetime, timezone, timedelta


@functools.lru_cache
def load_file(fileName):
    return rxr.open_rasterio(fileName, masked=True)


class Limao(object):
    def __init__(
        self, fileNameDSM, fileNameDTM, locLL, size, surfHeight=0.0, distHack=True
    ):

        self.fileNameDSM = fileNameDSM
        self.fileNameDTM = fileNameDTM
        self.locLL = locLL
        self.size = size

        self.distHack = distHack

        self.surfHeight = surfHeight

        self.startDate = datetime(2020, 1, 1, 0, 0, tzinfo=timezone.utc)

        self.locIdxX = self.size // 2
        self.locIdxY = self.size // 2

        # self.data = rxr.open_rasterio(self.fileNameDSM, masked=True)
        # self.dataDTM = rxr.open_rasterio(self.fileNameDTM, masked=True)

        self.data = load_file(self.fileNameDSM)
        self.dataDTM = load_file(self.fileNameDTM)

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

    def selectArea(self, data, loc, size=150):

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

        if self.distHack:

            # A hack so that the building at the location doesn't occlude itself
            dist[dist < 10] = 25

        return dist

    def observedSizes(self, region, regionDTM, surfHeight=0.0):

        size = region.shape[0]

        dist = self.distanceMatrix(size)

        # Only consider relative heights, using the terrain map for the 'ground'
        # floor, adjusting for height above the surface
        mags = (region - regionDTM[self.locIdxX][self.locIdxY] - surfHeight) / dist

        return mags

    def loadRegion(self, data, locLL, size=100):

        # Convert lat/lon to UK OS coordinates for the DSM/DTM

        locOS = latlon_to_os(*locLL)

        region = self.selectArea(data, locOS, size=size)

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

    def yearlyIntensity(self, size, progress=True):

        # Get surface and terrain maps

        region = self.loadRegion(self.data, self.locLL, size=size)

        # plt.imshow(region)
        # plt.savefig("data.pdf")
        # plt.clf()

        # f = h5py.File("data.h5", "w")
        # f.create_dataset("data", data=region)
        # f.close()

        # json.dump(region.tolist(), open("data.json", "w"), indent=4)

        regionDTM = self.loadRegion(self.dataDTM, self.locLL, size=size)

        # plt.imshow(region)
        # plt.savefig("region.png")
        # plt.clf()
        #
        # plt.imshow(regionDTM)
        # plt.savefig("regionDTM.png")
        # plt.clf()

        # Calculate observed heights of objects

        mags = self.observedSizes(region, regionDTM, self.surfHeight)

        # plt.imshow(mags)
        # plt.savefig("mags.png")
        # plt.clf()

        hoursInYear = 24 * 7 * 52

        altitudes = np.zeros(hoursInYear)
        azimuths = np.zeros(hoursInYear)
        intensities = np.zeros(hoursInYear)
        occluded = np.zeros(hoursInYear)
        dates = []

        date = self.startDate

        # For each hour over a year

        for i in tqdm(range(hoursInYear), disable=not progress):

            # Get sun properties for this date and location

            altitude = self.getAltitude(date, self.locLL)
            azimuth = self.getAzmuth(date, self.locLL)

            radiation = self.getRadiation(date, altitude)

            # Calculate the threshold below which `mags` occludes the sun

            threshold = self.getThreshold(altitude)

            # Draw a line (2D) from the sun to the test location

            ir, ic = line(
                self.locIdxX, self.locIdxY, *self.getPixelBoundary(region, azimuth)
            )

            # If any pixel is above the threshold, the direction is occluded

            if np.all(mags[ir, ic] < threshold):
                occluded[i] = False
            else:
                occluded[i] = True

            altitudes[i] = altitude
            azimuths[i] = azimuth
            intensities[i] = radiation
            dates.append(date)

            date = date + timedelta(hours=1)

        return altitudes, azimuths, intensities, occluded, np.array(dates)

    def yearlyIntensityTable(self, progress=True):

        altitudes, azimuths, intensities, occluded, dates = self.yearlyIntensity(
            self.size, progress
        )

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

        intensityTable["day"] = intensityTable["date"].map(
            lambda x: x.timetuple().tm_yday
        )
        intensityTable["week"] = intensityTable["date"].map(
            lambda x: x.isocalendar()[1]
        )

        return intensityTable

    def intensityOnElevation(self, intensityTable, elevationOrientation=0):

        # Split the overall intensity into two faces, each side of an elevation

        # North, elevationOrientation = 0, going clockwise
        # How to specify the 'normal'?
        left_az = np.mod(90 + elevationOrientation, 360)
        right_az = np.mod(270 + elevationOrientation, 360)

        az_north = (intensityTable["azimuth"] < left_az) | (
            intensityTable["azimuth"] > right_az
        )

        return (
            az_north,
            intensityTable[az_north]["intensity_passed"],
            intensityTable[~az_north]["intensity_passed"],
        )
