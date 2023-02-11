# Limão
Projected direct sunlight intensities over time for arbitrary building elevations. See how direct sunlight varies throughout the year for windows at different orientations, and heights above ground level, or compare sunlight intensities for different locations.

Primarily this depends on the the building orientation, the latitude, and the time of the year, but also on occluding obstacles. Limão incorporates raw time and location dependent sunlight intensity projections, and digital surface and terrain maps, to determine an estimated overall intensity as a function of time.

(Currently this assumes UK DSM/DTMs as the lat/lon is converted into the UK OS map coordinates.)

# Usage

Requires digital surface models (DSM) and digital terrain models (DTM) for the region of interest, in addition to the latitude and longitude. For the UK these are available from [DEFRA](https://www.data.gov.uk/dataset/f0db0249-f17b-4036-9e65-309148c97ce4/national-lidar-programme) ([direct link](https://environment.data.gov.uk/DefraDataDownload/?Mode=survey)), and are derived from satellite LIDAR measurements. Latitude and longitude can be obtained from Google Maps.

At the moment this ignores occlusions within a configurable distance from the indicated location to avoid 'self' intersections, but this will probably be improved later.

```bash

usage: limao.py [-h] [-s FILENAMEDSM] [-t FILENAMEDTM] [-lat LAT] [-lon LON] [--size SIZE]

optional arguments:
  -h, --help      show this help message and exit
  -s FILENAMEDSM  DSM input file (.tif).
  -t FILENAMEDTM  DTM input file (.tif).
  -lat LAT        Latitude.
  -lon LON        Longitude.
  --size SIZE     Region Size around location.

```
