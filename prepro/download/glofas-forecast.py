#!/usr/bin/env python

import cdsapi
import datetime
import pathlib

start_date = datetime.datetime(2025,4,6)
# max lat, min lon, min lat, max lon
area = [ 47.5, -10, 30, 42 ]

dataset = "cems-glofas-forecast"
request = {
    "system_version": ["operational"],
    "hydrological_model": ["lisflood"],
    "product_type": ["control_forecast"],
    "variable": "river_discharge_in_the_last_24_hours",
    "year": [f"{start_date:%Y}"],
    "month": [f"{start_date:%m}"],
    "day": [f"{start_date:%d}"],
    "leadtime_hour": [
        "24",
        "48",
        "72",
        "96",
        "120",
        "144",
        "168",
        "192",
        "216",
        "240",
        "264",
        "288",
        "312",
        "336",
        "360",
        "384",
        "408",
        "432",
        "456",
        "480",
        "504",
        "528",
        "552",
        "576",
        "600",
        "624",
        "648",
        "672",
        "696",
        "720"
    ],
    "area" : area,
    "data_format": "netcdf",
    "download_format": "unarchived"
}

client = cdsapi.Client(url = "https://ewds.climate.copernicus.eu/api", key = "778b31b7-7f4a-475b-a46b-476e50fc17a7")

outpathname = pathlib.Path(f"../../../GLOFAS/glofas_{start_date:%Y%m%d}.nc")

outpathname.parent.mkdir(parents=True, exist_ok=True)

client.retrieve(dataset, request).download(outpathname)


