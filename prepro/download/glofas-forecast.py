#!/usr/bin/env python

import cdsapi
import datetime
import pathlib
import csv
import pandas
import xarray

start_date = datetime.datetime(2025,4,6)
# max lat, min lon, min lat, max lon
area = [ 47.5, -10, 30, 42 ]
rivers = {
"rhone": (4.8499526,	43.3499475),
"ebre":	(0.8999684,	40.7166247),
"po":	(12.5499218,	44.9666077),
"nil":	(30.3665172,	31.5166615),
"marmara1":	(28.9831894,	41.0166235),
"marmara2":	(28.999856,	41.0166235),
}

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

outpathname = pathlib.Path(f"../../../GLOFAS/glofas_{start_date:%Y%m%d}.nc")
outpathname.parent.mkdir(parents=True, exist_ok=True)

#client = cdsapi.Client()
#client.retrieve(dataset, request).download(outpathname)

ds = xarray.open_dataset(outpathname)
ds = ds.squeeze("forecast_reference_time")
ds = ds.drop_vars(["number", "surface","forecast_reference_time", "forecast_period"])
ds = ds.rename({"forecast_period":"valid_time"})
ds = ds.set_xindex(["valid_time"]).rename({"valid_time":"time", "dis24": "Qbar"})
for name, pos in rivers.items():
     ilon = ds.indexes["longitude"].get_indexer([pos[0]],method="nearest")[0]
     ilat = ds.indexes["latitude"].get_indexer([pos[1]],method="nearest")[0]
     print(name, ilon,ilat)
     ds_data = ds.isel(longitude=slice(ilon-2,ilon+3), latitude=slice(ilat-2,ilat+3)).max(["longitude", "latitude"])
     df_data = ds_data.to_dataframe()
     df_data.index = df_data.index - pandas.Timedelta(12, unit="h")
     outpname = pathlib.Path(f"../../../GLOFAS/{name}_{start_date:%Y%m%d}.csv")
     outpname = pathlib.Path(f"../../../GLOFAS/{name}.dat")
     outpname.parent.mkdir(parents=True, exist_ok=True)
     with open(outpname, mode='w') as f:
         f.write("#\n")
         f.write("#\n")
         f.write("#\n")
         f.write("#\n")
     df_data.to_csv(outpname, mode="a", sep=" ", header=False, quotechar=" ", date_format="%d/%m/%Y %H:%M:%S")

