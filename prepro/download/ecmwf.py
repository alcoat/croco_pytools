#!/usr/bin/env python

import herbie
import numpy
import pandas
import cfgrib
import xarray
import croco_class
import pathlib
import metpy.calc

grd = croco_class.CROCO_grd("../../croco_gibrtwo_inno_energy_grd.nc")

minimum_longitude = grd.lonmin() - 0.5
maximum_longitude = grd.lonmax() + 0.5
minimum_latitude = grd.latmin() - 0.5
maximum_latitude = grd.latmax() + 0.5

start_date=pandas.Timestamp("2025-03-17 00:00")

fxx = []
fxx += list(range(0,144, 3))
fxx += list(range(144, 240+6, 6))

searchString='(?:10[uv]:|:2t:|:2d:|:msl:|:ssr:|:str:|:tp:)'

outdir = pathlib.Path("../../../ECMWF")

while True:
    FH = herbie.FastHerbie([start_date], model="ecmwf", product="oper", fxx=fxx, save_dir=outdir)

    files = FH.download(searchString)

    allok = True
    allds = []
    for fpath in files:
        ds = cfgrib.open_datasets(fpath)
        ds = xarray.merge(ds, compat="override").sortby('latitude')
        if len(ds.data_vars)<8:
            print("bad file", fpath, len(ds.data_vars))
            print(ds)
            fpath.unlink()
            alltok = False
            continue
        ds = ds.sel(latitude=slice(minimum_latitude, maximum_latitude), longitude=slice(minimum_longitude,maximum_longitude))
        allds.append(ds)

    if allok is True and len(allds)>0:
        break

ds = xarray.concat(allds, dim="valid_time")
ds = ds.drop('time').rename({'valid_time': 'time'})
ds.time.encoding['units'] = "days since 1900-01-01T00:00:00Z"

# RH: =100*(EXP((17.625*TD)/(243.04+TD))/EXP((17.625*T)/(243.04+T)))
# q :  equation 4.24, Pg 96 Practical Meteorolgy (Roland Stull) https://github.com/Unidata/MetPy/issues/791
#ds["rh2m"] = numpy.exp(17.625*ds["d2m"] / (243.04 + ds["d2m"])) / numpy.exp(17.625*ds["t2m"]/(243.04+ds["t2m"]))
#ds["rh2m"] = xarray.where(ds['rh2m']<1, ds['rh2m'], 1)
ds["hu2m"] = metpy.calc.relative_humidity_from_dewpoint(ds["t2m"], ds["d2m"])
#ds["hu2m"] = xarray.where(ds['hu2m']<1, ds['hu2m'], 1)
ds = ds.drop("d2m")
ds['ssr'] = ds['ssr'].differentiate("time", datetime_unit="s")
ds['ssr'] = xarray.where(ds['ssr']>0, ds['ssr'], 0)
ds['str'] = - ds['str'].differentiate("time", datetime_unit="s")
ds['tp'] = ds['tp'].differentiate("time", datetime_unit="s")

ds = ds.rename({"u10": "u10m", "v10": "v10m", "msl": "pmer", "ssr": "flsolaire", "str": "fltherm", "tp": "prate"})
ds = ds.rename({"latitude": "lat", "longitude": "lon" })
encoding =  { var : { "zlib": True, "shuffle": True} for var in ds.data_vars }
outname = f"ecmwf_{start_date:%Y%m%d%H}.nc"
outpath = outdir / outname
ds.to_netcdf(outpath, encoding=encoding)

