#!/usr/bin/env python

import datetime
import pathlib
import pandas
import xarray
import copernicusmarine
import croco_class

start_date = "2025-03-17 00:00"
pd_start_date = pandas.Timestamp(start_date)
pd_end_date = pd_start_date + datetime.timedelta(days=9,hours=18)

grdpathname = "../../croco_gibrtwo_inno_energy_grd.nc"
destpath = "/home/coat/dev/CROCO/MERCATOR"

datasets = [
   { "id": "cmems_mod_glo_phy-cur_anfc_0.083deg_PT6H-i", "vars": [ "uo", "vo" ] },
   { "id": "cmems_mod_glo_phy-thetao_anfc_0.083deg_PT6H-i", "vars": [ "thetao" ] },
   { "id": "cmems_mod_glo_phy-so_anfc_0.083deg_PT6H-i", "vars": [ "so" ] },
   { "id": "cmems_mod_glo_phy_anfc_merged-sl_PT1H-i", "vars": [ "sea_surface_height" ] },
   ]

grd = croco_class.CROCO_grd(grdpathname)

buffer_zone = 0.5

cm_request = {
    "username" : "acoat",
    "password" : "JAjyteva", 
    "start_datetime": pd_start_date,
    "end_datetime": pd_end_date,
    "minimum_longitude": grd.lonmin() - buffer_zone,
    "maximum_longitude": grd.lonmax() + buffer_zone,
    "minimum_latitude": grd.latmin() - buffer_zone,
    "maximum_latitude": grd.latmax() + buffer_zone,
    "disable_progress_bar": True,
    "no_metadata_cache": True,
    "overwrite_metadata_cache": False,
}

time_ax = None
allds = {}
for dataset in datasets:
     req = cm_request.copy()
     req["dataset_id"] = dataset["id"]
     req["variables"] = dataset["vars"]
     ds = copernicusmarine.open_dataset(**req)
     if "sea_surface_height" in dataset["vars"]:
         ds = ds.rename({"sea_surface_height": "zos"}).squeeze("depth")
     else:
         time_ax = ds.time
     allds[dataset["id"]] = ds

selds = []
for did, ds in allds.items():
    if time_ax is not None:
        ds = ds.sel(time=time_ax)
    encoding =  { var : { "zlib": True, "shuffle": True} for var in ds.data_vars }
    outname = f"{did}_{pd_start_date:%Y%m%d%H}.nc"
    outpath = pathlib.Path(destpath) / outname
    outpath.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(outpath, encoding=encoding)

