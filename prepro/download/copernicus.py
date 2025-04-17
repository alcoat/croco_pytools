#!/usr/bin/env python

import argparse
import datetime
import pathlib
import sys
import logging

import copernicusmarine
import croco_class
import pandas
import xarray

logger = logging.getLogger(__name__)

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s %(levelname)s] [%(filename)s:%(lineno)s - %(funcName)s() ] %(message)s",
)

# configname = "croco_med"
# grdpathname = "/home/shom_simuref/CROCO/fournitures_22AC05/MED_1.8/CROCO_FILES/test2.nc"
# destpath = "/home6/datawork/acoat/CROCO/MERCATOR"


def download(start_date, end_date, work_dir, grdpathname, cfg):
    datasets = [
        {"id": "cmems_mod_glo_phy-cur_anfc_0.083deg_PT6H-i", "vars": ["uo", "vo"]},
        {"id": "cmems_mod_glo_phy-thetao_anfc_0.083deg_PT6H-i", "vars": ["thetao"]},
        {"id": "cmems_mod_glo_phy-so_anfc_0.083deg_PT6H-i", "vars": ["so"]},
        {
            "id": "cmems_mod_glo_phy_anfc_merged-sl_PT1H-i",
            "vars": ["sea_surface_height"],
        },
    ]

    grd = croco_class.CROCO_grd(grdpathname)

    buffer_zone = 0.5

    cm_request = {
        "username": "acoat",
        "password": "JAjyteva",
        "start_datetime": start_date,
        "end_datetime": end_date,
        "minimum_longitude": grd.lonmin() - buffer_zone,
        "maximum_longitude": grd.lonmax() + buffer_zone,
        "minimum_latitude": grd.latmin() - buffer_zone,
        "maximum_latitude": grd.latmax() + buffer_zone,
        #"disable_progress_bar": True,
        #"no_metadata_cache": True,
        #"overwrite_metadata_cache": False,
        "chunk_size_limit":0
    }

    time_ax = None
    allds = {}
    for dataset in datasets:
        logger.info("Request dataset %s", dataset)
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
        logger.info("processing %s", did)
        if (ds.time == time_ax).all() is not True:
            if time_ax is not None:
                ds = ds.sel(time=time_ax)
        encoding = {var: {"zlib": True, "shuffle": True} for var in ds.data_vars}
        outname = f"{cfg}_{did}_{start_date:%Y%m%dT%H%M}_{end_date:%Y%m%dT%H%M}.nc"
        outpath = work_dir / outname
        outpath.parent.mkdir(parents=True, exist_ok=True)
        logger.info("Writing %s", outpath)
        ds.to_netcdf(outpath, encoding=encoding)


def get_args():
    parser = argparse.ArgumentParser(
        description="Copernicus Mercator Download",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # add_logging_parser_arguments(parser)
    parser.add_argument(
        "--begindate", type=pandas.Timestamp, help="first date for the croco run"
    )
    parser.add_argument(
        "--enddate", type=pandas.Timestamp, help="end date for the croco run"
    )
    parser.add_argument(
        "--workdir",
        type=pathlib.Path,
        default=pathlib.Path("/home/opsys/data/CROCO/GTT_CROCO/MERCATOR"),
        help="working directory path",
    )
    parser.add_argument(
        "--cfgname",
        type=str,
        default="gibr",
        help="configuration name",
    )
    parser.add_argument(
        "--grd",
        type=pathlib.Path,
        default=pathlib.Path(
            "/home/opsys/data/CROCO/GTT_CROCO/CROCO_FILES/croco_gibrtwo_inno_energy_grd.nc"
        ),
        help="grid pathname",
    )

    args = parser.parse_args()

    return args


def run_main():
    args = get_args()

    if args.begindate is not None:
        begindate = args.begindate
    else:
        begindate = pandas.Timestamp.utcnow().replace(
            hour=0, minute=0, second=0, microsecond=0, nanosecond=0, tzinfo=None
        ) - pandas.Timedelta(1, unit="day")
    if args.enddate is not None:
        enddate = args.enddate
    else:
        enddate = (
            begindate
            + pandas.Timedelta(5, unit="day")
            #+ pandas.Timedelta(10, unit="day")
            #- pandas.Timedelta(6, unit="hour")
        )

    download(begindate, enddate, args.workdir, args.grd, args.cfgname)

    return 0


if __name__ == "__main__":
    sys.exit(run_main())
