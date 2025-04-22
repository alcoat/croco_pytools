#!/usr/bin/env python
"""
Download ECMWF
"""
import logging
import pathlib
import sys
import argparse

import cfgrib
import croco_class
import herbie
import metpy.calc
import pandas
import xarray

logger = logging.getLogger(__name__)

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s %(levelname)s] [%(filename)s:%(lineno)s - %(funcName)s() ] %(message)s",
)


def download(start_date, end_date, work_dir, grdpathname, outpath):
    """ download data on grid """
    grd = croco_class.CROCO_grd(grdpathname)

    minimum_longitude = grd.lonmin() - 0.5
    maximum_longitude = grd.lonmax() + 0.5
    minimum_latitude = grd.latmin() - 0.5
    maximum_latitude = grd.latmax() + 0.5

    fxx = []
    fxx += list(range(0, 144, 3))
    fxx += list(range(144, 240 + 6, 6))

    for i, f in enumerate(fxx):
        if start_date + pandas.Timedelta(f, unit="h") >= end_date:
            break
    if i+1<len(fxx):
        fxx = fxx[:i+2]
    searchString = "(?:10[uv]:|:2t:|:2d:|:msl:|:ssr:|:str:|:tp:)"

    while True:
        FH = herbie.FastHerbie(
            [start_date], model="ecmwf", product="oper", fxx=fxx, save_dir=work_dir
        )

        files = FH.download(searchString)

        allok = True
        allds = []
        for fpath in files:
            ds = cfgrib.open_datasets(fpath)
            ds = xarray.merge(ds, compat="override").sortby("latitude")
            if len(ds.data_vars) < 8:
                print("bad file", fpath, len(ds.data_vars))
                print(ds)
                fpath.unlink()
                alltok = False
                continue
            ds = ds.sel(
                latitude=slice(minimum_latitude, maximum_latitude),
                longitude=slice(minimum_longitude, maximum_longitude),
            )
            allds.append(ds)

        if allok is True and len(allds) > 0:
            break

    ds = xarray.concat(allds, dim="valid_time")
    ds = ds.drop("time").rename({"valid_time": "time"})
    ds.time.encoding["units"] = "days since 1900-01-01T00:00:00Z"

    # RH: =100*(EXP((17.625*TD)/(243.04+TD))/EXP((17.625*T)/(243.04+T)))
    # q :  equation 4.24, Pg 96 Practical Meteorolgy (Roland Stull)
    # https://github.com/Unidata/MetPy/issues/791
    # ds["rh2m"] = numpy.exp(17.625*ds["d2m"] / (243.04 + ds["d2m"])) /
    # numpy.exp(17.625*ds["t2m"]/(243.04+ds["t2m"]))
    # ds["rh2m"] = xarray.where(ds['rh2m']<1, ds['rh2m'], 1)
    ds["hu2m"] = metpy.calc.relative_humidity_from_dewpoint(ds["t2m"], ds["d2m"])
    # ds["hu2m"] = xarray.where(ds['hu2m']<1, ds['hu2m'], 1)
    ds = ds.drop("d2m")
    ds["ssr"] = ds["ssr"].differentiate("time", datetime_unit="s")
    ds["ssr"] = xarray.where(ds["ssr"] > 0, ds["ssr"], 0)
    ds["str"] = -ds["str"].differentiate("time", datetime_unit="s")
    ds["tp"] = ds["tp"].differentiate("time", datetime_unit="s")

    ds = ds.rename(
        {
            "u10": "u10m",
            "v10": "v10m",
            "msl": "pmer",
            "ssr": "flsolaire",
            "str": "fltherm",
            "tp": "prate",
        }
    )
    ds = ds.rename({"latitude": "lat", "longitude": "lon"})
    encoding = {var: {"zlib": True, "shuffle": True} for var in ds.data_vars}
    ds.to_netcdf(outpath, encoding=encoding)


def get_args():
    """ get args """
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
        default=pathlib.Path("/home/opsys/data/CROCO/GIBRTWO/ECMWF"),
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
            "/home/opsys/data/CROCO/GIBRTWO/CROCO_FILES/croco_gibrtwo_inno_energy_grd.nc"
        ),
        help="grid pathname",
    )

    args = parser.parse_args()

    return args


def run_main():
    """ main """
    args = get_args()

    today = pandas.Timestamp.utcnow().replace(
        hour=0, minute=0, second=0, microsecond=0, nanosecond=0, tzinfo=None
    )
    if args.begindate is not None:
        begindate = min(args.begindate, today)
    else:
        begindate = today - pandas.Timedelta(1, unit="day")
    if args.enddate is not None:
        enddate = args.enddate
    else:
        enddate = (
            begindate
            + pandas.Timedelta(5, unit="day")
            # + pandas.Timedelta(10, unit="day")
            # - pandas.Timedelta(6, unit="hour")
        )

    outname = f"{args.cfgname}_ecmwf_{args.begindate:%Y%m%d%H}.nc"
    outpath = args.workdir / outname
    download(begindate, enddate, args.workdir, args.grd, outpath)

    return 0


if __name__ == "__main__":
    sys.exit(run_main())
