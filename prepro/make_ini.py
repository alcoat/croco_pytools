#!/usr/bin/env python
"""
===========================================================================
Further Information:
  http://www.croco-ocean.org

This file is part of CROCOTOOLS

Create a CROCO initial file
In the current state the script can handle:
    - mercator (glorys)
    - soda
    - eccov4  (V4.3)

To add a new dataset you just have to go in Readers/ibc_readers.py and
create a dico with correspondance between croco_var(ssh,u,v,temp,salt)
and dataset variable names.

The periodicity of a dataset (if it is -180/180 or 0/360) is handle by the
script. Just download the dataset and it should be fine

The script works as follow:
    - reads croco_grd
    - reads input data and restricts the area of spatial coverage
    - creates croco_ini.nc
    - Loop on var with:
        * horizontal interpolation
        * vertical interpolation
    - Writes data in netcdf
===========================================================================
"""
__author__ = "Mathieu Le Corre"
__email__ = "mathieu.le.corre@shom.fr"
__date__ = "2022-09"
__license__ = "GPL3"

# --- Dependencies ---------------------------------------------------------

import logging
import pathlib
import sys
import argparse

import netCDF4
import numpy as np
import pandas
import pylab as plt
import yaml

sys.path.append("./Modules/")
sys.path.append("./Readers/")

import Cgrid_transformation_tools as grd_tools
import croco_class
import ibc_class
import interp_tools
import sigmagrid_tools as sig_tools

logger = logging.getLogger()
handler = logging.StreamHandler(sys.stderr)
handler.setFormatter(
    logging.Formatter(
        fmt=(
            "[%(asctime)s %(levelname)s] "
            "[%(filename)s:%(lineno)s - %(funcName)s() ] "
            "%(message)s"
        ),
        datefmt="%H:%M:%S",
    )
)
logger.addHandler(handler)
logger.setLevel(logging.INFO)


def get_args():
    """ get args """
    parser = argparse.ArgumentParser(
        description="make initial condition",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--config_file",
        type=pathlib.Path,
        default=pathlib.Path("make_ini_def.yml"),
        help="Config file",
    )
    args = parser.parse_args()

    return args


def run_make_ini():
    """ main """
    args = get_args()
    with open(args.config_file, encoding="utf8") as infile:
        ini_def = yaml.safe_load(infile)

    inputdata_type = ini_def["inputdata"]

    # edit ini_filename to add starting date
    begindate = pandas.Timestamp(ini_def["begindate"])
    ini_filename = ini_def["ini_filename"].format(date_b=begindate)

    # Load croco_grd
    croco_dir = pathlib.Path(ini_def["croco_dir"])
    grid_pathname = croco_dir / ini_def["croco_grd"]
    crocogrd = croco_class.CROCO_grd(grid_pathname, ini_def["sigma_params"])

    # --- Load input (restricted to croco_grd) ----------------------------

    inpdat = ibc_class.getdata(
        inputdata_type,
        ini_def["input_file"],
        crocogrd,
        ini_def["multi_files"],
        ini_def["tracers"],
    )

    ini_filepath = croco_dir / ini_filename
    print(" ")
    print(f" Making initial file: {ini_filepath}")
    print(" ")

    # --- Create the empty initial file -----------------------------------

    croco_class.CROCO.create_ini_nc(
        ini_filepath,
        crocogrd,
        tracers=ini_def["tracers"],
    )

    # --- Handle initial time ---------------------------------------------
    ini_date_num = begindate.to_pydatetime()
    ini_date_num = plt.date2num(ini_date_num)

    origindate = pandas.Timestamp(ini_def["origindate"])
    day_zero_num = origindate.to_pydatetime()
    day_zero_num = plt.date2num(day_zero_num)

    tstart = 0

    if ini_date_num != day_zero_num:
        tstart = ini_date_num - day_zero_num  # days

    # scrumt = tstart * 3600 * 24  # convert in second
    oceant = tstart * 3600 * 24
    scrumt = tstart
    # oceant = tstart
    tend = 0.0

    #  --- Compute and save variables on CROCO grid ---------------

    for l_vars in ["ssh", "tracers", "velocity"]:
        print(f"\nProcessing *{l_vars}*")

        nc = netCDF4.Dataset(ini_filepath, "a")
        if l_vars == "ssh":
            (zeta, NzGood) = interp_tools.interp_tracers(
                inpdat, l_vars, -1, crocogrd, ini_def["tndx"], ini_def["tndx"]
            )
            nc.variables["zeta"][0, :, :] = zeta * crocogrd.maskr
            nc.Input_data_type = inputdata_type
            nc.variables["ocean_time"][:] = oceant
            nc.variables["scrum_time"][:] = scrumt
            nc.variables["scrum_time"].units = (
                "days since {origindate:%Y-%m-%d %H:%M:%S}"
            )
            nc.variables["tstart"][:] = tstart
            nc.variables["tend"][:] = tend

            z_rho = crocogrd.scoord2z_r(zeta=zeta)
            z_w = crocogrd.scoord2z_w(zeta=zeta)

        elif l_vars == "tracers":
            for tra in ini_def["tracers"]:
                print(f"\nIn tracers processing {tra}")
                trac_3d = interp_tools.interp(
                    inpdat,
                    tra,
                    ini_def["Nzgoodmin"],
                    z_rho,
                    crocogrd,
                    ini_def["tndx"],
                    ini_def["tndx"],
                )
                nc.variables[tra][0, :, :, :] = trac_3d * crocogrd.mask3d()

        elif l_vars == "velocity":

            cosa = np.cos(crocogrd.angle)
            sina = np.sin(crocogrd.angle)

            [ubar_ogcm, vbar_ogcm, ubar, vbar] = interp_tools.compute_uvbar_ogcm(
                inpdat, cosa, sina, crocogrd, ini_def["tndx"], ini_def["tndx"]
            )

            [u, v] = interp_tools.interp_uv(
                inpdat,
                ini_def["Nzgoodmin"],
                z_rho,
                cosa,
                sina,
                crocogrd,
                ini_def["tndx"],
                ini_def["tndx"],
            )

            if ini_def["conserv"] == 1:
                (ubar_croco, h0) = sig_tools.vintegr(
                    u, grd_tools.rho2u(z_w), grd_tools.rho2u(z_rho), np.nan, np.nan
                ) / grd_tools.rho2u(crocogrd.h)
                (vbar_croco, h0) = sig_tools.vintegr(
                    v, grd_tools.rho2v(z_w), grd_tools.rho2v(z_rho), np.nan, np.nan
                ) / grd_tools.rho2v(crocogrd.h)

                u = u - ubar_croco
                u = u + np.tile(ubar, (z_rho.shape[0], 1, 1))
                v = v - vbar_croco
                v = v + np.tile(vbar, (z_rho.shape[0], 1, 1))

            nc.variables["u"][0, :, :, :] = u * crocogrd.umask3d()
            nc.variables["v"][0, :, :, :] = v * crocogrd.vmask3d()
            nc.variables["ubar"][0, :, :] = ubar * crocogrd.umask
            nc.variables["vbar"][0, :, :] = vbar * crocogrd.vmask

        nc.close()
    return 0

if __name__ == "__main__":
    sys.exit(run_make_ini())
