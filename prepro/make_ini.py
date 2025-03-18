__author__ = "Mathieu Le Corre"
__email__ = "mathieu.le.corre@shom.fr"
__date__ = "2022-09"
__license__ = "GPL3"
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

# --- Dependencies ---------------------------------------------------------

import glob
import sys
from datetime import datetime

import netCDF4 as netcdf
import numpy as np
import pylab as plt

sys.path.append("./Modules/")
sys.path.append("./Readers/")
import Cgrid_transformation_tools as grd_tools
import croco_class
import ibc_class
import interp_tools
import sigmagrid_tools as sig_tools

import logging

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

# --- USER CHANGES ---------------------------------------------------------
ini_def = {
    # Dates
    # starting date
    "Yini": "2025",  # Month and days need to be 2-digits format
    "Mini": "03",  # Month and days need to be 2-digits format
    "Dini": "17",  # Month and days need to be 2-digits format
    "Hini": "00",  # Month and days need to be 2-digits format
    # reference time (default = ini time)
    "Yorig": "2000",  # Month and days need to be 2-digits format
    "Morig": "01",  # Month and days need to be 2-digits format
    "Dorig": "01",  # Month and days need to be 2-digits format
    # Input data information and formating
    #"inputdata": "mercator_croco",  # Input data dictionnary as defined in the Readers/ibc_reader.py
    "inputdata": "mercator",  # Input data dictionnary as defined in the Readers/ibc_reader.py
    "input_dir": "../../MERCATOR_GLOB_2013/",
    "input_prefix": "mercator_",
    #"multi_files": False,  # If variables are in different netcdf
    #"input_file": "../../MERCATOR_GLOB_2013/mercator_Y2013M01.cdf",
    "multi_files": True,  # If variables are in different netcdf
    "input_file" : { 'ssh'  : '../../MERCATOR/glo12_rg_6h-i_*-2D-zos_fcst_R20250317.nc',\
                      'temp' : '../../MERCATOR/glo12_rg_6h-i_*-3D-thetao_fcst_R20250317.nc',\
                       'salt' : '../../MERCATOR/glo12_rg_6h-i_*-3D-so_fcst_R20250317.nc',\
                       'u'    : '../../MERCATOR/glo12_rg_6h-i_*-3D-uovo_fcst_R20250317.nc',\
                       'v'    : '../../MERCATOR/glo12_rg_6h-i_*-3D-uovo_fcst_R20250317.nc'\
                    },
    # time index to use in the file
    "tndx": 0,
    # default value to consider a z-level fine to be used
    "Nzgoodmin": 4,
    # tracers
    "tracers": ["temp", "salt"],
    # CROCO grid informations
    #"croco_dir": "../../CROCO_FILES/",
    #"croco_grd": "croco_grd.nc",
    "croco_dir": "../",
    "croco_grd": "croco_gibrtwo_inno_energy_grd.nc",
    "sigma_params": {
        "theta_s": 7,
        "theta_b": 2,
        "N": 32,
        "hc": 200,
    },  # Vertical streching, sig_surf/sig_bot/ nb level/critical depth
    # Ini file informations
    "ini_filename": "croco_ini.nc",  # output will be put in croco_dir by default
    # Conserv OGCM transport option
    "conserv": 1,  # Correct the horizontal transport i.e. remove
    # the integrated tranport and add the OGCM transport
}

# --- END USER CHANGES -----------------------------------------------------

if __name__ == "__main__":
    # edit ini_filename to add starting date
    ini_filename = ini_def["ini_filename"].replace(
        ".nc",
        "_%s_Y%sM%s.nc" % (ini_def["inputdata"], ini_def["Yini"], ini_def["Mini"]),
    )

    # Load croco_grd
    crocogrd = croco_class.CROCO_grd(
        "".join((ini_def["croco_dir"], ini_def["croco_grd"])), ini_def["sigma_params"]
    )

    # --- Load input (restricted to croco_grd) ----------------------------

    inpdat = ibc_class.getdata(
        ini_def["inputdata"],
        ini_def["input_file"],
        crocogrd,
        ini_def["multi_files"],
        ini_def["tracers"],
    )

    print(" ")
    print(" Making initial file: " + ini_filename)
    print(" ")

    # --- Create the initial file -----------------------------------------

    croco_class.CROCO.create_ini_nc(
        None,
        "".join((ini_def["croco_dir"] + ini_filename)),
        crocogrd,
        tracers=ini_def["tracers"],
    )

    # --- Handle initial time ---------------------------------------------
    ini_date_num = datetime(
        int(ini_def["Yini"]),
        int(ini_def["Mini"]),
        int(ini_def["Dini"]),
        int(ini_def["Hini"]),
    )
    ini_date_num = plt.date2num(ini_date_num)

    day_zero_num = datetime(
        int(ini_def["Yorig"]), int(ini_def["Morig"]), int(ini_def["Dorig"])
    )
    day_zero_num = plt.date2num(day_zero_num)

    tstart = 0

    if ini_date_num != day_zero_num:
        tstart = ini_date_num - day_zero_num  # days

    scrumt = tstart * 3600 * 24  # convert in second
    oceant = tstart * 3600 * 24
    tend = 0.0

    #  --- Compute and save variables on CROCO grid ---------------

    for l_vars in ["ssh", "tracers", "velocity"]:
        print(f"\nProcessing *{l_vars}*")
        nc = netcdf.Dataset(ini_def["croco_dir"] + ini_filename, "a")
        if l_vars == "ssh":
            (zeta, NzGood) = interp_tools.interp_tracers(
                inpdat, l_vars, -1, crocogrd, ini_def["tndx"], ini_def["tndx"]
            )
            nc.variables["zeta"][0, :, :] = zeta * crocogrd.maskr
            nc.Input_data_type = ini_def["inputdata"]
            nc.variables["ocean_time"][:] = oceant
            nc.variables["scrum_time"][:] = scrumt
            nc.variables["scrum_time"].units = "seconds since %s-%s-%s 00:00:00" % (
                ini_def["Yorig"],
                ini_def["Morig"],
                ini_def["Dorig"],
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
