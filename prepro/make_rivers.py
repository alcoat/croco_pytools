__author__ = "Mathieu Le Corre"
__email__ = "mathieu.le.corre@shom.fr"
__date__ = "2023-01"
__license__ = "GPL3"

"""
===========================================================================
Further Information:  
  http://www.croco-ocean.org
  
This file is part of CROCOTOOLS

Create a CROCO river file
In the current state the script can handle:
    - 2D NETCDF reanalysis (https://cds.climate.copernicus.eu/cdsapp!/dataset/cems-glofas-historical )
    - Multiple/single station NETCDF (http://www.marineinsitu.eu/dashboard/)
    - .txt or .dat file for each rivers


To use this script you need a .txt or .dat file following the format:
    river.txt     lon   lat  # river data to put at lon/lat
    river.nc      lon   lat  # read river netcdf on put it at lon/lat, works for 1 river
    river.nc       X     X   # read rivers in netcdf and look for lon/lat in it
    river_2d.nc    x    Qmin # Read a 2D river map and will automatically place rivers
                             # with flow>=Qmin 

If river data is in format .txt or .dat it needs to follow the standart:
%d/%m/%Y %H:%M:%S     lon     lat   
or
%Y/%m/%d %H:%M:%S     lon     lat

In the current state the script can not handle temperature time series for rivers and will put a constant value for each rivers. Feel free to change this value in for_croco_in.txt output file.

The script is only working for interranual data but you can create a cycle to build your own climatology.

The script works ad follow:
    - reads input text file 
    - reads river data and restrictis to temporal coverage
    - interpolates data on the same time axis
    - reads croco_grd
    - positions rivers on croco_grd (with eventually manual edition)
    - creates croco_runoff and fills it 
    - creates text file to fill croco.in
    - summaries river positions

Thanks to Pierrick Penven and Guillaume Charria for their help
"""

# --- Dependencies ---------------------------------------------------------

import glob
import sys

import cftime
import matplotlib.pyplot as py
import netCDF4 as netcdf
import numpy as np
import pylab as plt
from dateutil.relativedelta import relativedelta

sys.path.append("./Modules/")
import croco_class as Croco
import tools_make_river as riv_tools

# --- USER CHANGES ---------------------------------------------------------

rivers_def = {
    # Dates
    "Yorig": 2013,  # year origin of time : days since Yorig-01-01
    "Ystart": 2013,
    "Mstart": 1,  # Starting month
    "Yend": 2013,  # Ending month
    "Mend": 3,  # Ending month
    # Input data informations
    "input_dir": "./",  # where .txt file can be found
    "input_file": "river_list.txt",  # txt file name
    "Crange": 3,  # only used in 2D reanalysis, number of iterations to broaden the coast mask
    # CROCO path and filename informations
    # "croco_dir": "../../CROCO_FILES/",
    "croco_dir": "../",
    "croco_grd": "croco_grd.nc",
    # Rivers file informations
    "river_filename": "croco_runoff.nc",  # output will be put in croco_dir by default
    "river_file_format": "MONTHLY",  # How outputs are split (MONTHLY,YEARLY,FULL)
    "rivers_cyl": 0.0,  # if cycle is needed
    "rivers_output_frequency": "DAILY",  # output frequency for rivers (HOURLY,DAILY,MONTHLY depending on your data)
}

# --- END USER CHANGES -----------------------------------------------------

# --- START MAIN SCRIPT ----------------------------------------------------

add_ts = False  # Script do not handle Temperature and salinity yet

# Handle time
time_units = "days since %i-01-01" % rivers_def["Yorig"]

# Put start and end date to the right format
rstr = plt.datetime.datetime(rivers_def["Ystart"], rivers_def["Mstart"], 1, 0, 0, 0)
rend = plt.datetime.datetime(
    rivers_def["Yend"], rivers_def["Mend"], 1, 12, 0, 0
) + relativedelta(
    months=1, days=-1
)  # Last day of the ending month


# --- Load croco_grd --------------------------------------------------

crocogrd = Croco.CROCO_grd("".join((rivers_def["croco_dir"], rivers_def["croco_grd"])))
geolim = [crocogrd.lonmin(), crocogrd.lonmax(), crocogrd.latmin(), crocogrd.latmax()]
# --- Read input_data file --------------------------------------------
data = np.genfromtxt(rivers_def["input_file"], dtype=str, comments="#")
if data.ndim == 1:
    data = data.reshape((1, 3))

if data.ndim == 0 or data.shape[1] != 3:
    print("Input data files should follow the following format:")
    print("    river_file lon_mouth lat_mouth\nor")
    print(" river_file.nc    X         X\nor")
    print("    2d_file.nc    X        Qmin")


list_river_files = np.array([x[0] for x in data])
lon_river = np.array([])
lat_river = np.array([])
for x in data:
    try:
        lon_river = np.append(lon_river, float(x[1]))
    except:
        try:
            lon_river = np.append(lon_river, int(x[1]))
        except:
            lon_river = np.append(lon_river, np.nan)
    #######
    try:
        lat_river = np.append(lat_river, float(x[2]))
    except:
        try:
            lat_river = np.append(lat_river, int(x[2]))
        except:
            lat_river = np.append(lat_river, np.nan)

# --- Read rivers data ------------------------------------------------
river_obs = riv_tools.read_river(
    list_river_files,
    lon_river,
    lat_river,
    rstr,
    rend,
    time_units,
    geolim=geolim,
    Crange=rivers_def["Crange"],
)

# --- Get river default indexes on crocogrd ---------------------------
river_obs = riv_tools.get_river_index(river_obs, crocogrd)

# --- Put rivers data on same time axis --------------------------------
river_ext = riv_tools.fill_period(
    river_obs, rstr, rend, time_units, rivers_def["rivers_output_frequency"]
)

# --- Manual river positionning ---------------------------------------
Question = input(
    "Do you want to manually edit river positions using interactive GUI ?: y,[n] "
)

if Question.lower() == ("y") or Question.lower() == ("yes"):
    river_ext = riv_tools.correc_qtriver(
        river_ext, grd=crocogrd
    )  # load_previous_corrections=False)
river_ext = riv_tools.locateji_croco(river_ext, grd=crocogrd, graph=False)


# --- Create and write netcdf runoff file -----------------------------
startloc = plt.datetime.datetime(rivers_def["Ystart"], rivers_def["Mstart"], 1)
if rivers_def["river_file_format"].upper() == "MONTHLY":
    endloc = startloc + relativedelta(months=1)
elif rivers_def["river_file_format"].upper() == "YEARLY":
    if plt.datetime.datetime.strptime(str(startloc), "%Y-%m-%d %H:%M:%S").year == int(
        rivers_def["Yend"]
    ):
        endloc = rend.replace(tzinfo=None)
    else:
        endloc = plt.datetime.datetime(
            int(rivers_def["start_date"][:4]) + 1, 1, 1, 0
        )  # BUG no start_date
elif rivers_def["river_file_format"].upper() == "FULL":
    endloc = rend.replace(tzinfo=None)
else:
    print(
        '\n Output file format "%s" is not setup. Please change it to MONTHLY, YEARLY or FULL'
    )
    sys.exit()
    # --- Loop on monthly/yearly/full data ----------------------------
while plt.date2num(endloc) <= plt.date2num(rend):
    loc_time = river_ext[list(river_ext.keys())[0]]["time"]
    # find index for the time range
    ind = np.where(
        (loc_time >= cftime.date2num(startloc, time_units))
        & (loc_time <= cftime.date2num(endloc, time_units))
    )

    if len(ind[0]) == 0:
        #            print('\nData is missing for range %s to %s' % (startloc ,endloc))
        sys.exit()

    [dtmin, dtmax] = np.min(ind), np.max(ind)
    # create monthly file
    tmp_date = plt.datetime.datetime.strptime(str(startloc), "%Y-%m-%d %H:%M:%S")
    # file name depending on format chosen
    if rivers_def["river_file_format"].upper() == "MONTHLY":
        river_outname = rivers_def["croco_dir"] + rivers_def["river_filename"].replace(
            ".nc", "_Y%sM%02i.nc" % (tmp_date.year, tmp_date.month)
        )
    elif rivers_def["river_file_format"].upper() == "YEARLY":
        river_outname = rivers_def["croco_dir"] + rivers_def["river_filename"].replace(
            ".nc", "_Y%s.nc" % (tmp_date.year)
        )
    elif rivers_def["river_file_format"].upper() == "FULL":
        river_outname = rivers_def["croco_dir"] + rivers_def["river_filename"]
    Croco.CROCO.create_river_nc(
        None, river_outname, crocogrd, len(river_ext.keys()), add_ts
    )

    # --- Check if data availabla for the surrounded months -----------
    dtmin = dtmin - 1  # create overlap before ( done inside fill_period)
    dtmax = dtmax + 1  # create overlap after  ( done inside fill_period)

    # --- Write data in netcdf ----------------------------------------
    nc = netcdf.Dataset(river_outname, "a")
    nc.variables["qbar_time"][:] = river_ext[list(river_ext.keys())[0]]["time"][
        dtmin : dtmax + 1
    ]
    nc.variables["qbar_time"].units = time_units
    if add_ts:
        nc.variables["temp_time"][:] = river_ext[list(river_ext.keys())[0]]["time"][
            dtmin : dtmax + 1
        ]
        nc.variables["temp_time"].units = time_units
        nc.variables["salt_time"][:] = river_ext[list(river_ext.keys())[0]]["time"][
            dtmin : dtmax + 1
        ]
        nc.variables["salt_time"].units = time_units

    for ir, r in enumerate(river_ext):
        nc.variables["runoff_name"][ir] = r
        nc["Qbar"][ir] = river_ext[r]["flow"][dtmin : dtmax + 1]
        if add_ts:
            nc["temp_src"][ir] = river_ext[r]["temp"][dtmin : dtmax + 1]
            nc["salt_src"][ir] = river_ext[r]["salt"][dtmin : dtmax + 1]

    nc.close()

    # --- Preparing time for next loop ------------------------------------
    startloc = endloc
    if rivers_def["river_file_format"].upper() == "MONTHLY":
        endloc = startloc + relativedelta(months=1)
    elif rivers_def["river_file_format"].upper() == "YEARLY":
        yearloc = plt.datetime.datetime.strptime(str(startloc), "%Y-%m-%d %H:%M:%S")
        if startloc == rend:
            endloc = startloc + relativedelta(days=1)
        elif yearloc.year == int(rivers_def["Yend"]):
            endloc = rend.replace(tzinfo=None)
        else:
            endloc = plt.datetime.datetime(int(yearloc.year) + 1, 1, 1, 0)
    elif rivers_def["river_file_format"].upper() == "FULL":
        endloc = startloc + relativedelta(days=1)
    print("end")

# --- Write runoff file for croco.in ----------------------------------
riv_tools.write_croco_in(
    river_ext, rivers_def["croco_dir"], rivers_def["river_filename"]
)


# --- Plot a summary of the river positions ---------------------------
Question = input("Do you want a summary of the position of the rivers ?: y,[n] ")

if Question.lower() == ("y") or Question.lower() == ("yes"):
    py.figure()
    py.title("Summary of the rivers in the domain")
    ax = py.pcolormesh(
        crocogrd.lon, crocogrd.lat, np.ma.masked_where(crocogrd.maskr == 0, crocogrd.h)
    )
    py.colorbar(ax, label="Bathymetry [m]")
    for ll in river_ext.keys():
        if river_ext[ll]["dsrc"] == 0:
            py.plot(
                crocogrd.lonu[river_ext[ll]["jj"] - 1, river_ext[ll]["ii"] - 1],
                crocogrd.latu[river_ext[ll]["jj"] - 1, river_ext[ll]["ii"] - 1],
                "+r",
            )
            py.quiver(
                crocogrd.lonu[river_ext[ll]["jj"] - 1, river_ext[ll]["ii"] - 1],
                crocogrd.latu[river_ext[ll]["jj"] - 1, river_ext[ll]["ii"] - 1],
                [river_ext[ll]["qbardir"]],
                [0],
            )
        else:
            py.plot(
                crocogrd.lonv[river_ext[ll]["jj"] - 1, river_ext[ll]["ii"] - 1],
                crocogrd.latv[river_ext[ll]["jj"] - 1, river_ext[ll]["ii"] - 1],
                "+r",
            )
            py.quiver(
                crocogrd.lonv[river_ext[ll]["jj"] - 1, river_ext[ll]["ii"] - 1],
                crocogrd.latv[river_ext[ll]["jj"] - 1, river_ext[ll]["ii"] - 1],
                [0],
                [river_ext[ll]["qbardir"]],
            )
    py.show()
