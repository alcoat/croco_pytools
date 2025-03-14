__author__ = "Mathieu Le Corre"
__email__ = "mathieu.le.corre@shom.fr"
__date__ = "2022-09"
__license__ = "GPL3"
"""
===========================================================================
Further Information:  
  http://www.croco-ocean.org
  
This file is part of CROCOTOOLS

Create a CROCO grid file
In the current state the script can handle:
    - etopo (5,2,1)
    - srtm30
    - gebco

To add a new dataset you just have to go into Readers/topo_readers.py and
create a dictionnary with the names for lon,lat,topo in the dataset.
At this time it only handles 1D lon/lat 

The script makes a grid with mercator projection centred at the equator
before rotating the sphere to put the middle of the grid in tra_lon/tra_lat
position.

Then, it reads topo's dataset and apply the desired smoothing
The mask is generated using a shapefile (.shp. default is gshhs dataset)

Smoothing uses fortran routines for better performance
"""

# --- Dependencies ---------------------------------------------------------

import os
import sys

import numpy as np

os.environ["ETS_TOOLKIT"] = "wx"
sys.path.append("./Modules/")
sys.path.append("./Readers/")
sys.path.append("./Modules/graphicUI_tools/")
import tools_make_grid
from main_window import *

# --- USER CHANGES ---------------------------------------------------------

grid_def = {
    # Grid center [degree]
    "tra_lon": 15,  # Longitude of the grid center
    "tra_lat": -32,  # Latitude of the grid center
    # Grid size [km]
    "size_x": 1556,
    "size_y": 1334,
    # Grid number of points
    # Note: grid resolution is grid size / number of points
    "nx": 39,
    "ny": 40,
    # Grid rotation [degree]
    "rot": 0,
    # Smoothing parameters
    # (see online documentation for more details)
    "hmin": 50,  # Minimum depth [m]
    "hmax": 6000,  # Maximum depth [m]
    "interp_rad": 2,  # Interpolation radius in number of points (usually between 2 and 8)
    "rfact": 0.2,  # Maximum r-fact to reach (the lower it is, the smoother it will be)
    "smooth_meth": "lsmooth",  # Smoothing method ('smooth', 'lsmooth', 'lsmooth_legacy', 'lsmooth2', 'lsmooth1', 'cond_rx0_topo')
    # Topo/Bathy file
    "topofile": "/home/opsys/DATA/BATHYMETRIE/ETOPO/ETOPO2v2c_f4.nc",
    # Coastline file (for the mask)
    "shp_file": "/home/opsys/DATA/SHORELINE/GSHHG/GSHHS_shp-2.3.7/i/GSHHS_i_L1.shp",
    # Single Connect [Mask water not connected to the main water body]
    "sgl_connect": [
        False,
        20,
        20,
    ],  # True or False, point indices inside the main water body
    # Output grid file
    "output_file": "../croco_grd.nc",
}

# --- END USER CHANGES -----------------------------------------------------


# CST
r_earth = 6371315  # Mean earth radius in meter


if __name__ == "__main__":

    Question = input(
        "Do you want to use interactive grid maker ? \
                      \n (e.g., for grid rotation or parameter adjustments) : y,[n] "
    )

    if Question.lower() == ("y") or Question.lower() == ("yes"):

        # --- Building grid with graphicUI -------------------------------------

        print("In interactive mode")
        MainWindow().configure_traits()
    elif (
        Question.lower() == ("n")
        or Question.lower() == ("no")
        or Question.lower() == ("")
    ):

        # --- Building grid without graphicUI ----------------------------------

        print("In normal mode")

        from croco_class import CROCO
        from main_window import Outputs
        from tools_make_grid import (EasyGrid, GetMask, GetTopo, inputs,
                                     inputs_smth)

        # --- Create inputs and outputs class -----------------------------

        inputs = inputs(
            grid_def["tra_lon"],
            grid_def["tra_lat"],
            grid_def["size_x"],
            grid_def["size_y"],
            grid_def["nx"],
            grid_def["ny"],
            grid_def["rot"],
        )
        inputs_smth = inputs_smth(
            grid_def["hmin"],
            grid_def["hmax"],
            grid_def["interp_rad"],
            grid_def["rfact"],
            grid_def["smooth_meth"],
        )
        outputs = Outputs()

        # --- Create lon/lat grid -----------------------------------------

        EasyGrid.easygrid(None, inputs, outputs)

        # --- Build mask and topo -----------------------------------------
        GetTopo.topo(
            None,
            outputs,
            grid_def["topofile"],
            grid_def["shp_file"],
            smooth=inputs_smth,
            sgl_connect=grid_def["sgl_connect"],
        )

        # --- Save netcdf -------------------------------------------------

        print("Writing Topography")
        CROCO.create_grid_nc(None, grid_def["output_file"], inputs, outputs)
