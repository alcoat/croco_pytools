#--- Dependencies ---------------------------------------------------------

import sys
sys.path.append("./graphicUI_tools/")
from tools_make_grid import inputs,inputs_smth,EasyGrid,GetMask,GetTopo
from main_window import Outputs





#---F U N C T I ON S -------------------------------------------------------

def create_classic_grid(tra_lon, tra_lat, size_x, size_y, nx, ny, rot,
                  hmin, hmax, interp_rad, rfact, smooth_meth, topofile, shp_file, sgl_connect, output_file):
    """
    Processes data by creating grids, building masks, and saving the output in a NetCDF file.

    Parameters:
    tra_lon (float): Longitude of the grid center.
    tra_lat (float): Latitude of the grid center.
    size_x (int): Grid size in x direction [km].
    size_y (int): Grid size in y direction [km].
    nx (int): Number of points in the x direction.
    ny (int): Number of points in the y direction.
    rot (float): Grid rotation [degree].
    hmin (float): Minimum depth [m].
    hmax (float): Maximum depth [m].
    interp_rad (int): Interpolation radius in number of points.
    rfact (float): Maximum r-fact to reach.
    smooth_meth (str): Smoothing method.
    topofile (str): Path to the topography file.
    shp_file (str): Path to the shapefile for the coastline.
    sgl_connect (list): Single Connect mask water parameters.
    output_file (str): Path to the output NetCDF file.

    Returns:
    outputs: Processed outputs variable.
    """
    
    # --- Create inputs and outputs class -----------------------------
    inputs_ = inputs(tra_lon, tra_lat, size_x, size_y, nx, ny, rot)
    inputs_smth_ = inputs_smth(hmin, hmax, interp_rad, rfact, smooth_meth)
    outputs = Outputs()

    # --- Create lon/lat grid -----------------------------------------
    EasyGrid.easygrid(None, inputs_, outputs)

    # --- Build mask and topo -----------------------------------------
    GetTopo.topo(None, outputs, topofile, shp_file, smooth=inputs_smth_, sgl_connect=sgl_connect)
    
    # Return the outputs variable
    return outputs
