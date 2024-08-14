#--- Dependencies ---------------------------------------------------------

import sys
#sys.path.append("./graphicUI_tools/")
from tools_make_grid import inputs,inputs_smth,EasyGrid,GetMask,GetTopo, topo_prt, Inputs_smth_c2c, Inputs_c2c
import numpy as np
#from main_window import Outputs

#--- Classes --------------------------------------------------------------
class Outputs:
    """
    Outputs object for handling geographical or grid data.
    """
    def __init__(self,
                 lon_rho=None,
                 lat_rho=None,
                 lon_u=None,
                 lat_u=None,
                 lon_v=None,
                 lat_v=None,
                 h=None,
                 hraw=None,
                 pm=None,
                 pn=None,
                 angle=None,
                 f=None,
                 mask_rho=None):
        # Initialize with provided values or empty arrays if not provided
        self.lon_rho = np.array(lon_rho, dtype=np.float64) if lon_rho is not None else np.array([])
        self.lat_rho = np.array(lat_rho, dtype=np.float64) if lat_rho is not None else np.array([])
        self.lon_u = np.array(lon_u, dtype=np.float64) if lon_u is not None else np.array([])
        self.lat_u = np.array(lat_u, dtype=np.float64) if lat_u is not None else np.array([])
        self.lon_v = np.array(lon_v, dtype=np.float64) if lon_v is not None else np.array([])
        self.lat_v = np.array(lat_v, dtype=np.float64) if lat_v is not None else np.array([])
        self.h = np.array(h, dtype=np.float64) if h is not None else np.array([])
        self.hraw = np.array(hraw, dtype=np.float64) if hraw is not None else np.array([])
        self.pm = np.array(pm, dtype=np.float64) if pm is not None else np.array([])
        self.pn = np.array(pn, dtype=np.float64) if pn is not None else np.array([])
        self.angle = np.array(angle, dtype=np.float64) if angle is not None else np.array([])
        self.f = np.array(f, dtype=np.float64) if f is not None else np.array([])
        self.mask_rho = np.array(mask_rho, dtype=np.float64) if mask_rho is not None else np.array([])





#---F U N C T I ON S -------------------------------------------------------

def create_classic_grid(tra_lon, tra_lat, size_x, size_y, nx, ny, rot,
                  hmin, hmax, interp_rad, rfact, smooth_meth, topofile, shp_file, sgl_connect):
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
    return outputs, inputs_

###########################################

def create_child_grid_offline(tra_lon, tra_lat, size_x, size_y, nx, ny, rot,
                  hmin, hmax, interp_rad, rfact, smooth_meth, topofile, shp_file, sgl_connect, sgl_index1, sgl_index2, open_boundary, parent_grid):
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
    open_boundary (dict): Contains checkboxes booleans for NESW directions, merging_area and 'WESN' which contains direction in a text format

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
    sgl_connect_grouped= [sgl_connect, sgl_index1, sgl_index2]
    prt= topo_prt(parent_grid)
    GetTopo.topo(None, outputs, topofile, shp_file, smooth=inputs_smth_, sgl_connect=sgl_connect_grouped, prt_grd=prt)

    # --- Match parent/child topo -------------------------------------
    obc= open_boundary
    GetTopo.match_topo(None,prt,outputs,[obc['dirs'],obc['MERGING_AREA']])
    
    # Return the outputs variable
    return outputs, inputs_, prt

###########################################

def create_child_grid_AGRIF(coef, imin, imax, jmin, jmax, hmin, hmax, interp_rad, rfact, smooth_meth, topofile, shp_file, sgl_connect, sgl_index1, sgl_index2, open_boundary,parent_grid):

    """
    Processes data by creating grids, building masks, and saving the output in a NetCDF file.

    Parameters:
    coef (int): Refinement coefficient.
    imin (int): Parent grid left limit's indice matching with first child grid column.
    imax (int): Parent grid right limit's indice matching with last child grid column.
    jmin (int): Parent grid lower limit's indice matching with first child grid row.
    jmax (int): Parent grid upper limit's indice matching with last child grid row.
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
    open_boundary (dict): Contains checkboxes booleans for NESW directions, merging_area and 'WESN' which contains direction in a text format

    Returns:
    outputs: Processed outputs variable.
    """
    
    # --- Create inputs and outputs class -----------------------------
    inputs= Inputs_c2c(coef,imin,imax,jmin,jmax)
    inputs_smth= Inputs_smth_c2c(interp_rad,rfact,smooth_meth)
    outputs = Outputs()

    # --- Create lon/lat grid -----------------------------------------
    prt= topo_prt(parent_grid)
    EasyGrid.AGRIFgrid(None,prt,inputs,outputs)

    # --- Build mask and topo -----------------------------------------
    sgl_connect_grouped= [sgl_connect, sgl_index1, sgl_index2]
    GetTopo.topo(None,outputs, topofile, shp_file, smooth= inputs_smth, hmin=np.nanmin(prt.h),
                  hmax=np.nanmax(prt.h), sgl_connect=sgl_connect_grouped, prt_grd=prt, coef= inputs.coef)

    # --- Match parent/child topo -------------------------------------
    obc= open_boundary
    GetTopo.match_topo(None,prt,outputs,[obc['dirs'],obc['MERGING_AREA']])
    
    # Return the outputs variable
    return outputs, inputs, prt



