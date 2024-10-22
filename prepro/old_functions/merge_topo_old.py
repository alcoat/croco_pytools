#--- Dependencies ---------------------------------------------------------
import os,sys

import numpy as np
from scipy.interpolate import griddata
import xarray as xr
import matplotlib.pyplot as plt

import scipy.ndimage as ndimage

from scipy.interpolate import RectBivariateSpline
from scipy.ndimage import binary_dilation, binary_erosion

#--- Functions ------------------------------------------------------------

# ============================
# MAIN FUNCTION FOR MERGING
# ============================

def merge_topo(high_res, low_res, smoothing_method, merging_area, output_file, gisbase=None, EPSG_num=None):   
    
    """
    Merges two topography grids (the high-resolution grid must be entirely within the extent of the low-resolution grid): supports NetCDF and .grd file formats.
    
    1) The low-resolution grid is interpolated to match the resolution of the high-resolution grid over the same extent -> new_grid.
    2) The high-resolution grid is interpolated onto the corresponding area of new_grid to ensure alignment of coordinates.
    3) Any NaN (hole) values in the high-resolution grid are replaced with interpolated data from the low-resolution grid.
    4) Bathymetry smoothing is applied in a transition zone between the high-resolution grid boundary and a virtual boundary defined by a user-specified number of cells:
       "merging_area". The merging area is part of the interpolated low-resolution grid. Three smoothing options are available:
       - A: Simple linear weighting.
       - B: Weighted averaging between the linear interpolation and the existing low-resolution data in the merging area.
       - C: Use GRASS GIS add-on "r.mblend" (in this case, the "merging_area" parameter is used to define the "FAR_EDGE" option).
    5) The final output is saved as a NetCDF file.

    Parameters:
    ----------
    high_res (str): Path to the high resolution topography file.
    low_res (str): Path to the low resolution topography file.
    smoothing_method (str): Smoothing method.(A,B,C)
    output_file (str): Path to the output NetCDF file.
    gisbas (str): Path to grass if using r.mblend smoothing method ('C'). ex: "/usr/lib/grass78"
    EPSG_num (str): EPSG (referential) number where grid will be projected. ex:'EPSG:4326'

    """
    if smoothing_method in ('A', 'B'):
        
        merge_smooth(high_res, low_res, smoothing_method, merging_area, output_file)
        
    if smoothing_method=='C':

        """
        TRAVAUX EN COURS -> CF PLUS BAS FONCTION R.MBLEND VERSION PYTHON (CONDITION: AVOIR GRASS + GSCRIPT ETC..)
        """
        mblend_smoothing(high_res, low_res, far_edge, output_file, gisbase, EPSG_num)



# ============================
# SUB-FUNCTION FOR R.MBLEND
# ============================


def mblend_smoothing(high_res, low_res, far_edge, output_file, gisbase, EPSG_num):
    """
    Python command lines for creating a temporary grass dataset and running r.mblend addon over high/low resolution grids
    More info: "https://grass.osgeo.org/grass84/manuals/addons/r.mblend.html"

    Mblend: DEM merging method proposed by Leitão et al. (2016). It deals with cases where a study area is only partially covered by a high resolution DEM,
    with a coarser DEM available for the remainder (as in the case shown below). r.mblend merges the two DEMs, producing a smooth transition from the high resolution DEM to the     low resolution DEM.
        
    Parameters:
    ----------
    high_res (str): Path to the high resolution topography file.
    low_res (str): Path to the low resolution topography file.
    output_file (str): Path to the output NetCDF file.
    gisbas (str): Path to the GRASS GIS installation. ex: "/usr/lib/grass78"
    EPSG_num (str): EPSG (referential) number where grid will be projected. ex:'EPSG:4326'
    far_edge (int): Percentage of distance to high resolution raster used to determine far edge. Number between 0 and 100.
                    When the blending occurs along a single edge a number closer to 100 tends to produce more even results. With more blending edges (e.g. high resolution DEM                       sits on the middle of the low resolution DEM) a lower number may produce a more regular blend.
    
    """
    
    import subprocess
    import shutil
    import binascii
    import tempfile
    import grass.script as gscript
    import grass.script.setup as gsetup
    
    # Configuration paths for GRASS GIS
    gisdbase = tempfile.mkdtemp()  # Create a temporary directory for the GRASS GIS database
    
    # INPUTS
    input1 = high_res  # Path to the high resolution raster
    input2 = low_res  # Path to the low resolution raster
    output = "fusion_grilles"  # Name for the output raster
    
    def reprojeter_raster(input_file, output_file):
        """
        Reprojects a raster file to the specified coordinate reference system (EPSG:4326).
        
        Parameters:
        - input_file: Path to the input raster file.
        - output_file: Path to the output reprojected raster file.
        """
        try:
            subprocess.run([
                'gdalwarp', 
                '-t_srs', 'EPSG:4326', 
                input_file, 
                output_file
            ], check=True)
            print(f"Raster successfully reprojected: {output_file}")
        except subprocess.CalledProcessError as e:
            print(f"Error during reprojection: {e}")
    
    # Reproject the rasters to EPSG:4326
    output_file_2 = input2[:-4] + '_reproj_WGS84.tif'
    output_file_1 = input1[:-4] + '_reproj_WGS84.tif'
    reprojeter_raster(input1, output_file_1)
    reprojeter_raster(input2, output_file_2)
    
    # Prepare raster names for GRASS GIS
    raster1 = 'raster1_final'  # Name for the high resolution raster in GRASS
    raster2 = 'raster2_final'  # Name for the low resolution raster in GRASS
    
    # Create a new GRASS GIS location and mapset
    location = binascii.hexlify(os.urandom(16)).decode()  # Generate a random name for the location
    mapset = 'PERMANENT'  # Default mapset name
    location_path = os.path.join(gisdbase, location)  # Path to the new location
    startcmd = f'grass78 -c epsg:4326 -e {location_path}'  # Command to create the location
    
    print(f"Creating location with command: {startcmd}")
    subprocess.run(startcmd, shell=True, check=True)  # Execute the command to create the location
    
    # Initialize the GRASS GIS environment
    gsetup.init(gisbase, gisdbase, location, mapset)
    
    # Import rasters into GRASS GIS
    gscript.run_command('r.in.gdal', input=output_file_1, output=raster1)
    gscript.run_command('r.in.gdal', input=output_file_2, output=raster2)
    
    # Define the region for processing based on the input rasters
    gscript.run_command('g.region', raster=f"{raster1},{raster2}")
    
    # Run the r.mblend command to blend the rasters
    gscript.run_command('r.mblend',
                        high=raster1,
                        low=raster2,
                        output=output)
    
    # Export the resulting raster to a GeoTIFF file
    #gscript.run_command('r.out.gdal', input=output, output=os.path.join(rep, output + '.tif'), format='GTiff')

    #Export as NETCDF file
    gscript.run_command('r.out.gdal', input=output, output= output_file, format='NCDF')
    
    # Clean up temporary files
    shutil.rmtree(location_path)  # Remove the temporary location directory
    print(f"Temporary location {location_path} removed.")
    
    print("Process completed successfully.")



# ================================
# SUB-FUNCTION FOR CLASSIC MERGING
# ================================

def merge_smooth(high_res, low_res, smoothing_method, merging_area, output_file):
    """
    This function aims to interpolate the low resolution gird to the highresolution's one, to replace low resolution data by the high one available, and to smooth their
    frontier, using an artifical buffer (mergin are). Holes in the high resolution datasets are replaced by low resolutionsmoothed data.

    Inputs have to be in Lat /lon decimal degrees format, regular grids, with no rotation. The high resolution grid has to be inluded within the low resolution grid extent.
    Files extent allowed: '.nc', '.tiff', '.grd'
    
    Parameters:
    ----------
    high_res (str): Path to the high resolution topography file.
    low_res (str): Path to the low resolution topography file.
    smoothing_method (str): Smoothing method.(A,B) A-> Linear ponderation in buffer area, B-> Mean of A's ponderation and low resolution data within the buffer area.
    output_file (str): Path to the output NetCDF file.
    
    """
    
# --- Loading of the data------------------------------------------
    #Opening of the topography grids as xarray datasets
    ds1 = xr.open_dataset(high_res)
    ds2 = xr.open_dataset(low_res)

    #If tiff files: drop the useless "band" variable-> 2D data 
    if input1[-4:]=='tiff':
        ds1= ds1.sel(band=1, drop=True)
    if input2[-4:]=='tiff':
        ds2= ds2.sel(band=1, drop=True)
    
    # Identify the variable containing the bathymetry data
    # We assume that the bathymetry variable is the only one with numerical values
    
    for var_name in ds1.data_vars:
        var_data = ds1[var_name]
        # Check if the variable is 2D and contains numerical values
        if var_data.ndim == 2 and var_data.dtype.kind in {'f', 'i'}:
             z_ds1 = var_name
             break
    else:
        print("No bathymetry variable found in ds1")
    
    for var_name in ds2.data_vars:
        var_data = ds2[var_name]
        # Check if the variable is 2D and contains numerical values
        if var_data.ndim == 2 and var_data.dtype.kind in {'f', 'i'}:
             z_ds2 = var_name
             break
    else:
        print("No bathymetry variable found in ds2")

    # Identify the names of the lat/lon coordinates in the files

    for coord_name in ds1.coords:
        coord_data = ds1[coord_name]
        # We assume that lon and lat are typically coordinates with geographical values
        if coord_name.lower() in {'lon', 'longitude'}:
            lon_coord_ds1 = coord_name
        elif coord_name.lower() in {'lat', 'latitude'}:
            lat_coord_ds1 = coord_name
        elif coord_data.dims == ('x',) or coord_data.dims == ('y',):
            if 'x' in coord_data.dims:
                lon_coord_ds1 = coord_name
            if 'y' in coord_data.dims:
                lat_coord_ds1 = coord_name

    for coord_name in ds2.coords:
        coord_data = ds2[coord_name]
        # We assume that lon and lat are typically coordinates with geographical values
        if coord_name.lower() in {'lon', 'longitude'}:
            lon_coord_ds2 = coord_name
        elif coord_name.lower() in {'lat', 'latitude'}:
            lat_coord_ds2 = coord_name
        elif coord_data.dims == ('x',) or coord_data.dims == ('y',):
            if 'x' is in coord_data.dims:
                lon_coord_ds2 = coord_name
            if 'y' is in coord_data.dims:
                lat_coord_ds2 = coord_name

    # --- I N T E R P O L A T I O N - G R I D 2 ------------------------
    
    # Switching the resolution of grid 2 to match the high-resolution grid (grid 1)
    
    # Calculating the resolution in the latitude and longitude directions
    resolution_lon= np.diff(ds1[lon_coord_ds1]).mean()
    resolution_lat= np.diff(ds1[lat_coord_ds1]).mean()
    
    lon2= ds2[lon_coord_ds2].values
    lat2= ds2[lat_coord_ds2].values
    z2= ds2[z_ds2].values
    
    # Creating the grid 2 corresponding meshgrid
    lon_grid_2, lat_grid_2 = np.meshgrid(lon2, lat2)
    
    # Get 1D  longitudes/latitudes from meshgrid 2
    lon_flat_2 = lon_grid_2.ravel()
    lat_flat_2 = lat_grid_2.ravel()
    z_flat_2 = z2.ravel()
    
    # Create a new grid for ds2 with the specified resolution
    new_lon_2 = np.arange(lon2.min(), lon2.max(), resolution_lon)
    new_lat_2 = np.arange(lat2.min(), lat2.max(), resolution_lat)
    new_lon_grid_2, new_lat_grid_2 = np.meshgrid(new_lon_2, new_lat_2)
    
    # Interpolate ds2 data onto the new grid
    z2_interp = griddata((lon_flat_2, lat_flat_2), z_flat_2, (new_lon_grid_2, new_lat_grid_2), method='nearest')
    

    # --- I N T E R P O L A T I O N - G R I D 2 ------------------------
    
    # Extract lon, lat, and z data from the high-resolution grid (ds1)
    lon1 = ds1[lon_coord_ds1].values
    lat1 = ds1[lat_coord_ds1].values
    z1 = ds1[z_ds1].values
    
    # Retrieve the original grid from ds1 (high-resolution grid)
    lon_grid_1, lat_grid_1 = np.meshgrid(lon1, lat1)
    
    # Flatten the 2D grids and their corresponding values
    lon_flat_1 = lon_grid_1.ravel()
    lat_flat_1 = lat_grid_1.ravel()
    z_flat_1 = z1.ravel()
    
    # Identify the overlapping region between the two grids
    lon_min_1, lon_max_1 = lon1.min(), lon1.max()  # Longitude extent of grid 1
    lat_min_1, lat_max_1 = lat1.min(), lat1.max()  # Latitude extent of grid 1
    
    # Create a mask on the resampled grid of ds2 that only covers the area of ds1
    mask_overlap = (new_lon_grid_2 >= lon_min_1) & (new_lon_grid_2 <= lon_max_1) & \
                   (new_lat_grid_2 >= lat_min_1) & (new_lat_grid_2 <= lat_max_1)
    
    # Interpolate the high-resolution data (ds1) onto the overlapping area of ds2
    z1_interp_on_z2_overlap = griddata(
        (lon_flat_1, lat_flat_1),  # Source points (grid 1)
        z_flat_1,                  # Source values (grid 1)
        (new_lon_grid_2[mask_overlap], new_lat_grid_2[mask_overlap]),  # Target points (overlapping region in grid 2)
        method='nearest'           # Interpolation method (nearest neighbor)
    )

    # Save the original interpolated ds2 grid for later use
    z2_save = z2_interp.copy()
    
    # Replace the values in z2_interp (low-res grid) with those from grid 1 in the overlapping region
    z2_interp[mask_overlap] = z1_interp_on_z2_overlap
    
    
# --- SMOOTHING ---------------------------------------------------
    
    # --- FILLING THE GAPS IN HIGH RESOLUTION GRID + SMOOTHING WITH SPLINE METHOD ----------------------
    
    #Identify which holes in grid 1 could be filled with low resolution grid available bathymetry
    mask_nan = np.isnan(z2_interp) & ~np.isnan(z2_save)
    
    #Create a virtual grid for smoothing operations: data is filled, smoothed
    z2_combined= z2_interp.copy()
    z2_combined[mask_nan]= z2_save[mask_nan]
    
    # Creating an interpolator tool for the studied 2D grid
    spline = RectBivariateSpline(new_lat_grid_2[:, 0], new_lon_grid_2[0, :], z2_combined)
    
    # Interpolation with the spline method to smooth the combined datasets
    z2_smoothed = spline(new_lat_grid_2[:, 0], new_lon_grid_2[0, :])
    
    # Repllce ONLY original NaNs by the processed data blended and smoothed z2_smoothed
    z2_interp[mask_nan]= z2_smoothed[mask_nan]

    # --- CREATING A MERGING AREA (BUFFER) FOR SMOOTH TRANSITION BETWEEN GRIDS --------------------------

    # Define the width of the buffer zone in number of cells
    buffer_width = merging_area
    
    # Edges of the high-resolution grid within the final grid
    lon_min = new_lon_grid_2[mask_overlap].min()
    lon_max = new_lon_grid_2[mask_overlap].max()
    lat_min = new_lat_grid_2[mask_overlap].min()
    lat_max = new_lat_grid_2[mask_overlap].max()
    
    # Calculate the outer limits of the buffer zone
    lon_buffer_min = lon_min - buffer_width * resolution_lon
    lon_buffer_max = lon_max + buffer_width * resolution_lon
    lat_buffer_min = lat_min - buffer_width * resolution_lat
    lat_buffer_max = lat_max + buffer_width * resolution_lat
    
    # Mask for the buffer zone + high-resolution grid
    mask_buff_extremum = (new_lon_grid_2 >= lon_buffer_min) & (new_lon_grid_2 <= lon_buffer_max) & \
                         (new_lat_grid_2 >= lat_buffer_min) & (new_lat_grid_2 <= lat_buffer_max)
    
    # Buffer zone excluding the high-resolution grid
    mask_buffer = mask_buff_extremum & ~mask_overlap
    

    # --- CATCHING THE USEFULL CELLS (NEAREST FROM BUFFER) FOR LINEAR PONDERATION TO COME --------------------------

    # Binary dilatation to add external neighbors meshes
    dilated_mask_buffer_ext = binary_dilation(mask_buff_extremum)
    
    # External borders of the buffer from low resolution grid only
    border_buffer_ext = dilated_mask_buffer_ext & ~mask_buff_extremum

    #Binary erosion to get mask_overlap without the frontier
    eroded_mask = binary_erosion(mask_overlap)
    
    # Get inner borders only
    border_mask_inner = mask_overlap & ~eroded_mask
    
    #LON/LAT of these borders masks
    lon_inner = new_lon_grid_2[border_mask_inner]
    lat_inner = new_lat_grid_2[border_mask_inner]
    lon_outer = new_lon_grid_2[border_buffer_ext]
    lat_outer = new_lat_grid_2[border_buffer_ext]

    #And the z2_interp corresponding values
    values_inner = z2_interp[border_mask_inner]
    values_outer = z2_interp[border_buffer_ext]

    # Buffer's cells coordinates
    buffer_lon = new_lon_grid_2[mask_buffer]
    buffer_lat = new_lat_grid_2[mask_buffer]

    # For each buffer's cell, compute the linear ponderation between inner and outer buffers's borders
    # /!\ HERE THE DISTANCES ARE CALCULATED WITH DECIMAL DEGREES
    # IT'S NOT AS ACCURATE THAN THE HAVERSINE DISTANCE BUT IT'S QUICKER
    #def haversine(lon1, lat1, lon2, lat2):
    #    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
    #    dlon = lon2 - lon1
    #    dlat = lat2 - lat1
    #    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    #    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    #    return R * c
    
    for i, (lon, lat) in enumerate(zip(buffer_lon, buffer_lat)):
        # Calculate distances to all cells in border_mask_inner and border_mask_ext

        # dist_inner = haversine(lon, lat, lon_inner, lat_inner)
        # dist_outer = haversine(lon, lat, lon_outer, lat_outer)
        
        #'Decimal degrees' distance
        dist_inner = np.sqrt((lon_inner - lon)**2 + (lat_inner - lat)**2)
        dist_outer = np.sqrt((lon_outer - lon)**2 + (lat_outer - lat)**2)
                
        # Find the indices of the nearest borders
                nearest_inner_idx = np.argmin(dist_inner)
                nearest_outer_idx = np.argmin(dist_outer)
                
        # Extract the nearest values
                value_inner = values_inner[nearest_inner_idx]
                value_outer = values_outer[nearest_outer_idx]
                
        # Calculate the total distance for weighting
                total_distance = dist_inner[nearest_inner_idx] + dist_outer[nearest_outer_idx]
                
        # Calculate weights for interpolation
                weight_inner = dist_outer[nearest_outer_idx] / total_distance
                weight_outer = dist_inner[nearest_inner_idx] / total_distance
        
        if smoothing_method=='A':        
            # Compute the interpolated value with linear ponderation only
            interpolated_values[i] = weight_inner * value_inner + weight_outer * value_outer
        elif smoothing_method=='N':
            # Compute the interpolated value with linear ponderation + meaning with low res'
            interpolated_values[i] = weight_inner * value_inner + weight_outer * value_outer
            

    #replace low resolution buffer by interpolated values
    z2_interp[mask_buffer] = interpolated_values
    

    # --- SAVE DATASET AS NETCDF --------------------------------------
    
   # Obtenez les coordonnées uniques
    new_lon_unique = np.unique(new_lon_grid_2)
    new_lat_unique = np.unique(new_lat_grid_2)
    
    # Créer un DataArray pour les données interpolées
    ds_interpolated = xr.DataArray(
        z2_interp, 
        coords=[('lat', new_lat_unique), ('lon', new_lon_unique)],  # Ajouter les coordonnées
        dims=['lat', 'lon']  # Spécifier les dimensions
    )
    
    # Créer un Dataset pour organiser les variables
    ds_to_save = xr.Dataset({
        'z': ds_interpolated  # Ajouter la variable interpolée
    })
    
    # Sauvegarder dans un fichier NetCDF
    ds_to_save.to_netcdf(output_file)

    
    print("DEM have been merged and the generated grid has been saved in: " + output_file)