#--- Dependencies ---------------------------------------------------------
import os,sys
import numpy as np
from scipy.interpolate import griddata
import xarray as xr
import matplotlib.pyplot as plt
from scipy.ndimage import binary_erosion, binary_dilation
from scipy.spatial import cKDTree
import pyinterp
import pyinterp.backends.xarray
# Module that handles the filling of undefined values.
import pyinterp.fill
import gc
import rasterio
import rioxarray
#Configure import paths
sys.path.extend(["./../Readers/"])
import topo_reader
from pyproj import Transformer, CRS
#Dependencies for interpolatig .txt file
import pandas as pd
import geopandas as gpd
import matplotlib.tri as tri
from shapely import vectorized
from itertools import product
import regionmask


# /!\ The dependencies used for r.mblend method (2nd function) are imported inside the function

#--- Functions ---------------------------------------------------------

# ================================================
# FUNCTION FOR CREATING BATHY GRID WITH .TXT FILE
# ================================================

def create_depth_grid_from_points(data_file_path, shapefile_path, lat_bounds, lon_bounds, epsg_int, grid_resolution, output_file=None):
    """
    Create a 2D grid of interpolated depth values from depth data points and geographical boundaries specified by a shapefile.
    
    This function reads depth data from a specified text file, filters it based on provided latitude and longitude bounds,
    and performs cubic triangulation interpolation to generate a grid. The function also masks the grid to exclude areas 
    not covered by specified geographical boundaries.
    
    Parameters:
    - data_file_path (str): Path to the text file containing depth data in three columns (latitude, longitude, depth).
    - shapefile_path (str): Path to the shapefile containing the geographical boundaries to mask the grid.
    - lat_bounds (list): A list containing the minimum and maximum latitude [min_lat, max_lat].
    - lon_bounds (list): A list containing the minimum and maximum longitude [min_lon, max_lon].
    - epsg_int (int): EPSG code for the coordinate reference system (CRS) of the depth data.
    - grid_resolution (float): Desired resolution for the interpolated grid, expressed in the units of the EPSG.
    - output_file (str, optional): Path to save the output NetCDF file containing the depth grid.
    
    Returns:
    - xarray.Dataset: A dataset containing the interpolated depth grid with dimensions of latitude and longitude.
                      If `output_file` is specified, the dataset is also saved to a NetCDF file.
    
    Raises:
    - ValueError: If the shapefile's CRS is not defined or if no depth data falls within the specified bounds.
    """

# ──────────────────────────────────────────────────────────
# ▶ SECTION 1: Reading bathy text file + filtering◀
# ──────────────────────────────────────────────────────────
    
    print("🗃️ Reading depth data from the text file...")
    
    # Read the data from the text file (no header)
    df = pd.read_csv(data_file_path, sep='\s+', header=None, names=['lat', 'lon', 'depth'])
    
    # Remove duplicates and sort the data
    df = df.drop_duplicates().sort_values(by=['lat', 'lon'])
    
    print("🔍 Filtering data within specified bounds...")
    # Filter data within the specified latitude and longitude bounds
    filtered_data = df[(df['lat'] >= lat_bounds[0]) & (df['lat'] <= lat_bounds[1]) &
                       (df['lon'] >= lon_bounds[0]) & (df['lon'] <= lon_bounds[1])]
    
    # Extract longitude, latitude, and depth from the filtered data
    lon = filtered_data['lon'].values
    lat = filtered_data['lat'].values
    depth = filtered_data['depth'].values

# ──────────────────────────────────────────────────────────
# ▶ SECTION 2: Shapefile data local extraction◀
# ──────────────────────────────────────────────────────────
    
    print("📜 Reading the shapefile...")
    # Load shapefile metadata to get its CRS (loading only one row to retrieve CRS)
    try:
        shapefile_meta = gpd.read_file(shapefile_path, rows=1)
        shapefile_crs = shapefile_meta.crs
        
        # Check if the shapefile has a defined CRS
        if shapefile_crs is None:
            raise ValueError("❌ No CRS found in the shapefile")
            
        # Reproject lon_bounds and lat_bounds to the shapefile's CRS (epsg_int -> shapefile CRS)
        transformer = Transformer.from_crs(f"EPSG:{epsg_int}", shapefile_crs, always_xy=True)
        min_lon, min_lat = transformer.transform(lon_bounds[0], lat_bounds[0])
        max_lon, max_lat = transformer.transform(lon_bounds[1], lat_bounds[1])
        
    except Exception as e:
        #In case the shoreline file is .json/kml ... and not a shapefile
        print(f"⚠️ The shoreline file seems to not be a shapefile: {e}")
        print("Assuming the shoreline file is in the same reference system as the bathymetry...")


    # Use the reprojected bounding box to load only the necessary data
    gdf = gpd.read_file(shapefile_path, bbox=(min_lon, min_lat, max_lon, max_lat))

# ───────────────────────────────────────────────────────────
# ▶ SECTION 3: Shapefile CRS conversion to match bathy's one◀
# ───────────────────────────────────────────────────────────
        
    # Check and ensure that the CRS of the shapefile is defined
    if gdf.crs is None:
        print("❌ Error: Shapefile CRS is not defined.")
        return None
    
    # If the CRS of the shapefile and the text data are different, we need to reproject the shapefile
    if gdf.crs.to_epsg() != epsg_int:
        print("🔄 Reprojecting shapefile to match text data EPSG...")
        gdf = gdf.to_crs(epsg=epsg_int)  # Reproject the shapefile to match text data EPSG
    
    # Filter the shapefile to keep only geometries within the specified bounds
    #gdf = gdf.cx[lon_bounds[0]:lon_bounds[1], lat_bounds[0]:lat_bounds[1]]

# ──────────────────────────────────────────────────────────────────
# ▶ SECTION 4: Output grid creation + interpolation of bathy points◀
# ──────────────────────────────────────────────────────────────────
    
    print("📏 Defining grid for interpolation...")
    # Define the grid for interpolation using the specified resolution in degrees
    num_points_lon = int((lon_bounds[1] - lon_bounds[0]) / grid_resolution)
    num_points_lat = int((lat_bounds[1] - lat_bounds[0]) / grid_resolution)
    grid_lon, grid_lat = np.mgrid[lon_bounds[0]:lon_bounds[1]:num_points_lon*1j, 
                                   lat_bounds[0]:lat_bounds[1]:num_points_lat*1j]
    
    print("⚙️ Performing interpolation using cubic triangulation...")
    # Interpolation using Matplotlib's Triangulation
    triang = tri.Triangulation(lon, lat)
    grid_depth = tri.CubicTriInterpolator(triang, depth)(grid_lon, grid_lat)

# ──────────────────────────────────────────────────────────────
# ▶ SECTION 5: Erasing terrestrial zones from interpolated data◀
# ──────────────────────────────────────────────────────────────
    
    print("🏝️ Creating mask for grid points covering terrestrial areas...")
    # Create a mask for the grid points based on the polygons in the shapefile
    mask = np.zeros_like(grid_depth, dtype=bool)

    ### OLD WAY TO DETRMINE COORDS INSIDE POLYGONS ##################################################################################
    # Extract x and y coordinates from the grid
    #grid_lon_flat = grid_lon.ravel()
    #grid_lat_flat = grid_lat.ravel()
    
    # Apply vectorized contains for each polygon in the GeoDataFrame
    #for poly in gdf.geometry:
    #    mask |= vectorized.contains(poly, grid_lon_flat, grid_lat_flat).reshape(grid_depth.shape)
    #################################################################################################################################

    ny, nx = grid_lon.shape
    mask = np.zeros_like(grid_lon, dtype=bool)  # Initialize the mask as a boolean array with the same shape as grid_lon
    n2max = 50000  # Maximum number of points per chunk
    
    # Determine if the data should be processed in chunks based on the size of gdf
    if gdf.shape[0] > 3000:
        # Calculate the number of chunks based on the grid size and the maximum number of points per chunk
        nchunk = int(np.max([np.sqrt((ny)*(nx)/n2max), 1]))  # At least 1 chunk
    else:
        nchunk = 1  # No chunking if gdf is small
    
    # If chunking is applied, print the chunk format
    if nchunk > 1:
        print(f"Chunk format (y,x):({nchunk},{nchunk})")

    
    # Loop over the chunks in both x and y directions
    for i, j in product(list(range(nchunk)), list(range(nchunk))):
        if nchunk > 1:
            print(f"🧩 Processing chunk ({i+1}/{nchunk}, {j+1}/{nchunk})")

    
        # Overlap between chunks to avoid edge effects
        dx1 = 2; dx2 = 2; dy1 = 2; dy2 = 2
        if i == 0: dx1 = 0  # No overlap on the left-most chunk
        if i == nchunk - 1: dx2 = 0  # No overlap on the right-most chunk
        if j == 0: dy1 = 0  # No overlap on the top-most chunk
        if j == nchunk - 1: dy2 = 0  # No overlap on the bottom-most chunk
    
        # Define the indices for the current chunk in the x and y dimensions
        nx1i = int(i * (nx) / nchunk - 2 * dx1)  # Start index in x
        nx2i = int((i + 1) * (nx) / nchunk + 2 * dx2)  # End index in x
        ny1i = int(j * (ny) / nchunk - 2 * dy1)  # Start index in y
        ny2i = int((j + 1) * (ny) / nchunk + 2 * dy2)  # End index in y
    
        # Find the geographic boundaries of the current chunk in lon/lat
        llcrnrlon = np.nanmin(grid_lon[ny1i:ny2i, nx1i:nx2i])  # Lower-left corner longitude
        urcrnrlon = np.nanmax(grid_lon[ny1i:ny2i, nx1i:nx2i])  # Upper-right corner longitude
        llcrnrlat = np.nanmin(grid_lat[ny1i:ny2i, nx1i:nx2i])  # Lower-left corner latitude
        urcrnrlat = np.nanmax(grid_lat[ny1i:ny2i, nx1i:nx2i])  # Upper-right corner latitude
    
        # Clip the geodataframe to the chunk bounding box
        gs = gdf.clip((llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat))
    
        try:  # Try to apply region masking if polygons exist in the chunk
            # Create a mask for the chunk using regionmask and the grid coordinates (lon/lat)
            rmask = regionmask.mask_geopandas(
                      gs.geometry, grid_lon[ny1i:ny2i, nx1i:nx2i], grid_lat[ny1i:ny2i, nx1i:nx2i])
    
            # Update the main mask: set areas where regionmask is NaN to 1 (i.e., inside polygons)
            mask[ny1i + dy1:ny2i - dy2, nx1i + dx1:nx2i - dx2]\
                [np.isnan(rmask[dy1:ny2i - ny1i - dy2, dx1:nx2i - nx1i - dx2])] = 1
            
            print(f"✅ Chunk ({i+1}/{nchunk}, {j+1}/{nchunk}) processed successfully!")
            
        except:  # If there are no polygons in the chunk, set the mask to 1 (marking the whole chunk)
            print(f"ℹ️ No polygons found in chunk ({i+1}/{nchunk}, {j+1}/{nchunk})")
            mask[ny1i + dy1:ny2i - dy2, nx1i + dx1:nx2i - dx2] = 1
            continue  # Continue to the next chunk
    
    # Apply the mask to set depths inside polygons to NaN
    grid_depth_masked = np.where(mask==False, np.nan, grid_depth)

# ──────────────────────────────────────────────────────────
# ▶ SECTION 6: Dataset creation for export◀
# ──────────────────────────────────────────────────────────
    print("Shape of grid_depth_masked:", grid_depth_masked.shape)
    print("Length of grid_lon:", len(grid_lon[:, 0]))
    print("Length of grid_lat:", len(grid_lat[0, :]))
    #ds = xr.Dataset({"depth": (["lon", "lat"], grid_depth_masked)}, coords={"lon": (["lon"], grid_lon[:, 0]), "lat": (["lat"], grid_lat[0, :])})
    ds = xr.Dataset(
        {
            "depth": (["lat", "lon"], np.transpose(grid_depth_masked))
        },
        coords={
            "lon": (["lon"], grid_lon[:, 0]),
            "lat": (["lat"], grid_lat[0, :])
        }
    )
    
    # Add spatial reference attribute (optional)
    ds.rio.write_crs(CRS.from_epsg(epsg_int), inplace=True)

    print("✅ Depth grid created successfully from textfile!")

# ──────────────────────────────────────────────────────────
# ▶ SECTION 7: Optionnal saving◀
# ──────────────────────────────────────────────────────────
    
    # Save the grid_depth to a NetCDF file if output_file is specified
    if output_file:
        print("💾 Saving the depth grid to NetCDF file...")
        ds.to_netcdf(output_file)
        print("✅ Depth grid saved successfully!")
    
    return ds


# ================================================
# FUNCTION CHUNKING INTERPOLATION
# ================================================

def interpolate_large_grid(latitudes, longitudes, data_array, chunk_size=1000, overlap=25, method='nearest'):
    """
    Interpolates a large grid using chunking to avoid memory overload.

    Parameters:
    -----------
    latitudes : np.ndarray
        2D array of latitude values for the grid to interpolate over.
    longitudes : np.ndarray
        2D array of longitude values for the grid to interpolate over.
    data_array : xr.DataArray
        DataArray containing the data to interpolate, with dimensions 
        corresponding to latitude and longitude.
    chunk_size : int, optional
        Size of the chunks to process at a time (default: 1000).
    overlap : int, optional
        Size of overlap between chunks to ensure smooth transitions (default: 25).
    method : str, optional
        Interpolation method, default is 'nearest'. Other options include 'linear'.

    Returns:
    --------
    interpolated_result : np.ndarray
        2D array of interpolated values matching the shape of the input grid.
    """

    # Initialize the grid interpolator using the provided DataArray
    grid_interpolator = pyinterp.backends.xarray.RegularGridInterpolator(data_array, geodetic=False)
    
    grid_shape = longitudes.shape  # Get the shape of the grid to interpolate
    
    interpolated_result = np.full(grid_shape, np.nan)  # Initialize the result grid with NaN

    # Check if the grid is large enough to require chunking
    if grid_shape[0] * grid_shape[1] > 3000000:
        # Loop through the grid in chunks to interpolate
        for row_start in range(0, grid_shape[0], chunk_size - overlap):
            for col_start in range(0, grid_shape[1], chunk_size - overlap):
                # Calculate the end of the current chunk, not exceeding the grid boundaries
                row_end = min(row_start + chunk_size, grid_shape[0])
                col_end = min(col_start + chunk_size, grid_shape[1])

                # Adjust the chunk to include overlap, without exceeding grid boundaries
                row_start_overlap = max(row_start - overlap, 0)
                col_start_overlap = max(col_start - overlap, 0)

                # Extract longitude and latitude values for this chunk                    
                lon_chunk = longitudes[row_start_overlap:row_end, col_start_overlap:col_end].ravel() #.flatten()
                lat_chunk = latitudes[row_start_overlap:row_end, col_start_overlap:col_end].ravel() #.flatten()

                # Prepare the coordinate pairs for the interpolator
                coords_chunk = {
                    'longitude': lon_chunk,
                    'latitude': lat_chunk
                }

                # Perform interpolation for this chunk
                interpolated_chunk = grid_interpolator(coords_chunk, method=method).reshape(row_end - row_start_overlap, col_end - col_start_overlap)

                # Update the result grid, replacing NaN values only
                interpolated_result[row_start_overlap:row_end, col_start_overlap:col_end] = np.where(
                    np.isnan(interpolated_result[row_start_overlap:row_end, col_start_overlap:col_end]), 
                    interpolated_chunk, 
                    interpolated_result[row_start_overlap:row_end, col_start_overlap:col_end]
                )
    else:
        # For smaller grids, perform interpolation directly
        flat_coords = {
            'longitude': longitudes.ravel(), #.flatten(),
            'latitude': latitudes.ravel() #.flatten()
        }
        # Interpolate directly across the entire grid
        interpolated_result = grid_interpolator(flat_coords, method=method).reshape(longitudes.shape)

    return interpolated_result


# ================================================
# FUNCTION FOR FINDING COORDINATES'NAMES (NETCDF)
# ================================================

def find_coords(ds):
    """
    Find the coordinate names for longitude and latitude in the given dataset (ds).
    
    It checks if common variations of longitude and latitude ('lon', 'lat', 'longitude', 
    'latitude', 'x', 'y', etc.) are contained in the coordinate names, even if they have 
    additional characters (e.g., 'longitude025').
    
    Parameters:
    ds (xarray.Dataset): The dataset in which to find longitude and latitude coordinates.
    
    Returns:
    lon_coord, lat_coord (str, str): The names of the longitude and latitude coordinates.
    """
    possible_lon_names = {'lon', 'longitude', 'x'}
    possible_lat_names = {'lat', 'latitude', 'y'}
    
    lon_coord = None
    lat_coord = None

    # Loop through the coordinates and check for known substrings of longitude and latitude names
    for coord_name in ds.coords:
        coord_lower = coord_name.lower()  # Convert the coordinate name to lowercase for comparison
        
        # Check if any known longitude or latitude substrings are contained in the coordinate name
        if any(lon in coord_lower for lon in possible_lon_names):
            lon_coord = coord_name
        if any(lat in coord_lower for lat in possible_lat_names):
            lat_coord = coord_name
        
        # If both lon and lat are found, no need to continue
        if lon_coord and lat_coord:
            break

    # Final fallback: if no match found, raise an error
    if lon_coord is None or lat_coord is None:
        raise ValueError(f"❌ Could not find valid longitude/latitude coordinates in the dataset.")

    return lon_coord, lat_coord

# ================================================
# FUNCTION FOR MERGING - LINEAR PONDERATION METHOD
# ================================================

def merge_smooth(high_res, low_res, buffer_width, output_file, target_epsg='EPSG:4326', coarsen_factor=None, downscale_bounds=None, params=None):

    """
    Process high and low resolution grids by interpolating missing values and creating a blended output.

    This function performs the following steps:
    1. Sets up the coordinate grid for interpolation based on the high resolution grid (`high_res`).
    2. Generates a buffer around the high resolution grid, and identifies borders to mame a quickest ponderation with direct neighbors.
    3. Interpolates values in the buffer area where the data is missing (NaN) using the low resolution grid (`input2`).
    4. Fills the NaN in high_res with data interpolated from low_res + smoothing around NaN-filled areas
    5. Fills the last Nan with pyinterp loess method
    6. Saves the processed grid to an output file in NetCDF format.

    Parameters:
    ----------
    - high_res : (file for xarray.DataArray) The high resolution grid containing the primary data to be processed.
    - low_res : (file for xarray.DataArray) The low resolution grid used for interpolating values in the buffer area of `input1`.
    - buffer_width : (int) The width of the buffer to be applied around the high resolution grid. This determines the extent of the dilation. 
    - output_file : (str) The path to the output NetCDF file where the processed grid will be saved.
    - target_epsg : (str, optional) The EPSG code of the target coordinate reference system (CRS) to which the datasets will be reprojeted.
                    If both datasets have undefined CRS, they will be assumed to be in the specified `target_epsg` CRS. Defaults to 'EPSG:4326'.
    - coarsen_factor: (int) Resolution reduction factor for high-resolution grid. OPTION
    - downscale_bounds: (List) [lon_min, lon_max, lat_min, lat_max] for downscaling the low-resolution grid. OPTION
    - params (dict, optional): 
        A dictionary containing parameters for `create_depth_grid_from_points` when high_res is a .txt file.
        The dictionary must include the following keys:
        - shapefile_path (str): Path to the shapefile containing geographical boundaries.
        - lat_bounds (list): Latitude bounds [min_lat, max_lat].
        - lon_bounds (list): Longitude bounds [min_lon, max_lon].
        - epsg_int (int): EPSG code for the coordinate reference system.
        - grid_resolution (float): Desired grid resolution.
        - output_file (str, optional): Path to save the output NetCDF file. If not provided, will default to None.
    

    Returns:
    -------
    None
        The function saves the processed grid to the specified output file and does not return any value.

    Notes:
    ------
    - The buffer area is created around the high resolution grid with a 'mergin_area' extent which can be lowered if high_res borders too close from the border of low_res 
    ones
    - The function performs linear ponderation within this buffer using the high/low resolution grids.
    - /!\ Ensure that the input grids are compatible in terms of coordinate systems and dimensions.
    - The output file will be in NetCDF format with the processed data.
    
    """


# ──────────────────────────────────────────────────────────
# ▶ SECTION 1: Data recovery and conversion◀
# ──────────────────────────────────────────────────────────

    # Open the datasets
    if high_res.endswith('.txt'):
        if params is None:
            raise ValueError("params must be provided when high_res is a .txt file.")
        print("🔄 Loading high resolution data from text file...")
        
        ds1 = create_depth_grid_from_points(
            high_res,
            params['shapefile_path'],
            params['lat_bounds'],
            params['lon_bounds'],
            params['epsg_int'],
            params['grid_resolution'],
            params.get('output_file', None)
        )

    else:
        ds1 = xr.open_dataset(high_res)
        
    ds2 = xr.open_dataset(low_res)
    
   # If high_res is a TIFF file, select the first band
    if high_res.endswith('.tiff') or high_res.endswith('.tif'):
        ds1 = ds1.sel(band=1, drop=True)

    # If low_res is a TIFF file, select the first band
    if low_res.endswith('.tiff') or low_res.endswith('.tif'):
        ds2 = ds2.sel(band=1, drop=True)

    # Check for CRS
    if 'spatial_ref' in ds1.coords:
        proj_high_res = ds1.rio.crs
    else:
        proj_high_res = None
    
    if 'spatial_ref' in ds2.coords:
        proj_low_res = ds2.rio.crs
    else:
        proj_low_res = None

    # Compare projections if both are defined
    # Reproject using pyproj if CRS is defined and different
    if proj_high_res and proj_low_res:
        if proj_high_res != proj_low_res:
            print(f"⚠️ The datasets are in different projections. Checking for necessary reprojections to {target_epsg}...")
        
            # Set the target CRS (EPSG code provided or WGS84 by default)
            target_crs = CRS.from_user_input(target_epsg)
        
            # Reproject high resolution dataset only if not already in the target CRS
            if proj_high_res != target_crs:
                print(f"🔄 Reprojecting high resolution dataset to {target_epsg}...")
                transformer_high_to_target = Transformer.from_crs(proj_high_res, target_crs, always_xy=True)
                high_res_coords = transformer_high_to_target.transform(ds1['lon'].values, ds1['lat'].values)
                ds1['lon'], ds1['lat'] = high_res_coords
            else:
                print(f"✅ High resolution dataset is already in {target_epsg}, no reprojection needed.")
        
            # Reproject low resolution dataset only if not already in the target CRS
            if proj_low_res != target_crs:
                print(f"🔄 Reprojecting low resolution dataset to {target_epsg}...")
                transformer_low_to_target = Transformer.from_crs(proj_low_res, target_crs, always_xy=True)
                low_res_coords = transformer_low_to_target.transform(ds2['lon'].values, ds2['lat'].values)
                ds2['lon'], ds2['lat'] = low_res_coords
            else:
                print(f"✅ Low resolution dataset is already in {target_epsg}, no reprojection needed.")
        
            print(f"✅ Both datasets are now aligned to {target_epsg}.")
        else:
            print("✅ Both datasets are already in the same projection.")
    
    # Handle cases where one or both datasets lack a CRS
    else:
        if proj_high_res is None and proj_low_res is None:
            print("⚠️ Warning: Both high and low resolution datasets lack a defined CRS.")
            print(f"➡️ Assuming both datasets are in CRS {target_epsg}.")
        elif proj_high_res is None:
            print("⚠️ Warning: High resolution dataset lacks a defined CRS.")
            print(f"➡️ Assuming it matches the low resolution dataset CRS and is in {target_epsg}.")
        elif proj_low_res is None:
            print("⚠️ Warning: Low resolution dataset lacks a defined CRS.")
            print(f"➡️ Assuming it matches the high resolution dataset CRS and is in {target_epsg}.")

    # Identify the variable containing bathymetry data in ds1
    # Assuming the bathymetry variable is the only one with numerical values
    for var_name in ds1.data_vars:
        var_data = ds1[var_name]
        # Check if the variable is 2D and contains numerical values
        if var_data.ndim == 2 and var_data.dtype.kind in {'f', 'i'}:
            z_ds1 = var_name
            break
    else:
        raise ValueError("❌ No valid bathymetry variable found in the high resolution dataset (ds1)")
    
    # Identify the variable containing bathymetry data in ds2
    # Assuming the bathymetry variable is the only one with numerical values
    for var_name in ds2.data_vars:
        var_data = ds2[var_name]
        # Check if the variable is 2D and contains numerical values
        if var_data.ndim == 2 and var_data.dtype.kind in {'f', 'i'}:
            z_ds2 = var_name
            break
    else:
        raise ValueError("❌ No valid bathymetry variable found in the low resolution dataset (ds2)")
    
    # Identify the coordinate names for latitude/longitude in ds1
    lon_coord_ds1, lat_coord_ds1 = find_coords(ds1)
    lon_coord_ds2, lat_coord_ds2 = find_coords(ds2)

    # Apply resolution reduction with coarsen if specified
    if coarsen_factor:
        ds1 = ds1.coarsen({lat_coord_ds1: coarsen_factor, lon_coord_ds1: coarsen_factor}, boundary='trim').mean()
        print("✅ Coarsening completed. The high-resolution grid has been successfully downsampled.")

    # Apply downsample if limits are specified in the list
    # Validation des bornes de downscale_bounds
    if downscale_bounds:
        if len(downscale_bounds) != 4:
            raise ValueError("❌ downscale_bounds must be a list of four elements [lon_min, lon_max, lat_min, lat_max]")
        
        lon_min, lon_max, lat_min, lat_max = downscale_bounds
        if lon_min >= lon_max or lat_min >= lat_max:
            raise ValueError("❌ Longitude/Latitude bounds are not valid: ensure lon_min < lon_max and lat_min < lat_max")

        #ds2 = ds2.sel({lat_coord_ds2: slice(lat_min, lat_max), lon_coord_ds2: slice(lon_min, lon_max)})

        # Reverse latitude bounds if necessary (for descending latitudes)
        if ds2[lat_coord_ds2].values[0] > ds2[lat_coord_ds2].values[-1]:
            ds2 = ds2.sel({lat_coord_ds2: slice(lat_max, lat_min)})
        else:
            ds2 = ds2.sel({lat_coord_ds2: slice(lat_min, lat_max)})
        
        # Reverse longitude bounds if necessary (for descending longitudes)
        if ds2[lon_coord_ds2].values[0] > ds2[lon_coord_ds2].values[-1]:
            ds2 = ds2.sel({lon_coord_ds2: slice(lon_max, lon_min)})
        else:
            ds2 = ds2.sel({lon_coord_ds2: slice(lon_min, lon_max)})

        
        print("✅ Downscaling completed. The low-resolution grid has been successfully cropped to the specified bounds.")

# ────────────────────────────────────────────────────────────────────
# ▶ SECTION 2: INTERPOLATION OF GRID 2 TO THE RESOLUTION OF GRID 1  ◀
# ────────────────────────────────────────────────────────────────────

    # Switching the resolution of grid 2 to match the high-resolution grid (grid 1)
    
    # Calculating the resolution in the latitude and longitude directions
    resolution_lon= np.abs(np.diff(ds1[lon_coord_ds1]).mean())
    resolution_lat= np.abs(np.diff(ds1[lat_coord_ds1]).mean())
    
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

    # Create a DataArray for z2 with unique 1D latitude and longitude coordinates
    z2_dataarray = xr.DataArray(
        z2,
        coords={
            'latitude': lat2,  # 1D array of latitudes
            'longitude': lon2    # 1D array of longitudes
        },
        dims=['latitude', 'longitude']  # Dimensions should match the DataArray shape
    )

    '''
    # Initialize the grid interpolator for z2
    interpolator_2 = pyinterp.backends.xarray.RegularGridInterpolator(z2_dataarray, geodetic=False)
    

    # Get the shape of the interpolation grid
    shape = new_lon_grid_2.shape  
    
    # Check if the grid size exceeds 3,000,000 elements
    if shape[0] * shape[1] > 3000000:
        
        # Parameters for chunking
        chunk_size = 1000  # Maximum size of a chunk
        overlap = 25       # Size of the overlap
        
        # Initialize the interpolation grid, filled with NaN
        z2_interp = np.full(shape, np.nan)  
        
        # Loop through the grid in chunks
        for i in range(0, shape[0], chunk_size - overlap):
            for j in range(0, shape[1], chunk_size - overlap):
                # Determine the size of the chunk without exceeding the limits
                end_i = min(i + chunk_size, shape[0])
                end_j = min(j + chunk_size, shape[1])
                
                # Ensure the overlap does not exceed the limits
                start_i = max(i - overlap, 0)
                start_j = max(j - overlap, 0)
    
                # Create coordinates for this chunk
                lon_coords = new_lon_grid_2[start_i:end_i, start_j:end_j].flatten()  # Flattened longitude coordinates for the chunk
                lat_coords = new_lat_grid_2[start_i:end_i, start_j:end_j].flatten()  # Flattened latitude coordinates for the chunk
                
                # Prepare the coordinate pairs for the interpolator
                coords_subset = {
                    'longitude': lon_coords,
                    'latitude': lat_coords
                }
    
                # Use the interpolator on this chunk
                z2_interp_subset = interpolator_2(coords_subset, method='nearest').reshape(end_i - start_i, end_j - start_j)
    
                # Update the interpolation grid with results
                z2_interp[start_i:end_i, start_j:end_j] = np.where(
                    np.isnan(z2_interp[start_i:end_i, start_j:end_j]), 
                    z2_interp_subset,  # Use the interpolated values
                    z2_interp[start_i:end_i, start_j:end_j]  # Keep existing values
                )
        
        print("✅ Low-resolution grid interpolation completed with chunks' method.")
    
    else:
        # Prepare new coordinates for interpolation if the grid is smaller
        new_coords_2 = {
            'longitude': new_lon_grid_2.flatten(),  # Flattened array of new longitudes
            'latitude': new_lat_grid_2.flatten()     # Flattened array of new latitudes
        }
        
        # Perform interpolation using the nearest method
        z2_interp = interpolator_2(new_coords_2, method='nearest').reshape(new_lon_grid_2.shape)
        
        print("✅ Low-resolution grid interpolation completed (No chunks).")

    # Prepare new coordinates for interpolation
        #new_coords_2 = {
         #   'longitude': new_lon_grid_2.flatten(),  # Flattened array of new longitudes
          #  'latitude': new_lat_grid_2.flatten()     # Flattened array of new latitudes
       # }
    # Perform interpolation using the nearest method
    #z2_interp = interpolator_2(new_coords_2, method='nearest').reshape(new_lon_grid_2.shape)

    #print("✅ Low-resolution grid interpolation completed.")
    '''
    z2_interp = interpolate_large_grid(new_lat_grid_2, new_lon_grid_2, z2_dataarray)

    print("✅ Low-resolution grid interpolation completed.")


# ───────────────────────────────────────────────────────────────
# ▶ SECTION 3: HIGH RESOLUTION DATA INTEGRATION IN NEW GRID◀
# ──────────────────────────────────────────────────────────

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

    # Create a DataArray for z1
    z1_dataarray = xr.DataArray(
        z1,  # The values to be interpolated
        coords={
            'latitude': lat1,  # 1D array of latitudes corresponding to z1
            'longitude': lon1   # 1D array of longitudes corresponding to z1
        },
        dims=['latitude', 'longitude']  # Define the dimensions of the DataArray
    )

    '''
    # Create the interpolator for z1
    interpolator = pyinterp.backends.xarray.RegularGridInterpolator(z1_dataarray, geodetic=False)
    
    # Prepare new coordinates for interpolation within the overlap region
    new_coords_1 = {
        'longitude': new_lon_grid_2[mask_overlap].flatten(),  # Flattened longitudes from the overlap region
        'latitude': new_lat_grid_2[mask_overlap].flatten()    # Flattened latitudes from the overlap region
    }
    
    # Apply the interpolator to the new coordinates
    z1_interp_on_z2_overlap = interpolator(new_coords_1, method='nearest')
    
    # Reshape the result to match the shape of the grid 2
    z1_interp_on_z2_overlap = z1_interp_on_z2_overlap.reshape(new_lon_grid_2[mask_overlap].shape)
    '''
    #Interpolates grid 1 over all the grid 2 extent (=> NAN on grid 1's sides)
    z1_interp= interpolate_large_grid(new_lat_grid_2, new_lon_grid_2, z1_dataarray, chunk_size=500, overlap=25, method='nearest')
    #z1_interp= z1_interp[mask_overlap]
    
    #--->Clean memory:
    gc.collect()

    # Save the original interpolated ds2 grid for later use
    #z2_save = z2_interp.copy()
    
    # Replace the values in z2_interp (low-res grid) with those from grid 1 in the overlapping region
    #z2_interp[mask_overlap] = z1_interp_on_z2_overlap
    #z2_interp[mask_overlap] = z1_interp

    print("✅ High-resolution grid interpolation completed.")

    #---> Clean memory:
    #del z1_interp_on_z2_overlap, z2_dataarray, z1_dataarray, lon_grid_1, lat_grid_1, ds1, ds2
    #del z1_interp, z2_dataarray, z1_dataarray, lon_grid_1, lat_grid_1, ds1, ds2
    del z2_dataarray, z1_dataarray, lon_grid_1, lat_grid_1, ds1, ds2
    gc.collect()

# ──────────────────────────────────────────────────────────
# ▶ SECTION 4: REPLACE NAN IN THE HIGH-RESOLUTION ZONE WITH DATA FROM THE LOW-RESOLUTION GRID + SMOOTHING ◀
# ──────────────────────────────────────────────────────────

    #Function used to catch NaNs in high_res wich can be filled and smooth thanks to low_res data, and do that
    #----------------
    def nan_buffer_linear_pond(blended_data, low_res_interp, buffer_width):
        # Create a mask for NaN values in blended_data where low_res_interp has values
        mask_nan = np.isnan(blended_data) & ~np.isnan(low_res_interp)
        new_z = blended_data
        
        # Create a structuring element of size buffer_width x buffer_width
        structure = np.ones((3, 3), dtype=bool)  # Structuring element for erosion and dilation
        mask_eroded = binary_erosion(mask_nan, structure=structure)

        for i in range(0,buffer_width-1):
            mask_eroded = binary_erosion(mask_eroded, structure=structure)
        
        # Apply binary dilation to the eroded mask
        mask_dilated = binary_dilation(mask_eroded, structure=structure)
        
        for i in range(0,buffer_width-1):
            mask_eroded = binary_erosion(mask_eroded, structure=structure)
        
        
        # Compute the transition mask as the difference between dilated and eroded masks
        mask_transition = mask_dilated & ~mask_eroded
        
        # Compute the interpolation mask for areas not covered by the dilated mask
        mask_interpolation = mask_nan & ~mask_dilated
        
        # Apply binary dilation to the transition mask to expand it
        transition_dilated = binary_dilation(mask_transition, structure=np.ones((3, 3)))
        # Compute transition borders as the difference between expanded and transition masks
        transition_borders = transition_dilated & ~mask_transition
        
        # Assign low resolution data to the areas covered by the eroded mask
        new_z[mask_eroded] = low_res_interp[mask_eroded]  # Low-res data in the eroded zone
        
        # Get coordinates of transition and transition border areas
        coords_mask_transition = np.array(np.where(mask_transition)).T
        coords_transition_borders_int = np.array(np.where(transition_borders & mask_eroded)).T
        coords_transition_borders_ext = np.array(np.where(transition_borders & ~mask_eroded)).T
        
        # Build KDTree for internal and external transition borders
        tree_transition_borders_int = cKDTree(coords_transition_borders_int)
        tree_transition_borders_ext = cKDTree(coords_transition_borders_ext)
        
        # Compute minimum distances from mask transition coordinates to internal and external borders
        distances_int, indices_int = tree_transition_borders_int.query(coords_mask_transition)
        distances_ext, indices_ext = tree_transition_borders_ext.query(coords_mask_transition)
        
        def get_values_from_indices(indices, coords_borders):
            # Extract values from z2_combined based on indices
            values = np.full(len(indices), np.nan)
            for i, idx in enumerate(indices):
                if idx < len(coords_borders):
                    x, y = coords_borders[idx]
                    values[i] = new_z[x, y]
            return values
        
        # Get values from internal and external border coordinates
        values_int = get_values_from_indices(indices_int, coords_transition_borders_int)
        values_ext = get_values_from_indices(indices_ext, coords_transition_borders_ext)
        
        # Compute weighted average for transition values
        tot_dist = distances_int + distances_ext
        weight_out = distances_int / tot_dist
        weight_in = distances_ext / tot_dist
        ponderated_transition = weight_in * values_int + weight_out * values_ext
        
        def apply_ponderated_values(new_z, mask_transition, ponderated_values, coords_mask_transition):
            # Apply weighted values to z2_combined based on transition mask
            for i, (coord, value) in enumerate(zip(coords_mask_transition, ponderated_values)):
                x, y = coord
                if mask_transition[x, y]:
                    new_z[x, y] = value
            return new_z
        
        # Apply weighted values to the combined matrix
        new_z = apply_ponderated_values(new_z, mask_transition, ponderated_transition, coords_mask_transition)
    
        return new_z
    #----------------

    rows, cols = np.where(mask_overlap)

    # Définir les limites du rectangle autour de la région de True
    min_row, max_row = rows.min(), rows.max()
    min_col, max_col = cols.min(), cols.max()
    
    # Extraire la sous-section rectangulaire du tableau
    zoom_z1 = z1_interp[min_row:max_row+1, min_col:max_col+1]
    zoom_z2 = z2_interp[min_row:max_row+1, min_col:max_col+1]

    # List of buffer widths to apply
    buffer_widths = [10, 5, 3, 2]

    # Initialize the variable with the initial high-resolution grid
    #z_less_nan = z2_interp
    z_tempo = zoom_z1

    # Apply buffers of decreasing widths
    for width in buffer_widths:
        print(f'⚙️ Filling NaN areas with lower-resolution data and smoothing the buffer of {width} cells...')
        #z_less_nan = nan_buffer_linear_pond(z_less_nan, z2_save, buffer_width=width)
        z_tempo = nan_buffer_linear_pond(z_tempo, zoom_z2, buffer_width=width)
        print(f'✅ Successfully created a {width}-cell buffer, enhancing data quality around NaN areas.')
        #---> Clean memory:
        gc.collect()
    
    print('🎉 Groups of NaNs in the high-resolution grid have been filled with lower-resolution grid data and smoothed.')

    # Replace the values in z2_interp (low-res grid) with those from grid 1 in the overlapping region
    z2_interp[min_row:max_row+1, min_col:max_col+1] = z_tempo
    z_less_nan=  z2_interp
    
    #---> Clean memory:
    #del z2_save, z2_interp
    del z1_interp, z_tempo
    gc.collect()

# ──────────────────────────────────────────────────────────
# ▶ SECTION 5: FILLING THE REMAINING NANs◀
# ──────────────────────────────────────────────────────────
    '''
    # Create axes for the coordinates
    lon_axis = pyinterp.Axis(np.ravel(new_lon_2))
    lat_axis = pyinterp.Axis(np.ravel(new_lat_2))
    
    # Create a Grid2D object with the axes
    data = pyinterp.Grid2D(lat_axis, lon_axis, z_less_nan)

    # Apply LOESS interpolation with a 3x3 grid
    print('⚙️ Applying LOESS interpolation to fill punctual NaN ...')
    z = pyinterp.fill.loess(data, nx=3, ny=3)

    print('✅ The ponctual Nans have been filled with loess interpolation method')

    #---> Clean memory:
    del data, lon_axis, lat_axis
    gc.collect()
    '''
    # Size threshold to decide the processing method
    size_threshold = 1000  # Define this threshold according to your available memory
    chunk_size = 1000  # Size of each chunk for processing
    overlap = 25  # Define the overlap size

    # Create axes for the coordinates
    lon_axis = pyinterp.Axis(np.ravel(new_lon_2))
    lat_axis = pyinterp.Axis(np.ravel(new_lat_2))
    
    
    # Check the size of new_lon_2 to decide the method to use
    if new_lon_2.size > size_threshold:
        print('⚙️ The grid is large; using chunk processing with overlap ...')
        
        # Loop over the grid with overlap
        for i in range(0, z_less_nan.shape[0], chunk_size - overlap):
            for j in range(0, z_less_nan.shape[1], chunk_size - overlap):
                i_end = min(i + chunk_size, z_less_nan.shape[0])
                j_end = min(j + chunk_size, z_less_nan.shape[1])
                
                # Adjust indices for the axes based on the current chunk
                lat_chunk_axis = lat_axis[i:i_end]
                lon_chunk_axis = lon_axis[j:j_end]
                
                # Extract the chunk
                chunk_data = z_less_nan[i:i_end, j:j_end]
                data_chunk = pyinterp.Grid2D(lat_chunk_axis, lon_chunk_axis, chunk_data)
                
                # Apply LOESS on the chunk
                z_chunk = pyinterp.fill.loess(data_chunk, nx=3, ny=3)
                    
                # Integrate the result back into z_less_nan
                z_less_nan[i:i_end, j:j_end] = z_chunk
                
                # Clean up
                del data_chunk, chunk_data
                gc.collect()
    
    else:
        print('⚙️ The grid is small; applying LOESS interpolation directly ...')
        # Initialize the data grid
        data = pyinterp.Grid2D(lat_axis, lon_axis, z_less_nan)
        # Apply LOESS interpolation on the entire grid
        z = pyinterp.fill.loess(data, nx=3, ny=3)
    
    # Indication of the end of interpolation
    print('✅ The punctual NaNs have been filled using the LOESS interpolation method')
    
    #---> Clean memory:
    del data, lon_axis, lat_axis
    gc.collect()

# ──────────────────────────────────────────────────────────
# ▶ SECTION 6: BUFFER ZONE SMOOTHING ◀
# ──────────────────────────────────────────────────────────

    # CREATE THE BUFFER MASK
    
    # Grid dimensions
    grid_shape = new_lon_grid_2.shape
    
    # Get the indices where the overlap mask is True
    true_indices = np.array(np.where(mask_overlap))
    
    # Get the coordinates of the four edges of the mask_overlap
    top = true_indices[0].min()   # Coordinate of the first row with True (top)
    bottom = true_indices[0].max()  # Coordinate of the last row with True (bottom)
    left = true_indices[1].min()  # Coordinate of the first column with True (left)
    right = true_indices[1].max()  # Coordinate of the last column with True (right)
    
    # Buffer width
    min_buffer_width = 3  # Minimum buffer width
    
    # Initialize an empty buffer mask
    mask_buffer = np.zeros(grid_shape, dtype=bool)
    
    # Function to create a buffer if the distance from the edge is sufficient
    def create_buffer_on_side(start_idx, end_idx, buffer_limit, buffer_width, max_limit, side):
        if buffer_limit <= min_buffer_width:  # Not enough space for a proper buffer
            print(f"⚠️ Warning: Not enough space for a buffer on the {side} side. No smoothing applied.")
            return 0
        elif buffer_limit < buffer_width:  # Reduce buffer size
            print(f"⚠️ Warning: The buffer on the {side} side has been reduced to {buffer_limit} cells.")
            buffer_width = buffer_limit
        return buffer_width
    
    # Function to check for valid data outside the buffer
    def has_valid_data_outside(start_row, end_row, start_col, end_col):
        return np.any(new_lon_grid_2[start_row:end_row, start_col:end_col])
    
    # Create the buffer only on sides where the mask_overlap is not too close to the grid edges
    # and ensure that there is valid data to interpolate.
    
    # If the top of the mask_overlap is not at the top edge of the grid
    buffer_top = create_buffer_on_side(top, bottom, top, buffer_width, grid_shape[0], "top")
    if buffer_top > 0 and has_valid_data_outside(top - buffer_top, top, left, right + 1): 
        mask_buffer[top - buffer_top:top, left:right + 1] = True
    
    # If the bottom of the mask_overlap is not at the bottom edge of the grid
    buffer_bottom = create_buffer_on_side(bottom, bottom + buffer_width, grid_shape[0] - bottom - 1, buffer_width, grid_shape[0], "bottom")
    if buffer_bottom > 0 and has_valid_data_outside(bottom + 1, bottom + 1 + buffer_bottom, left, right + 1): 
        mask_buffer[bottom + 1:bottom + 1 + buffer_bottom, left:right + 1] = True
    
    # If the left of the mask_overlap is not at the left edge of the grid
    buffer_left = create_buffer_on_side(left, right, left, buffer_width, grid_shape[1], "left")
    if buffer_left > 0 and has_valid_data_outside(top, bottom + 1, left - buffer_left, left): 
        mask_buffer[top:bottom + 1, left - buffer_left:left] = True
    
    # If the right of the mask_overlap is not at the right edge of the grid
    buffer_right = create_buffer_on_side(right, right + buffer_width, grid_shape[1] - right - 1, buffer_width, grid_shape[1], "right")
    if buffer_right > 0 and has_valid_data_outside(top, bottom + 1, right + 1, right + 1 + buffer_right): 
        mask_buffer[top:bottom + 1, right + 1:right + 1 + buffer_right] = True
    
    # Include corners if necessary
    def include_corners(mask_buffer, top, bottom, left, right, buffer_width):
        # Top-left corner
        if top - buffer_width >= 0 and left - buffer_width >= 0:
            mask_buffer[top - buffer_width:top, left - buffer_width:left] = True
        # Top-right corner
        if top - buffer_width >= 0 and right + buffer_width < grid_shape[1]:
            mask_buffer[top - buffer_width:top, right + 1:right + 1 + buffer_width] = True
        # Bottom-left corner
        if bottom + buffer_width < grid_shape[0] and left - buffer_width >= 0:
            mask_buffer[bottom + 1:bottom + 1 + buffer_width, left - buffer_width:left] = True
        # Bottom-right corner
        if bottom + buffer_width < grid_shape[0] and right + buffer_width < grid_shape[1]:
            mask_buffer[bottom + 1:bottom + 1 + buffer_width, right + 1:right + 1 + buffer_width] = True
    
    include_corners(mask_buffer, top, bottom, left, right, buffer_width)
    
    # Exclude the overlap area from the buffer
    mask_buffer[mask_overlap] = False

    print('✅ The buffer mask around high resolution grid has been created')

    # BUFFER's BORDERS MASKS
    dilated_mask_buffer= binary_dilation(mask_buffer)
    border_buffer_ext = dilated_mask_buffer & ~mask_buffer & ~mask_overlap
    border_buffer_inner = mask_overlap & dilated_mask_buffer
    
    # LINEAR PONDERATION ON BUFFER WITH BORDERS' DATA
    # Get coordinates of borders and buffer
    coords_border_inner = np.array(np.where(border_buffer_inner)).T
    coords_border_outer = np.array(np.where(border_buffer_ext)).T
    coords_buffer = np.array(np.where(mask_buffer)).T
    
    # Get values from z2_interp at the border points
    values_inner = z[border_buffer_inner]
    values_outer = z[border_buffer_ext]
    
    # Build KDTree for inner and outer borders
    tree_border_inner = cKDTree(coords_border_inner)
    tree_border_outer = cKDTree(coords_border_outer)
    
    # Initialize array for interpolated values
    interpolated_values_buffer = np.zeros(len(coords_buffer))
    
    # For each point in the buffer, find nearest border points
    for i, coord in enumerate(coords_buffer):
        # Query the KDTree for nearest points
        dist_inner, idx_inner = tree_border_inner.query(coord)
        dist_outer, idx_outer = tree_border_outer.query(coord)
        
        # Extract values of nearest border points
        value_inner = values_inner[idx_inner]
        value_outer = values_outer[idx_outer]
        
        # Calculate total distance and weights
        total_distance = dist_inner + dist_outer
        if total_distance == 0:
            # Handle the case where total distance is zero
            print(f"⚠️ Warning: Total distance is zero for point {coord}")
            interpolated_values_buffer[i] = np.nan
        else:
            weight_inner = dist_outer / total_distance
            weight_outer = dist_inner / total_distance
            
            # Calculate interpolated value
            interpolated_values_buffer[i] = weight_inner * value_inner + weight_outer * value_outer
            
    z[mask_buffer] = interpolated_values_buffer
    
    print('✅ The mask_buffer has been smoothed with linear ponderation')

# ──────────────────────────────────────────────────────────
# ▶ SECTION 7: SAVING TO NETCDF FORMAT ◀
# ──────────────────────────────────────────────────────────

    # Get unique coordinates
    new_lon_unique = np.unique(new_lon_grid_2)
    new_lat_unique = np.unique(new_lat_grid_2)

    ### OLD WAY TO SAVE DATA AS DOUBLE ##########################################################################
    # Create a DataArray for the interpolated data
    #ds_interpolated = xr.DataArray(
    #    z, 
    #    coords=[('lat', new_lat_unique), ('lon', new_lon_unique)],  # Add coordinates
    #    dims=['lat', 'lon']  # Specify dimensions
    #)
    
    # Create a Dataset to organize the variables
    #ds_to_save = xr.Dataset({
    #    'topo': ds_interpolated  # Add the interpolated variable
    #})
    
    # SAVE WITH CHUNKS TO LIMIT MEMORY USAGE
    #ds_to_save.to_netcdf(output_file, encoding={'topo': {'chunksizes': (100, 100)}})
    #############################################################################################################

    # Checking source of grid 2
    try:
        topo_type= topo_reader.lookvar(low_res)

    except:
        topo_type= { 'lon':'lon',\
                 'lat':'lat',\
                 'topo':'topo',\
                 'zaxis':'up'\
               }

    # Create a DataArray for the interpolated data and cast it to float32
    ds_interpolated = xr.DataArray(
        z.astype('float32'),  # Cast to float32
        coords={topo_type['lat']: new_lat_unique, topo_type['lon']: new_lon_unique},  # Add coordinates dynamically
        dims=[topo_type['lat'], topo_type['lon']]  # Specify dimensions dynamically
    )
    
    # Create a Dataset to organize the variables
    ds_to_save = xr.Dataset({
        topo_type['topo']: ds_interpolated  # Add the interpolated variable
    })
    
    # SAVE WITH CHUNKS TO LIMIT MEMORY USAGE
    ds_to_save.to_netcdf(output_file, encoding={topo_type['topo']: {'chunksizes': (100, 100), 'dtype': 'float32'}})


    print(f"📂 The merged grid has been saved as: {output_file}")
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------

# =========================================
# FUNCTION FOR R.MBLENDMETHOD (GRASS TOOLS)
# =========================================



def mblend_smoothing(high_res, low_res, far_edge, output_file, gisbase, EPSG_num):
    """
    Python command lines for creating a temporary grass dataset and running r.mblend addon over high/low resolution grids
    More info: "https://grass.osgeo.org/grass84/manuals/addons/r.mblend.html"

    Mblend: DEM merging method proposed by Leitão et al. (2016). It deals with cases where a study area is only partially covered by a high resolution DEM,
    with a coarser DEM available for the remainder (as in the case shown below). r.mblend merges the two DEMs, producing a smooth transition from the high resolution DEM to the   
    low resolution DEM.
        
    Parameters:
    ----------
    high_res (str): Path to the high resolution topography file.
    low_res (str): Path to the low resolution topography file.
    output_file (str): Path to the output NetCDF file.
    gisbas (str): Path to the GRASS GIS installation. ex: "/usr/lib/grass78"
    EPSG_num (str): EPSG (referential) number where grid will be projected. ex:'EPSG:4326'
    far_edge (int): Percentage of distance to high resolution raster used to determine far edge. Number between 0 and 100.
                    When the blending occurs along a single edge a number closer to 100 tends to produce more even results. With more blending edges (e.g. high resolution DEM 
                    sits on the middle of the low resolution DEM) a lower number may produce a more regular blend.
    
    """

    try:
        import subprocess
        import shutil
        import binascii
        import tempfile
        import grass.script as gscript
        import grass.script.setup as gsetup
    except ImportError as e:
        print(f"Required module is missing: {e}. Please install GRASS GIS and related dependencies.")

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
    
    
    
    
    
    
    
    
    
    
    
    



















