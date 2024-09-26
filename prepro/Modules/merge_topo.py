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

# /!\ The dependencies used for r.mblend method (2nd function) are imported inside the function

#--- Functions ---------------------------------------------------------

# ================================================
# FUNCTION FOR MERGING - LINEAR PONDERATION METHOD
# ================================================

def merge_smooth(high_res, low_res, buffer_width, output_file, coarsen_factor=None, downscale_bounds=None):

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
    - coarsen_factor: (int) Resolution reduction factor for high-resolution grid. OPTION
    - downscale_bounds: (List) [lon_min, lon_max, lat_min, lat_max] for downscaling the low-resolution grid. OPTION

    Returns:
    -------
    None
        The function saves the processed grid to the specified output file and does not return any value.

    Notes:
    ------
    - The buffer area is created around the high resolution grid with a 'mergin_area' extent which can be lowered if high_res borders too close from the border of low_res ones
    - The function performs linear ponderation within this buffer using the high/low resolution grids.
    - /!\ Ensure that the input grids are compatible in terms of coordinate systems and dimensions.
    - The output file will be in NetCDF format with the processed data.
    
    """


# ──────────────────────────────────────────────────────────
# ▶ SECTION 1: Data recovery and conversion◀
# ──────────────────────────────────────────────────────────

    # Open the datasets
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
    if proj_high_res and proj_low_res:
        if proj_high_res != proj_low_res:
            raise ValueError("❌ The high resolution and low resolution datasets are in different projections.")
    else:
        if proj_high_res is None:
            print("⚠️ Warning: High resolution dataset does not have a defined CRS.")
        if proj_low_res is None:
            print("⚠️ Warning: Low resolution dataset does not have a defined CRS.")
    
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
    for coord_name in ds1.coords:
        coord_data = ds1[coord_name]
        # Assume lon and lat are typically coordinates with geographic values
        if coord_name.lower() in {'lon', 'longitude'}:
            lon_coord_ds1 = coord_name
        elif coord_name.lower() in {'lat', 'latitude'}:
            lat_coord_ds1 = coord_name
        elif coord_data.dims == ('x',) or coord_data.dims == ('y',):
            if 'x' in coord_data.dims:
                lon_coord_ds1 = coord_name
            if 'y' in coord_data.dims:
                lat_coord_ds1 = coord_name
    
    # Identify the coordinate names for latitude/longitude in ds2
    for coord_name in ds2.coords:
        coord_data = ds2[coord_name]
        # Assume lon and lat are typically coordinates with geographic values
        if coord_name.lower() in {'lon', 'longitude'}:
            lon_coord_ds2 = coord_name
        elif coord_name.lower() in {'lat', 'latitude'}:
            lat_coord_ds2 = coord_name
        elif coord_data.dims == ('x',) or coord_data.dims == ('y',):
            if 'x' in coord_data.dims:
                lon_coord_ds2 = coord_name
            if 'y' in coord_data.dims:
                lat_coord_ds2 = coord_name

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

    ### OLD WAY TO INTERPOLATE ##################################################################################
    # Interpolate ds2 data onto the new grid
    #z2_interp = griddata((lon_flat_2, lat_flat_2), z_flat_2, (new_lon_grid_2, new_lat_grid_2), method='nearest')
    #############################################################################################################

    ### NEW WAY TO INTERPOLATE ##################################################################################
    # Create a DataArray for z2 with unique 1D latitude and longitude coordinates
    z2_dataarray = xr.DataArray(
        z2,
        coords={
            'latitude': lat2,  # 1D array of latitudes
            'longitude': lon2    # 1D array of longitudes
        },
        dims=['latitude', 'longitude']  # Dimensions should match the DataArray shape
    )
    
    # Initialize the grid interpolator for z2
    interpolator_2 = pyinterp.backends.xarray.RegularGridInterpolator(z2_dataarray, geodetic=False)
    
    # Prepare new coordinates for interpolation
    new_coords_2 = {
        'longitude': new_lon_grid_2.flatten(),  # Flattened array of new longitudes
        'latitude': new_lat_grid_2.flatten()     # Flattened array of new latitudes
    }
    
    # Perform interpolation using the nearest method
    z2_interp = interpolator_2(new_coords_2, method='nearest').reshape(new_lon_grid_2.shape)

    print("✅ Low-resolution grid interpolation completed.")
    #############################################################################################################


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

    ### OLD WAY TO INTERPOLATE ##################################################################################
    # Interpolate the high-resolution data (ds1) onto the overlapping area of ds2
    #z1_interp_on_z2_overlap = griddata(
        #(lon_flat_1, lat_flat_1),  # Source points (grid 1)
        #z_flat_1,                  # Source values (grid 1)
        #(new_lon_grid_2[mask_overlap], new_lat_grid_2[mask_overlap]),  # Target points (overlapping region in grid 2)
        #method='nearest'           # Interpolation method (nearest neighbor)
    #)
    #############################################################################################################

    ### NEW WAY TO INTERPOLATE ##################################################################################
    # Create a DataArray for z1
    z1_dataarray = xr.DataArray(
        z1,  # The values to be interpolated
        coords={
            'lat': lat1,  # 1D array of latitudes corresponding to z1
            'lon': lon1   # 1D array of longitudes corresponding to z1
        },
        dims=['lat', 'lon']  # Define the dimensions of the DataArray
    )
    
    # Create the interpolator for z1
    interpolator = pyinterp.backends.xarray.RegularGridInterpolator(z1_dataarray, geodetic=False)
    
    # Prepare new coordinates for interpolation within the overlap region
    new_coords_1 = {
        'lon': new_lon_grid_2[mask_overlap].flatten(),  # Flattened longitudes from the overlap region
        'lat': new_lat_grid_2[mask_overlap].flatten()    # Flattened latitudes from the overlap region
    }
    
    # Apply the interpolator to the new coordinates
    z1_interp_on_z2_overlap = interpolator(new_coords_1, method='nearest')
    
    # Reshape the result to match the shape of the grid 2
    z1_interp_on_z2_overlap = z1_interp_on_z2_overlap.reshape(new_lon_grid_2[mask_overlap].shape)
    #############################################################################################################

    # Save the original interpolated ds2 grid for later use
    z2_save = z2_interp.copy()
    
    # Replace the values in z2_interp (low-res grid) with those from grid 1 in the overlapping region
    z2_interp[mask_overlap] = z1_interp_on_z2_overlap

    print("✅ High-resolution grid interpolation completed.")

    #---> Clean memory:
    del z1_interp_on_z2_overlap, z2_dataarray, z1_dataarray, lon_grid_1, lat_grid_1, ds1, ds2
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

    # Apply buffer with width 10
    z_less_nan = nan_buffer_linear_pond(z2_interp, z2_save, buffer_width=10)
    
    # Apply buffer with width 5
    z_less_nan = nan_buffer_linear_pond(z_less_nan, z2_save, buffer_width=5)
    
    # Apply buffer with width 3
    z_less_nan = nan_buffer_linear_pond(z_less_nan, z2_save, buffer_width=3)
    
    # Apply buffer with width 2
    z_less_nan = nan_buffer_linear_pond(z_less_nan, z2_save, buffer_width=2)

    print('✅ Groups of Nans in the high resolution grid have been filled with lower resolution grid data and smoothed')

    #---> Clean memory:
    del z2_save, z2_interp
    gc.collect()

# ──────────────────────────────────────────────────────────
# ▶ SECTION 5: FILLING THE REMAINING NANs◀
# ──────────────────────────────────────────────────────────

    # Create axes for the coordinates
    lon_axis = pyinterp.Axis(np.ravel(new_lon_2))
    lat_axis = pyinterp.Axis(np.ravel(new_lat_2))
    
    # Create a Grid2D object with the axes
    data = pyinterp.Grid2D(lat_axis, lon_axis, z_less_nan)
    
    # Apply LOESS interpolation with a 3x3 grid
    z = pyinterp.fill.loess(data, nx=3, ny=3)

    print('✅ The ponctual Nans have been filled with loess interpolation method')

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
    
    # Create a DataArray for the interpolated data
    ds_interpolated = xr.DataArray(
        z, 
        coords=[('lat', new_lat_unique), ('lon', new_lon_unique)],  # Add coordinates
        dims=['lat', 'lon']  # Specify dimensions
    )
    
    # Create a Dataset to organize the variables
    ds_to_save = xr.Dataset({
        'topo': ds_interpolated  # Add the interpolated variable
    })
    
    # SAVE WITH CHUNKS TO LIMIT MEMORY USAGE
    ds_to_save.to_netcdf(output_file, encoding={'topo': {'chunksizes': (100, 100)}})

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
    
    
    
    
    
    
    
    
    
    
    
    



















