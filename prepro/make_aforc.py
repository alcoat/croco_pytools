#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 14:47:46 2024

@author: annelou
"""
#--- Dependencies ---------------------------------------------------------
import xarray as xr
import pylab as plt
import numpy as np
import glob as glob
from dateutil.relativedelta import relativedelta
import json
import time
import os
# Readers
import sys
sys.path.append("./Readers")
from aforc_reader import lookvar
from aforc_class import aforc_class, create_class
sys.path.append("./Modules")
from interp_tools import make_xarray
from aforc_netcdf import *
from data_consistency import consistency
from aforc_transformation import *
import scipy.interpolate as itp
import pyinterp.backends.xarray

# *******************************************************************************
#                         U S E R  *  O P T I O N S
# *******************************************************************************

# -------------------------------------------------
# INPUT :
# -------------------------------------------------
data_origin = 'era_dataref' # era_dataref, era_ecmwf, cfsr
input_dir = '/path/in/'
input_prefix = 'era_5-copernicus__*' # For multifiles, if the name of the file begin with the variable name, just write '*'
#input_prefix = '*' # For multifiles, if the name of the file begin with the variable name, just write '*'
multi_files = False # If one file per variable in input : True

# -------------------------------------------------
# OUTPUT :
# -------------------------------------------------
output_dir = '/path/out/'
output_file_format = "MONTHLY" # How output files are split (MONTHLY,DAILY)

# -------------------------------------------------
# Grid size : 
ownArea = 0 # 0 if area from croco_grd.nc +/- 5°
            # 1 if own area
if ownArea == 0:
    croco_grd = '/pathin/to/your/croco/grid/croco_grd.nc'
else:
    lon_min,lon_max,lat_min,lat_max = 4,27,-40,-24

# Dates limits
Yorig = 1950                 # year defining the origin of time as: days since Yorig-01-01
Ystart, Mstart = 1980,3   # Starting month
Yend, Mend  = 1981,1  # Ending month

# -------------------------------------------------
# OPTIONS :
# -------------------------------------------------
# To convert the atmospheric pressure : True
READ_PATM = False

# If there is no STRD variable in raw data, it will be calculated with STR and SST. 
# SST may or may not be extrapolated to the coast. 
# If it is not, this may result in temperature spike at the coast at certain points,
# however, note that extrapolation increases pre-processing time.
# If STRD is in raw data, extrapolation_sst will no be considered.
extrapolation_sst = True # /!\ if sst = ts (surf temp) no need

# *****************************************************************************
#                      E N D     U S E R  *  O P T I O N S
# *****************************************************************************

# -------------------------------------------------
# Read variables
# -------------------------------------------------
variables = create_class(lookvar(data_origin),multi_files)

# -------------------------------------------------
# Setting processed output directory
# -------------------------------------------------
# Get the current directory
os.makedirs(output_dir,exist_ok=True)

# -------------------------------------------------
# python Dictionary from JSON file
# -------------------------------------------------
# with open('ERA5_variables.json', 'r') as jf:
with open('croco_variables.json', 'r') as jf:
    croco_variables = json.load(jf)

# -------------------------------------------------
# Input files : find paths
# Warning : only for path type /data_origin/file.nc or /data_origin/Y/M/file.nc
# -------------------------------------------------
input_file = find_input(variables,input_dir,input_prefix,Ystart,Mstart,Yend,Mend,multi_files,READ_PATM)
print(input_file)

# -------------------------------------------------
# Read croco grid to find emprise
# -------------------------------------------------
if ownArea == 0:
    grid = xr.open_dataset(croco_grd)
    lon_min,lon_max,lat_min,lat_max = np.floor(grid['lon_rho'].min())-5,np.ceil(grid['lon_rho'].max())+5,np.floor(grid['lat_rho'].min())-5,np.ceil(grid['lat_rho'].max())+5

# *****************************************************************************
#                                MAIN FUNCTION
# *****************************************************************************
start_time = time.time()

if __name__ == "__main__":   # if multi files, 'input_file' not used 
# -----------------------------------
# READ THE METADATA
# -----------------------------------
    if multi_files: # Need to merge variable input files
        for var in range(len(input_file)): # Loop on variables
            if var == 0:
                dataxr = xr.open_mfdataset(input_file[var],combine='nested',concat_dim='time', engine='netcdf4', chunks={'time': 'auto',  variables.get_var('lon'): 'auto',  variables.get_var('lat'): 'auto'})
            else:
                dataxr = xr.merge([dataxr,xr.open_mfdataset(input_file[var],combine='nested',concat_dim='time', engine='netcdf4', chunks={'time': 'auto',  variables.get_var('lon'): 'auto',  variables.get_var('lat'): 'auto'})])    
    else: # If not multi_files :
        dataxr = xr.open_mfdataset(input_file,combine='nested',concat_dim='time', engine='netcdf4', chunks={'time': 'auto',  variables.get_var('lon'): 'auto',  variables.get_var('lat'): 'auto'})
   
    # Rename lon/lat
    if variables.get_var('lon') != 'lon': 
        dataxr = dataxr.rename({variables.get_var('lon') : 'lon'})
        dataxr = dataxr.rename({variables.get_var('lat') : 'lat'})
    
    # Il y aura toujours lon et lat en variable mais pas forcement en dim (peut etre x,y notamment pour grille irreg et coord 2D)
    lon_dim = dataxr['lon'].dims
    lat_dim = dataxr['lat'].dims
    
    
    # Script to know if irregular/regular grid
    if len(lon_dim) == 2 and len(lat_dim) == 2:
        irreg = 1
        # Find the indices to cut longitudes and latitudes
        # Only need the first time
        var_temp = list(variables.raw_name.items())[2]
        if var_temp[0] != 'lon' and var_temp[0] != 'lat':
            ix_min,ix_max,iy_min,iy_max = ind_irreg_grid(dataxr.isel(time=0),var_temp[1],lon_min,lon_max,lat_min,lat_max)
            
        else:
            print('Please put the longitude and the latitude in the 1st and 2nd place in Readers/aforc_reader.py')
            sys.exit()
    else:
        
        irreg = 0


        
        
    # loop on years
    for year_inprocess in range(Ystart,Yend+1): 
        
        # Put start and end date to the right format
        if year_inprocess == Ystart: 
            start_date = plt.datetime.datetime(year_inprocess,Mstart,1,0)
        else:
            start_date = plt.datetime.datetime(year_inprocess,1,1,0)
        
        if year_inprocess == Yend:
            end_date = plt.datetime.datetime(year_inprocess,Mend,1,0) + relativedelta(months=1,hours=-1) # Last day of the ending month
        else:
            end_date = plt.datetime.datetime(year_inprocess,12,1,0) + relativedelta(months=1,hours=-1) # Last day of the ending month
        
# ------------------------------------
# SELECT AREA AND READ DATA OVER IT
# ------------------------------------

        if irreg == 0: # Regular grid : Coords = 1D array
        # The selection is made by values
        # First, will check if latitudes are inversed
        
        
            if dataxr['lat'][1].values > dataxr['lat'][0].values: # from lat_min to lat_max
                dataxr_inprocess = dataxr.sel( time = slice(start_date,end_date) , \
                                    lon = slice(lon_min,lon_max) ,\
                                    lat = slice(lat_min,lat_max) )
            else:
                dataxr_inprocess = dataxr.sel( time = slice(start_date,end_date) , \
                                    lon = slice(lon_min,lon_max) ,\
                                    lat = slice(lat_max,lat_min) )                        
        else: # Irregular grid : Coords = 2D array
        # The selection is made by indices
        # For now, suppose latitudes are from lat_min to lat_max
            sel_args = {lon_dim[0]: slice(iy_min,iy_max),lon_dim[1]:slice(ix_min,ix_max)}
            dataxr_inprocess = dataxr.sel(time=slice(start_date,end_date)).isel(**sel_args)
            
        
        if output_file_format == 'MONTHLY':
            data_grouped = dataxr_inprocess.groupby('time.month')
            
        elif output_file_format == 'DAILY':
            data_grouped = dataxr_inprocess.groupby('time.dayofyear')
        
        # Creation of a group index
        labels = list(data_grouped.groups.keys())
        index_to_label = {i: label for i, label in enumerate(labels)}
        
        # loop on groups (months or days)
        for ii in range(len(labels)): 
            i = index_to_label[ii]

            print('\n-----------------------------------')
            if output_file_format.upper() == "DAILY":
                print(' Processing Year %s - Month %02i - Day %02i' %(data_grouped[i]['time'].dt.year[0].item(),data_grouped[i]['time'].dt.month[0].item(),data_grouped[i]['time'].dt.day[0].item()))
            elif output_file_format.upper() == "MONTHLY":
                print(' Processing Year %s - Month %02i' %(data_grouped[i]['time'].dt.year[0].item(),data_grouped[i]['time'].dt.month[0].item()))
            print('-----------------------------------')
            
            # Data consistency : if not enough data, will stop
            consistency(data_grouped[i])

            # loop on variables
            for var in variables.raw_name: 
                if var == 'msl' and not READ_PATM:
                    continue
                elif var == 'lon' or var == 'lat' or var == 'dswrf' or var == 'sst':
                    continue
                print('  Processing variable: var/vname :'+var+'/'+variables.get_var(var))

# -----------------------------------
# CONVERT DATA IF NECESSARY
# -----------------------------------
                data = flip_data(data_grouped[i][variables.get_var(var)],data_origin)
                data = unit_conversion(data,var,variables)
                if var == 'str':
                    sst = flip_data(data_grouped[i][variables.get_var('sst')],data_origin)
                    if extrapolation_sst:
                        sst = extrapolation(sst.values,sst.lon.values,sst.lat.values)
                    else:
                        sst = sst.values
                    data = strd_calculation(data,sst,variables,croco_variables)
                    var = 'strd'
                elif var == 'q':
                    data = r_calculation(data,data_grouped[i][variables.get_var('t2m')],croco_variables)
                    var = 'r'
                elif var == 'uswrf':
                    data = ssr_calculation(data,data_grouped[i][variables.get_var('dswrf')],croco_variables)
                    var = 'ssr'
                else:
                    data = attr(data,var,variables,croco_variables)
                data = time_origin(data,Yorig)
                data = metadata(data)
                data,encoding = missing_data(data,var)
                
# -----------------------------------
# PUT DATA IN OUTPUT FILE
# -----------------------------------
                data = create_netcdf(data,output_dir,output_file_format,encoding)


end_time = time.time()
time_taken = end_time - start_time
print("Computation time:", time_taken, "sec")













