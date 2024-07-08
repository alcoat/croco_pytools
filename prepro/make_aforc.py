#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 01:21:22 2024

@author: menkes
"""
#--- Dependencies ---------------------------------------------------------
import xarray as xr
import dask
import pylab as plt
import numpy as np
import glob as glob
from dateutil.relativedelta import relativedelta
import json
import time
import os
# Readers
import sys
sys.path.append("/Readers")
from aforc_reader import lookvar
from aforc_class import aforc_class, create_class
# Extrapolation sst
sys.path.append("/Modules")
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
input_prefix = '*' # For multifiles, if the name of the file begin with the variable name, just write '*'
multi_files = False # If one file per variable in input : True

# -------------------------------------------------
# OUTPUT :
# -------------------------------------------------
output_dir = '/path/out/'
output_file_format = "DAILY" # How outputs are split (MONTHLY,DAILY)

# -------------------------------------------------
# Grid size : 
ownArea = 1 # 0 if area from croco_grd.nc +/- 5°
            # 1 if own area
if ownArea == 0:
    croco_grd = '../../CONFIGS/your_run_folder/CROCO_FILES/croco_grd.nc'
else:
    lon_min,lon_max,lat_min,lat_max = 6,24,-40,-24

# Dates limits
Yorig = 2000                 # year defining the origin of time as: days since Yorig-01-01
Ystart, Mstart = 2005,1   # Starting month
Yend, Mend  = 2005,1    # Ending month

# -------------------------------------------------
# OPTIONS :
# -------------------------------------------------
# To convert the atmospheric pressure : True
READ_PATM = False

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
with open('/Readers/Aforc_variables.json', 'r') as jf:
    croco_variables = json.load(jf)

# -------------------------------------------------
# Input files : find paths
# Warning : only for path type /data_origin/file.nc or /data_origin/Y/M/file.nc
# -------------------------------------------------
input_file = find_input(variables,input_dir,input_prefix,Ystart,Mstart,Yend,Mend,multi_files)
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
        ind = 0
        for var in variables.raw_name: # Loop on variables
            if var == 'lon' or var == 'lat':
                continue
            elif variables.get_filename('u10m') == variables.get_filename('v10m') and var == 'v10m':
                continue
            elif var == 'msl' and not READ_PATM:
                continue 
            #-------------------------
            if ind == 0:
                dataxr = xr.open_mfdataset(eval(var),combine='nested',concat_dim='time', engine='netcdf4', chunks={'time': 'auto',  variables.get_var('lon'): 'auto',  variables.get_var('lat'): 'auto'})
            else:
                dataxr = xr.merge([dataxr,xr.open_mfdataset(eval(var),combine='nested',concat_dim='time', engine='netcdf4', chunks={'time': 'auto',  variables.get_var('lon'): 'auto',  variables.get_var('lat'): 'auto'})])
            ind +=1
    else: # If not multi_files :
        dataxr = xr.open_mfdataset(input_file,combine='nested',concat_dim='time', engine='netcdf4', chunks={'time': 'auto',  variables.get_var('lon'): 'auto',  variables.get_var('lat'): 'auto'})
   
    # Rename lon/lat
    if variables.get_var('lon') != 'lon': 
        dataxr = dataxr.rename({variables.get_var('lon') : 'lon'})
        dataxr = dataxr.rename({variables.get_var('lat') : 'lat'})
    
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
        if data_origin == 'era_ecmwf' or data_origin == 'cfsr': # latitudes are from lat_max to lat_min
            dataxr_inprocess = dataxr.sel( time = slice(start_date,end_date) , \
                                lon = slice(lon_min,lon_max) ,\
                                lat = slice(lat_max,lat_min) )
        else:
            dataxr_inprocess = dataxr.sel( time = slice(start_date,end_date) , \
                                lon = slice(lon_min,lon_max) ,\
                                lat = slice(lat_min,lat_max) )
        
        
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
                if var == 't2m':
                    data = kelvin_2_celsius(data)
                if var == 'str':
                    sst = extrapolation(data_grouped[i][variables.get_var('sst')].values,data_grouped[i].lon.values,data_grouped[i].lat.values)
                    data = strd_calculation(data,sst,variables,croco_variables)
                elif var == 'q':
                    data = r_calculation(data,data_grouped[i][variables.get_var('t2m')],croco_variables)
                elif var == 'uswrf':
                    data = ssr_calculation(data,data_grouped[i][variables.get_var('dswrf')],croco_variables)
                else:
                    data = attr(data,var,variables,croco_variables)
                data = time_origin(data,Yorig)
                data = metadata(data)
                data = missing_data(data)
                
# -----------------------------------
# PUT DATA IN OUTPUT FILE
# -----------------------------------
                data = create_netcdf(data,output_dir,output_file_format)


end_time = time.time()
time_taken = end_time - start_time
print("Computation time:", time_taken, "sec")













