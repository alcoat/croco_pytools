#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 17:22:13 2025

@author: annelou
"""
import numpy as np
import os
import sys
import glob
from netCDF4 import num2date
import logging

# ---------------------------------------------------
# FUNCTIONS USED BY make_aforc.py TO FIND INPUTFILES
# ---------------------------------------------------
def find_type_path(input_dir,year_inprocess,month_inprocess):
# Find the final path of data
# ---------------------------
# input_dir : path to the data (if the path is '/data_folder/Year/Month/', just write '/data_folder/')
# year_inprocess, month_inprocess : int, year and month in process
# -----------------------
# Warning : data must follow these kinds of paths :
#               - /input_dir/*.nc (or grb or grib)
#               - /input_dir/YYYY/MM/*.nc
#               - /input_dir/YYYY/*.nc
#               - /input_dir/MM/*.nc
    if os.path.exists(input_dir + str(year_inprocess)):
        pathin = input_dir + str(year_inprocess) + '/'
        if os.path.exists(input_dir + str(year_inprocess) + '/' + str(month_inprocess).zfill(2)):
            pathin += str(month_inprocess).zfill(2) + '/'
            
    elif os.path.exists(input_dir + str(month_inprocess).zfill(2)):
        pathin = input_dir + str(month_inprocess).zfill(2) + '/'
        
    else:
        pathin = input_dir
        
    return pathin

def find_input(variables,pathin,input_prefix,year_inprocess,month_inprocess,multi_files,READ_PATM,var=None):
# To find file paths and names
# ----------------------------
# variables : class of variables with the name of the variable, his factor of conversion to attempt the needed unity and the nomenclature of the variable file if multifile
# pathin : data path
# input_prefix : prefix of file names
# year_inprocess, month_inprocess : int, year and month in process
# multi_files : bool, if variables are separated in different files
# READ_PATM : to find data files for atmospheric pressure
# var : None if not multi_files ; name of the variable in process if multi_files
    if multi_files == True and var == None:
        print('You need to put a variable name in find_input fonction')
        sys.exit()
        
    input_file = []
    if multi_files:    
        if var == 'lon' or var == 'lat':
            pass
        elif var == 'msl' and not READ_PATM:
            pass
        else:
            if input_prefix[0] == '*':
                input_file=sorted(glob.glob(pathin + variables.get_filename(var) + input_prefix[1:] + str(year_inprocess) + str(month_inprocess).zfill(2)+'*'))
            else:
                input_file=sorted(glob.glob(pathin + input_prefix + variables.get_filename(var) + '*' + str(year_inprocess) + str(month_inprocess).zfill(2)+'*'))
    else:
        input_file=sorted(glob.glob(pathin + input_prefix[:-1] + str(year_inprocess) + str(month_inprocess).zfill(2) + '*' ))

    if not input_file:
        raise FileNotFoundError('Please follow these kinds of path : \n /input_dir/files  \n /input_dir/YYYY/MM/files \n /input_dir/YYYY/files \n /input_dir/MM/files')
         
    return input_file

def find_engine(input_file):
# To find if the file are netcdf or grib
# --------------------------------------
# input_file : one file name just to see its extension
    if input_file[-2:] == 'nc':
        enginein = 'netcdf4'
    elif input_file[-4:] == 'grib' or input_file[-3:] == 'grb':
        enginein = 'cfgrib'
        logging.getLogger('cfgrib').setLevel(logging.CRITICAL) # Cause grib returns many warnings
    return enginein

def remove_grb_indx(pathin):
    grib_idx = glob.glob(pathin + '*.idx')
    for iii in range(len(grib_idx)):
        os.remove(grib_idx[iii])
    return

# --------------------------------------------------------
# FUNCTIONS USED BY make_aforc.py TO FIND LON/LAT INDICES
# --------------------------------------------------------
# REGULARD GRID
#--------------
def find_nearest(array,value):
# array : 1D
# value : find the indice of 'value' in 'array'
    idx = (np.abs(array-value)).argmin()
    return idx

# IRREGULAR GRID
#---------------
def ind_irreg_grid(dataxr,var,lon_min,lon_max,lat_min,lat_max):
# dataxr : xarray.DataArray of one frame
# var : any variable
# lonmin,lonmax,latmin,latmax : values of wanted grid limits
    data_mask = dataxr.where((dataxr.lat < lat_max) & (dataxr.lat > lat_min) & (dataxr.lon > lon_min) & (dataxr.lon < lon_max))
    a = np.squeeze(data_mask[var].values)
    # Find longitude min :
    for i in range(len(a[0,:])):
        if np.where(~np.isnan(a[:,i]))[0].size > 0:
            xmin = i
            break
    # Find longitude max :
    for i in range(len(a[0,:])-1,-1,-1):
        if np.where(~np.isnan(a[:,i]))[0].size > 0:
            xmax = i
            break
    # Find latitude min :
    for i in range(len(a[:,0])):
        if np.where(~np.isnan(a[i,:]))[0].size > 0:
            ymin = i
            break
    # Find latitude max :
    for i in range(len(a[:,0])-1,-1,-1):
        if np.where(~np.isnan(a[i,:]))[0].size > 0:
            ymax = i
            break
    try:
        xmin,xmax,ymin,ymax
    except NameError:
        print("The croco_grid does not appear to be included in the atmospheric forcing grid.")
        sys.exit()
    return xmin,xmax,ymin,ymax


# ------------------------------------------------------
# FUNCTIONS USED BY make_aforc.py TO PREPARE NEW NETCDF
# ------------------------------------------------------
def attr(data,var,variables,croco_variables):
# data : xarray.DataArray of one variable
# var : name of the variable
# variables : class of variables with the name of the variable, \
# his factor of conversion to attempt the needed unity and the nomenclature of the variable file if multifile
# croco_variables : dictionary with long name and unity of the variable \
# (see ../Readers/croco_variables.json)
    data.attrs = { 'units': croco_variables[var][1], 'long_name': croco_variables[var][0] }
    data.name = var.upper()
    return data

def time_origin(data,Yorig):
# Time origin changement
# ---------------------
# data : xarray.DataArray of one variable
# Yorigin : Year to the time origin
    dt_origins = (np.datetime64(str(Yorig)+'-01-01T00:00:00.000000000') - np.datetime64('1900-01-01T00:00:00.000000000')) 
    dt_old = ( data['time'] - np.datetime64('1900-01-01T00:00:00.000000000')) 
    time_num = ((dt_old - dt_origins) / np.timedelta64(1, 's')) / (24 * 3600) # en jrs, precision sec
    data['time']=time_num
    data['time'].encoding['units'] = 'days since ' + str(Yorig) + '-1-1'
    data['time'].attrs['units'] = 'days since ' + str(Yorig) + '-1-1'
    return data

def metadata(data):
# data : xarray.DataArray of one variable
    data['lon'].attrs['long_name'] = 'longitude of RHO-points'
    data['lat'].attrs['long_name'] = 'latitude of RHO-points'
    data['time'].attrs['long_name'] = 'Time'
    data['lon'].attrs['units'] = 'degree_east'
    data['lat'].attrs['units'] = 'degree_north'
    return data

def missing_data(data,var):
# data : xarray.DataArray of one variable
# var : name of the variable
    encoding = {'lat': {'_FillValue': None},
            'lon': {'_FillValue': None},
            'time': {'_FillValue': None},
            var.upper(): {'_FillValue': 9999.0,
                          'missing_value': 9999.0}
            }
    
    return data,encoding

# ---------------------------------------------------
# FUNCTIONS USED BY make_aforc.py TO CREATE NETCDF
# ---------------------------------------------------
def output_name(data,output_dir,output_file_format):
# data : xarray.DataArray of one variable
# output_dir : path to output files
# output_file_format : DAILY or MONTHLY : indicate how output will be separate
    time = num2date(data.time[0].values, data.time.encoding['units'])
    if output_file_format.upper() == "MONTHLY":
        aforc_filename = output_dir + f'{data.name.upper()}_Y%sM%02i.nc' %(time.year,time.month)
    elif output_file_format.upper() == "DAILY":
        aforc_filename = output_dir + f"{data.name.upper()}_Y%sM%02iD%02i.nc"%(time.year,time.month,time.day)
    return aforc_filename

def create_netcdf(data,output_dir,output_file_format,encoding):
# data : xarray.DataArray of one variable
# output_dir : path to output files
# output_file_format : DAILY or MONTHLY : indicate how output will be separate
# encoding : give information to missing and fill values
    filename_out = output_name(data,output_dir,output_file_format)
    data = data.astype(np.float32)
    data = data.to_netcdf(filename_out,engine='netcdf4',encoding = encoding, compute='False')
    return data  















