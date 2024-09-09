#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  6 01:53:59 2024

@author: menkes
"""
import numpy as np
import os
import sys
import glob
from netCDF4 import num2date

# ---------------------------------------------------
# FUNCTIONS USED BY make_aforc.py TO FIND INPUTFILES
# ---------------------------------------------------
def find_input(variables,input_dir,input_prefix,Ystart,Mstart,Yend,Mend,multi_files,READ_PATM):
# Finding input_file paths
# -----------------------
# variables : class of variables with the name of the variable, his factor of conversion to attempt the needed unity and the nomenclature of the variable file if multifile
# input_dir : path to the data (if the path is '/data_folder/Year/Month/', just write '/data_folder/')
# Ystart,Mstart,Yend,Mend : Year and month of start and end
# multi_files : bool, if variables are separated in different files
# -----------------------
# Warning : only for path type /data_origin/file.nc or /data_origin/Y/M/file.nc
    if os.path.exists(input_dir + str(Ystart) + '/' + str(Mstart).zfill(2)):
        type_path = '/origin/Y/M/'
    else:
        type_path = '/origin/'
    
    if multi_files:
        varlist = []
        for var in variables.raw_name:
            if var == 'lon' or var == 'lat':
                continue
            globals()[var]=[]
            varlist.append(var)
    
    input_file = []
    for j in range(Ystart, Yend + 1):
        if j == Ystart: i = Mstart
        else: i = 1
        while i <= 12 and (j < Yend or (j == Yend and i <= Mend)):
            
            if multi_files:    
                for var in variables.raw_name:
                    if var == 'lon' or var == 'lat':
                        continue
                    elif variables.get_filename('u10m') == variables.get_filename('v10m') and var == 'v10m':
                        continue
                    elif var == 'msl' and not READ_PATM:
                        continue
                    else:
                        if type_path == '/origin/Y/M/':
                            eval(var).extend(glob.glob(input_dir + str(j) + '/' + str(i).zfill(2) + '/'+input_prefix +variables.get_filename(var)+'*.nc'))
                            input_file.extend(glob.glob(input_dir + str(j) + '/' + str(i).zfill(2) + '/'+input_prefix +variables.get_filename(var)+'*.nc'))
                        elif type_path == '/origin/':
                            eval(var).extend(sorted(glob.glob(input_dir + input_prefix +variables.get_filename(var)+'*.nc')))
                            input_file.append(sorted(glob.glob(input_dir + input_prefix +variables.get_filename(var)+'*.nc')))
                            i = 100 # end of while loop
                        else:
                            print('Please put your data either in a single common file or with a path like DATA_TYPE/YEAR/MONTH/file.nc')
                            sys.exit()
      
            else:
                if type_path == '/origin/Y/M/':
                    input_file.extend(glob.glob(input_dir + str(j) + '/' + str(i).zfill(2) + '/' + input_prefix))
                elif type_path == '/origin/':
                    input_file.extend(glob.glob(input_dir + input_prefix[:-1] + str(j) + '*' ))
                    i = 100 # end of while loop
                else:
                    print('Please put your data either in a single common file or with a path like DATA_TYPE/YEAR/MONTH/file.nc')
                    sys.exit()
            i += 1
    if not multi_files:
        input_file = sorted(input_file)  
    return input_file


# -----------------------------------------------------------------------------
# FUNCTIONS USED BY make_aforc.py TO FIND LON/LAT INDICES TO CUT IRREGULAR GRID
# -----------------------------------------------------------------------------
def ind_irreg_grid(dataxr,var,lonmin,lonmax,latmin,latmax):
# dataxr : xarray.DataArray of one frame
# var : any variable
# lonmin,lonmax,latmin,latmax : values of wanted grid limits
    data_mask = dataxr.where((dataxr.lat < latmax) & (dataxr.lat > latmin) & (dataxr.lon > lonmin) & (dataxr.lon < lonmax))
    a = np.squeeze(data_mask[var].values)
    # Find longitude min :
    for i in range(len(a[0,:])):
        if np.where(~np.isnan(a[:,i]))[0].size > 0:
            print(str(i),'not empty')
            xmin = i
            if xmin != 0:
                xmin = i -1
            break
    # Find longitude max :
    for i in range(len(a[0,:])-1,-1,-1):
        if np.where(~np.isnan(a[:,i]))[0].size > 0:
            print(str(i),'not empty')
            xmax = i +1
            break
    # Find latitude min :
    for i in range(len(a[:,0])):
        if np.where(~np.isnan(a[i,:]))[0].size > 0:
            print(str(i),'not empty')
            ymin = i
            if ymin != 0:
                ymin = i -1
            break
    # Find latitude max :
    for i in range(len(a[:,0])-1,-1,-1):
        if np.where(~np.isnan(a[i,:]))[0].size > 0:
            print(str(i),'not empty')
            ymax = i +1
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















