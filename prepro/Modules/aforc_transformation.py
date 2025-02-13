#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 14:52:10 2025

@author: annelou
"""

import pylab as plt
import numpy as np
import sys
sys.path.append("/Modules")
import scipy.interpolate as itp
import time
import datetime as dt
import xarray as xr

# -------------------------------------------------
# FUNCTIONS USED BY make_aforc.py
# -------------------------------------------------
def flip_data(data):
# The latitude of some data sources are not increasing, so the script will flip them
# ---------------------
# data : xarray.DataArray of one variable
    if len(data['lat'].dims) == 1: # Regular grid
        if data['lat'][1].values < data['lat'][0].values:
            data = data.isel(lat=slice(None,None,-1))
        else: data = data.isel(lat=slice(None,None,1))
    return data

def unit_conversion(data,var,variables):
# The script will multiply data by a factor to change the unit
# (see ../Readers/croco_variables.json to find the wanted unit)
# ---------------------
# data : xarray.DataArray of one variable
# var : name of the variable
# variables : class of variables with the name of the variable, \
# his factor of conversion to attempt the needed unity and the nomenclature of the variable file if multifile
    if variables.get_conv_cff(var) == 273.15:
        data1 = data - variables.get_conv_cff(var)
    else:
        data1 = data * variables.get_conv_cff(var)
    return data1

def extrapolation(data,lon,lat):
# Extrapolate data with the nearest neighbor method
# ---------------------
# data : xarray.DataArray of one variable
# lon : longitude 1D
# lat : latitude 1D
    longitude, latitude = np.meshgrid(lon,lat)
    data[data==0]=np.nan
    crocolon = longitude
    Vout = np.zeros([data.shape[0],crocolon.shape[0],crocolon.shape[1]])

    for tt in range(data.shape[0]):
        igood = np.where(np.isnan(data[tt,:])==False)
        ibad  = np.where(np.isnan(data[tt,:]))
 
        spline = itp.NearestNDInterpolator((longitude[igood].ravel(),latitude[igood].ravel()),data[tt,igood[0],igood[1]])
        Vinfilled =np.copy(np.squeeze(data[tt,:]))
        Vinfilled[ibad] = spline(longitude[ibad[0],ibad[1]],latitude[ibad[0],ibad[1]])
        Vout[tt,:] = Vinfilled
    return Vout

def strd_calculation(str_var,sst,variables,croco_variables):
# Calculation of the surface thermal radiation downwards
# ---------------------
# str_var : surface net thermal radiation data xarray
# sst : sea surface temperature data xarray
# variables : class of variables with the name of the variable, \
# his factor of conversion to attempt the needed unity and the nomenclature of the variable file if multifile
# croco_variables : dictionary with long name and unity of the variable \
# (see ../Readers/croco_variables.json)
    emiss = 0.985 ; sigma = 5.6697e-8
    strd = str_var + emiss*sigma*sst**4
    strd.name = 'strd'.upper()
    strd.attrs = { 'units': croco_variables['strd'][1], 'long_name': croco_variables['strd'][0] }
    return strd
    
def r_calculation(q,t2m,croco_variables):
# Calculation of the relative humidity (0-1)
# ---------------------
# q : specific humidity
# t2m : temperature at 2m in Celsius
# croco_variables : dictionary with long name and unity of the variable \
# (see ../Readers/croco_variables.json)
    Pref= 1020
    ew = 6.1121*(1.0007+3.46e-6*Pref)* np.exp((17.502*(t2m))/(240.97+(t2m)))
    Qsat = 0.62197*(ew/(Pref-0.378*ew))
    r = q / Qsat
    r.name = 'r'.upper()
    r.attrs = { 'units': croco_variables['r'][1], 'long_name': croco_variables['r'][0] } 
    return r
    
def ssr_calculation(uswrf,dswrf,croco_variables):
# Calculation of surface net solar radiation
# ---------------------
# uswrf : upward long wave rad flux surface
# dswrf : downward long wave rad flux surface
# croco_variables : dictionary with long name and unity of the variable \
# (see ../Readers/croco_variables.json)
    ssr = dswrf - uswrf
    ssr.name = 'ssr'.upper()
    ssr.attrs = { 'units': croco_variables['ssr'][1], 'long_name': croco_variables['ssr'][0] }
    return ssr

def remove_cumul(var_cumul,cumul_step):
# Remove accumulation for data with 'cumul' in the reader
# ---------------------
# var_cumul : accumulated data
# cumul_step : accumulation period in hours
    # Find indices where data is not accumulated :
    c_step = int(dt.timedelta(hours=cumul_step).total_seconds())
    mask_resetcumul = (var_cumul.time - var_cumul.time[0]) % np.timedelta64(c_step,'s')
    mask_resetcumul = mask_resetcumul == np.timedelta64(0)
    # Compute var[t=i] - var[t=i-1] :
    var_dif = var_cumul.diff(dim='time',n=1)
    
    # Create a frame of zeros to add to var_dif and represent t=0 :
    if len(var_cumul.lon.dims) == 2: # 2D lon/lat
        zero_frame = xr.DataArray(
                np.zeros((1, var_cumul.sizes[var_cumul.lon.dims[0]], var_cumul.sizes[var_cumul.lon.dims[1]])),  
                dims=['time', var_cumul.lon.dims[0], var_cumul.lon.dims[1]], 
                coords={'time': [var_cumul.time.values[0]], 
                      var_cumul.lon.dims[0] : var_dif[var_cumul.lon.dims[0]], var_cumul.lon.dims[1]: var_dif[var_cumul.lon.dims[1]]}  
            )
    else: # 1D lon/lat
        zero_frame = xr.DataArray(
            np.zeros((1, var_cumul.sizes[var_cumul.lat.dims[0]], var_cumul.sizes[var_cumul.lon.dims[0]])),  
            dims=['time', var_cumul.lat.dims[0], var_cumul.lon.dims[0]], 
            coords={'time': [var_cumul.time.values[0]], 
                    var_cumul.lat.dims[0] : var_dif[var_cumul.lat.dims[0]], var_cumul.lon.dims[0]: var_dif[var_cumul.lon.dims[0]]}  
        )
    var_dif = xr.concat([zero_frame, var_dif], dim='time')
    
    # For indices where data is not accumulated : var_cumul
    # For others                                : var_dif
    var_tmp = mask_resetcumul * var_cumul + (1 - mask_resetcumul) * var_dif

    # To have [var] = unit per hour :
    step_time = var_cumul[var_cumul.dims[0]][1] - var_cumul[var_cumul.dims[0]][0]
    step_time = float(step_time / np.timedelta64(1, 'h'))
    var_nocumul = var_tmp / step_time
    
    return var_nocumul








