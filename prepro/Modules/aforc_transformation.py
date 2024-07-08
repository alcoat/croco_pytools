#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  6 01:50:45 2024

@author: menkes
"""
import pylab as plt
import numpy as np
import sys
# Extrapolation sst
sys.path.append("/home/menkes/Documents/annelou/gitlab_release/croco_pytools/prepro/Modules")
import scipy.interpolate as itp

# -------------------------------------------------
# FUNCTIONS USED BY make_aforc.py
# -------------------------------------------------


def flip_data(subset_data,data_origin):
    if data_origin == 'era_ecmwf' or data_origin == 'cfsr':
        subset_data0 = subset_data.isel(lat=slice(None,None,-1))
    else: subset_data0 = subset_data.isel(lat=slice(None,None,1))
    return subset_data0

def unit_conversion(data,var,variables): # CFF To change unit
    data1 = data * variables.get_conv_cff(var)
    return data1

def kelvin_2_celsius(data):
    data = data - 273.15
    return data

def extrapolation(Vin,lon,lat):
    longitude, latitude = np.meshgrid(lon,lat)
    Vin[Vin==0]=np.nan
    crocolon = longitude
    Vout = np.zeros([Vin.shape[0],crocolon.shape[0],crocolon.shape[1]])

    for tt in range(Vin.shape[0]):
        igood = np.where(np.isnan(Vin[tt,:])==False)
        ibad  = np.where(np.isnan(Vin[tt,:]))
 
        spline = itp.NearestNDInterpolator((longitude[igood].ravel(),latitude[igood].ravel()),Vin[tt,igood[0],igood[1]])
        Vinfilled =np.copy(np.squeeze(Vin[tt,:]))
        Vinfilled[ibad] = spline(longitude[ibad[0],ibad[1]],latitude[ibad[0],ibad[1]])
        Vout[tt,:] = Vinfilled
    return Vout

def strd_calculation(str_var,sst,variables,croco_variables):
    emiss = 0.985 ; sigma = 5.6697e-8
    strd = str_var + emiss*sigma*sst**4
    strd.name = 'strd'.upper()
    strd.attrs = { 'units': croco_variables['strd'][1], 'long_name': croco_variables['strd'][0] }
    return strd
    
def r_calculation(q,t2m,variables,croco_variables):
    Pref= 1020
    ew = 6.1121*(1.0007+3.46e-6*Pref)* np.exp((17.502*(t2m-273.15))/(240.97+(t2m-273.15)))
    Qsat = 0.62197*(ew/(Pref-0.378*ew))
    r = q / Qsat
    r.name = 'r'.upper()
    r.attrs = { 'units': croco_variables['r'][1], 'long_name': croco_variables['r'][0] } 
    return r
    
def ssr_calculation(uswrf,dswrf,variables,croco_variables):
    ssr = dswrf - uswrf
    ssr.name = 'ssr'.upper()
    ssr.attrs = { 'units': croco_variables['ssr'][1], 'long_name': croco_variables['ssr'][0] }
    return ssr

def attr(data,var,variables,croco_variables):
    data.attrs = { 'units': croco_variables[variables.get_var(var)][1], 'long_name': croco_variables[variables.get_var(var)][0] }
    data.name = var.upper()
    return data

def time_origin(data,Yorig): # Time origin changement
    dt_origins = (np.datetime64(str(Yorig)+'-01-01T00:00:00.000000000') - np.datetime64('1900-01-01T00:00:00.000000000')) 
    dt_old = ( data['time'] - np.datetime64('1900-01-01T00:00:00.000000000')) 
    time_num = ((dt_old - dt_origins) / np.timedelta64(1, 's')) / (24 * 3600) # en jrs, precision sec
    data['time']=time_num
    data['time'].encoding['units'] = 'days since ' + str(Yorig) + '-1-1'
    data['time'].attrs['units'] = 'days since ' + str(Yorig) + '-1-1'
    return data

def metadata(data):
    data['lon'].attrs['long_name'] = 'longitude of RHO-points'
    data['lat'].attrs['long_name'] = 'latitude of RHO-points'
    data['time'].attrs['long_name'] = 'Time'
    data['lon'].attrs['units'] = 'degree_east'
    data['lat'].attrs['units'] = 'degree_north'
    return data

def missing_data(data):
    data = data.fillna(9999.)
    data.encoding['missing_value']=9999.
    data.encoding['_FillValue']=9999.
    return data

