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
import scipy.interpolate as itp
import pyinterp.backends.xarray

# *******************************************************************************
#                         U S E R  *  O P T I O N S
# *******************************************************************************

#############
# INPUT :
#############
data_origin = 'era_dataref' # era_dataref, era_ecmwf, cfsr
input_dir = '/path/to/data/'
# For multifiles, if the name of the file begin with the variable name, just write '*'
input_prefix = 'era_5-copernicus__*' # Please use * to include all files
multi_files = False # If one file per variable in input : True

#############
# OUTPUT :
#############
output_dir = '/path/to/output/data/'
bry_filename = 'croco_aforc.nc'
output_file_format = "DAILY" # How outputs are split (MONTHLY,DAILY)

# Grid size : 
ownArea = 1 # 0 if area from croco_grd.nc -> +/- 5°
            # 1 if own area
if ownArea == 0:
    croco_grd = '/../../CONFIGS/your_config/CROCO_FILES/croco_grd.nc'
else:
    lon_min,lon_max,lat_min,lat_max = 6,24,-40,-24 # in degrees

# Dates limits
Yorig = 2000                 # year defining the origin of time as: days since Yorig-01-01
Ystart, Mstart = 2005,1   # Starting month
Yend, Mend  = 2005,1    # Ending month

#############
# OPTIONS :
#############
# To convert the atmospheric pressure : True
READ_PATM = False

# *******************************************************************************
#                         E N D     U S E R  *  O P T I O N S
# *******************************************************************************

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
# Loading ERA5 variables's information as 
# python Dictionary from JSON file
# -------------------------------------------------
# with open('ERA5_variables.json', 'r') as jf:
with open('/Readers/Aforc_variables.json', 'r') as jf:
    varlong = json.load(jf)

# -------------------------------------------------
# Test que pour multi_files = False
# Warning : only for path type /data_origin/file.nc or /data_origin/Y/M/file.nc
# -------------------------------------------------
# Construction of input_file
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
                else:
                    if type_path == '/origin/Y/M/':
                        eval(var).extend(glob.glob(input_dir + str(j) + '/' + str(i).zfill(2) + '/'+input_prefix +variables.get_filename(var)+'*.nc'))
                        input_file.extend(glob.glob(input_dir + str(j) + '/' + str(i).zfill(2) + '/'+input_prefix +variables.get_filename(var)+'*.nc'))
                    elif type_path == '/origin/':
                        eval(var).extend(sorted(glob.glob(input_dir + input_prefix +variables.get_filename(var)+'*.nc')))
                        input_file.extend(sorted(glob.glob(input_dir + input_prefix +variables.get_filename(var)+'*.nc')))
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
print(input_file)
###############################################################################
# -------------------------------------------------
# Read croco grid to find emprise
# -------------------------------------------------
if ownArea == 0:
    grid = xr.open_dataset(croco_grd)
    lon_min,lon_max,lat_min,lat_max = np.floor(grid['lon_rho'].min())-5,np.ceil(grid['lon_rho'].max())+5,np.floor(grid['lat_rho'].min())-5,np.ceil(grid['lat_rho'].max())+5

# -------------------------------------------------
# Functions
# -------------------------------------------------
def output_name(subset_data):
    if output_file_format.upper() == "MONTHLY":
        bdy_filename = output_dir+bry_filename.replace('croco_aforc.nc', 'Y%sM%02i.nc' %(subset_data['time'].dt.year[0].item(),subset_data['time'].dt.month[0].item()))
    elif output_file_format.upper() == "DAILY":
        bdy_filename = output_dir+bry_filename.replace('croco_aforc.nc', 'Y%sM%02iD%02i.nc' %(subset_data['time'].dt.year[0].item(),subset_data['time'].dt.month[0].item(),subset_data['time'].dt.day[0].item()))
    
    print('\n-----------------------------------')
    if output_file_format.upper() == "DAILY":
        print(' Processing Year %s - Month %02i - Day %02i' %(subset_data['time'].dt.year[0].item(),subset_data['time'].dt.month[0].item(),subset_data['time'].dt.day[0].item()))
    elif output_file_format.upper() == "MONTHLY":
        print(' Processing Year %s - Month %02i' %(subset_data['time'].dt.year[0].item(),subset_data['time'].dt.month[0].item()))
    print('-----------------------------------')
    
    if subset_data.dims['time']==0 :
        print('\nData is missing - No Data for the subset')
        sys.exit()
    
    if np.nanvar(np.gradient(plt.date2num(subset_data['time'].values))) >=5: # Abnormal distribution of days
        Question = input( "Abnormal distribution of days (variance to high) \
                \nThis may be due to the use of different temproral resolution dataset.\
                \n Do you want to proceed?: y,[n] ") or 'no'
        if Question.lower() == ("n") or Question.lower() == ("no"):
            print('Aborting')
            sys.exit()
     
    return bdy_filename

def flip_data(subset_data):
    if data_origin == 'era_ecmwf' or data_origin == 'cfsr':
        subset_data0 = subset_data.isel(lat=slice(None,None,-1))
    else: subset_data0 = subset_data.isel(lat=slice(None,None,1))
    return subset_data0

def unit_conversion(data,var): # CFF To change unit
    data1 = data * variables.get_conv_cff(var)
    if var == 'tp':
        data1 = data1 * 8640
    elif variables.get_units(var) == '%':
        data1 = data1 / 100
    return data1

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

def strd_calculation(str_var,sst):
    emiss = 0.985 ; sigma = 5.6697e-8
    strd = str_var + emiss*sigma*sst**4
    strd.name = 'strd'.upper()
    strd.attrs = { 'units': variables.get_units('str'), 'long_name': 'surface_net_thermal_radiation' }
    return strd
    
def r_calculation(q,t2m):
    Pref= 1020
    ew = 6.1121*(1.0007+3.46e-6*Pref)* np.exp((17.502*(t2m-273.15))/(240.97+(t2m-273.15)))
    Qsat = 0.62197*(ew/(Pref-0.378*ew))
    r = q / Qsat
    r.name = 'r'.upper()
    r.attrs = { 'units': '%', 'long_name': 'relative_humidity' } 
    return r
    
def ssr_calculation(uswrf,dswrf):
    ssr = dswrf - uswrf
    ssr.name = 'ssr'.upper()
    ssr.attrs = { 'units': variables.get_units('uswrf'), 'long_name': 'Net Short Wave' }
    return ssr

def attr(data,var,vlong):
    if var == 't2m':
        data.attrs = { 'units': 'C', 'long_name': vlong }
    elif var == 'tp':
        data.attrs = { 'units': 'cm/day', 'long_name': vlong }
    elif variables.get_units(var) == '%':
        data.attrs = { 'units': '(0-1)', 'long_name': vlong }
    else:
        data.attrs = { 'units': variables.get_units(var), 'long_name': vlong }
    data.name = var.upper()
    return data

def time_origin(data): # Time origin changement
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

def create_netcdf(data,bdy_filename):
    if output_file_format == 'MONTHLY':
        filename_out = bdy_filename[:-11]+data.name+'_'+bdy_filename[-11:]
    elif output_file_format == 'DAILY':
        filename_out = bdy_filename[:-14]+data.name+'_'+bdy_filename[-14:]
    data = data.astype(np.float32)
    data = data.to_netcdf(filename_out,engine='netcdf4',compute='False')
    return data    

# -------------------------------------------------
# Main function
# -------------------------------------------------
def formatage(input_file):   # if multi files, 'input_file' not used 

    if multi_files:
        ind = 0
        for var in variables.raw_name:
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
    else:
        dataxr = xr.open_mfdataset(input_file,combine='nested',concat_dim='time', engine='netcdf4', chunks={'time': 'auto',  variables.get_var('lon'): 'auto',  variables.get_var('lat'): 'auto'})
   
    # Rename lon/lat
    if variables.get_var('lon') != 'lon':
        dataxr = dataxr.rename({variables.get_var('lon') : 'lon'})
        dataxr = dataxr.rename({variables.get_var('lat') : 'lat'})
    
    for year_inprocess in range(Ystart,Yend+1): # loop on years
        # -------------------------------------------------
        # Put start and end date to the right format
        # -------------------------------------------------
        if year_inprocess == Ystart:
            start_date = plt.datetime.datetime(year_inprocess,Mstart,1,0)
        else:
            start_date = plt.datetime.datetime(year_inprocess,1,1,0)
        
        if year_inprocess == Yend:
            end_date = plt.datetime.datetime(year_inprocess,Mend,1,0) + relativedelta(months=1,hours=-1) # Last day of the ending month
        else:
            end_date = plt.datetime.datetime(year_inprocess,12,1,0) + relativedelta(months=1,hours=-1) # Last day of the ending month
    
        # Cut of the grid
        if data_origin == 'era_ecmwf' or data_origin == 'cfsr': # Inverse latitudes
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
        
        for ii in range(len(labels)): # loop on groups (months or days)
            i = index_to_label[ii]
            bdy_filename = output_name(data_grouped[i])
            for var in variables.raw_name: # loop on var
                if var == 'msl' and not READ_PATM:
                    continue
                elif var == 'lon' or var == 'lat' or var == 'dswrf' or var == 'sst':
                    continue
                print('  Processing variable: var/vname :'+var+'/'+variables.get_var(var))
                vlong = varlong[variables.get_var(var)][1]
                
                if vlong != variables.get_units(var):
                    if var == 'r' and variables.get_units(var) == '%': 
                        print('Relative humidity will be convert in ratio (0-1)')
                    else:
                        print("Please be sure that the units are correct : \
                          \n You must not change the units in aforc_reader.py but give the conversion to access it (from your data unit)")
                        sys.exit()
                
                data = flip_data(data_grouped[i][variables.get_var(var)])
                data = unit_conversion(data,var)
                if var == 't2m':
                    data = data - 273.15
                if var == 'str':
                    sst = extrapolation(data_grouped[i][variables.get_var('sst')].values,data_grouped[i].lon.values,data_grouped[i].lat.values)
                    data = strd_calculation(data,sst)
                elif var == 'q':
                    data = r_calculation(data,data_grouped[i][variables.get_var('t2m')])
                elif var == 'uswrf':
                    data = ssr_calculation(data,data_grouped[i][variables.get_var('dswrf')])
                else:
                    data = attr(data,var,vlong)
                data = time_origin(data)
                data = metadata(data)
                data = missing_data(data)
                data = create_netcdf(data,bdy_filename)

start_time = time.time()
res = formatage(input_file)
end_time = time.time()
time_taken = end_time - start_time
print("Computation time:", time_taken, "sec")













