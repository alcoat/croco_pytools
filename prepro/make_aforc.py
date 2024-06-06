#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 10:43:38 2024

@author: annelou
"""
#--- Dependencies ---------------------------------------------------------
import xarray as xr
from datetime import date
import dask
import pylab as plt
import numpy as np
import glob as glob
from dateutil.relativedelta import relativedelta
import json
import time
import os
# Reader
import sys
sys.path.append("./Readers/")
from aforc_reader import lookvar
from aforc_class import aforc_class, create_class

# *******************************************************************************
#                         U S E R  *  O P T I O N S
# *******************************************************************************

#############
# INPUT :
#############
data_origin = 'era_dataref' # era_dataref, era_ecmwf
input_dir = '/path/to/data/'
input_prefix = 'prefix__*'  # Please use * to include all files
multi_files = False # If one file per variable in input : True

#############
# OUTPUT :
#############
output_dir = '/output/path/'
bry_filename = 'croco_aforc.nc'
output_file_format = "MONTHLY" # How outputs are split (MONTHLY,DAILY)
output_file_dt = "HOURLY" # Time step of the output : Daily or hourly

# Grid size : Cut
lon_grd = [90,289.95] # min, max
lat_grd = [-30,30] # min, max

# Dates limits
Yorig = 1950                    # year defining the origin of time as: days since Yorig-01-01
Ystart, Mstart = 1980, 3   # Starting month
Yend, Mend  = 1980,12      # Ending month

#############
# OPTIONS :
#############
# To convert the atmospheric pressure : True
READ_PATM = False

# To add wave variables : True
wave_extract = False

# *******************************************************************************
#                         E N D     U S E R  *  O P T I O N S
# *******************************************************************************

# -------------------------------------------------
# Read variables
# -------------------------------------------------
variables = create_class(lookvar(data_origin))

# -------------------------------------------------
# Add waves variables if True
# -------------------------------------------------
if wave_extract:
    wave_var=['swh', 'mwd', 'pp1d' ,'cdww'];variables.extend(wave_var)
    wave_conv_cff=[1.,  1., 1. , 1.];  conv_cff.extend(wave_conv_cff)
    wave_units=['m','Degrees true','s', 'dimensionless']; units.extend(wave_units)

# -------------------------------------------------
# Setting processed output directory
# -------------------------------------------------
# Get the current directory
os.makedirs(output_dir,exist_ok=True)

# -------------------------------------------------
# Loading ERA5 variables's information as 
# python Dictionary from JSON file
# -------------------------------------------------
with open('ERA5_variables.json', 'r') as jf:
    era5 = json.load(jf)

# -------------------------------------------------
# Test que pour multi_files = False
# Attention que pour des dossiers du type /data/Y/M/file.nc
# -------------------------------------------------
if multi_files: etan=[], theta=[], salt=[], evel=[], nvel=[]
else: input_file = []
for j in range(Ystart, Yend + 1):
    if j == Ystart: i = Mstart
    else: i = 1
    while i <= 12 and (j < Yend or (j == Yend and i <= Mend)):
        if multi_files:
            etan.extend(glob.glob(input_dir + str(j) + '/' + str(i).zfill(2) + '/' + input_prefix+'ETAN.*.nc'))
            temp.extend(glob.glob(input_dir + str(j) + '/' + str(i).zfill(2) + '/' + input_prefix+'THETA.*.nc'))           
            salt.extend(glob.glob(input_dir + str(j) + '/' + str(i).zfill(2) + '/' + input_prefix+'SALT.*.nc'))
            evel.extend(glob.glob(input_dir + str(j) + '/' + str(i).zfill(2) + '/' + input_prefix+'EVEL.*.nc')) 
            nvel.extend(glob.glob(input_dir + str(j) + '/' + str(i).zfill(2) + '/' + input_prefix+'NVEL.*.nc')) 
        else:
            input_file.extend(glob.glob(input_dir + str(j) + '/' + str(i).zfill(2) + '/' + input_prefix))
        i += 1

if multi_files:
    input_file = {'ssh':sorted(etan),\
                  'temp':sorted(temp),\
                  'salt':sorted(salt),\
                  'u':sorted(evel),\
                  'v':sorted(nvel)\
                }
else:    
    input_file = sorted(input_file)
    

#--------------------------------------------------------------------------------------
print('Go into __main__')
if __name__ == '__main__':
    
    start_time_dask = time.time()

    # Put start and end date to the right format
    if output_file_dt.upper() == "DAILY":
        start_date = str(Ystart)+str(Mstart).zfill(2)+'01'+'12'  # defaut start day is 1st
        dtenddt = plt.datetime.datetime(Yend,Mend,1,12) \
            + relativedelta(months=1,days=-1) # Last day of the ending month    

    elif output_file_dt.upper() == "HOURLY":
        start_date = str(Ystart)+str(Mstart).zfill(2)+'01'+'01'  # defaut start day is 1st
        dtenddt = plt.datetime.datetime(Yend,Mend,1,0) \
            + relativedelta(months=1,hours=-1) # Last day of the ending month
    
    startloc=plt.datetime.datetime(int(start_date[:4]),
                                    int(start_date[4:6]),
                                    int(start_date[6:8]),
                                    int(start_date[8:]))
    print(input_file)
    # --- Initialize input data class -------------------------------------
    if data_origin == 'era_dataref':
        dataxr = xr.open_mfdataset(input_file,combine='nested',concat_dim='time',\
                             engine='netcdf4', chunks={'time': 'auto', 'latitude025': 'auto', 'longitude025': 'auto'})
        # Decoupe selon start/end + emprise de la grille + marge
        dataxr = dataxr.sel( time = slice(startloc,dtenddt) , \
                            longitude025 = slice(lon_grd[0]-5,lon_grd[1]+5) ,\
                                latitude025 = slice(lat_grd[0]-5,lat_grd[1]+5) )
    else:
        dataxr = xr.open_mfdataset(input_file,combine='nested',concat_dim='time',\
                                 engine='netcdf4', chunks={'time': 'auto', 'latitude': 'auto', 'longitude': 'auto'})
        # Decoupe selon start/end + emprise de la grille + marge
        dataxr = dataxr.sel( time = slice(startloc,dtenddt) , \
                            longitude = slice(lon_grd[0]-5,lon_grd[1]+5) ,\
                                latitude = slice(lat_grd[1]+5,lat_grd[0]-5) )


    if not 'latitude' in dataxr.dims and 'latitude025' in dataxr.dims:
        dataxr = dataxr.rename({'latitude025' : 'latitude'})
    if not 'longitude' in dataxr.dims and 'longitude025' in dataxr.dims:
        dataxr = dataxr.rename({'longitude025' : 'longitude'})

 
    if output_file_dt.upper() == "DAILY":
        dataxr = dataxr.sel(time=dataxr['time.hour'] == 12)
    
    print('Open data with mfdataset check')

    def formatting(subset_data):

        if output_file_format.upper() == "MONTHLY":
            bdy_filename = output_dir+bry_filename.replace('.nc', '_%s_Y%sM%02i.nc' %(data_origin,subset_data['time'].dt.year[0].item(),subset_data['time'].dt.month[0].item()))
        elif output_file_format.upper() == "DAILY":
            bdy_filename = output_dir+bry_filename.replace('.nc', '_%s_Y%sM%02iD%02i.nc' %(data_origin,subset_data['time'].dt.year[0].item(),subset_data['time'].dt.month[0].item(),subset_data['time'].dt.day[0].item()))

        
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
         
            
            
        for var in variables.raw_name:
            vlong = era5[variables.get_var(var)][0]
            print('  Processing variable: var/vname :'+var+'/'+variables.get_var(var))
            
            if var == 'msl' and not READ_PATM:
                continue 
            elif var == 'dswrf' or var == 'sst':
                continue
            
            
            
            ## flip
            if data_origin == 'era_dataref':
                subset_data0 = subset_data.isel(latitude=slice(None,None,-1))
            else: subset_data0 = subset_data.isel(latitude=slice(None,None,1))

            ## IMPORT DATA
            data = subset_data0[variables.get_var(var)]
            # data = subset_data[var]

            ## MISSING VALUES
            mvalue = data.encoding.get('_FillValue',np.nan)
            
            # Indices of fillvalue, for after calculations
            data = data.where(data != mvalue)
           
            # CFF To change unit
            data = data * variables.get_conv_cff(var)
            
            # Calcul de nouvelles variables si besoin :
            if var=='str':
                sst = subset_data0[variables.get_var('sst')]
                emiss = 0.985 ; sigma = 5.6697e-8
                data = data + emiss*sigma*sst**4
                data.name = 'strd'.upper()
                data.attrs = { 'units': variables.get_units(var), 'long_name': 'surface_net_thermal_radiation' }
            elif var == 'q':
                t2m = subset_data0[variables.get_var('t2m')]
                Pref= 1020
                ew = 6.1121*(1.0007+3.46e-6*Pref)* np.exp((17.502*(t2m-273.15))/(240.97+(t2m-273.15)))
                Qsat = 0.62197*(ew/(Pref-0.378*ew))
                data = data / Qsat
                data.name = 'r'.upper()
                data.attrs = { 'units': variables.get_units(var), 'long_name': 'relative_humidity' } 
            elif var == 'uswrf':
                dsw = subset_data0[variables.get_var('dswrf')]
                data = dsw - data
                data.name = 'ssr'.upper()
                data.attrs = { 'units': variables.get_units(var), 'long_name': 'Net Short Wave' }
            else:
                data.attrs = { 'units': variables.get_units(var), 'long_name': vlong }
                data.name = var.upper()
            
            # Changement de l'origine du temps
            dt_origins = (np.datetime64(str(Yorig)+'-01-01T00:00:00.000000000') - np.datetime64('1900-01-01T00:00:00.000000000')) 
            dt_old = ( data['time'] - np.datetime64('1900-01-01T00:00:00.000000000')) 
            time_num = ((dt_old - dt_origins) / np.timedelta64(1, 's')) / (24 * 3600) # en jrs, precision sec
            data['time']=time_num
            data['time'].encoding['units'] = 'days since ' + str(Yorig) + '-1-1'
            data['time'].attrs['units'] = 'days since ' + str(Yorig) + '-1-1'
    
            data['longitude'].attrs['long_name'] = 'longitude of RHO-points'
            data['latitude'].attrs['long_name'] = 'latitude of RHO-points'
            data['time'].attrs['long_name'] = 'Time'
            data['longitude'].attrs['units'] = 'degree_east'
            data['latitude'].attrs['units'] = 'degree_north'
    
            data.compute()
            print('Computed....')
    
            data = data.fillna(9999.)
            
            data.encoding['missing_value']=9999.
            data.encoding['_FillValue']=9999.
            
            data.to_netcdf(bdy_filename[:-3]+'_'+data.name+'.nc')
            print('Created.....')
            
        return subset_data
            
            
    if output_file_format == 'MONTHLY':
        dataxr.groupby('time.month').map(formatting)
        

    elif output_file_format == 'DAILY':
        dataxr.groupby('time.dayofyear').map(formatting)


    end_time_dask = time.time()
    time_taken = end_time_dask - start_time_dask
    print("Temps de calcul:", time_taken, "secondes")





















