#!/usr/bin/env python
#
# ERA5_convert.py
#
# Script to convert ERA5 store files in Datarmor /dataref/ref19/intranet/ERA5/
# a format and using units which can be used by the online interpolation of CROCO
#
#
# -------------------------------------------------
# Getting libraries and utilities
# -------------------------------------------------

from datetime import date
import json
from netCDF4 import Dataset as netcdf
import numpy as np
import os
from calendar import monthrange

# -------------------------------------------------
# Import my crocotools_param_python file
from era5_crocotools_param import *
# -------------------------------------------------

# -------------------------------------------------
# Setting processed output directory
# -------------------------------------------------
# Get the current directory
os.makedirs(era5_dir_processed,exist_ok=True)

# -------------------------------------------------
# Loading ERA5 variables's information as 
# python Dictionary from JSON file
# -------------------------------------------------

with open('ERA5_variables.json', 'r') as jf:
    era5 = json.load(jf)

#
# -------------------------------------------------
# Loop on Years and Months
# -------------------------------------------------
# 

for iyear in range(year_start,year_end+1):
    for imonth in range(month_start,month_end+1):
#
# -------------------------------------------------
# Loop on variables names
# -------------------------------------------------
# 
        for var in variables.raw_name:
#
# Variable's name, long-name and level-type
#
            vname = variables.get_var(var)
            vlong = era5[vname][0]

            print('  Processing variable: '+var)

            if var=='sst' or var=='tp' or var=='ssr' or var=='t2m' or var=='u10m' or var=='v10m' or var=='strd':
                continue

            if not READ_PATM and var == 'msl':
                print(' MSL is not formatted because READ_PATM == FALSE')
                continue

            for iday in range(1,monthrange(iyear,imonth)[1]+1):
#
# Read input filedate.toordinal(date(Yorig,1,1))
#
                if data_origin == 'era_dataref':
                    fname_in = era5_dir_raw + '/' + str(iyear) + '/' + str(imonth).rjust(2,"0") +  '/era_5-copernicus__' + str(iyear) + str(imonth).rjust(2,"0") + str(iday).rjust(2,"0") + '.nc'
                    latitude = 'latitude025' ; longitude = 'longitude025'
                    onefile = True # all variables in on file

                elif data_origin == 'era_ecmwf':
                    fname_in = era5_dir_raw + '/ERA5_ecmwf_' + vname.upper() + '_Y' + str(iyear) + 'M' + str(imonth).zfill(2) + '.nc'
                    latitude = 'latitude' ; longitude = 'longitude'
                    onefile = False # one file per variable
     
                else:
                    print('ERROR : Data origin unknown')
                    sys.exit()

                # Lecture du fichier avec netcdf4
                nc = netcdf(fname_in,'r',format='NETCDF4')
                time_day = nc.variables['time'][:]
                lat = nc.variables[latitude][:]
                lon = nc.variables[longitude][:]
                data_day = nc.variables[vname][:,:,:]

                # Lecture du fichier avec xarray
                #nc = xr.open_dataset(fname_in)
                #time_day = nc['time'].values
                #lat = nc[latitude].values
                #lon = nc[longitude].values
                #data_day = nc[vname].values

                # LOAD Data for variables needed for variable transformation/calculation
                if var=='str':
                    if onefile: sst_day = nc.variables[variables.get_var('sst')][:,:,:]
                    else:
                        fname_sst = era5_dir_raw + '/ERA5_ecmwf_SST_Y' + str(iyear) + 'M' + str(imonth).zfill(2) + '.nc'
                        nc_sst = netcdf(fname_sst,'r',format='NETCDF4')
                        sst_day = nc_sst.variables[variables.get_var('sst')][:,:,:]

                    if iday==1: sst = sst_day
                    else: sst = np.concatenate((sst,sst_day),axis=0)

                if var=='q':
                    if onefile: t2m_day = nc.variables[variables.get_var('t2m')][:,:,:]
                    else: 
                        fname_t2m = era5_dir_raw + '/ERA5_ecmwf_T2M_Y' + str(iyear) + 'M' + str(imonth).zfill(2) + '.nc'
                        nc_t2m = netcdf(fname_t2m,'r',format='NETCDF4')
                        t2m_day = nc_t2m.variables[variables.get_var('t2m')][:,:,:]

                    if iday==1: t2m = t2m_day
                    else: t2m = np.concatenate((t2m,t2m_day),axis=0)

                if iday==1: data=data_day ; time=time_day
                else: data = np.concatenate((data,data_day),axis=0) ; time=np.concatenate((time,time_day))
                
                nc.close()
    
                if data_origin == 'era_ecmwf': # Input = monthly file, no loop on days
                    break # Will be change when I will put the time reading into the script
#
# Flip latitudes (to have increasing latitudes...)
# Only for ERA5 ?

            lat = np.flip(lat, axis=0)
            data = np.flip(data, axis=1)

#
# Missing values and multiply by cff to change unit
#

            try:
                mvalue=data.fill_value
            except AttributeError:
                print ('No fill value.. use nan')
                mvalue=np.nan

            data=np.array(data)
            data=variables.get_conv_cff(var) * data
            fillvalue_ind = np.where(data==mvalue)
            name = var # For the netcdf creation

#
# Calculate the downward long-wave radiation (strd) from the net long-wave radiation (str) and the sst.
#
            if var=='str':
                emiss = 0.985 ; sigma = 5.6697e-8
                data = data + emiss*sigma*sst**4
                var = 'strd'
                vlong = 'surface_net_thermal_radiation' 

#
# Calculate the relative humidity from the specific humidity
#
            if var == 'q':
                Pref= 1020
                ew = 6.1121*(1.0007+3.46e-6*Pref)* np.exp((17.502*(t2m-273.15))/(240.97+(t2m-273.15)))
                Qsat = 0.62197*(ew/(Pref-0.378*ew))
                data = data / Qsat
                var = 'r'
                vlong = 'relative_humidity'
#
# Missing values after calculations
#
            data[fillvalue_ind]=9999.

#
# Convert time from hours since 1900-1-1 0:0:0 into days since Yorig-1-1 0:0:0
#
            time = time / 24.
            time = time - date.toordinal(date(Yorig,1,1)) \
	            + date.toordinal(date(1900,1,1))

#
# Create and write output netcdf file
#
            fname_out = era5_dir_processed + '/' + var.upper() + '_Y' + str(iyear) + 'M' + str(imonth).rjust(2,"0") + '.nc'

            nw = netcdf(fname_out,mode='w',format='NETCDF4')

            dimlon  = nw.createDimension('lon',  len(lon))
            dimlat  = nw.createDimension('lat',  len(lat))
            dimtime = nw.createDimension('time', None)

            varlon = nw.createVariable('lon', 'f4',('lon',))
            varlat = nw.createVariable('lat', 'f4',('lat',))
            vartime = nw.createVariable('time', 'f4',('time',))
            vardata = nw.createVariable(var.upper(), 'f4',('time','lat','lon'))
            varlon.long_name = 'longitude of RHO-points'
            varlat.long_name = 'latitude of RHO-points'
            vartime.long_name = 'Time'
            varlon.units = 'degree_east'
            varlat.units = 'degree_north'
            vartime.units = 'days since '+str(Yorig)+'-1-1'
            vardata.missing_value = 9999.
            print(var,name)
            vardata.units = variables.get_units(name) #Same as var
            vardata.long_name = vlong
	
            varlon[:]=lon
            varlat[:]=lat
            vartime[:]=time
            vardata[:]=data
	
            nw.close()

# Print last message on screen
print(' ')
print(' ERA5 files conversion done ')
print(' ')




