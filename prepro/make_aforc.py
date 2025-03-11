#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 17:51:49 2025

@author: annelou
"""
#--- Dependencies ---------------------------------------------------------
import xarray as xr
import pylab as plt
import numpy as np
import numpy.ma as ma
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
data_origin = 'era_dataref'
input_dir = '/path/in/'
input_prefix = 'era_5-copernicus__*' # For multifiles, if the name of the file begin with the variable name, write '*' before sufix

multi_files = False # If one file per variable in input : True

# -------------------------------------------------
# OUTPUT :
# -------------------------------------------------
output_dir = '/path/out/'
output_file_format = "MONTHLY" # How output files are split (MONTHLY,DAILY)

# -------------------------------------------------
# Grid size : 
ownArea = 0 # 0 if area from croco_grd.nc
            # 1 if own area
if ownArea == 0:
    croco_grd = '/pathin/to/your/croco/grid/croco_grd.nc'
else:
    lon_min,lon_max,lat_min,lat_max = 4,27,-40,-24

# Dates limits
Yorig = 1950                 # year defining the origin of time as: days since Yorig-01-01
Ystart, Mstart = 2005,1   # Starting month
Yend, Mend  = 2005,2  # Ending month

# -------------------------------------------------
# OPTIONS :
# -------------------------------------------------
# If cumul, indicate the step in hours :
cumul_step = 1 # in hour

# To convert the atmospheric pressure : True
READ_PATM = False

# If there is no STRD variable in raw data, it will be calculated with STR and SST. 
# SST may or may not be extrapolated to the coast. 
# If it is not, this may result in temperature spike at the coast at certain points,
# however, note that extrapolation increases pre-processing time.
# If STRD is in raw data, extrapolation_sst will no be considered.
extrapolation_sst = True # /!\ if sst = ts (surf temp) no need
if extrapolation_sst: 
    data_with_mask = 'sst'
    sst_interp = 'nearest' # 'nearest' or 'cloughtocher' (= cubic, not recommended)

# Extrapolation of all atmospheric data to the coast
# Cut data from the atmospheric land/sea mask
# Need to have a masked variable (for example sst data)
drowning = False
if drowning: 
    data_with_mask = 'sst' # Name from variable dictionnary
    drowning_interp = 'nearest' # 'nearest' or 'cloughtocher' (= cubic, not recommended)
    drowning_plot = True # Save a plot to show the interpolation
    if drowning_plot:
        drowning_plot_folder = output_dir
    
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
with open('./Readers/croco_variables.json', 'r') as jf:
    croco_variables = json.load(jf)

# -------------------------------------------------
# Read croco grid to find emprise
# -------------------------------------------------
if ownArea == 0:
    grid = xr.open_dataset(croco_grd)
    lon_min,lon_max,lat_min,lat_max = grid['lon_rho'].min(),grid['lon_rho'].max(),grid['lat_rho'].min(),grid['lat_rho'].max()

# *****************************************************************************
#                                MAIN FUNCTION
# *****************************************************************************
start_time = time.time() 
    
if __name__ == "__main__":
# -----------------------------------
# IF DROWNING, WILL EXTRACT THE MASK
# -----------------------------------
    if extrapolation_sst or drowning:
        pathin = find_type_path(input_dir,Ystart,Mstart)
        if multi_files: input_file = find_input(pathin,input_prefix,Ystart,Mstart,multi_files,READ_PATM,variables.get_var(data_with_mask))
        else: input_file = find_input(pathin,input_prefix,Ystart,Mstart,multi_files,READ_PATM)
        enginein = find_engine(input_file[0])
        land_mask = xr.open_dataset(input_file[0],engine=enginein)[variables.get_var(data_with_mask)].to_masked_array()
        land_mask = ma.getmask(land_mask[0,:,:])

# -----------------------------------
# LOOP ON YEARS AND MONTHS
# -----------------------------------
    for year_inprocess in range(Ystart,Yend+1): 
        if year_inprocess == Ystart: month_inprocess = Mstart
        else: month_inprocess = 1
        
        while month_inprocess <= 12 and (year_inprocess < Yend or (year_inprocess == Yend and month_inprocess <= Mend)):
            
            # Find the path of data (subfolders in years and months or not)
            pathin = find_type_path(input_dir,year_inprocess,month_inprocess)
            
            # Remove index grib files if necessary
            remove_grb_indx(pathin)
                
            # Put start and end date to the right format 
            start_date = plt.datetime.datetime(year_inprocess,month_inprocess,1,0)
            end_date = plt.datetime.datetime(year_inprocess,month_inprocess,1,0) + relativedelta(months=1,hours=-1) # Last day of the ending month

            # Flag not to read files more than once if they are not multi_files
            # If multi_files, will remains True
            flag_frst = True
            
# -----------------------------------
# LOOP ON VARIABLES
# -----------------------------------
            for var in variables.raw_name: 
                if var == 'msl' and not READ_PATM:
                    continue
                elif var == 'lon' or var == 'lat' or var == 'dswrf' or var == 'sst':
                    continue

# -----------------------------------
# READ THE DATA
# -----------------------------------

# If multi_files : 
# ----------------
                if multi_files:
                    # Find file names :
                    input_file = find_input(pathin,input_prefix,year_inprocess,month_inprocess,multi_files,READ_PATM,variables.get_var(var))
                    if input_file==[]:
                        continue
                    else:
                        # Find if netcdf or grib files
                        enginein = find_engine(input_file[0])
                        # Open data
                        dataxr = xr.open_mfdataset(input_file,combine='nested',concat_dim='time', engine=enginein, chunks={'time': 'auto',  variables.get_var('lon'): 'auto',  variables.get_var('lat'): 'auto'})
                        
                    # Some variables need other for calculation :
                    if var == 'str' or var == 'q' or var == 'uswrf':  
                        # Remove index grib files if necessary
                        if enginein == 'cfgrib': remove_grb_indx(pathin)
                        # Find the needeed variable :
                        var2 = [i for i in range(len(['str','q','uswrf'])) if ['str','q','uswrf'][i] == var]
                        var2 = ['sst', 't2m', 'dswrf'][var2[0]]
                        input_file = find_input(pathin,input_prefix,year_inprocess,month_inprocess,multi_files,READ_PATM,variables.get_var(var2))
                        if input_file==[]:
                            print('/!\ ', var2, 'files cannot be read and are necessary for calculations relative to', var, '/!\ ')
                            sys.exit()
                            
                        # Merge data in the same Dataset
                        dataxr = xr.merge([dataxr,xr.open_mfdataset(input_file,combine='nested',concat_dim='time', engine=enginein, chunks={'time': 'auto',  variables.get_var('lon'): 'auto',  variables.get_var('lat'): 'auto'})])    
                        
# If not multi_files and not already read :
# -----------------------------------------
                elif flag_frst==True: 
                    # Find file names :
                    input_file = find_input(pathin,input_prefix,year_inprocess,month_inprocess,multi_files,READ_PATM)
                    # Find if netcdf or grib files
                    enginein = find_engine(input_file[0])
                    # Open data
                    dataxr = xr.open_mfdataset(input_file,combine='nested',concat_dim='time', engine=enginein, chunks={'time': 'auto',  variables.get_var('lon'): 'auto',  variables.get_var('lat'): 'auto'})
                    
                # Remove index grib files if necessary
                if enginein == 'cfgrib': remove_grb_indx(pathin)
                
# -----------------------------------
# IF MULTI_FILES OR FIRST VARIABLE :
# -----------------------------------       
                if flag_frst==True:
                    if not multi_files:
                        flag_frst = False
                        
# Rename dimension lon/lat :
# --------------------------
                    if variables.get_var('lon') != 'lon': 
                        dataxr = dataxr.rename({variables.get_var('lon') : 'lon'})
                        dataxr = dataxr.rename({variables.get_var('lat') : 'lat'})

# Drop dimensions and coordinates that are not in ['time', 'lon', 'lat', 'x', 'y'] :
# ----------------------------------------------------------------------------------
                    # Remove dimension with len = 1 (if any)
                    dataxr = np.squeeze(dataxr)
                    drop_dims = [dim for dim in dataxr[variables.get_var(var)].dims if dim not in ['time', 'lon', 'lat', 'x', 'y']]
                    drop_coords = [coords for coords in dataxr.coords if coords not in ['time', 'lon', 'lat', 'x', 'y']]
                    dataxr = dataxr.drop_vars(drop_coords, errors="ignore")
                    isel_dim = {dim: 0 for dim in drop_dims}
                    indvar2 = [i for i in range(len(['str','q','uswrf'])) if ['str','q','uswrf'][i] == var]
                    dataxr = dataxr.isel(**isel_dim)
                    if drop_dims != []:
                        print("WARNING : the variable has dimension other than ['time', 'lon', 'lat', 'x', 'y'] \n ",\
                              "the index 0 has been taken :",drop_dims)

# ------------------------------------
# SELECT AREA AND TIME :
# ------------------------------------

                    # Check if lon/lat are 1D or 2D :
                    lon_dim = dataxr['lon'].dims
                    lat_dim = dataxr['lat'].dims

# Irregular grid (2D lon/lat) :
# -----------------------------
                    if len(lon_dim) == 2 and len(lat_dim) == 2:
                        irreg = 1
                        # Find lon/lat indices to cut the grid :
                        ix_min,ix_max,iy_latmin,iy_latmax = ind_irreg_grid(dataxr.isel(time=0),variables.get_var(var),lon_min,lon_max,lat_min,lat_max)
                        
                        # Latitudes from min to max :
                        if iy_latmin < iy_latmax: iy_min = iy_latmin  ; iy_max = iy_latmax
                        # Reversed latitudes (max to min) :
                        else: iy_min = iy_latmax  ; iy_max = iy_latmin

                        # Take margin :
                        iy_min -= 4 ; ix_min -= 4
                        if iy_min < 0: iy_min = 0
                        if ix_min < 0: ix_min = 0
                        iy_max += 4 ; ix_max += 4
                        if iy_max // len(dataxr[lon_dim[0]]) == 1: iy_max = len(dataxr[lon_dim[0]])
                        if ix_max // len(dataxr[lon_dim[1]]) == 1: ix_max = len(dataxr[lon_dim[1]])
                        
                        sel_args = {lon_dim[0]: slice(iy_min,iy_max),lon_dim[1]:slice(ix_min,ix_max)}
              
# Regular grid (1D lon/lat) :
# ---------------------------             
                    else: 
                        irreg = 0
                        # Find lon/lat indices to cut the grid :
                        ix_min = int(find_nearest(dataxr['lon'],lon_min))
                        ix_max = int(find_nearest(dataxr['lon'],lon_max))
                        iy_latmin = int(find_nearest(dataxr['lat'],lat_min))
                        iy_latmax = int(find_nearest(dataxr['lat'],lat_max))

                        # Latitudes from min to max :
                        if iy_latmin < iy_latmax: iy_min = iy_latmin  ; iy_max = iy_latmax
                        # Reversed latitudes (max to min) :
                        else: iy_min = iy_latmax  ; iy_max = iy_latmin

                        # Take margin :
                        iy_min -= 4 ; ix_min -= 4
                        if iy_min < 0: iy_min = 0
                        if ix_min < 0: ix_min = 0
                        iy_max += 4 ; ix_max += 4
                        if iy_max // len(dataxr['lat']) == 1: iy_max = len(dataxr['lat'])
                        if ix_max // len(dataxr['lon']) == 1: ix_max = len(dataxr['lon'])

                        # Selection :
                        sel_args = {lat_dim[0]: slice(iy_min,iy_max),lon_dim[0]:slice(ix_min,ix_max)}
                    
                    # Reduced dataset in time, lon and lat :
                    dataxr_inprocess = dataxr.sel(time=slice(start_date,end_date)).isel(**sel_args)
                    
                    # Reduced land_mask in lon and lat :
                    if extrapolation_sst or drowning:
                        land_mask_cut = land_mask[iy_min:iy_max,ix_min:ix_max]

# -----------------------------------------------------
# GROUP BY MONTHS OR DAYS DEPEND ON THE WANTED OUTPUT :
# -----------------------------------------------------
                    
                    if output_file_format == 'MONTHLY':
                        data_grouped = dataxr_inprocess.groupby('time.month')
                        
                    elif output_file_format == 'DAILY':
                        data_grouped = dataxr_inprocess.groupby('time.dayofyear')
            
                    # Creation of a group index :
                    labels = list(data_grouped.groups.keys())
                    index_to_label = {i: label for i, label in enumerate(labels)}

# -----------------------------------
# LOOP ON GROUPS (MONTHS OR DAYS) :
# -----------------------------------
                for ii in range(len(labels)): 
                    i = index_to_label[ii]
        
                    print('\n-----------------------------------')
                    if output_file_format.upper() == "DAILY":
                        print(' Processing Year %s - Month %02i - Day %02i' %(data_grouped[i]['time'].dt.year[0].item(),data_grouped[i]['time'].dt.month[0].item(),data_grouped[i]['time'].dt.day[0].item()))
                    elif output_file_format.upper() == "MONTHLY":
                        print(' Processing Year %s - Month %02i' %(data_grouped[i]['time'].dt.year[0].item(),data_grouped[i]['time'].dt.month[0].item()))
                    print('-----------------------------------')
                    
                    print('Processing variable : var/vname :'+var+'/'+variables.get_var(var))
    
# Check data consistency : 
# ------------------------
                    consistency(data_grouped[i])
                    
# -----------------------------------
# CONVERT DATA IF NECESSARY : 
# -----------------------------------

# Flip data : 
# -----------
                    # If reversed latitudes, will flip them :
                    data = flip_data(data_grouped[i][variables.get_var(var)])

# Remove accumulation period : 
# ----------------------------
                    # Depend if the variable is flag with 'cumul' (cf aforc_reader.py) :
                    if variables.get_iscumul(var) == 'cumul':
                        data = remove_cumul(data,cumul_step)
                        
# Put data in the wanted unit : 
# -----------------------------
                    data = unit_conversion(data,var,variables)
                    
# If the variable is not the final one and calculations are required : 
# --------------------------------------------------------------------
                    # To access to strd : 
                    if var == 'str':
                        sst = flip_data(data_grouped[i][variables.get_var('sst')])
                        sst = unit_conversion(sst,'sst',variables)
                        if extrapolation_sst:
                            sst = extrapolation(sst.values,sst.lon.values,sst.lat.values,land_mask_cut,sst_interp)
                        else:
                            sst = sst.values
                        data = strd_calculation(data,sst,variables,croco_variables)
                        var = 'strd'
                    
                    # To access to r :
                    elif var == 'q':
                        t2m = flip_data(data_grouped[i][variables.get_var('t2m')])
                        t2m = unit_conversion(t2m,'t2m',variables)
                        data = r_calculation(data,t2m,croco_variables)
                        var = 'r'
                        
                    # To access to ssr :
                    elif var == 'uswrf':
                        dswrf = flip_data(data_grouped[i][variables.get_var('dswrf')])
                        if variables.get_iscumul('dswrf') == 'cumul':
                            dswrf = remove_cumul(dswrf,cumul_step)
                        dswrf = unit_conversion(dswrf,'dswrf',variables)
                        data = ssr_calculation(data,dswrf,croco_variables)
                        var = 'ssr'
                        
# If the variable is the final one : 
# ----------------------------------
                    else:
                        data = attr(data,var,variables,croco_variables)
                        
# If drowning, will extrapolate the sea data to the land data :
# -------------------------------------------------------------
                    if drowning:
                        data_extra = extrapolation(data.values,data.lon.values,data.lat.values,land_mask_cut,drowning_interp)
                        
                        if drowning_plot and year_inprocess == Ystart and month_inprocess == Mstart:
                            extrapolation_plot(data,data_extra,land_mask_cut,drowning_plot_folder)
                            
                        data = xr.DataArray(data_extra, dims = data.dims, coords = data.coords, name = data.name, attrs = data.attrs) 
                        
                        
# Put time in the right format :
# ------------------------------
                    data = time_origin(data,Yorig)
                    
# Put the good metadata :
# -----------------------
                    data = metadata(data)
                    
# Put the good encoding : 
# -----------------------
                    data,encoding = missing_data(data,var)
                    
# -----------------------------------
# CREATE THE NETCDF :
# -----------------------------------
                    data = create_netcdf(data,output_dir,output_file_format,encoding)

# Will go on next month :
# -----------------------
            month_inprocess += 1
            
end_time = time.time()
time_taken = end_time - start_time
print("Computation time:", time_taken, "sec")

