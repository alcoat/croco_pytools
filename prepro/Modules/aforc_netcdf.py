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

# ---------------------------------------------------
# FUNCTIONS USED BY make_aforc.py TO FIND INPUTFILES
# ---------------------------------------------------
def find_input(variables,input_dir,input_prefix,Ystart,Mstart,Yend,Mend,multi_files):
# -------------------------------------------------
# Warning : only for path type /data_origin/file.nc or /data_origin/Y/M/file.nc
# -------------------------------------------------
# Finding input_file paths
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
    return input_file




# ---------------------------------------------------
# FUNCTIONS USED BY make_aforc.py TO CREATE NETCDF
# ---------------------------------------------------

# subset_data : 
# output_dir : 
# output_filename : 
# output_file_format : 

def output_name(subset_data,output_dir,output_filename,output_file_format):
    if output_file_format.upper() == "MONTHLY":
        aforc_filename = output_dir+output_filename.replace('croco_aforc.nc', 'Y%sM%02i.nc' %(subset_data['time'].dt.year[0].item(),subset_data['time'].dt.month[0].item()))
    elif output_file_format.upper() == "DAILY":
        aforc_filename = output_dir+output_filename.replace('croco_aforc.nc', 'Y%sM%02iD%02i.nc' %(subset_data['time'].dt.year[0].item(),subset_data['time'].dt.month[0].item(),subset_data['time'].dt.day[0].item()))  
    return aforc_filename

# data :
# aforc_filename : 
# output_file_format :

def create_netcdf(data,aforc_filename,output_file_format):
    if output_file_format == 'MONTHLY':
        filename_out = aforc_filename[:-11]+data.name+'_'+aforc_filename[-11:]
    elif output_file_format == 'DAILY':
        filename_out = aforc_filename[:-14]+data.name+'_'+aforc_filename[-14:]
    data = data.astype(np.float32)
    data = data.to_netcdf(filename_out,engine='netcdf4',compute='False')
    return data  