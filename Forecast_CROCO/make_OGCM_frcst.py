######################################################################
#
# Create and fill CROCO clim and bry files with OGCM data.
# for a forecast run

############################# LYBRARIES ##############################
from datetime import datetime, timedelta
import numpy as np
from netCDF4 import Dataset
import pandas as pd
import Forecast_tools as ft
import Preprocessing_tools as ppt

##################### USERS DEFINED VARIABLES ########################
#
# Common parameters
from croco_tools_params import *
# Get date
#
rundate_str=datetime.today()
rundate=rundate_str-datetime(Yorig,1,1)
rundate=rundate.days


# Set generic OGCM file name
#
FRCST_prefix=OGCM+'_'
OGCM_name=FRCST_dir+FRCST_prefix+str(rundate)+'.cdf'

#  MERCATOR : see get_file_python_mercator.m
#
raw_mercator_name=FRCST_dir+'raw_motu_mercator_'+str(rundate)+'.nc'

################################################################
# end of user input  parameters
################################################################
#
if level==0:
    nc_suffix='.nc'
else:
    nc_suffix='.nc.'+str(level)
    grdname=grdname+'.'+str(level)
#
# Get the model grid
#
grd=Dataset(grdname, 'r',set_auto_maskandscale=False)   
lon=grd['lon_rho'][:].data
lat=grd['lat_rho'][:].data
angle=grd['angle'][:].data
h=grd['h'][:].data
pm=grd['pm'][:].data
pn=grd['pn'][:].data
rmask=grd['mask_rho'][:].data
grd.close()

#---------------------------------------------------------------
# Extract data from the Internet
#---------------------------------------------------------------
if Download_data == 1:
    #
    # Get model limits
    #
    lonmin=np.min(lon, axis=(0,1)) 
    lonmax=np.max(lon, axis=(0,1))
    latmin=np.min(lat, axis=(0,1))
    latmax=np.max(lat, axis=(0,1))

    # Use Motu python
    #
    print('Download data...')
    dwn = ft.download_mercator(pathMotu,user,password,mercator_type,
    motu_url_fcst,service_id_fcst,product_id_fcst,
    hdays,fdays,
    lonmin,lonmax,latmin,latmax,hmax,
    FRCST_dir,FRCST_prefix,raw_mercator_name,Yorig)

#---------------------------------------------------------------
# Get OGCM grid 
#---------------------------------------------------------------
ogcm_grd=Dataset(OGCM_name, 'r', set_auto_maskandscale=False)
lonT=ogcm_grd['lonT'][:].data
latT=ogcm_grd['latT'][:].data
lonU=ogcm_grd['lonU'][:].data
latU=ogcm_grd['latU'][:].data
lonV=ogcm_grd['lonV'][:].data
latV=ogcm_grd['latV'][:].data
Z=-ogcm_grd['depth'][:].data
NZ=len(Z)
NZ=NZ-rmdepth
Z=Z[0:NZ]

#---------------------------------------------------------------
# Get time array 
#---------------------------------------------------------------

OGCM_time=ogcm_grd['time'][:].data
time_cycle=0
delta=1 # >1 if subsampling
trange=np.arange(0,len(OGCM_time),delta)
time = np.zeros_like(OGCM_time)
for i in range(len(trange)):
    time[i] = OGCM_time[trange[i]]

#---------------------------------------------------------------
# Initial file 
#---------------------------------------------------------------
if makeini==1:
    ininame=ini_prefix+str(rundate)+nc_suffix
    print('Create an initial file for '+str(rundate))
    nc_ini = ppt.create_inifile(ininame,grdname,CROCO_title,
    theta_s,theta_b,hc,N,
    rundate-hdays,'clobber',vtransform); # starts at 00:00
    
    nc_ini = Dataset(ininame, 'r+')
    nc_ini, nc_ept = ft.interp_OGCM_frcst(OGCM_name,Roa,interp_method,
    lonU,latU,lonV,latV,lonT,latT,Z,0,
    nc_ini,[],lon,lat,angle,h,pm,pn,rmask,
    0,vtransform,obc)

    nc_ini.close()

  #eval(['!cp ',ininame,' ',ini_prefix,'hct',nc_suffix])


#---------------------------------------------------------------
# Clim and Bry files
#---------------------------------------------------------------
if makeclim==1 or makebry==1:
    if makebry==1:
        bryname=bry_prefix+str(rundate)+nc_suffix
        nc_bry = ppt.create_bryfile(bryname,grdname,CROCO_title,[1, 1, 1, 1],
        theta_s,theta_b,hc,N,
        time,time_cycle,'clobber',vtransform)
        
        nc_bry = Dataset(bryname, 'r+')
    else:
        nc_bry=[]
    
    if makeclim==1:
        clmname=clm_prefix+str(rundate)+nc_suffix
        nc_clm = ppt.create_climfile(clmname,grdname,CROCO_title,
        theta_s,theta_b,hc,N,
        time,time_cycle,'clobber',vtransform)
        
        nc_clm = Dataset(clmname,'r+')
    else:
        nc_clm=[]

#---------------------------------------------------------------
# Perform interpolations for all selected records
#---------------------------------------------------------------
for tndx in range(0,len(time),1):
    print(' Time step : '+str(tndx+1)+' of '+str(len(time))+' :')
    nc_clm, nc_bry = ft.interp_OGCM_frcst(OGCM_name,Roa,interp_method,
    lonU,latU,lonV,latV,lonT,latT,Z,trange[tndx],
    nc_clm,nc_bry,lon,lat,angle,h,pm,pn,rmask,
    tndx,vtransform,obc)

nc_bry.close()
nc_clm.close()

print('TERMINADO')