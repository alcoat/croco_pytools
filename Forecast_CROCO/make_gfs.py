######################################################################
# Create and fill frc and bulk files with GFS data.
# for a forecast run
#
# The on-line reference to GFS is at
# http://nomad3.ncep.noaa.gov/
################################################################
#
# Common parameters
#
from datetime import datetime, timedelta
import numpy as np
from netCDF4 import Dataset
import pandas as pd
import Forecast_tools as ft
import Preprocessing_tools as ppt
from scipy.interpolate import interp2d
from croco_tools_params import *
#
it=2
#
frc_prefix=frc_prefix+'_GFS_'
blk_prefix=blk_prefix+'_GFS_'
#
################################################################
# end of user input  parameters
################################################################
#
# time (in matlab time)
#
rundate_str=datetime.today()
rundate=rundate_str-datetime(Yorig,1,1)
rundate=rundate.days
#
# GFS data name
#
gfs_name=FRCST_dir+'GFS_'+str(rundate)+'.nc'
#
#
if level==0:
  nc_suffix='.nc'
else:
  nc_suffix='.nc.'+str(level)
  grdname=grdname+'.'+str(level)
#
# Get the model grid
#
nc=Dataset(grdname, 'r',set_auto_maskandscale=False)
lon=nc['lon_rho'][:]
lat=nc['lat_rho'][:]
angle=nc['angle'][:]
h=nc['h'][:]
nc.close()
cosa=np.cos(angle)
sina=np.sin(angle)
#
# Extract data over the internet
#
#Download_data=0
if Download_data==1:
  #
  # Get model limits
  #
  lonmin=np.min(lon, axis=(0,1)) 
  lonmax=np.max(lon, axis=(0,1))
  latmin=np.min(lat, axis=(0,1))
  latmax=np.max(lat, axis=(0,1))
  #
  # Download data with DODS (the download matlab routine depends on the OGCM)
  # 
  print('Download data...')
  ft.download_GFS(rundate_str,lonmin,lonmax,latmin,latmax,FRCST_dir,Yorig,it)
#
#end
#
# Get the GFS grid 
# 
nc=Dataset(gfs_name)
lon1=nc['lon'][:]
lat1=nc['lat'][:]
time=nc['time'][:]
mask=nc['mask'][:]
tlen=len(time)
#
# bulk and forcing files
#
blkname=blk_prefix+str(rundate)+nc_suffix
print('Create a new bulk file: '+blkname)
ppt.create_bulk(blkname,grdname,CROCO_title,time,0)
nc_blk=Dataset(blkname,'r+')
frcname=frc_prefix+str(rundate)+nc_suffix
print('Create a new forcing file: '+frcname)
ppt.create_forcing(frcname,grdname,CROCO_title,time,0,0,0,0,0,0,0,0,0,0,0)
nc_frc=Dataset(frcname,'r+')

#
# Loop on time
#
missval=float('nan')
default=float('nan')

for l in range(tlen):
  print('time index: ',str(l+1),' of total: '+str(tlen))
  var=nc['tair'][l,:,:]
  if np.mean(np.isnan(var)!=1):
    var=ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    f = interp2d(lon1,lat1,var,interp_method)
    nc_blk['tair'][l,:,:]= f(lon, lat) 
  else:
    var=nc['tair'][l-1,:,:] 
    var=ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    f = interp2d(lon1,lat1,var,interp_method)
    nc_blk['tair'][l,:,:]= f(lon, lat)
  #end

  var=nc['rhum'][l,:,:]
  if np.mean(np.isnan(var)!=1):
    var=ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    f = interp2d(lon1,lat1,var,interp_method)
    nc_blk['rhum'][l,:,:]= f(lon, lat)
  else:
    var=nc['rhum'][l-1,:,:] 
    var=ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    f = interp2d(lon1,lat1,var,interp_method)
    nc_blk['rhum'][l,:,:]= f(lon, lat)
  #end
  
  var=nc['prate'][l,:,:]
  if np.mean(np.isnan(var)!=1):
    var=ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    f = interp2d(lon1,lat1,var,interp_method)
    nc_blk['prate'][l,:,:]= f(lon, lat) 
  else:
    var=nc['prate'][l-1,:,:] 
    var=ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    f = interp2d(lon1,lat1,var,interp_method)
    nc_blk['prate'][l,:,:]= f(lon, lat)
  #end
  
  var=nc['wspd'][l,:,:]
  if np.mean(np.isnan(var)!=1):
    var=ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    f = interp2d(lon1,lat1,var,interp_method)
    nc_blk['wspd'][l,:,:]= f(lon, lat) 
  else:
    var=nc['wspd'][l-1,:,:] 
    var=ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    f = interp2d(lon1,lat1,var,interp_method)
    nc_blk['wspd'][l,:,:]= f(lon, lat)
  #end
  
  #Zonal wind speed
  var=nc['uwnd'][l,:,:]
  if np.mean(np.isnan(var)!=1):
    uwnd=ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    f = interp2d(lon1,lat1,uwnd,interp_method)
    uwnd= f(lon, lat) 
  else:
    var=nc['uwnd'][l-1,:,:] 
    uwnd=ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    f = interp2d(lon1,lat1,uwnd,interp_method)
    uwnd = f(lon, lat)
  #end
  
  #Meridian wind speed
  var=nc['vwnd'][l,:,:]
  if np.mean(np.isnan(var)!=1):
    vwnd=ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    f = interp2d(lon1,lat1,vwnd,interp_method)
    vwnd = f(lon, lat) 
  else:
    var=nc['vwnd'][l-1,:,:] 
    vwnd=ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    f = interp2d(lon1,lat1,vwnd,interp_method)
    vwnd= f(lon, lat)
  #end
  
  nc_frc['uwnd'][l,:,:]=ppt.rho2u_2d(uwnd*cosa+vwnd*sina)
  nc_frc['vwnd'][l,:,:]=ppt.rho2v_2d(vwnd*cosa-uwnd*sina)
  
  nc_blk['uwnd'][l,:,:]=ppt.rho2u_2d(uwnd*cosa+vwnd*sina)
  nc_blk['vwnd'][l,:,:]=ppt.rho2v_2d(vwnd*cosa-uwnd*sina)
  
  #Net longwave flux
  var=nc['radlw'][l,:,:]
  if np.mean(np.isnan(var)!=1):
    var=ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    f = interp2d(lon1,lat1,var,interp_method)
    nc_blk['radlw'][l,:,:]= f(lon, lat) 
  else:
    var=nc['radlw'][l-1,:,:] 
    var=ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    f = interp2d(lon1,lat1,var,interp_method)
    nc_blk['radlw'][l,:,:]= f(lon, lat)
  #end
  
  #Downward longwave flux
  var=nc['radlw_in'][l,:,:]
  if np.mean(np.isnan(var)!=1):
    var=ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    f = interp2d(lon1,lat1,var,interp_method)
    nc_blk['radlw_in'][l,:,:]= f(lon, lat)
  else:
    var=nc['radlw_in'][l-1,:,:] 
    var=ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    f = interp2d(lon1,lat1,var,interp_method)
    nc_blk['radlw_in'][l,:,:]= f(lon, lat)
  #end
   
  #Net solar short wave radiation
  var=nc['radsw'][l,:,:]
  if np.mean(np.isnan(var)!=1):
    var=ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    f = interp2d(lon1,lat1,var,interp_method)
    nc_blk['radsw'][l,:,:]= f(lon, lat)
  else:
    var=nc['radsw'][l-1,:,:] 
    var=ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    f = interp2d(lon1,lat1,var,interp_method)
    nc_blk['radsw'][l,:,:]= f(lon, lat)
  #end
  
  var=nc['tx'][l,:,:]
  if np.mean(np.isnan(var)!=1):
    tx=ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    f = interp2d(lon1,lat1,tx,interp_method)
    tx = f(lon, lat)
  else:
    var=nc['tx'][l,:,:]
    tx=ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    f = interp2d(lon1,lat1,tx,interp_method)
    tx = f(lon, lat)
  #end
  
  var=nc['ty'][l,:,:]
  if np.mean(np.isnan(var)!=1):
    ty=ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    f=interp2d(lon1,lat1,ty,interp_method)
    ty = f(lon, lat)
  else:
    var=nc['ty'][l,:,:]
    ty=ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    f = interp2d(lon1,lat1,ty,interp_method)
    ty = f(lon, lat)
  #end
  
  nc_frc['sustr'][l,:,:]=ppt.rho2u_2d(tx*cosa+ty*sina)
  nc_frc['svstr'][l,:,:]=ppt.rho2v_2d(ty*cosa-tx*sina)
  
  nc_blk['sustr'][l,:,:]=ppt.rho2u_2d(tx*cosa+ty*sina)
  nc_blk['svstr'][l,:,:]=ppt.rho2v_2d(ty*cosa-tx*sina)
#end
# 
nc_frc.close()
nc_blk.close()
nc.close()