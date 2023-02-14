from math import dist
import numpy as np
import sys
from scipy.interpolate import interp2d
import Preprocessing_tools as ppt
from netCDF4 import Dataset

def ext_data_OGCM(nc,X,Y,vname,tndx,lon,lat,k,Roa,interp_method):
  #
  # extrapolation parameters
  #
  default=0
  if vname=='SAVE' or vname=='salt':
    default=34.6
  #
  # Get the CROCO grid extension + a little margin (~ 2 data grid points)
  #
  dx=max(np.abs(np.gradient(X)))
  dy=max(np.abs(np.gradient(Y)))
  dl=2*max([dx, dy])
  #
  lonmin=(np.min(lon, axis=(0,1)))-dl
  lonmax=(np.max(lon, axis=(0,1)))+dl
  latmin=(np.min(lat, axis=(0,1)))-dl
  latmax=(np.max(lat, axis=(0,1)))+dl
  #
  # Extract a data subgrid
  #
  j=np.where((Y>=latmin) & (Y<=latmax))
  i1=np.where((X-360>=lonmin) & (X-360<=lonmax))
  i2=np.where((X>=lonmin) & (X<=lonmax))
  i3=np.where((X+360>=lonmin) & (X+360<=lonmax))
  if len(i2)!=0:
    x=X[i2]
  else:
    x=[]
  #
  if len(i1)!=0:
    x=np.concatenate(((X[i1]-360),x))
  #
  if len(i3)!=0:
    x=np.concatenate((x,(X[i3]+360)))
  #
  y=Y[j]
  #
  #  Get dimensions
  #
  #vname
  ndims=len(nc[vname][:].shape)
  #
  # Get data (Horizontal 2D matrix)
  #
  if len(i2[0])!=0:
    if ndims==2:
      data=(nc[vname][j[0],i2[0]])
    elif ndims==3:
      data=(nc[vname][tndx,j[0],i2[0]])
    elif ndims==4:
      data=(nc[vname][tndx,k,j[0],i2[0]])
    else:
      sys.exit('Bad dimension number ',str(ndims))

    #
  else:
    data=[]
  #  
  if len(i1[0])!=0:
    if ndims==2:
      data=np.concatenate(((nc[vname][j[0],i1[0]]),data), axis=1)
    elif ndims==3:
      print(nc[vname][tndx,j[0]])
      print(i1)
      data=np.concatenate(((nc[vname][tndx,j[0],i1[0]]),data), axis=1)
    elif ndims==4:
      data=np.concatenate(((nc[vname][tndx,k,j[0],i1[0]]),data), axis=1)
    else:
      sys.exit('Bad dimension number ',str(ndims))

    #
  #
  if len(i3[0])!=0:
    if ndims==2:
      data=np.concatenate((data,(nc[vname][j[0],i3[0]])), axis=1)
    elif ndims==3:
      data=np.concatenate((data,(nc[vname][tndx,j[0],i3[0]])), axis=1)
    elif ndims==4:
      data=np.concatenate((data,(nc[vname][tndx,k,j[0],i3[0]])), axis=1)
    else:
      sys.exit('Bad dimension number ',str(ndims))

  #
  #
  # Perform the extrapolation
  #

  data,interp_flag=ppt.get_missing_val(x,y,data.data,float('nan'),Roa,default)
  #
  # Interpolation on the CROCO grid
  #
  if interp_flag==0:
    f = interp2d(x,y,data.data,kind='linear')
    data= f(lon[0,:], lat[:,0])     
  else:
    f = interp2d(x,y,data.data,kind=interp_method)
    data= f(lon[0,:], lat[:,0])
  #
  #
  return data


######################################################################
# Create the OGCM file
######################################################################
def create_OGCM(fname,lonT,latT,lonU,latU,lonV,latV,depth,time,temp,salt,u,v,ssh,Yorig):

  missval=float('nan')
  print('    Create the OGCM file')
  
  nc = Dataset(fname,  "w", format="NETCDF3_CLASSIC", clobber=True)

  nc.createDimension('lonT', len(lonT))
  nc.createDimension('latT', len(latT))
  nc.createDimension('lonU', len(lonU))
  nc.createDimension('latU', len(latU))
  nc.createDimension('lonV', len(lonV))
  nc.createDimension('latV', len(latV))
  nc.createDimension('depth', len(depth))
  nc.createDimension('time', len(time))
  #redef(nc)

  temp_nc = nc.createVariable('temp', 'f', ('time','depth','latT','lonT'))
  temp_nc.long_name='TEMPERATURE'
  temp_nc.units='deg. C'
  temp_nc.missing_value=missval

  salt_nc = nc.createVariable('salt', 'f', ('time','depth','latT','lonT'))
  salt_nc.long_name='SALINITY'
  salt_nc.units='ppt'
  salt_nc.missing_value=missval

  u_nc = nc.createVariable('u', 'f', ('time','depth','latU','lonU'))
  u_nc.long_name='ZONAL VELOCITY'
  u_nc.units='m/sec'
  u_nc.missing_value=missval

  v_nc = nc.createVariable('v', 'f', ('time','depth','latV','lonV'))
  v_nc.long_name='MERIDIONAL VELOCITY'
  v_nc.units='m/sec'
  v_nc.missing_value=missval

  ubar_nc = nc.createVariable('ubar', 'f', ('time','latU','lonU'))
  ubar_nc.long_name='ZONAL BAROTROPIC VELOCITY'
  ubar_nc.units='m/sec'
  ubar_nc.missing_value=missval

  vbar_nc = nc.createVariable('vbar', 'f', ('time','latV','lonV'))
  vbar_nc.long_name='MERIDIONAL BAROTROPIC VELOCITY'
  vbar_nc.units='m/sec'
  vbar_nc.missing_value=missval

  ssh_nc = nc.createVariable('ssh', 'f', ('time','latT','lonT'))
  ssh_nc.long_name='SEA LEVEL HEIGHT'
  ssh_nc.units='m'
  ssh_nc.missing_value=missval

  lonT_nc = nc.createVariable('lonT', 'd', ('lonT'))
  lonT_nc.units='degrees_east'

  latT_nc = nc.createVariable('latT', 'd', ('latT'))
  latT_nc.units='degrees_north'

  lonU_nc = nc.createVariable('lonU', 'd', ('lonU'))
  lonU_nc.units='degrees_east'

  latU_nc = nc.createVariable('latU', 'd', ('latU'))
  latU_nc.units='degrees_north'

  lonV_nc = nc.createVariable('lonV', 'd', ('lonV'))
  lonV_nc.units='degrees_east'

  latV_nc = nc.createVariable('latV', 'd', ('latV'))
  latV_nc.units='degrees_north'

  depth_nc = nc.createVariable('depth', 'd', ('depth'))
  depth_nc.units='meters'

  time_nc = nc.createVariable('time', 'd', ('time'))
  time_nc.units='days since 1-Jan-'+str(Yorig)+' 00:00:0.0'

  #
  # Fill the file
  #
  print('    Fill the OGCM file')
  depth_nc[:]=depth
  latT_nc[:]=latT
  lonT_nc[:]=lonT
  latU_nc[:]=latU
  lonU_nc[:]=lonU
  latV_nc[:]=latV
  lonV_nc[:]=lonV
  #

  nc['time'][:] = time
  ssh_nc[:,:,:]=ssh
  u1=u
  v1=v
  u_nc[:,:,:,:]=u1
  v_nc[:,:,:,:]=v1
  temp_nc[:,:,:,:]=temp
  salt_nc[:,:,:,:]=salt
  
  for tndx in range(len(time)):
    if len(time)==1:
      u1=u
      v1=v
    else:
      u1=u[tndx,:,:,:] #
      v1=v[tndx,:,:,:] #
    #
    # Compute the barotropic velocities
    #
    masku=np.ones(u1.shape)
    masku[u1.mask==True]=0
    maskv=np.ones(v1.shape)
    maskv[v1.mask==True]=0

    u1[u1.mask==True]=0
    v1[v1.mask==True]=0
    dz=np.gradient(depth)
    NZ=len(depth)

    du=0*u1[0,:,:]
    zu=du
    dv=0*v1[0,:,:]
    zv=dv
    for k in range(0,NZ):
      du=du+dz[k]*u1[k,:,:]
      zu=zu+dz[k]*masku[k,:,:]
      dv=dv+dz[k]*v1[k,:,:]
      zv=zv+dz[k]*maskv[k,:,:]
    #
    du[zu==0]=float('nan')
    dv[zv==0]=float('nan')
    zu[zu==0]=float('nan')
    zv[zv==0]=float('nan')
    ubar=du/zu
    vbar=dv/zv
    #
    ubar_nc[tndx,:,:]=ubar
    vbar_nc[tndx,:,:]=vbar
    #
  nc.close()
  #
  return
