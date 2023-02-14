######################################################################
#
# Extract a subgrid from mercator to get a CROCO forcing
# Store that into monthly files.
# Take care of the Greenwitch Meridian.

############################# LYBRARIES ##############################
from datetime import datetime, timedelta
import numpy as np
from netCDF4 import Dataset 
import pandas as pd
import get_file_python_mercator as gfm
from pydap.client import open_url
import pprint
import subprocess
import Preprocessing_tools as ppt
import oforc_OGCM as ofgm
import os
import sys
from croco_tools_params import *

######################################################################
#
# Extract a subset from Marcator using python motu client (cls)
# Write it in a local file (keeping the classic SODA netcdf format)
######################################################################
def write_mercator_frcst(FRCST_dir,FRCST_prefix,raw_mercator_name,mercator_type,vars,time,Yorig):
    #
    print('    writing MERCATOR file')
    #
    # Get grid and time frame
    #
    nc = Dataset(raw_mercator_name, set_auto_maskandscale=True)
    if mercator_type==1:
        lon = nc.variables['longitude'][:]
        lat = nc.variables['latitude'][:]
        depth = nc.variables['depth'][:]
        time = nc.variables['time'][:]
        diff_dates = (datetime(1950,1,1) - datetime(Yorig,1,1)).days
        #time = time.astype('float32') / 24. + diff_dates.astype('float32')
        time = time/24. + diff_dates
    else:
        lon = nc.variables['lon'][:]
        lat = nc.variables['lat'][:]
        depth = nc.variables['depth'][:]
        time = nc.variables['time'][:]
        diff_dates = (datetime(2014,1,9) - datetime(Yorig,1,1)).days
        #time = time.astype('float32') / 86400 + diff_dates.astype('float32')
        time = time/86400. + diff_dates
        
    #
    # Get SSH
    #
    print('    ...SSH')
    vname=vars[0]
    ncc=nc[vname]
    ssh=ncc[:,:,:]
    #
    # Get U
    #
    print('    ...U')
    vname=vars[1]
    ncc=nc[vname]
    u=ncc[:,:,:,:]
    #
    # Get V
    #
    print('    ...V')
    vname=vars[2]
    ncc=nc[vname]
    v=ncc[:,:,:,:]
    #
    # Get TEMP
    #
    print('    ...TEMP')
    vname=vars[3]
    ncc=nc[vname]
    temp=ncc[:,:,:,:]
    #
    # Get SALT
    #
    print('    ...SALT')
    vname=vars[4]
    ncc=nc[vname]
    salt=ncc[:,:,:,:]

    #
    # Create the Mercator file
    #
    rundate_str=datetime.today()
    rundate=rundate_str-datetime(Yorig,1,1)
    nc.close()


    ofgm.create_OGCM(FRCST_dir+FRCST_prefix+str(rundate.days)+'.cdf',lon,lat,lon,lat,lon,lat,depth,time,\
      (temp),(salt),(u),(v),(ssh),Yorig)

    #
    return


###############################################################################
def download_mercator(pathMotu,user,password,mercator_type,motu_url, service_id,product_id,lh,lf,lonmin,lonmax,latmin,latmax,zmax,FRCST_dir,FRCST_prefix,raw_mercator_name,Yorig):
  # pathMotu is a deprecated parameter !
  raw_exist = os.path.exists(raw_mercator_name)
  if not raw_exist: 
    download_raw_data=1
    print(' ')
    print('Downloading Raw Mercator File')
  else:
    download_raw_data=0
    print(' ')
    print('Raw Mercator file already downloaded, converting to croco_tools')
    print(' ')

  
  convert_raw2crocotools=1 # convert -> crocotools format data
  #
  # Set variable names according to mercator type data
  #
  if mercator_type==1:
    vars = ' --variable zos --variable uo --variable vo --variable thetao --variable so '
    var_names = ['zos', 'uo', 'vo', 'thetao', 'so']
  else:
    vars = ' --variable zos --variable uo --variable vo --variable thetao --variable so '
    var_names = ['zos', 'uo', 'vo', 'thetao', 'so']

  #
  # Get dates
  #
  rundate_str=datetime.today()
  rundate=rundate_str-datetime(Yorig,1,1)
  rundate=rundate.days

  time1 = []
  for i in range(1,lh+1):
    time1.append(rundate_str-timedelta(days=lh+2-i))

  time2=rundate_str
  time3 = []
  for j in range(1,lf+2):
    time3.append(rundate_str+timedelta(days=j))

  time=[*time1,time2,*time3]
  tiempo_inicial = time1[0]
  tiempo_final = time3[-1]

  if(lonmin > 180):
    lonmin = lonmin - 360
    
  if(lonmax > 180):
    lonmax = lonmax - 360

  print(' ')
  print('Get data for '+str(rundate_str))
  print('Minimum Longitude: '+str(lonmin))
  print('Maximum Longitude: '+str(lonmax))
  print('Minimum Latitude:  '+str(latmin))
  print('Maximum Latitude:  '+str(latmax))
  print(' ')

  if download_raw_data==1:
    #
    # Get data 
    #
    gfm.get_mercator(pathMotu,mercator_type,
    motu_url,service_id,product_id,vars,
    [str(lonmin-1), str(lonmax+1), str(latmin-1), str(latmax+1), str(0), str(zmax)],
    [str(tiempo_inicial), str(tiempo_final)],
    [user, password],
    raw_mercator_name)

  if convert_raw2crocotools==1: 
    #
    # Convert data format and write in a more CROCOTOOLS 
    # compatible input file 
    #
    print('Making output data directory '+FRCST_dir) # create directory
    frcst_ext = os.path.exists(FRCST_dir)
    if not frcst_ext:
      os.makedirs(FRCST_dir)
        #
    mercator_name=FRCST_dir+FRCST_prefix+str(rundate)+'.cdf'
    if mercator_name != None:
      print('Mercator file already exist => overwrite it')
    
    write_mercator_frcst(FRCST_dir,FRCST_prefix,raw_mercator_name,mercator_type,var_names,time,Yorig)

def interp_OGCM_frcst(OGCM_name,Roa,interp_method,lonU,latU,lonV,latV,lonT,latT,Z,tin,nc_clm,nc_bry,lon,lat,angle,h,pm,pn,rmask,tout,vtransform,obc):
  conserv=1 # same barotropic velocities as the OGCM
  #
  print(['  Horizontal interpolation: ',OGCM_name])
  #
  # CROCO grid angle
  #
  cosa=np.cos(angle)
  sina=np.sin(angle)
  #
  # Open the OGCM file
  #
  nc=Dataset(OGCM_name, 'r', set_auto_maskandscale=False)
  #
  # Interpole data on the OGCM Z grid and CROCO horizontal grid
  #
  #
  # Read and extrapole the 2D variables
  #
  zeta=ofgm.ext_data_OGCM(nc,lonT,latT,'ssh',tin,lon,lat,1,Roa,interp_method)
  u2d=ofgm.ext_data_OGCM(nc,lonU,latU,'ubar',tin,lon,lat,1,Roa,interp_method)
  v2d=ofgm.ext_data_OGCM(nc,lonV,latV,'vbar',tin,lon,lat,1,Roa,interp_method)
  ubar=ppt.rho2u_2d(u2d*cosa+v2d*sina)
  vbar=ppt.rho2v_2d(v2d*cosa-u2d*sina)
  #
  # Read and extrapole the 3D variables
  #
  NZ=len(Z)
  [M,L]=np.shape(lon)
  dz=np.gradient(Z)
  temp=np.zeros((NZ,M,L))
  salt=np.zeros((NZ,M,L))
  u=np.zeros((NZ,M,L-1))
  v=np.zeros((NZ,M-1,L))
  for k in range(NZ):
    if (k%10)==0:
      print('  Level ',str(k),' of ',str(NZ))
    #
    u2d=ofgm.ext_data_OGCM(nc,lonU,latU,'u',tin,lon,lat,k,Roa,interp_method)
    v2d=ofgm.ext_data_OGCM(nc,lonV,latV,'v',tin,lon,lat,k,Roa,interp_method)
    u[k]=ppt.rho2u_2d(u2d*cosa+v2d*sina)                                 
    v[k]=ppt.rho2v_2d(v2d*cosa-u2d*sina)
    temp[k]=ofgm.ext_data_OGCM(nc,lonT,latT,'temp',tin,lon,lat,k,Roa,interp_method)
    salt[k]=ofgm.ext_data_OGCM(nc,lonT,latT,'salt',tin,lon,lat,k,Roa,interp_method)
  #
  #
  # Close the OGCM file
  #
  nc.close()
  #
  # Get the CROCO vertical grid
  #
  print('  Vertical interpolations')
  if type(nc_clm)!=list:
    theta_s=nc_clm['theta_s'][:]
    theta_b=nc_clm['theta_b'][:]
    hc=nc_clm['hc'][:]
    N=nc_clm.dimensions['s_rho'].size
  #
  if type(nc_bry)!=list:
    theta_s=nc_bry['theta_s'][:]
    theta_b=nc_bry['theta_b'][:]
    hc=nc_bry['hc'][:]
    N=nc_bry.dimensions['s_rho'].size
  #
  zr=ppt.zlevs(h,zeta,theta_s,theta_b,hc,N,'r',vtransform)
  zu=ppt.rho2u_3d(zr)    
  zv=ppt.rho2v_3d(zr)
  zw=ppt.zlevs(h,zeta,theta_s,theta_b,hc,N,'w',vtransform)
  dzr=zw[1:,:,:]-zw[:-1,:,:]
  dzu=ppt.rho2u_3d(dzr)
  dzv=ppt.rho2v_3d(dzr)
  #
  # Add an extra bottom layer (-100000m) and an extra surface layer (+100m)
  # to prevent vertical extrapolations
  #
  pre_z = np.zeros(len(Z)+2)
  pre_z[0] = 100.
  pre_z[-1] = -100000.
  pre_z[1:-1]=Z
  Z=pre_z
  u=np.concatenate((np.reshape(u[0,:,:], (1,u.shape[1], u.shape[2])),u))
  u=np.concatenate((u,np.reshape(u[-1,:,:], (1,u.shape[1], u.shape[2]))))
  v=np.concatenate((np.reshape(v[0,:,:], (1,v.shape[1], v.shape[2])),v))
  v=np.concatenate((v,np.reshape(v[-1,:,:], (1,v.shape[1], v.shape[2]))))
  temp=np.concatenate((np.reshape(temp[0,:,:], (1, temp.shape[1], temp.shape[2])),temp))
  temp=np.concatenate((temp,np.reshape(temp[-1,:,:], (1,temp.shape[1], temp.shape[2]))))
  salt=np.concatenate((salt,np.reshape(salt[-1,:,:], (1,salt.shape[1], salt.shape[2]))))
  salt=np.concatenate((np.reshape(salt[0,:,:], (1,salt.shape[1], salt.shape[2])),salt))
  # 
  # Perform the vertical interpolations 
  #
  temp=ppt.ztosigma(np.flip(temp,0),zr,np.flipud(Z))
  salt=ppt.ztosigma(np.flip(salt,0),zr,np.flipud(Z))
  u=ppt.ztosigma(np.flip(u,0),zu,np.flipud(Z))
  v=ppt.ztosigma(np.flip(v,0),zv,np.flipud(Z))
  #
  # Correct the horizontal transport 
  # i.e. remove the interpolated tranport and add 
  #      the OGCM transport
  #
  if conserv==1:
    u=u-ppt.tridim((np.sum(u*dzu, axis=0)/np.sum(dzu, axis=0)),N)
    v=v-ppt.tridim((np.sum(v*dzv, axis=0)/np.sum(dzv, axis=0)),N)
    u=u+ppt.tridim(ubar,N)
    v=v+ppt.tridim(vbar,N)
  #end
  #
  # Barotropic velocities
  #
  ubar=(np.sum(u*dzu, axis=0)/np.sum(dzu, axis=0))
  vbar=(np.sum(v*dzv, axis=0)/np.sum(dzv, axis=0))
  #
  #  fill the files
  #
  if type(nc_clm)!=list: 
    nc_clm['zeta'][tout,:,:]=zeta
    if 'SSH' in list(nc_clm.variables):
      nc_clm['SSH'][tout,:,:]=zeta
    nc_clm['temp'][tout,:,:,:]=temp
    nc_clm['salt'][tout,:,:,:]=salt
    nc_clm['u'][tout,:,:,:]=u
    nc_clm['v'][tout,:,:,:]=v
    nc_clm['ubar'][tout,:,:]=ubar
    nc_clm['vbar'][tout,:,:]=vbar
  #end
  if type(nc_bry)!=list:
    for obcndx in range(0,4):
      if obcndx==0:
        nc_bry['zeta_south'][tout,:]=zeta[0,:]
        nc_bry['temp_south'][tout,:,:]=temp[:,0,:]
        nc_bry['salt_south'][tout,:,:]=salt[:,0,:]
        nc_bry['u_south'][tout,:,:]=u[:,0,:]
        nc_bry['v_south'][tout,:,:]=v[:,0,:]
        nc_bry['ubar_south'][tout,:]=ubar[0,:]
        nc_bry['vbar_south'][tout,:]=vbar[0,:]
      elif obcndx==1:
        nc_bry['zeta_east'][tout,:]=zeta[:,-1]
        nc_bry['temp_east'][tout,:,:]=temp[:,:,-1]
        nc_bry['salt_east'][tout,:,:]=salt[:,:,-1]
        nc_bry['u_east'][tout,:,:]=u[:,:,-1]
        nc_bry['v_east'][tout,:,:]=v[:,:,-1]
        nc_bry['ubar_east'][tout,:]=ubar[:,-1]
        nc_bry['vbar_east'][tout,:]=vbar[:,-1]
      elif obcndx==2:
        nc_bry['zeta_north'][tout,:]=zeta[-1,:]
        nc_bry['temp_north'][tout,:,:]=temp[:,-1,:]
        nc_bry['salt_north'][tout,:,:]=salt[:,-1,:]
        nc_bry['u_north'][tout,:,:]=u[:,-1,:]
        nc_bry['v_north'][tout,:,:]=v[:,-1,:]
        nc_bry['ubar_north'][tout,:]=ubar[-1,:]
        nc_bry['vbar_north'][tout,:]=vbar[-1,:]
      elif obcndx==3:
        nc_bry['zeta_west'][tout,:]=zeta[:,0]
        nc_bry['temp_west'][tout,:,:]=temp[:,:,0]
        nc_bry['salt_west'][tout,:,:]=salt[:,:,0]
        nc_bry['u_west'][tout,:,:]=u[:,:,0]
        nc_bry['v_west'][tout,:,:]=v[:,:,0]
        nc_bry['ubar_west'][tout,:]=ubar[:,0]
        nc_bry['vbar_west'][tout,:]=vbar[:,0]
      #end
    #end
  return nc_clm, nc_bry

######################################################################
#  Reproduce the old readattribute behavior when using loaddap library
#  But uses Built-in Support for OPeNDAP from Matlab >= 2012a 
#
#  Get the attribute of an OPENDAP dataset
######################################################################
def loaddap(url):
  x={}
  #
  try:
    ncid = open_url(url)
    
    nvar=len(ncid.keys())
    for ii in range(nvar):
      varname = list(ncid.keys())[ii]
      nattvar = len(ncid[varname].attributes)
      varname2=varname
      varname2=varname.replace('-','_2d')
      for jj in range(nattvar):
        attname = list(ncid[varname].attributes)[jj]
        attdict = ncid[varname].attributes
        
        if attname.find('missing_value')!=-1 or attname.find('FillValue')!=-1:
          attval = attdict[attname]
          #x[varname2].setncattr('missing_value',attval)
          #x[varname2].setncattr('ml__FillValue',attval)
          x['missing_value'] = attval
          x['ml__FillValue'] = attval
          #end
        if attname.find('scale_factor')!=-1:
          attval = attdict[attname]
          #x[varname2].setncattr('scale_factor',attval)
          x['scale_factor'] = attval
          #end
        if attname.find('add_offset')!=-1:
          attval = attdict[attname]
          #x[varname2].setncattr('add_offset',attval)
          x['add_offset'] = attval
          #end
        if attname.find('units')!=-1:
          attval = attdict[attname]
          #x[varname2].setncattr('units',attval)
          x['units'] = attval
          #end
        if attname.find('DODS_ML_Size')!=-1:
          attval = attdict[attname]
          #x[varname2].setncattr('DODS_ML_Size',attval)
          x['DODS_ML_Size'] = attval
  except:
    print('PB with loaddap')
    x=[]

  return x

def find_all(a_str, sub):
  start = 0
  while True:
    start = a_str.find(sub, start)
    if start == -1: return
    yield start
    start += len(sub)

######################################################################
#  Reproduce the old readdap behavior when using loaddap library
#  But uses Built-in Support for OPeNDAP from Matlab >= 2012a 
#  Retry (100 times) in case of network failure.
######################################################################
#
def readdap(url,varname,query):
  nmax=100
  data=[]
  ntry=0
  nargin = readdap.__code__.co_argcount
  #
  if nargin <2:
    print('not engough input argments')
  elif nargin <3 or len(query)==0:
    print('READDAP_New: Extract : '+ varname)
    while len(data)==0:
      if ntry>nmax:
        sys.exit('READDAP_New: repeated failures after '+str(nmax)+' queries')
      #end
      ntry=ntry+1
      ncid = open_url(url)
      if len(ncid)!=0:  
        data = ncid[varname][:].data
      else:
        data=[]
        print('READDAP_New: did not work at '+str(ntry)+' try: lets try again.')
      #end
    #end
      
  else:
    print('READDAP_New: Extract : ', varname, query)
    ind1 = list(find_all(query, '['))
    ind2 = list(find_all(query, ']'))

    nb_dims=len(ind1)
    start2=[]
    count2=[]

    for ii in range(nb_dims):
      str_tmp=query[ind1[ii]+1:ind2[ii]]
      if str_tmp.find(':')==-1:
        start2.append(int(str_tmp))
      else:
        start2.append(int(str_tmp[:str_tmp.find(':')]))
        count2.append(int(str_tmp[str_tmp.find(':')+1:]))#-start2[ii]+1
      #end
    #end
    #start=np.flip(start2)
    #count=np.flip(count2)
      #[startcount]
      
  while len(data)==0:
    if ntry>nmax:
      sys.exit('READDAP_New: repeated failures after '+str(nmax)+' queries')
      #end
    ntry=ntry+1
    try:
      ncid = open_url(url)

      if len(start2)==1:
        if start2[0]==count2[0]:
          data = ncid[varname][start2[0]].data
        else:
          data = ncid[varname][start2[0]:count2[0]].data
      elif len(start2)==2:
        data = ncid[varname].array[start2[0]:count2[0], start2[1]:count2[1]].data
        #data=np.transpose(data,(1,0))
      elif len(start2)==3:
        if start2[0]==count2[0]:
          data = ncid[varname].array[start2[0], start2[1]:count2[1], start2[2]:count2[2]].data
        else:
          data = ncid[varname].array[start2[0]:count2[0], start2[1]:count2[1], start2[2]:count2[2]].data
        #data=np.transpose(data,(2, 1, 0))               
      elif len(start2)==4:
        data = ncid[varname].array[start2[0]:count2[0], start2[1]:count2[1], start2[2]:count2[2], start2[3]:count2[3]].data
        #data=np.transpose(data,(3, 2, 1, 0))                 
        #end
    except:
      data=[]
      print('READDAP_New: did not work at '+str(ntry)+' try: lets try again.')
  #
  return data

######################################################################
#
#  var=getdap(path,fname,vname,trange,krange,jrange,...
#             i1min,i1max,i2min,i2max,i3min,i3max)
#
#  Download a data subsets from a OPENDAP server.
#
#  Take care of the greenwitch meridian
#  (i.e. get 3 subgrids defined by i1min,i1max,i2min,i2max,i3min,i3max
#  and concatenate them).

######################################################################
def getdap(path,fname,vname,trange,krange,jrange,i1min,i1max,i2min,i2max,i3min,i3max):
#
  url=path+fname
  #
  var=[]
  #
  if i2min:
    irange='['+str(i2min)+':'+str(i2max)+']'
    var=readdap(url,vname,trange+krange+jrange+irange)
  #end
  #
  if i1min:
    irange='['+str(i1min)+':'+str(i1max)+']'
    var0=readdap(url,vname,trange+krange+jrange+irange)
    if var:
      var=np.concatenate((var0,var), axis=-1)
    else:
      var = var0
  #end
  #
  if i3min:
    irange='[',str(i3min)+':'+str(i3max)+']'
    var0=readdap(url,vname,trange+krange+jrange+irange)
    if var:
      var=np.concatenate((var,var0), axis=-1)
    else:
      var=var0
  #end
  #
  return var


######################################################################
#  Give the GFS url for a given date.
######################################################################
def get_GFS_fname(time,gfs_run_time,gfstype):
  #
  # set URL
  #
  url='https://nomads.ncep.noaa.gov'
  #url='https://nomads-cprk.ncep.noaa.gov'

  # set file types
  if gfstype==0:
    gfsname ='fnl'
    gfsname1='fnlflx'     # 1/2 GDAS data
  else:
    gfsname ='gfs'
    gfsname1='gfs_0p25'   # 1/4 deg res GFS data
    #gfsname1='gfs_0p50'  # 1/2 deg res GFS data
  #end
  #
  # Get the date
  #
  str_time = time.strftime('%Y%m%d')
  stry=str_time[:4]
  strm=str_time[4:6]
  strd=str_time[6:]
  #end
  if gfs_run_time < 10:
    strh='_0'
  else:
    strh='_'
  #end
  #
  # Get the grid
  #
  if gfstype==0:
    gfsdir =url+'/dods/'+gfsname+'/'+gfsname+stry+strm+strd+'/'
    fname=gfsdir+gfsname1+strh+str(gfs_run_time)+'z'
  else:
    gfsdir =url+'/dods/'+gfsname1+'/'+gfsname+stry+strm+strd+'/'
    fname=gfsdir+gfsname1+strh+str(gfs_run_time)+'z'
  #end
  return fname

         
######################################################################
# Get the indices for a GFS subgrid 
######################################################################
def get_GFS_grid(fname,lonmin,lonmax,latmin,latmax):
  dl=1
  lonmin=lonmin-dl
  lonmax=lonmax+dl
  latmin=latmin-dl
  latmax=latmax+dl
  #
  # Get the global grid
  #
  nc=open_url(fname)
  lon = nc['lon'][:].data
  lat = nc['lat'][:].data
  #
  # Get a subgrid
  #
  # 1 Longitude: take care of greenwitch
  #
  i1=np.where((lon-360>=lonmin) & (lon-360<=lonmax))
  i2=np.where((lon>=lonmin) & (lon<=lonmax))
  i3=np.where((lon+360>=lonmin) & (lon+360<=lonmax))
  #
  lon=np.concatenate((lon[i1]-360,lon[i2],lon[i3]+360), axis=0)
  #
  if i1[0].shape[0]!=0:
    i1min=np.min(i1)-1
    i1max=np.max(i1)
  else:
    i1min=[]
    i1max=[]
  #end

  if i2[0].shape[0]!=0:
    i2min=np.min(i2)-1
    i2max=np.max(i2)
  else:
    i2min=[]
    i2max=[]
  #end

  if i3[0].shape[0]!=0:
    i3min=np.min(i3)-1
    i3max=np.max(i3)
  else:
    i3min=[]
    i3max=[]
  #end
  #
  # 2 Latitude
  #
  j=np.where((lat>=latmin) & (lat<=latmax))
  lat=lat[j]
  jmin=np.min(j)-1
  jmax=np.max(j)
  jrange='['+str(jmin)+':'+str(jmax)+']'
  #
  return [i1min,i1max,i2min,i2max,i3min,i3max,jrange,lon,lat]


         
######################################################################
# Download one full subset of GFS for CROCO bulk for 1 time step
# Put them in the CROCO units
######################################################################
def get_GFS(fname,mask,tndx,jrange,i1min,i1max,i2min,i2max,i3min,i3max,missvalue):
  if type(tndx)==int:
    trange='['+str(tndx)+':'+str(tndx)+']'
  else:
    trange='['+str(min(tndx))+':'+str(max(tndx))+']'
  #
  # Get GFS variables for 1 time step
  #
  print(' ')
  print('====================================================')

  t=readdap(fname,'time',trange)
  print('TRANGE='+str(trange))
  print('GFS raw time='+str(t))
  #t=t+365 # put it in "matlab" time
  print('GFS: ', datetime(1,1,1)+t*timedelta(days=1)-timedelta(days=2))
  print('====================================================')

  #print('u...')
  u=mask*getdap('',fname,'ugrd10m',trange,'',jrange,i1min,i1max,i2min,i2max,i3min,i3max)
  u[abs(u)>=missvalue]=float('nan')

  #print('v...')
  v=mask*getdap('',fname,'vgrd10m',trange,'',jrange,i1min,i1max,i2min,i2max,i3min,i3max)
  v[abs(v)>=missvalue]=float('nan')

  #print('ty...')
  ty=mask*getdap('',fname,'vflxsfc',trange,'',jrange,i1min,i1max,i2min,i2max,i3min,i3max)
  ty[abs(ty)>=missvalue]=float('nan')

  #print('tx...')
  tx=mask*getdap('',fname,'uflxsfc',trange,'',jrange,i1min,i1max,i2min,i2max,i3min,i3max)
  tx[abs(tx)>=missvalue]=float('nan')


  #print('skt...')
  #skt=mask.*getdap('',fname,'tmpsfc',trange,'',jrange,...
  #                  i1min,i1max,i2min,i2max,i3min,i3max)
  #skt(abs(skt)>=missvalue)=float('nan')

  #print('tair...')
  tair=mask*getdap('',fname,'tmp2m',trange,'',jrange,i1min,i1max,i2min,i2max,i3min,i3max)
  tair[abs(tair)>=missvalue]=float('nan')

  #print('rhum...')
  rhum=mask*getdap('',fname,'rh2m',trange,'',jrange,i1min,i1max,i2min,i2max,i3min,i3max)
  rhum[abs(rhum)>=missvalue]=float('nan')

  #print('prate...')
  prate=mask*getdap('',fname,'pratesfc',trange,'',jrange,i1min,i1max,i2min,i2max,i3min,i3max)
  prate[abs(prate)>=missvalue]=float('nan')

  #print('down radlw...')
  dradlw=mask*getdap('',fname,'dlwrfsfc',trange,'',jrange,i1min,i1max,i2min,i2max,i3min,i3max)
  dradlw[abs(dradlw)>=missvalue]=float('nan')

  #print('up radlw')
  uradlw=mask*getdap('',fname,'ulwrfsfc',trange,'',jrange,i1min,i1max,i2min,i2max,i3min,i3max)
  uradlw[abs(uradlw)>=missvalue]=float('nan')

  #print('down radsw')
  dradsw=mask*getdap('',fname,'dswrfsfc',trange,'',jrange,i1min,i1max,i2min,i2max,i3min,i3max)
  dradsw[abs(dradsw)>=missvalue]=float('nan')

  #print('up radsw')
  uradsw=mask*getdap('',fname,'uswrfsfc',trange,'',jrange,i1min,i1max,i2min,i2max,i3min,i3max)
  uradsw[abs(uradsw)>=missvalue]=float('nan')
  #
  # Transform the variables
  #
  #
  # 1: Air temperature: Convert from Kelvin to Celsius
  #
  tair=tair-273.15
  #
  # 2: Relative humidity: Convert from # to fraction
  #
  rhum=rhum/100
  #
  # 3: Precipitation rate: Convert from [kg/m^2/s] to cm/day
  #
  prate=prate*0.1*(24*60*60.0)
  prate[abs(prate)<1.e-4]=0
  #
  # 4: Net shortwave flux: [W/m^2]
  #    CROCO convention: positive downward: same as GFS
  # ?? albedo ??
  #
  radsw=dradsw - uradsw
  radsw[radsw<1.e-10]=0
  #
  # 5: Net outgoing Longwave flux:  [W/m^2]
  #    CROCO convention: positive upward (opposite to nswrs)
  #    GFS convention: positive downward --> * (-1)
  #    input: downward longwave rad. and
  #    skin temperature.
  #
  #skt=skt-273.15 
  #radlw=-lwhf(skt,radlw) 
  radlw = uradlw - dradlw
  radlw_in=dradlw

  #
  # 6: Wind speed
  #
  wspd=(u**2+v**2)**0.5
  # 7:  Wind vectors
  uwnd=u # rho point
  vwnd=v # rho point
  # #
  # # 7: Compute the stress following large and pond
  # #
  # [Cd,uu]=cdnlp(wspd,10.)
  # rhoa=air_dens(tair,rhum*100)
  # tx=Cd.*rhoa.*u.*wspd
  # ty=Cd.*rhoa.*v.*wspd
  #
  return [t,tx,ty,tair,rhum,prate,wspd,uwnd,vwnd,radlw,radlw_in,radsw]



######################################################################
#  function write_GFS(fname,Yorig,lon,lat,mask,time,tx,ty,...
#                     tair,rhum,prate,wspd,radlw,radsw)
#  Write into a GFS file
######################################################################
def write_GFS(fname,Yorig,lon,lat,mask,time,tx,ty,tair,rhum,prate,wspd,uwnd,vwnd,radlw,radlw_in,radsw):
  #
  print(['Create ',fname])
  nc=Dataset(fname,'w')
  #
  nc.createDimension('lon',len(lon))
  nc.createDimension('lat',len(lat))
  #
  # nc('latu') = len(lat)
  # nc('latv') = len(lat)-1
  # nc('lonu') = len(lon)-1
  # nc('lonv') = len(lon)
  #
  nc.createDimension('time', len(time))
  #
  nc_lon = nc.createVariable('lon','f',('lon'))
  nc_lon.long_name = 'longitude of RHO-points'
  nc_lon.units = 'degree_east'
  # 
  nc_lat = nc.createVariable('lat','f',('lat'))
  nc_lat.long_name = 'latitude of RHO-points'
  nc_lat.units = 'degree_north'
  #
  nc_time = nc.createVariable('time','f',('time'))
  nc_time.long_name = 'Time'
  nc_time.units = 'days since 1-Jan-'+str(Yorig)+' 00:00:0.0'
  #
  nc_mask = nc.createVariable('mask','f',('lat','lon'))
  nc_tx = nc.createVariable('tx','f',('time','lat','lon'))
  nc_ty = nc.createVariable('ty','f',('time','lat','lon'))
  nc_tair = nc.createVariable('tair','f',('time','lat','lon'))
  nc_rhum = nc.createVariable('rhum','f',('time','lat','lon'))
  nc_prate = nc.createVariable('prate','f',('time','lat','lon'))
  nc_wspd = nc.createVariable('wspd','f',('time','lat','lon'))
  nc_radlw = nc.createVariable('radlw','f',('time','lat','lon'))
  nc_radsw = nc.createVariable('radsw','f',('time','lat','lon'))
  nc_radlw_in = nc.createVariable('radlw_in','f',('time','lat','lon'))
  nc_uwnd = nc.createVariable('uwnd','f',('time','lat','lon'))
  nc_vwnd = nc.createVariable('vwnd','f',('time','lat','lon'))
  #
  #endef(nc)
  #
  
  nc_lon[:]=lon.data
  nc_lat[:]=lat.data
  nc_time[:]=time.data
  nc_mask[:]=mask[0].data
  nc_tx[:]=tx.data
  nc_ty[:]=ty.data
  nc_tair[:]=tair.data
  nc_rhum[:]=rhum.data
  nc_prate[:]=prate.data
  nc_wspd[:]=wspd.data
  nc_uwnd[:]=uwnd.data
  nc_vwnd[:]=vwnd.data
  nc_radlw[:]=radlw.data
  nc_radlw_in[:]=radlw_in.data
  nc_radsw[:]=radsw.data
  #
  nc.close()


######################################################################
#  download_GFS(today,lonmin,lonmax,latmin,latmax,FRCST_dir,Yorig)
#  Extract a subgrid from GFS to get a CROCO forcing
#  Store that into monthly files (to limit the problems
#  of bandwith...).
#  Take care of the Greenwitch Meridian.
#  Transform variables in the CROCO format.
######################################################################
def download_GFS(today,lonmin,lonmax,latmin,latmax,FRCST_dir,Yorig,it):
  #
  # Put the date in 'Yorig' time
  #
  #rundate=datenum(today)-datenum(Yorig,1,1)
  gfs_date0 = datetime.today()-timedelta(days=8)
  rundate_str=datetime.today()
  rundate=rundate_str-datetime(Yorig,1,1)
  rundate=rundate.days
  #
  # GFS output name
  #
  gfsftype=1  # GFS files
  gfs_name=FRCST_dir+'GFS_'+str(rundate)+'.nc'
  #
  # start
  #
  print(' ')
  print('Get GFS data for ', str(today))
  print('Minimum Longitude: ',str(lonmin))
  print('Maximum Longitude: ',str(lonmax))
  print('Minimum Latitude:  ',str(latmin))
  print('Maximum Latitude:  ',str(latmax))
  print(' ')
  #
  # Create directory if needed
  #
  print('Making output data directory '+FRCST_dir) # create directory
  frcst_ext = os.path.exists(FRCST_dir)
  if not frcst_ext:
    os.makedirs(FRCST_dir)
  #
  # Get GFS file name (first check if a forecast is available)
  #
  gfs_run_time=0
  gfs_date=today-timedelta(hours=12)
  foundfile=0
  while foundfile==0:
    fname=get_GFS_fname(gfs_date,gfs_run_time,gfsftype)
    x=loaddap(fname)
    #print(x)
    if len(x)!=0:
      foundfile=1
    else:
      foundfile=0
    #end
    if foundfile==1 & len(x)!=0:
      print('  File found')
    else:
      foundfile=0
      print('  GFS : did not found ',fname)
      gfs_run_time=gfs_run_time-6
      if gfs_run_time<0:
        gfs_date=gfs_date-timedelta(days=1)
        gfs_run_time=18
        if gfs_date<gfs_date0:
          sys.exit('  GFS: did not found anything')
  #
  # Get subgrid for GFS extraction
  #
  gfs_run_time_GFS=gfs_run_time 
  gfs_date_GFS=gfs_date
  fname=get_GFS_fname(gfs_date_GFS,gfs_run_time_GFS,gfsftype)
  [i1min,i1max,i2min,i2max,i3min,i3max,jrange,lon,lat] = get_GFS_grid(fname,lonmin,lonmax,latmin,latmax)
  #
  # Get mask on this grid
  #
  mask=getdap('',fname,'landsfc','[0:0]','',jrange,i1min,i1max,i2min,i2max,i3min,i3max)
  mask[mask==1]=float('nan')
  mask[np.isfinite(mask)]=1
  #
  # Initialize arrays with subgrid dimensions
  #
  Nrec=(hdays+fdays+1)*8/it
  [M,L]=mask[0].shape
  tx=np.zeros((int(Nrec),M,L))
  ty=tx
  tair=tx
  rhum=tx
  prate=tx
  wspd=tx
  radlw=tx
  radsw=tx
  radlw_in=tx
  uwnd=tx
  vwnd=tx
  #
  n=0
  gfstime=0*np.arange(0,Nrec)
  #
  #==================================================
  # 1: Get Hindcast data from past GFS forecast files
  #==================================================
  #
  print(' ... READ HINDCAST FILES ... ')
  #
  # Starting date : (hdays+1) days ago 06Z (12h will be added)
  #
  gfs_run_time=6
  gfs_date=today-(hdays+1)*timedelta(days=1) # (hdays+1) days ago 6Z
  #
  # Loop on past GFS forecast files 
  #     - Starts (hdays+1) days ago 18Z
  #     - Increment every 6h
  #     - Ends yesterday 12Z
  #
  for ihind in range(4*hdays):    # loop on files until time=yesterday 12Z 
    gfs_run_time=gfs_run_time+6  # add 6 hours
    if gfs_run_time>18:
      gfs_date=gfs_date+1*timedelta(days=1)
      gfs_run_time=0
    #end
  #
    tndxdap=2  # take 3rd rec. of each file => add another 6 hours
  #
  # 1.1: check if GFS file is available.
  #
    fname=get_GFS_fname(gfs_date,gfs_run_time,gfsftype)
    #warning off
    x=loaddap(fname)
    if len(x)!=0:
      foundfile=1
    else:
      foundfile=0
    #end
    if foundfile==1:
      print('  File found')
      gfs_run_time1=gfs_run_time  # Increment Forecast 
      gfs_date1=gfs_date          # Start time
      missvalue=x['missing_value']
    else:
      foundfile=0
      print('  GFS: file not found, try next file')
    #end
    #warning on
  #
  # 1.2: read file 
  #
    if foundfile==1: 
      n=n+1
      print('N_hindcast=',str(n))
      print(fname)

      [gfstime[n],tx[n,:,:],ty[n,:,:],tair[n,:,:],rhum[n,:,:],
      prate[n,:,:],wspd[n,:,:],uwnd[n,:,:],vwnd[n,:,:],
      radlw[n,:,:],radlw_in[n,:,:],radsw[n,:,:]]= get_GFS(fname,mask,tndxdap,jrange,i1min,i1max,i2min,i2max,i3min,i3max,missvalue)
    #end 
  #end # ihind (nb of hindcast record)
  #
  #==================================================================
  # 2: Get Forecast data, starting yesterday 18Z
  #==================================================================
  #
  print(' ... READ FORECAST FILES ... ')
  #
  # Get starting GFS date following last hindcast date
  #
  gfs_run_time=gfs_run_time1+6 # add 6 hours => 12Z
  gfs_date=gfs_date1
  if gfs_run_time>18:
    gfs_date=gfs_date+1
    gfs_run_time=0
  #end
  #
  # 2.1: check if GFS forecast is available for this time.
  #      if not: take previous one (but increment time index)
  #
  t1=2
  foundfile=0
  while foundfile==0:
    fname=get_GFS_fname(gfs_date,gfs_run_time,gfsftype)
    #warning off
    x=loaddap(fname)
    if len(x)!=0:
      foundfile=1
    else:
      foundfile=0
    #end
    if foundfile==1:
      print('  File found')
    else:
      foundfile=0
      print('  GFS : did not found ',fname)
      t1=t1+2
      gfs_run_time=gfs_run_time-6
      if gfs_run_time<0:
        gfs_date=gfs_date-1*timedelta(days=1)
        gfs_run_time=18
        if gfs_date<gfs_date0:
          sys.exit('  GFS: did not found anything')

  #
  # 2.2: loop on records in same forecast file (dt = 3h * it)
  #
  tend=(fdays+1)*8
  for tndx in range(t1+1,tend, it):
    tndxdap=tndx-1
    n=n+1
    print('N_forecast=',str(n))
    [gfstime[n],tx[n,:,:],ty[n,:,:],tair[n,:,:],rhum[n,:,:],
    prate[n,:,:],wspd[n,:,:],uwnd[n,:,:],vwnd[n,:,:],
    radlw[n,:,:],radlw_in[n,:,:],radsw[n,:,:]]= get_GFS(fname,mask,tndxdap,jrange,i1min,i1max,i2min,i2max, i3min,i3max,missvalue)
  #end
  #
  #==================================================================
  # 3: Finalize data processing and write in file
  #==================================================================
  #
  # Reduce the matrices
  #
  gfstime=gfstime[:n]
  tx=tx[:n,:,:]
  ty=ty[:n,:,:]
  tair=tair[:n,:,:]
  rhum=rhum[:n,:,:]
  prate=prate[:n,:,:]
  wspd=wspd[:n,:,:]
  uwnd=uwnd[:n,:,:]
  vwnd=vwnd[:n,:,:]
  radlw=radlw[:n,:,:]
  radlw_in=radlw_in[:n,:,:]
  radsw=radsw[:n,:,:]
  #
  # Set time in Yorig time
  #
  gfstime=gfstime-float((datetime(1,1,1)-datetime(Yorig,1,1)).days)
  #
  # Create the GFS output file and write everything down
  #
  mask[np.isnan(mask)]=0
  write_GFS(gfs_name,Yorig,lon,lat,mask,gfstime,tx,ty,tair,rhum,prate,wspd,uwnd,vwnd,radlw,radlw_in,radsw)
  #
  print('Download GFS: done')
  #
  return
