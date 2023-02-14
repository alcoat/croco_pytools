from croco_tools_params import *
from netCDF4 import Dataset 
import sys
from datetime import datetime
import numpy as np
from scipy.interpolate import griddata
import dill 

################################################################
# function var_u=rho2u_2d(var_rho)
# interpole a field at rho points to a field at u points
################################################################
def rho2u_2d(var_rho):
    [Mp,Lp]=var_rho.shape
    L=Lp-1
    var_u=0.5*(var_rho[:,:L]+var_rho[:,1:])
    return var_u

################################################################
# function var_v=rho2v_2d(var_rho)
# interpole a field at rho points to a field at v points
################################################################
def rho2v_2d(var_rho):
    [Mp,Lp]=var_rho.shape
    M=Mp-1
    var_v=0.5*(var_rho[:M,:]+var_rho[1:,:])
    return var_v

################################################################
# function var_u=rho2u_3d(var_rho)
# interpole a field at rho points to a field at u points
################################################################
def rho2u_3d(var_rho):
    [N,Mp,Lp]=var_rho.shape
    L=Lp-1
    var_u=0.5*(var_rho[:,:,:L]+var_rho[:,:,1:])
    return var_u

################################################################
# function var_v=rho2v_3d(var_rho)
# interpole a field at rho points to a field at v points
################################################################
def rho2v_3d(var_rho):
    [N,Mp,Lp]=var_rho.shape
    M=Mp-1
    var_v=0.5*(var_rho[:,:M,:]+var_rho[:,1:,:])
    return var_v

################################################################
# function  vnew = ztosigma(var,z,depth)
# This function transform a variable from z to sigma coordinates
#    warning: the variable must be in the form: var(k,j,i)
################################################################
def ztosigma(var,z,depth):
  [Ns,Mp,Lp]=z.shape
  Nz=len(depth)
  vnew = np.zeros_like(z)
  #
  # Find the grid position of the nearest vertical levels
  #
  for ks in range(Ns):
    sigmalev=(z[ks,:,:])
    thezlevs=0*sigmalev
    for kz in range(Nz):
      thezlevs[sigmalev>depth[kz]]=thezlevs[sigmalev>depth[kz]]+1
      #end
    if np.max(thezlevs)>=Nz or np.min(thezlevs)<=0:
      print('min sigma level = ',str(min(min(min(z)))),' - min z level = ',str(min(depth)))
      print('max sigma level = ',str(max(max(max(z)))),' - max z level = ',str(max(depth)))
    #end
    thezlevs = thezlevs.astype('int')
    
    [imat,jmat]=np.meshgrid(np.arange(1, Lp+1),np.arange(1,Mp+1))
    pos=(Nz*Mp*(imat-1)+Nz*(jmat-1)+thezlevs)-1
    pos=pos.astype('int') 
    
    z1=depth[thezlevs-1]
    z2=depth[thezlevs]

    v1=var[np.unravel_index(pos, var.shape, 'F')]
    v2=var[np.unravel_index((pos+1), var.shape, 'F')]

    vnew[ks,:,:]=(((v1-v2)*sigmalev+v2*z1-v1*z2)/(z1-z2))
  #end
  return vnew


################################################################
# function var3d=tridim(var2d,N)
#  Put a 2D matrix in 3D (reproduce it N times).
################################################################
def tridim(var2d,N):
    [M,L]=np.shape(var2d)
    var3d=np.reshape(var2d,(1,M,L))
    var3d=np.repeat(var3d, N, axis=0)
    return var3d

################################################################
def csf(sc, theta_s,theta_b):
    if theta_s > 0:
        csrf=(1-np.cosh(sc*theta_s))/(np.cosh(theta_s)-1)
    else:
        csrf=-sc^2
    #
    if theta_b > 0:
        h = (np.exp(theta_b*csrf)-1)/(1-np.exp(-theta_b))
    else:
        h  = csrf
    return h


################################################################
def scoordinate(theta_s,theta_b,N,hc,vtransform,*args):
  narg = scoordinate.__code__.co_argcount
  if narg<4:
    vtransform = 1  #Old vtransform
    hc=[]
  elif narg<5:
    vtransform = 1  #Old vtransform

  # Set S-Curves in domain [-1 < sc < 0] at vertical W- and RHO-points.
  sc_r=np.zeros((N))
  Cs_r=np.zeros((N))
  sc_w=np.zeros((N+1))
  Cs_w=np.zeros((N+1))

  if (vtransform == 2):
    ds=1./N
    sc_r= ds*(np.arange(1,N+1)-N-0.5)
    Cs_r=csf(sc_r, theta_s,theta_b)
    #
    sc_w[0] = -1.0
    sc_w[-1] =  0
    Cs_w[0] = -1.0
    Cs_w[-1] =  0
    sc_w[1:N] = ds*(np.arange(1,N)-N)
    Cs_w=csf(sc_w, theta_s,theta_b)
  else:
    cff1=1./np.sinh(theta_s)
    cff2=0.5/np.tanh(0.5*theta_s)
    sc_w=(np.arange(0,N+1)-N)/N
    Cs_w=(1.-theta_b)*cff1*np.sinh(theta_s*sc_w)+theta_b*(cff2*np.tanh(theta_s*(sc_w+0.5))-0.5)
    #
    sc_r=(np.arange(1,N+1)-N-0.5)/N
    Cs_r=(1.-theta_b)*cff1*np.sinh(theta_s*sc_r)+theta_b*(cff2*np.tanh(theta_s*(sc_r+0.5))-0.5)
  return [sc_r,Cs_r,sc_w,Cs_w]

################################################################
#  function z = zlevs(h,zeta,theta_s,theta_b,hc,N,type,vtransform);
#  this function compute the depth of rho or w points for CROCO
################################################################
def zlevs(h,zeta,theta_s,theta_b,hc,N,type,vtransform,*args):
  [M,L]=h.shape
  narg = zlevs.__code__.co_argcount
  if narg < 8:
      print('WARNING no vtransform defined')
      vtransform = 1  #old vtranform = 1
      print('Default S-coordinate system use : Vtransform=1 (old one)')
  #end
  # Set S-Curves in domain [-1 < sc < 0] at vertical W- and RHO-points.


  sc_r=np.zeros((N))
  Cs_r=np.zeros((N))
  sc_w=np.zeros((N+1))
  Cs_w=np.zeros((N+1))

  if vtransform == 2:
    ds=1./N
    if type=='w':
      sc_w[0] = -1.0
      sc_w[-1] = 0
      Cs_w[0] = -1.0
      Cs_w[-1] =  0
        
      sc_w[1:N] = ds*(np.arange(1,N)-N)
      Cs_w=csf(sc_w, theta_s,theta_b)
      N=N+1
    else:
      sc= ds*(np.arange(1,N+1)-N-0.5)    
      Cs_r=csf(sc, theta_s,theta_b)
      sc_r=sc


  else:
    cff1=1/np.sinh(theta_s)
    cff2=0.5/np.tanh(0.5*theta_s)
    if type=='w':
      sc=(np.arange(0,N+1)-N)/N
      N=N+1
    else:
      sc=(np.arange(1,N+1)-N-0.5)/N
      #end
      Cs=(1.-theta_b)*cff1*np.sinh(theta_s*sc)+theta_b*(cff2*np.tanh(theta_s*(sc+0.5))-0.5)
  #
  # Create S-coordinate system: based on model topography h(i,j),
  # fast-time-averaged free-surface field and vertical coordinate
  # transformation metrics compute evolving depths of of the three-
  # dimensional model grid. Also adjust zeta for dry cells.
  #  
  h[h==0]=1.e-2
  Dcrit=0.01   # min water depth in dry cells
  zeta[zeta<(Dcrit-h)]=Dcrit-h[zeta<(Dcrit-h)]
  #
  hinv=1./h
  z=np.zeros((N,M,L))
  if (vtransform == 2):
    if type=='w':
      cff1=Cs_w
      cff2=sc_w+1
      sc=sc_w
    else:
      cff1=Cs_r
      cff2=sc_r+1
      sc=sc_r
      #end
    h2=(h+hc)
    cff=hc*sc
    h2inv=1./h2
    for k in range(N):
      z0=cff[k]+cff1[k]*h
      z[k,:,:]=z0*h/[h2] + zeta*(1.+z0*h2inv)
      #end
  else:
    cff1=Cs
    cff2=sc+1
    cff=hc*(sc-Cs)
    cff2=sc+1
    for k in range(N):
      z0=cff[k]+cff1[k]*h
      z[k,:,:]=z0+zeta*(1.+z0*hinv)

  return z

#####################################################################
#dist=spheric_dist(lat1,lat2,lon1,lon2)
#####################################################################
def spheric_dist(lat1,lat2,lon1,lon2):
    #
  # Earth radius
  #
  R=6367442.76
  #
  #  Determine proper longitudinal shift.
  #
  l=abs(lon2-lon1)
  l[l>=180]=360-l[l>=180]
  #                  
  #  Convert Decimal degrees to radians.
  #
  deg2rad=np.pi/180
  lat1=lat1*deg2rad
  lat2=lat2*deg2rad
  l=l*deg2rad
  #
  #  Compute the distances
  #
  dist=R*np.arcsin(np.sqrt(((np.sin(l)*np.cos(lat2))^2)+(((np.sin(lat2)*np.cos(lat1))\
    -(np.sin(lat1)*np.cos(lat2)*np.cos(l)))^2)))
  return dist


#####################################################################
#  extrfield = oainterp(londata,latdata,data,lon,lat,ro,savefile)
#####################################################################
def oainterp(londata,latdata,data,lon,lat,ro,savefile,*args):
  narg = oainterp.__code__.co_argcount
  if narg < 5:
    sys.exit('Not enough input arguments')
  elif narg < 6:
    print('using default decorrelation scale:  ro = 500 km')
    ro = 5e5
  #
  if narg < 7: 
    savefile=2
  #
  mdata=np.mean(data)
  data=data-mdata
  #
  if savefile !=0:
    invro=1/ro
    i=np.arrange(1,len(londata)+1)
    j=np.arange(1,len(lon)+1)
    I,J=np.meshgrid(i,i)
    r1=spheric_dist(latdata(I),latdata(J),londata(I),londata(J))
  #
    I,J=np.meshgrid(i,j)
    r2=spheric_dist(lat(J),latdata(I),lon(J),londata(I))
  #
    if savefile ==1:
      dill.dump_session('tmp.pkl')
    #
  elif savefile==0:
    dill.load_session('tmp.pkl')
  #
  #
  extrfield=mdata+(np.exp(-r2*invro)/np.exp(-r1*invro))*data
  #
  return extrfield

#####################################################################
#  function field=get_missing_val(lon,lat,field,missvalue,ro,default)
######################################################################
def get_missing_val(lon,lat,field,missvalue,ro,default,savefile=None,*args):
  interp_flag=1
  narg = get_missing_val.__code__.co_argcount
  if narg<4:
    oa_interp=0
    missvalue=float('nan')
    default=0
    ro=0
  elif narg<5:
    oa_interp=0
    default=0
    ro=0
  elif narg<6:
    default=0
  #
  #  get a masking matrix and the good data matrix
  #
  if np.isnan(missvalue):
    #ismask=np.isnan(field)
    ismask=np.zeros_like(field)
    ismask[np.isnan(field)] = 1
  else: 
    #ismask=(field==missvalue)
    ismask=np.zeros_like(field)
    ismask[ismask==missvalue] = 1
  #
  isdata=np.ones_like(ismask)-ismask

  M,L= field.shape
  #
  if sum(np.shape(lon))==(len(lon)+1):
    [lon,lat]=np.meshgrid(lon,lat)
  #
  # Check dimensions
  #
  if lon.shape[0] != field.shape[1]:
    sys.exit('GET-MISSING_VALUE: sizes do not correspond')
  #
  # test if there are any data
  #
  if ((np.sum(isdata))==0): 
    #  print('no data')
    print('GET_MISSING_VAL: No data point -> using default value:',str(default))
    field=np.zeros((M,L))+default
    interp_flag=0
  
  elif ((np.sum(isdata))<6): 
    #  default=min(field(isdata==1))
    default=np.mean(field[isdata==1])
    print('GET_MISSING_VAL: only ',str(sum(sum(isdata))),' data points: using the mean value:',str(default))
    field=np.zeros((M,L))+default
    interp_flag=0

  #
  if (np.sum(ismask))==0: 
    print('no mask')
    
  #---------------------------------------------------------------
  # Objective Analysis
  #---------------------------------------------------------------
  if ro>0:
    if ((np.sum(ismask)))==1: 
  #   print('1 mask')
      [j,i]=np.where(ismask==1)
      lat0=lat(j,i)
      lon0=lon(j,i)
      if j>1:
        od1=1./spheric_dist(lat0,lat(j-1,i),lon0,lon(j-1,i))
        f1=field(j-1,i)
      else:
        od1=0
        f1=0
      #
      if j<M:
        od2=1./spheric_dist(lat0,lat(j+1,i),lon0,lon(j+1,i))
        f2=field(j+1,i)
      else:
        od2=0
        f2=0
      #
      if i>1:
        od3=1./spheric_dist(lat0,lat(j,i-1),lon0,lon(j,i-1))
        f3=field(j,i-1)
      else:
        od3=0
        f3=0
      #
      if i<L:
        od4=1./spheric_dist(lat0,lat(j,i+1),lon0,lon(j,i+1))
        f4=field(j,i+1)
      else:
        od4=0
        f4=0
      #
      field[j,i]=(od1*f1+od2*f2+od3*f3+od4*f4)/(od1+od2+od3+od4)
      
    #
  #
  # perform the oa (if savefile is given, store the distances matrices 
  #                 (r1 and r2) in a tmp.mat file)
  #
    if narg == 7:
      field[ismask]=oainterp(lon[~ismask],lat[~ismask],field[~ismask],lon[ismask],lat[ismask],ro,savefile)
    else: 
      field[ismask]=oainterp(lon[~ismask],lat[~ismask],field[~ismask],lon[ismask],lat[ismask],ro,2)
    #
  else:
  #---------------------------------------------------------------
  # Extrapolation using nearest values
  #---------------------------------------------------------------
    lat_lon, coord = np.meshgrid(lon, lat)

    #field[ismask]=griddata((lon[~(ismask[0,:])].ravel(),lat[~(ismask[:,0])].ravel()),\
     #field[~ismask].ravel(),(lon[ismask[0,:]].ravel(), lat[ismask[:,0]].ravel()),'nearest')
    
    field[ismask==1]=griddata((lat_lon[ismask==0]),field[ismask==0]\
      ,(lat_lon[ismask==1]),'nearest')
  #
  return field, interp_flag

################################################################
#
#  function nc=create_inifile(inifile,gridfile,theta_s,...
#                  theta_b,hc,N,t([''],time)s([''],time)u([''],time)... 
#                  cycle,clobber)
#
#   This function create the header of a Netcdf climatology 
#   file.
#
#   Input: 
# 
#   inifile      Netcdf initial file name (character string).
#   gridfile     Netcdf grid file name (character string).
#   theta_s      S-coordinate surface control parameter.(Real)
#   theta_b      S-coordinate bottom control parameter.(Real)
#   hc           Width (m) of surface or bottom boundary layer
#                where higher vertical resolution is required 
#                during stretching.(Real)
#   N            Number of vertical levels.(Integer)  
#   time         Initial time.(Real) 
#   clobber      Switch to allow or not writing over an existing
#                file.(character string) 
#
#   Output
#
#   nc       Output netcdf object.
################################################################
def create_inifile(inifile,gridfile,title,theta_s,theta_b,hc,N,time,clobber,vtransform=None):
  print(' ')
  print(' Creating the file : ',inifile)
  if vtransform == None:
    print(' NO VTRANSFORM parameter found')
    print(' USE TRANSFORM default value vtransform = 1')
    vtransform = 1 
  print(' VTRANSFORM = ',str(vtransform))
  #
  #  Read the grid file
  #
  nc=Dataset(gridfile, set_auto_maskandscale=False)
  h=nc['h'][:].data  
  mask=nc['mask_rho'][:].data

  hmin=np.min(h[mask==1])
  if vtransform ==1:
    if hc > hmin:
      sys.exit(' hc ('+str(hc)+' m) > hmin ('+str(hmin)+' m)')
    

  Mp,Lp=h.shape[0],h.shape[1]
  L=Lp-1
  M=Mp-1
  Np=N+1
  nc.close()
    
  type = 'INITIAL file'  
  history = 'CROCO' 
  #
  #  Create the initial file
  #
  nc_ini = Dataset(inifile,  "w", format="NETCDF4", clobber=True)

  #
  #  Create dimensions
  #
  nc_ini.createDimension('xi_u', L)
  nc_ini.createDimension('xi_v', Lp)
  nc_ini.createDimension('xi_rho', Lp)
  nc_ini.createDimension('eta_u', Mp)
  nc_ini.createDimension('eta_v', M)
  nc_ini.createDimension('eta_rho', Mp)
  nc_ini.createDimension('s_rho', N)
  nc_ini.createDimension('s_w', Np)
  nc_ini.createDimension('tracer', 2)
  nc_ini.createDimension('time', 0)
  nc_ini.createDimension('one', 1)
  #
  # Create variables
  #
  spherical_nc = nc_ini.createVariable('spherical', 'S1', ('one')) 
  Vtransform_nc = nc_ini.createVariable('Vtransform', 'i', ('one'))
  Vstretching_nc = nc_ini.createVariable('Vstretching', 'i', ('one'))
  tstart_nc = nc_ini.createVariable('tstart', 'd', ('one'))
  tend_nc = nc_ini.createVariable('tend', 'd', ('one'))
  theta_s_nc = nc_ini.createVariable('theta_s', 'd', ('one'))
  theta_b_nc = nc_ini.createVariable('theta_b', 'd', ('one'))
  Tcline_nc = nc_ini.createVariable('Tcline', 'd', ('one')) 
  hc_nc = nc_ini.createVariable('hc', 'd', ('one')) 
  sc_r_nc = nc_ini.createVariable('sc_r', 'd', ('s_rho')) 
  Cs_r_nc = nc_ini.createVariable('Csr_r', 'd', ('s_rho')) 
  ocean_time_nc = nc_ini.createVariable('ocean_time', 'd', ('time')) 
  scrum_time_nc = nc_ini.createVariable('scrum_time', 'd', ('time')) 
  u_nc = nc_ini.createVariable('u', 'd', ('time','s_rho','eta_u','xi_u')) 
  v_nc = nc_ini.createVariable('v', 'd', ('time','s_rho','eta_v','xi_v')) 
  ubar_nc = nc_ini.createVariable('ubar', 'd', ('time','eta_u','xi_u')) 
  vbar_nc = nc_ini.createVariable('vbar', 'd', ('time','eta_v','xi_v'))  
  zeta_nc = nc_ini.createVariable('zeta', 'd', ('time','eta_rho','xi_rho')) 
  temp_nc = nc_ini.createVariable('temp', 'd', ('time','s_rho','eta_rho','xi_rho')) 
  salt_nc = nc_ini.createVariable('salt', 'd', ('time','s_rho','eta_rho','xi_rho'))
    #
    #  Attributes
    #
  Vtransform_nc.long_name = 'vertical terrain-following transformation equation'
    #
  Vstretching_nc.long_name = 'vertical terrain-following stretching func_inition'
    #
  tstart_nc.long_name = 'start processing day' 
  tstart_nc.units = 'day'
    #
  tend_nc.long_name = 'end processing day' 
  tend_nc.units= 'day'
    #
  theta_s_nc.long_name = 'S-coordinate surface control parameter' 
  theta_s_nc.units='nondimensional'
    #
  theta_b_nc.long_name = 'S-coordinate bottom control parameter' 
  theta_b_nc.units='nondimensional'
    #
  Tcline_nc.long_name = 'S-coordinate surface/bottom layer width' 
  Tcline_nc.units='meter'
    #
  hc_nc.long_name = 'S-coordinate parameter critical depth' 
  hc_nc.units='meter'
    #
  sc_r_nc.long_name = 'S-coordinate at RHO-points' 
  sc_r_nc.units='nondimensional'
  sc_r_nc.valid_min=-1 
  sc_r_nc.valid_max=0
    #
  Cs_r_nc.long_name = 'S-coordinate stretching curves at RHO-points' 
  Cs_r_nc.units='nondimensional'
  Cs_r_nc.valid_min=-1 
  Cs_r_nc.valid_max=0
    #
  ocean_time_nc.long_name = 'time sinc_inie initialization' 
  ocean_time_nc.units='second'
    #
  scrum_time_nc.long_name = 'time sinc_inie initialization' 
  scrum_time_nc.units='second'
    #
  u_nc.long_name = 'u-momentum component' 
  u_nc.units='meter second-1'
    #
  v_nc.long_name = 'v-momentum component' 
  v_nc.units='meter second-1'
    #
  ubar_nc.long_name = 'vertically integrated u-momentum component' 
  ubar_nc.units='meter second-1'
    #
  vbar_nc.long_name = 'vertically integrated v-momentum component' 
  vbar_nc.units='meter second-1'
    #
  zeta_nc.long_name = 'free-surface' 
  zeta_nc.units='meter'
    #
  temp_nc.long_name = 'potential temperature' 
  temp_nc.units='Celsius'
    #
  salt_nc.long_name = 'salinity' 
  salt_nc.units='PSU'
    #
    # Create global attributes
    #
  nc_ini.title = title
  nc_ini.date = datetime.today().strftime('%Y-%m-%d')
  nc_ini.clim_file = inifile 
  nc_ini.grd_file = gridfile, 
  nc_ini.type = type
  nc_ini.history = history
    #
    # Compute S coordinates
    #
  [sc_r,Cs_r,sc_w,Cs_w] = scoordinate(theta_s,theta_b,N,hc,vtransform)
    #
    # Write variables
    #
  spherical_nc[:]='T'
  Vtransform_nc[:]=vtransform
  Vstretching_nc[:]=1
  tstart_nc[:] =  time 
  tend_nc[:] =  time 
  theta_s_nc[:] =  theta_s 
  theta_b_nc[:] =  theta_b 
  Tcline_nc[:] =  hc 
  hc_nc[:] =  hc 
  sc_r_nc[:] =  sc_r 
  Cs_r_nc[:] =  Cs_r 
  scrum_time_nc[0] =  time*24*3600 
  ocean_time_nc[0] =  time*24*3600 
  u_nc[:] =  0 
  v_nc[:] =  0 
  zeta_nc[:] =  0 
  ubar_nc[:] =  0 
  vbar_nc[:] =  0 
  temp_nc[:] =  0 
  salt_nc[:] =  0 
  #
  # Leave define mode
  #
  nc_ini.close()
  print('inifile: -'+inifile+'- created')

  return


################################################################
#
# function create_bryfile(bryname,grdname,title,obc...
#                          theta_s,theta_b,hc,N,...
#                          ([''],time)cycle,clobber);
#
#   This function create the header of a Netcdf climatology 
#   file.
#
#   Input:
#
#   bryname      Netcdf climatology file name (character string).
#   grdname      Netcdf grid file name (character string).
#   obc          open boundaries flag (1=open , [S E N W]).
#   theta_s      S-coordinate surface control parameter.(Real)
#   theta_b      S-coordinate bottom control parameter.(Real)
#   hc           Width (m) of surface or bottom boundary layer 
#                where higher vertical resolution is required 
#                during stretching.(Real)
#   N            Number of vertical levels.(Integer)
#   time         time.(vector)
#   cycle        Length (days) for cycling the climatology.(Real)
#   clobber      Switch to allow or not writing over an existing
#                file.(character string)
################################################################
def create_bryfile(bryname,grdname,title,obc,theta_s,theta_b,hc,N,time,cycle,clobber,vtransform, *args):
  print(' ')
  print([' Creating the file : ',bryname])
  print(' ')
  narg = create_bryfile.__code__.co_argcount
  if narg < 12:
    print([' NO VTRANSFORM parameter found'])
    print([' USE TRANSFORM default value vtransform = 1'])
    vtransform = 1; 


  print([' VTRANSFORM = ',str(vtransform)])
  #
  #  Read the grid file and check the topography
  #
  nc=Dataset(grdname, set_auto_maskandscale=False)
  h=nc['h'][:].data  
  maskr=nc['mask_rho'][:].data

  Lp=nc.dimensions['xi_rho'].size
  Mp=nc.dimensions['eta_rho'].size
  nc.close()
  hmin=np.min(h[maskr==1])
    
  if vtransform ==1:
    if hc > hmin:
      sys.exit(' hc (',str(hc),' m) > hmin (',str(hmin),' m)')
        #end
  #end
  L=Lp-1
  M=Mp-1
  Np=N+1
  #
  #  Create the boundary file
  #
  type = 'BOUNDARY file' 
  history = 'CROCO'

  nc_bry = Dataset(bryname, 'w', format="NETCDF4", clobber=True)#inifile,clobber)
  #
  #  Create dimensions
  #
  nc_bry.createDimension('xi_u', L)
  nc_bry.createDimension('xi_v', Lp)
  nc_bry.createDimension('xi_rho', Lp)
  nc_bry.createDimension('eta_u', Mp)
  nc_bry.createDimension('eta_v', M)
  nc_bry.createDimension('eta_rho', Mp)
  nc_bry.createDimension('s_rho', N)
  nc_bry.createDimension('s_w', Np)
  nc_bry.createDimension('tracer', 2)
  nc_bry.createDimension('bry_time', len(time))
  nc_bry.createDimension('tclm_time', len(time))
  nc_bry.createDimension('temp_time', len(time))
  nc_bry.createDimension('sclm_time', len(time))
  nc_bry.createDimension('salt_time', len(time))
  nc_bry.createDimension('uclm_time', len(time))
  nc_bry.createDimension('vclm_time', len(time))
  nc_bry.createDimension('v2d_time', len(time))
  nc_bry.createDimension('v3d_time', len(time))
  nc_bry.createDimension('ssh_time', len(time))
  nc_bry.createDimension('zeta_time', len(time))
  nc_bry.createDimension('one', 1)
  #
  #  Create variables and attributes
  #
  spherical_nc = nc_bry.createVariable('spherical', 'c', ('one')) 
  spherical_nc.long_name = 'grid type logical switch'
  spherical_nc.flag_values = 'T, F'
  spherical_nc.flag_meanings = 'spherical Cartesian'

  Vtransform_nc = nc_bry.createVariable('Vtransform', 'i', ('one'))
  Vtransform_nc.long_name = 'vertical terrain-following transformation equation'

  Vstretching_nc = nc_bry.createVariable('Vstreching', 'i', ('one'))
  Vstretching_nc.long_name = 'vertical terrain-following stretching function'

  tstart_nc = nc_bry.createVariable('tstart', 'd', ('one'))
  tstart_nc.long_name = 'start processing day'
  tstart_nc.units = 'day'

  tend_nc = nc_bry.createVariable('tend', 'd', ('one'))
  tend_nc.long_name = 'end processing day'
  tend_nc.units = 'day'

  theta_s_nc = nc_bry.createVariable('theta_s', 'd', ('one'))
  theta_s_nc.long_name = 'S-coordinate surface control parameter'
  theta_s_nc.units = 'nondimensional'

  theta_b_nc = nc_bry.createVariable('theta_b', 'd', ('one'))
  theta_b_nc.long_name = 'S-coordinate bottom control parameter'
  theta_b_nc.units = 'nondimensional'

  Tcline_nc = nc_bry.createVariable('Tcline', 'd',('one'))
  Tcline_nc.long_name = 'S-coordinate surface/bottom layer width'
  Tcline_nc.units = 'meter'

  hc_nc = nc_bry.createVariable('hc', 'd',('one'))
  hc_nc.long_name = 'S-coordinate parameter, critical depth'
  hc_nc.unitd = 'meter' 

  sc_r_nc = nc_bry.createVariable('sc_r', 'd', ('s_rho'))
  sc_r_nc.long_name = 'S-coordinate at RHO-points'
  sc_r_nc.valid_min = -1.
  sc_r_nc.valid_max = 0.
  sc_r_nc.positive = 'up'
  if vtransform==1:
    sc_r_nc.standard_name = 'ocena_s_coordinate_g1'
  elif vtransform==2:
    sc_r_nc.standard_name = 'ocena_s_coordinate_g2'
  sc_r_nc.formula_terms = 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc'

  sc_w_nc = nc_bry.createVariable('sc_w', 'd', ('s_w'))
  sc_w_nc.long_name = 'S-coordinate at W-points'  
  sc_w_nc.valid_min = -1.   
  sc_w_nc.valid_max = 0.  
  sc_w_nc.positive = 'up'
  if vtransform == 1:
    sc_w_nc.standard_name = 'ocean_s_coordinate_g1'
  elif vtransform == 2:
    sc_w_nc.standard_name = 'ocean_s_coordinate_g2' 
  sc_w_nc.formula_terms = 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc'

  Cs_r_nc = nc_bry.createVariable('Cs_r', 'd', ('s_rho'))
  Cs_r_nc.long_name = 'S-coordinate stretching curves at RHO-points'
  Cs_r_nc.units = 'nondimensional'
  Cs_r_nc.valid_min = -1
  Cs_r_nc.valid_max = 0

  Cs_w_nc = nc_bry.createVariable('Cs_w', 'd', ('s_w'))
  Cs_w_nc.long_name = 'S-coordinate stretching curves at W-points'
  Cs_w_nc.units = 'nondimensional'
  Cs_w_nc.valid_min = -1
  Cs_w_nc.valid_max = 0

  bry_time_nc = nc_bry.createVariable('bry_time', 'd', ('bry_time')) 
  bry_time_nc.long_name = 'time for boundary climatology'
  bry_time_nc.units = 'day'
  bry_time_nc.calendar = '360.0 days in every year'
  bry_time_nc.cycle_length = cycle
  #
  tclm_time_nc = nc_bry.createVariable('tclm_time', 'd', ('tclm_time')) 
  tclm_time_nc.long_name = 'time for temperature climatology'
  tclm_time_nc.units = 'day'
  tclm_time_nc.calendar = '360.0 days in every year'
  tclm_time_nc.cycle_length = cycle
  #
  temp_time_nc = nc_bry.createVariable('temp_time', 'd', ('temp_time')) 
  temp_time_nc.long_name = 'time for temperature climatology'
  temp_time_nc.units = 'day'
  temp_time_nc.calendar = '360.0 days in every year'
  temp_time_nc.cycle_length = cycle
  #
  sclm_time_nc = nc_bry.createVariable('sclm_time', 'd', ('sclm_time')) 
  sclm_time_nc.long_name = 'time for salinity climatology'
  sclm_time_nc.units = 'day'
  sclm_time_nc.calendar = '360.0 days in every year'
  sclm_time_nc.cycle_length = cycle
  #
  salt_time_nc = nc_bry.createVariable('salt_time', 'd', ('salt_time')) 
  salt_time_nc.long_name = 'time for salinity climatology'
  salt_time_nc.units = 'day'
  salt_time_nc.calendar = '360.0 days in every year'
  salt_time_nc.cycle_length = cycle
  #
  uclm_time_nc = nc_bry.createVariable('uclm_time', 'd', ('uclm_time')) 
  uclm_time_nc.long_name = 'time climatological u'
  uclm_time_nc.units = 'day'
  uclm_time_nc.calendar = '360.0 days in every year'
  uclm_time_nc.cycle_length = cycle
  #
  vclm_time_nc = nc_bry.createVariable('vclm_time', 'd', 'vclm_time') 
  vclm_time_nc.long_name = 'time climatological v'
  vclm_time_nc.units = 'day'
  vclm_time_nc.calendar = '360.0 days in every year'
  vclm_time_nc.cycle_length = cycle
  #
  v2d_time_nc = nc_bry.createVariable('v2d_time', 'd', ('v2d_time')) 
  v2d_time_nc.long_name = 'time for 2D velocity climatology'
  v2d_time_nc.units = 'day'
  v2d_time_nc.calendar = '360.0 days in every year'
  v2d_time_nc.cycle_length = cycle
  #
  v3d_time_nc = nc_bry.createVariable('v3d_time', 'd', ('v3d_time')) 
  v3d_time_nc.long_name = 'time for 3D velocity climatology'
  v3d_time_nc.units = 'day'
  v3d_time_nc.calendar = '360.0 days in every year'
  v3d_time_nc.cycle_length = cycle
  #
  ssh_time_nc = nc_bry.createVariable('ssh_time', 'd', ('ssh_time')) 
  ssh_time_nc.long_name = 'time for sea surface height'
  ssh_time_nc.units = 'day'
  ssh_time_nc.calendar = '360.0 days in every year'
  ssh_time_nc.cycle_length = cycle
  #
  zeta_time_nc = nc_bry.createVariable('zeta_time', 'd', ('zeta_time')) 
  zeta_time_nc.long_name = 'time for sea surface height'
  zeta_time_nc.units = 'day'
  zeta_time_nc.calendar = '360.0 days in every year'
  zeta_time_nc.cycle_length = cycle      

  if obc[0]==1:
  #
  #   Southern boundary
  #
    temp_south_nc = nc_bry.createVariable('temp_south', 'd',('temp_time','s_rho','xi_rho')) 
    temp_south_nc.long_name = 'southern boundary potential temperature'
    temp_south_nc.units = 'Celsius'
    temp_south_nc.coordinates = 'lon_rho s_rho temp_time'
  #
    salt_south_nc = nc_bry.createVariable('salt_south', 'd',('salt_time','s_rho','xi_rho')) 
    salt_south_nc.long_name = 'southern boundary salinity'
    salt_south_nc.units = 'PSU'
    salt_south_nc.coordinates = 'lon_rho s_rho salt_time'
  #
    u_south_nc = nc_bry.createVariable('u_south', 'd',('v3d_time','s_rho','xi_u'))
    u_south_nc.long_name = 'southern boundary u-momentum component'
    u_south_nc.units = 'meter second-1'
    u_south_nc.coordinates = 'lon_u s_rho u_time'
  #
    v_south_nc = nc_bry.createVariable('v_south', 'd',('v3d_time','s_rho','xi_rho')) 
    v_south_nc.long_name = 'southern boundary v-momentum component'
    v_south_nc.units = 'meter second-1'
    v_south_nc.coordinates = 'lon_v s_rho vclm_time'
  #
    ubar_south_nc = nc_bry.createVariable('ubar_south', 'd',('v2d_time','xi_u')) 
    ubar_south_nc.long_name = 'southern boundary vertically integrated u-momentum component'
    ubar_south_nc.units = 'meter second-1'
    ubar_south_nc.coordinates = 'lon_u uclm_time'
  #
    vbar_south_nc = nc_bry.createVariable('vbar_south', 'd',('v2d_time','xi_rho')) 
    vbar_south_nc.long_name = 'southern boundary vertically integrated v-momentum component'
    vbar_south_nc.units = 'meter second-1'
    vbar_south_nc.coordinates = 'lon_v vclm_time'
  #
    zeta_south_nc = nc_bry.createVariable('zeta_south', 'd',('zeta_time','xi_rho')) 
    zeta_south_nc.long_name = 'southern boundary sea surface height'
    zeta_south_nc.units = 'meter'
    zeta_south_nc.coordinates = 'lon_rho zeta_time'

  if obc[1]==1:
  #
  #   Eastern boundary
  #
    temp_east_nc = nc_bry.createVariable('temp_east', 'd',('temp_time','s_rho','eta_rho')) 
    temp_east_nc.long_name = 'eastern boundary potential temperature'
    temp_east_nc.units = 'Celsius'
    temp_east_nc.coordinates = 'lat_rho s_rho temp_time'
  #
    salt_east_nc = nc_bry.createVariable('salt_east', 'd',('salt_time','s_rho','eta_rho')) 
    salt_east_nc.long_name = 'eastern boundary salinity'
    salt_east_nc.units = 'PSU'
    salt_east_nc.coordinates = 'lat_rho s_rho salt_time'
  #
    u_east_nc = nc_bry.createVariable('u_east', 'd',('v3d_time','s_rho','eta_rho')) 
    u_east_nc.long_name = 'eastern boundary u-momentum component'
    u_east_nc.units = 'meter second-1'
    u_east_nc.coordinates = 'lat_u s_rho u_time'
  #
    v_east_nc = nc_bry.createVariable('v_east', 'd',('v3d_time','s_rho','eta_v')) 
    v_east_nc.long_name = 'eastern boundary v-momentum component'
    v_east_nc.units = 'meter second-1'
    v_east_nc.coordinates = 'lat_v s_rho vclm_time'
  #
    ubar_east_nc = nc_bry.createVariable('ubar_east', 'd',('v2d_time','eta_rho')) 
    ubar_east_nc.long_name = 'eastern boundary vertically integrated u-momentum component'
    ubar_east_nc.units = 'meter second-1'
    ubar_east_nc.coordinates = 'lat_u uclm_time'
  #
    vbar_east_nc = nc_bry.createVariable('vbar_east', 'd',('v2d_time','eta_v')) 
    vbar_east_nc.long_name = 'eastern boundary vertically integrated v-momentum component'
    vbar_east_nc.units = 'meter second-1'
    vbar_east_nc.coordinates = 'lat_v vclm_time'
  #
    zeta_east_nc = nc_bry.createVariable('zeta_east', 'd',('zeta_time','eta_rho')) 
    zeta_east_nc.long_name = 'eastern boundary sea surface height'
    zeta_east_nc.units = 'meter'
    zeta_east_nc.coordinates = 'lat_rho zeta_time'

  if obc[2]==1:
  #
  #   Northern boundary
  #
    temp_north_nc = nc_bry.createVariable('temp_north', 'd',('temp_time','s_rho','xi_rho')) 
    temp_north_nc.long_name = 'northern boundary potential temperature'
    temp_north_nc.units = 'Celsius'
    temp_north_nc.coordinates = 'lon_rho s_rho temp_time'
  #
    salt_north_nc = nc_bry.createVariable('salt_north', 'd',('salt_time','s_rho','xi_rho')) 
    salt_north_nc.long_name = 'northern boundary salinity'
    salt_north_nc.units = 'PSU'
    salt_north_nc.coordinates = 'lon_rho s_rho salt_time'
  #
    u_north_nc = nc_bry.createVariable('u_north', 'd',('v3d_time','s_rho','xi_u')) 
    u_north_nc.long_name = 'northern boundary u-momentum component'
    u_north_nc.units = 'meter second-1'
    u_north_nc.coordinates = 'lon_u s_rho u_time'
  #
    v_north_nc = nc_bry.createVariable('v_north', 'd',('v3d_time','s_rho','xi_rho')) 
    v_north_nc.long_name = 'northern boundary v-momentum component'
    v_north_nc.units = 'meter second-1'
    v_north_nc.coordinates = 'lon_v s_rho vclm_time'
  #
    ubar_north_nc = nc_bry.createVariable('ubar_north', 'd',('v2d_time','xi_u')) 
    ubar_north_nc.long_name = 'northern boundary vertically integrated u-momentum component'
    ubar_north_nc.units = 'meter second-1'
    ubar_north_nc.coordinates = 'lon_u uclm_time'
  #
    vbar_north_nc = nc_bry.createVariable('vbar_north', 'd',('v2d_time','xi_rho')) 
    vbar_north_nc.long_name = 'northern boundary vertically integrated v-momentum component'
    vbar_north_nc.units = 'meter second-1'
    vbar_north_nc.coordinates = 'lon_v vclm_time'

    zeta_north_nc = nc_bry.createVariable('zeta_north', 'd',('zeta_time','xi_rho')) 
    zeta_north_nc.long_name = 'northern boundary sea surface height'
    zeta_north_nc.units = 'meter'
    zeta_north_nc.coordinates = 'lon_rho zeta_time'

  if obc[3]==1:
  #
  #   Western boundary
  #
    temp_west_nc = nc_bry.createVariable('temp_west', 'd',('temp_time','s_rho','eta_rho')) 
    temp_west_nc.long_name = 'western boundary potential temperature'
    temp_west_nc.units = 'Celsius'
    temp_west_nc.coordinates = 'lat_rho s_rho temp_time'
  #
    salt_west_nc = nc_bry.createVariable('salt_west', 'd',('salt_time','s_rho','eta_rho')) 
    salt_west_nc.long_name = 'western boundary salinity'
    salt_west_nc.units = 'PSU'
    salt_west_nc.coordinates = 'lat_rho s_rho salt_time'
  #
    u_west_nc = nc_bry.createVariable('u_west', 'd',('v3d_time','s_rho','eta_rho')) 
    u_west_nc.long_name = 'western boundary u-momentum component'
    u_west_nc.units = 'meter second-1'
    u_west_nc.coordinates = 'lat_u s_rho u_time'
  #
    v_west_nc = nc_bry.createVariable('v_west', 'd',('v3d_time','s_rho','eta_v')) 
    v_west_nc.long_name = 'western boundary v-momentum component'
    v_west_nc.units = 'meter second-1'
    v_west_nc.coordinates = 'lat_v s_rho vclm_time'
  #
    ubar_west_nc = nc_bry.createVariable('ubar_west', 'd',('v2d_time','eta_rho')) 
    ubar_west_nc.long_name = 'western boundary vertically integrated u-momentum component'
    ubar_west_nc.units = 'meter second-1'
    ubar_west_nc.coordinates = 'lat_u uclm_time'
  #
    vbar_west_nc = nc_bry.createVariable('vbar_west', 'd',('v2d_time','eta_v')) 
    vbar_west_nc.long_name = 'western boundary vertically integrated v-momentum component'
    vbar_west_nc.units = 'meter second-1'
    vbar_west_nc.coordinates = 'lat_v vclm_time'
  #
    zeta_west_nc = nc_bry.createVariable('zeta_west', 'd',('zeta_time','eta_rho')) 
    zeta_west_nc.long_name = 'western boundary sea surface height'
    zeta_west_nc.units = 'meter'
    zeta_west_nc.coordinates = 'lat_rho zeta_time'
      
  #
  # Create global attributes
  #
  nc_bry.title = title
  nc_bry.date = datetime.today().strftime('%Y-%m-%d') 
  nc_bry.clim_file = bryname 
  nc_bry.grd_file = grdname 
  nc_bry.type = type 
  nc_bry.history = history
  #
  # Compute S coordinates
  #
  [sc_r,Cs_r,sc_w,Cs_w] = scoordinate(theta_s,theta_b,N,hc,vtransform)

    #
  # Write variables
  #
  spherical_nc[:]='T'
  Vtransform_nc[:]=vtransform
  Vstretching_nc[:]=1
  tstart_nc[:] =  np.min(time)
  tend_nc[:] =  np.max(time)
  theta_s_nc[:] =  theta_s 
  theta_b_nc[:] =  theta_b 
  Tcline_nc[:] =  hc 
  hc_nc[:] =  hc 
  sc_r_nc[:] = sc_r
  sc_w_nc[:] = sc_w
  Cs_r_nc[:] = Cs_r  
  Cs_w_nc[:] = Cs_w
  tclm_time_nc[:] =  time 
  temp_time_nc[:] =  time 
  sclm_time_nc[:] =  time 
  salt_time_nc[:] =  time 
  uclm_time_nc[:] =  time 
  vclm_time_nc[:] =  time 
  v2d_time_nc[:] =   time 
  v3d_time_nc[:] =   time 
  ssh_time_nc[:] =   time
  zeta_time_nc[:] =  time
  bry_time_nc[:] =  time 
  if obc[0]==1:
    u_south_nc[:] =  0 
    v_south_nc[:] =  0 
    ubar_south_nc[:] =  0 
    vbar_south_nc[:] =  0 
    zeta_south_nc[:] =  0 
    temp_south_nc[:] =  0 
    salt_south_nc[:] =  0
  if obc[1]==0:
    u_east_nc[:] =  0 
    v_east_nc[:] =  0 
    ubar_east_nc[:] =  0 
    vbar_east_nc[:] =  0 
    zeta_east_nc[:] =  0 
    temp_east_nc[:] =  0 
    salt_east_nc[:] =  0
  if obc[2]==1:
    u_north_nc[:] =  0 
    v_north_nc[:] =  0 
    ubar_north_nc[:] =  0 
    vbar_north_nc[:] =  0 
    zeta_north_nc[:] =  0 
    temp_north_nc[:] =  0 
    salt_north_nc[:] =  0
  if obc[3]==1:
    u_west_nc[:] =  0 
    v_west_nc[:] =  0 
    ubar_west_nc[:] =  0 
    vbar_west_nc[:] =  0 
    zeta_west_nc[:] =  0 
    temp_west_nc[:] =  0 
    salt_west_nc[:] =  0
    
  #
  # Leave define mode
  #
  nc_bry.close()

  return


################################################################
#
# function create_climfile(clmname,grdname,title,...
#                          theta_s,theta_b,hc,N,...
#                          time,cycle,clobber);
#
#   This function create the header of a Netcdf climatology 
#   file.
#
#   Input:
#
#   clmname      Netcdf climatology file name (character string).
#   grdname      Netcdf grid file name (character string).
#   theta_s      S-coordinate surface control parameter.(Real) 
#   theta_b      S-coordinate bottom control parameter.(Real)
#   hc           Width (m) of surface or bottom boundary layer
#                where higher vertical resolution is required
#                during stretching.(Real) 
#   N            Number of vertical levels.(Integer) 
#   time        Temperature climatology time.(vector) 
#   time        Salinity climatology time.(vector)
#   time        Velocity climatology time.(vector)
#   cycle        Length (days) for cycling the climatology.(Real)
#   clobber      Switch to allow or not writing over an existing 
#                file.(character string)
################################################################
def create_climfile(clmname,grdname,title,theta_s,theta_b,hc,N,time,cycle,clobber,vtransform,*args):
  print(' ')
  print(' Creating the file : ',clmname)
  print(' ')
  narg = create_climfile.__code__.co_argcount
  if narg < 11:
    print([' NO VTRANSFORM parameter found'])
    print([' USE TRANSFORM default value vtransform = 1'])
    vtransform = 1 
    
  print([' VTRANSFORM = ', str(vtransform)])
  #
  #  Read the grid file
  #
  nc = Dataset(grdname, 'r', set_auto_maskandscale=False)
  h=nc['h'][:].data
  maskr=nc['mask_rho'][:].data
  Lp=nc.dimensions['xi_rho'].size
  Mp=nc.dimensions['eta_rho'].size

  nc.close()
    
  hmin=np.min(h[maskr==1])
  if vtransform ==1:
    if hc > hmin:
      sys.exit(' hc (',str(hc),' m) > hmin (',str(hmin),' m)')
    
  L=Lp-1
  M=Mp-1
  Np=N+1
  #
  #  Create the climatology file
  #
  type = 'CLIMATOLOGY file'  
  history = 'CROCO'

  nc_clm = Dataset(clmname, 'w', format="NETCDF4", clobber=True)

  #
  #  Create dimensions
  #
  nc_clm.createDimension('xi_u', L)
  nc_clm.createDimension('xi_v', Lp)
  nc_clm.createDimension('xi_rho', Lp)
  nc_clm.createDimension('eta_u', Mp)
  nc_clm.createDimension('eta_v', M)
  nc_clm.createDimension('eta_rho', Mp)
  nc_clm.createDimension('s_rho', N)
  nc_clm.createDimension('s_w', Np)
  nc_clm.createDimension('tracer', 2)
  nc_clm.createDimension('bry_time', len(time))
  nc_clm.createDimension('tclm_time', len(time))
  nc_clm.createDimension('temp_time', len(time))
  nc_clm.createDimension('sclm_time', len(time))
  nc_clm.createDimension('salt_time', len(time))
  nc_clm.createDimension('uclm_time', len(time))
  nc_clm.createDimension('vclm_time', len(time))
  nc_clm.createDimension('v2d_time', len(time))
  nc_clm.createDimension('v3d_time', len(time))
  nc_clm.createDimension('ssh_time', len(time))
  nc_clm.createDimension('zeta_time', len(time))
  nc_clm.createDimension('one', 1)

  #
  #  Create variables and attributes
  #
  spherical_nc = nc_clm.createVariable('spherical', 'c', ('one')) 
  spherical_nc.long_name = 'grid type logical switch'
  spherical_nc.flag_values = 'T, F'
  spherical_nc.flag_meanings = 'spherical Cartesian'

  Vtransform_nc = nc_clm.createVariable('Vtransform', 'i', ('one'))
  Vtransform_nc.long_name = 'vertical terrain-following transformation equation'

  Vstretching_nc = nc_clm.createVariable('Vstreching', 'i', ('one'))
  Vstretching_nc.long_name = 'vertical terrain-following stretching function'

  tstart_nc = nc_clm.createVariable('tstart', 'd', ('one'))
  tstart_nc.long_name = 'start processing day'
  tstart_nc.units = 'day'

  tend_nc = nc_clm.createVariable('tend', 'd', ('one'))
  tend_nc.long_name = 'end processing day'
  tend_nc.units = 'day'

  theta_s_nc = nc_clm.createVariable('theta_s', 'd', ('one'))
  theta_s_nc.long_name = 'S-coordinate surface control parameter'
  theta_s_nc.units = 'nondimensional'

  theta_b_nc = nc_clm.createVariable('theta_b', 'd', ('one'))
  theta_b_nc.long_name = 'S-coordinate bottom control parameter'
  theta_b_nc.units = 'nondimensional'

  Tcline_nc = nc_clm.createVariable('Tcline', 'd',('one'))
  Tcline_nc.long_name = 'S-coordinate surface/bottom layer width'
  Tcline_nc.units = 'meter'

  hc_nc = nc_clm.createVariable('hc', 'd',('one'))
  hc_nc.long_name = 'S-coordinate parameter, critical depth'
  hc_nc.unitd = 'meter' 

  sc_r_nc = nc_clm.createVariable('sc_r', 'd', ('s_rho'))
  sc_r_nc.long_name = 'S-coordinate at RHO-points'
  sc_r_nc.valid_min = -1.
  sc_r_nc.valid_max = 0.
  sc_r_nc.positive = 'up'
  if vtransform==1:
    sc_r_nc.standard_name = 'ocena_s_coordinate_g1'
  elif vtransform==2:
    sc_r_nc.standard_name = 'ocena_s_coordinate_g2'
  sc_r_nc.formula_terms = 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc'

  sc_w_nc = nc_clm.createVariable('sc_w', 'd', ('s_w'))
  sc_w_nc.long_name = 'S-coordinate at W-points'  
  sc_w_nc.valid_min = -1.   
  sc_w_nc.valid_max = 0.  
  sc_w_nc.positive = 'up'
  if vtransform == 1:
    sc_w_nc.standard_name = 'ocean_s_coordinate_g1'
  elif vtransform == 2:
    sc_w_nc.standard_name = 'ocean_s_coordinate_g2' 
  sc_w_nc.formula_terms = 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc'

  Cs_r_nc = nc_clm.createVariable('Cs_r', 'd', ('s_rho'))
  Cs_r_nc.long_name = 'S-coordinate stretching curves at RHO-points'
  Cs_r_nc.units = 'nondimensional'
  Cs_r_nc.valid_min = -1
  Cs_r_nc.valid_max = 0

  Cs_w_nc = nc_clm.createVariable('Cs_w', 'd', ('s_w'))
  Cs_w_nc.long_name = 'S-coordinate stretching curves at W-points'
  Cs_w_nc.units = 'nondimensional'
  Cs_w_nc.valid_min = -1
  Cs_w_nc.valid_max = 0

  bry_time_nc = nc_clm.createVariable('bry_time', 'd', ('bry_time')) 
  bry_time_nc.long_name = 'time for boundary climatology'
  bry_time_nc.units = 'day'
  bry_time_nc.calendar = '360.0 days in every year'
  bry_time_nc.cycle_length = cycle
  #
  tclm_time_nc = nc_clm.createVariable('tclm_time', 'd', ('tclm_time')) 
  tclm_time_nc.long_name = 'time for temperature climatology'
  tclm_time_nc.units = 'day'
  tclm_time_nc.calendar = '360.0 days in every year'
  tclm_time_nc.cycle_length = cycle
  #
  temp_time_nc = nc_clm.createVariable('temp_time', 'd', ('temp_time')) 
  temp_time_nc.long_name = 'time for temperature climatology'
  temp_time_nc.units = 'day'
  temp_time_nc.calendar = '360.0 days in every year'
  temp_time_nc.cycle_length = cycle
  #
  sclm_time_nc = nc_clm.createVariable('sclm_time', 'd', ('sclm_time')) 
  sclm_time_nc.long_name = 'time for salinity climatology'
  sclm_time_nc.units = 'day'
  sclm_time_nc.calendar = '360.0 days in every year'
  sclm_time_nc.cycle_length = cycle
  #
  salt_time_nc = nc_clm.createVariable('salt_time', 'd', ('salt_time')) 
  salt_time_nc.long_name = 'time for salinity climatology'
  salt_time_nc.units = 'day'
  salt_time_nc.calendar = '360.0 days in every year'
  salt_time_nc.cycle_length = cycle
  #
  uclm_time_nc = nc_clm.createVariable('uclm_time', 'd', ('uclm_time')) 
  uclm_time_nc.long_name = 'time climatological u'
  uclm_time_nc.units = 'day'
  uclm_time_nc.calendar = '360.0 days in every year'
  uclm_time_nc.cycle_length = cycle
  #
  vclm_time_nc = nc_clm.createVariable('vclm_time', 'd', 'vclm_time') 
  vclm_time_nc.long_name = 'time climatological v'
  vclm_time_nc.units = 'day'
  vclm_time_nc.calendar = '360.0 days in every year'
  vclm_time_nc.cycle_length = cycle
  #
  v2d_time_nc = nc_clm.createVariable('v2d_time', 'd', ('v2d_time')) 
  v2d_time_nc.long_name = 'time for 2D velocity climatology'
  v2d_time_nc.units = 'day'
  v2d_time_nc.calendar = '360.0 days in every year'
  v2d_time_nc.cycle_length = cycle
  #
  v3d_time_nc = nc_clm.createVariable('v3d_time', 'd', ('v3d_time')) 
  v3d_time_nc.long_name = 'time for 3D velocity climatology'
  v3d_time_nc.units = 'day'
  v3d_time_nc.calendar = '360.0 days in every year'
  v3d_time_nc.cycle_length = cycle
  #
  ssh_time_nc = nc_clm.createVariable('ssh_time', 'd', ('ssh_time')) 
  ssh_time_nc.long_name = 'time for sea surface height'
  ssh_time_nc.units = 'day'
  ssh_time_nc.calendar = '360.0 days in every year'
  ssh_time_nc.cycle_length = cycle
  #
  zeta_time_nc = nc_clm.createVariable('zeta_time', 'd', ('zeta_time')) 
  zeta_time_nc.long_name = 'time for sea surface height'
  zeta_time_nc.units = 'day'
  zeta_time_nc.calendar = '360.0 days in every year'
  zeta_time_nc.cycle_length = cycle

  temp_nc = nc_clm.createVariable('temp', 'd', ('tclm_time','s_rho','eta_rho','xi_rho')) 
  temp_nc.long_name = 'potential temperature'
  temp_nc.units = 'Celsius'
  temp_nc.time = 'temp_time'
  temp_nc.coordinates = 'lon_rho lat_rho s_rho temp_time'
  #
  salt_nc = nc_clm.createVariable('salt', 'd',('sclm_time','s_rho','eta_rho','xi_rho')) 
  salt_nc.long_name = 'salinity'
  salt_nc.units = 'PSU'
  salt_nc.time = 'salt_time'
  salt_nc.coordinates = 'lon_rho lat_rho s_rho salt_time'
  #
  u_nc = nc_clm.createVariable('u', 'd',('uclm_time','s_rho','eta_u','xi_u')) 
  u_nc.long_name = 'u-momentum component'
  u_nc.units = 'meter second-1'
  u_nc.time = 'uclm_time'
  u_nc.coordinates = 'lon_u lat_u s_rho u_time'
  #
  v_nc = nc_clm.createVariable('v', 'd',('vclm_time','s_rho','eta_v','xi_v')) 
  v_nc.long_name = 'v-momentum component'
  v_nc.units = 'meter second-1'
  v_nc.time = 'vclm_time'
  v_nc.coordinates = 'lon_v lat_v s_rho vclm_time'
  #
  ubar_nc = nc_clm.createVariable('ubar', 'd',('uclm_time','eta_u','xi_u')) 
  ubar_nc.long_name = 'vertically integrated u-momentum component'
  ubar_nc.units = 'meter second-1'
  ubar_nc.time = 'uclm_time'
  ubar_nc.coordinates = 'lon_v lat_u uclm_time'
  #
  vbar_nc = nc_clm.createVariable('vbar', 'd',('vclm_time','eta_v','xi_v')) 
  vbar_nc.long_name = 'vertically integrated v-momentum component'
  vbar_nc.units = 'meter second-1'
  vbar_nc.time = 'vclm_time'
  vbar_nc.coordinates = 'lon_v lat_v vclm_time'
  #
  SSH_nc = nc_clm.createVariable('SSH', 'd', ('ssh_time','eta_rho','xi_rho'))
  SSH_nc.long_name = 'sea surface height'
  SSH_nc.units = 'meter'
  SSH_nc.time = 'zeta_time'
  SSH_nc.coordinates = 'lon_rho lat_rho zeta_time'
  #
  zeta_nc = nc_clm.createVariable('zeta', 'd',('zeta_time','eta_rho', 'xi_rho')) 
  zeta_nc.long_name = 'sea surface height'
  zeta_nc.units = 'meter'
  zeta_nc.time = 'zeta_time'
  zeta_nc.coordinates = 'lon_rho lat_rho zeta_time'   

  nc_clm.title = title
  nc_clm.date = datetime.today().strftime('%Y-%m-%d') 
  nc_clm.clim_file = clmname 
  nc_clm.grd_file = grdname 
  nc_clm.type = type 
  nc_clm.history = history   
  
  [sc_r,Cs_r,sc_w,Cs_w] = scoordinate(theta_s,theta_b,N,hc,vtransform)
    #
    #  Write variables
    #
  spherical_nc[:]='T'
  Vtransform_nc[:]=vtransform
  Vstretching_nc[:]=1
  tstart_nc[:] =  np.min(time)
  tend_nc[:] =  np.max(time)
  theta_s_nc[:] =  theta_s 
  theta_b_nc[:] =  theta_b 
  Tcline_nc[:] =  hc 
  hc_nc[:] =  hc 
  sc_r_nc[:] = sc_r
  sc_w_nc[:] = sc_w
  Cs_r_nc[:] =  Cs_r 
  Cs_w_nc[:] = Cs_w
  tclm_time_nc[:] =  time 
  temp_time_nc[:] =  time 
  sclm_time_nc[:] =  time 
  salt_time_nc[:] =  time 
  uclm_time_nc[:] =  time 
  vclm_time_nc[:] =  time 
  v2d_time_nc[:] =   time 
  v3d_time_nc[:] =   time 
  ssh_time_nc[:] =   time
  zeta_time_nc[:] =  time
  u_nc[:] =  0 
  v_nc[:] =  0 
  ubar_nc[:] =  0 
  vbar_nc[:] =  0 
  SSH_nc[:] =  0 
  zeta_nc[:] =  0 
  temp_nc[:] =  0 
  salt_nc[:] =  0 

  nc_clm.close()
  
  return nc_clm

################################################################
# 	Create an empty netcdf forcing file
#       frcname: name of the forcing file
#       grdname: name of the grid file
#       title: title in the netcdf file  
################################################################
def create_forcing(frcname,grdname,title,smst,shft,swft,srft,sstt,ssst,smsc,shfc,swfc,srfc,sstc,sssc):
  nc=Dataset(grdname,'r')
  L=nc.dimensions['xi_psi'].size
  M=nc.dimensions['eta_psi'].size
  nc.close()
  Lp=L+1
  Mp=M+1

  nw = Dataset(frcname, 'w', format='NETCDF4')
  #result = redef(nw)
  #
  #  Create dimensions
  #
  nw.createDimension('xi_u', L)
  nw.createDimension('eta_u', Mp)
  nw.createDimension('xi_v', Lp)
  nw.createDimension('eta_v', M)
  nw.createDimension('xi_rho', Lp)
  nw.createDimension('eta_rho', Mp)
  nw.createDimension('xi_psi', L)
  nw.createDimension('eta_psi', M)
  nw.createDimension('sms_time', len(smst))
  nw.createDimension('shf_time', shft)
  nw.createDimension('swf_time', swft)
  nw.createDimension('sst_time', sstt)
  nw.createDimension('srf_time', srft)
  nw.createDimension('sss_time', ssst)
  nw.createDimension('wwv_time', len(smst))
  #
  #  Create variables and attributes
  #
  nw_sms_time = nw.createVariable('sms_time','d',('sms_time'))
  nw_sms_time.long_name = 'surface momentum stress time'
  nw_sms_time.units = 'days'
  nw_sms_time.cycle_length = smsc

  nw_shf_time = nw.createVariable('shf_time', 'd',('shf_time'))
  nw_shf_time.long_name = 'surface heat flux time'
  nw_shf_time.units = 'days'
  nw_shf_time.cycle_length =shfc 

  nw_swf_time = nw.createVariable('swf_time', 'd',('swf_time'))
  nw_swf_time.long_name = 'surface freshwater flux time'
  nw_swf_time.units = 'days'
  nw_swf_time.cycle_length = swfc

  nw_sst_time = nw.createVariable('sst_time', 'd',('sst_time'))
  nw_sst_time.long_name = 'sea surface temperature time'
  nw_sst_time.units = 'days'
  nw_sst_time.cycle_length = sstc

  nw_sss_time = nw.createVariable('sss_time', 'd',('sss_time'))
  nw_sss_time.long_name = 'sea surface salinity time'
  nw_sss_time.units = 'days'
  nw_sss_time.cycle_length = sssc

  nw_srf_time = nw.createVariable('srf_time','d',('srf_time'))
  nw_srf_time.long_name = 'solar shortwave radiation time'
  nw_srf_time.units = 'days'
  nw_srf_time.cycle_length = srfc

  nw_wwv_time = nw.createVariable('wwv_time','d',('wwv_time'))
  nw_wwv_time.long_name = 'surface wave fields time'
  nw_wwv_time.units = 'days'
  nw_wwv_time.cycle_length = smsc


  nw_sustr = nw.createVariable('sustr','d',('sms_time', 'eta_u', 'xi_u'))
  nw_sustr.long_name = 'surface u-momentum stress'
  nw_sustr.units = 'Newton meter-2'

  nw_svstr = nw.createVariable('svstr','d',('sms_time', 'eta_v', 'xi_v'))
  nw_svstr.long_name = 'surface v-momentum stress'
  nw_svstr.units = 'Newton meter-2'

  nw_shflux = nw.createVariable('shflux', 'd',('shf_time', 'eta_rho', 'xi_rho'))
  nw_shflux.long_name = 'surface net heat flux'
  nw_shflux.units = 'Watts meter-2'

  nw_swflux = nw.createVariable('swflux','d',('swf_time', 'eta_rho', 'xi_rho'))
  nw_swflux.long_name = 'surface freshwater flux (E-P)'
  nw_swflux.units = 'centimeter day-1'
  nw_swflux.positive = 'net evaporation'
  nw_swflux.negative = 'net precipitation'

  nw_SST = nw.createVariable('SST','d',('sst_time', 'eta_rho', 'xi_rho'))
  nw_SST.long_name = 'sea surface temperature'
  nw_SST.units = 'Celsius'

  nw_SSS = nw.createVariable('SSS','d',('sss_time', 'eta_rho', 'xi_rho'))
  nw_SSS.long_name = 'sea surface salinity'
  nw_SSS.units = 'PSU'

  nw_dQdSST = nw.createVariable('dQdSST','d',('sst_time', 'eta_rho', 'xi_rho'))
  nw_dQdSST.long_name = 'surface net heat flux sensitivity to SST'
  nw_dQdSST.units = 'Watts meter-2 Celsius-1'

  nw_swrad = nw.createVariable('swrad','d',('srf_time', 'eta_rho', 'xi_rho'))
  nw_swrad.long_name = 'solar shortwave radiation'
  nw_swrad.units = 'Watts meter-2'
  nw_swrad.positive = 'downward flux, heating'
  nw_swrad.negative = 'upward flux, cooling'

  nw_Awave = nw.createVariable('Awave','d',('wwv_time', 'eta_rho', 'xi_rho'))
  nw_Awave.long_name = 'wind induced wave amplitude'
  nw_Awave.units = 'm'

  nw_Dwave = nw.createVariable('Dwave','d',('wwv_time', 'eta_rho', 'xi_rho'))
  nw_Dwave.long_name = 'wind induced wave direction'
  nw_Dwave.units = 'degree'

  nw_Pwave = nw.createVariable('Pwave','d',('wwv_time', 'eta_rho', 'xi_rho'))
  nw_Pwave.long_name = 'wind induced wave period'
  nw_Pwave.units = 'second'

  #result = endef(nw)

  #
  # Create global attributes
  #
  nw.title = title
  nw.date = datetime.today().strftime('%Y-%m-%d') 
  nw.grd_file = grdname
  nw.type = 'CROCO forcing file'

  #
  # Write time variables
  #

  nw['sms_time'][:] = smst
  nw['shf_time'][:] = shft
  nw['swf_time'][:] = swft
  nw['sst_time'][:] = sstt
  nw['srf_time'][:] = srft
  nw['sss_time'][:] = ssst
  nw['wwv_time'][:] = smst

  nw.close()
  return


#####################################################################
# 	Create an empty netcdf heat flux bulk forcing file
#       frcname: name of the forcing file
#       grdname: name of the grid file
#       title: title in the netcdf file  
####################################################################
def create_bulk(frcname,grdname,title,bulkt,bulkc):
  nc=Dataset(grdname,'r')
  L=nc.dimensions['xi_psi'].size
  M=nc.dimensions['eta_psi'].size
  nc.close()
  Lp=L+1
  Mp=M+1

  nw = Dataset(frcname, 'w', format='NETCDF4')
  #result = redef(nw)

  #
  #  Create dimensions
  #
  nw.createDimension('xi_rho', Lp)
  nw.createDimension('eta_rho', Mp)
  nw.createDimension('xi_psi', L)
  nw.createDimension('eta_psi', M)
  nw.createDimension('xi_u', L)
  nw.createDimension('eta_u', Mp)
  nw.createDimension('xi_v', Lp)
  nw.createDimension('eta_v', M)
  nw.createDimension('bulk_time', len(bulkt))
  #
  #  Create variables and attributes
  #
  nw_bulk_time = nw.createVariable('bulk_time', 'd',('bulk_time'))
  nw_bulk_time.long_name    = 'bulk formulation execution time'
  nw_bulk_time.units        = 'days'
  nw_bulk_time.cycle_length = bulkc

  nw_tair             = nw.createVariable('tair', 'd', ('bulk_time', 'eta_rho', 'xi_rho'))
  nw_tair.long_name   = 'surface air temperature'
  nw_tair.units       = 'Celsius'

  nw_rhum             = nw.createVariable('rhum','d',('bulk_time', 'eta_rho', 'xi_rho'))
  nw_rhum.long_name   = 'relative humidity'
  nw_rhum.units       = 'fraction'

  nw_prate            = nw.createVariable('prate','d',('bulk_time', 'eta_rho', 'xi_rho'))
  nw_prate.long_name  = 'precipitation rate'
  nw_prate.units      = 'cm day-1'

  nw_wspd             = nw.createVariable('wspd','d',('bulk_time', 'eta_rho', 'xi_rho'))
  nw_wspd.long_name   = 'wind speed 10m'
  nw_wspd.units       = 'm s-1'

  nw_radlw            = nw.createVariable('radlw','d',('bulk_time', 'eta_rho', 'xi_rho'))
  nw_radlw.long_name  = 'net outgoing longwave radiation'
  nw_radlw.units      = 'Watts meter-2'
  nw_radlw.positive   = 'upward flux, cooling water'

  nw_radlw_in            = nw.createVariable('radlw_in','d',('bulk_time', 'eta_rho', 'xi_rho'))
  nw_radlw_in.long_name  = 'downward longwave radiation'
  nw_radlw_in.units      = 'Watts meter-2'
  nw_radlw_in.positive   = 'downward flux, warming water'

  nw_radsw            = nw.createVariable('radsw','d',('bulk_time', 'eta_rho', 'xi_rho'))
  nw_radsw.long_name  = 'shortwave radiation'
  nw_radsw.units      = 'Watts meter-2'
  nw_radsw.positive   = 'downward flux, heating water'

  nw_sustr = nw.createVariable('sustr','d',('bulk_time', 'eta_u', 'xi_u'))
  nw_sustr.long_name = 'surface u-momentum stress'
  nw_sustr.units = 'Newton meter-2'

  nw_svstr = nw.createVariable('svstr','d',('bulk_time', 'eta_v', 'xi_v'))
  nw_svstr.long_name = 'surface v-momentum stress'
  nw_svstr.units = 'Newton meter-2'

  nw_uwnd = nw.createVariable('uwnd','d',('bulk_time', 'eta_u', 'xi_u'))
  nw_uwnd.long_name = 'u-wind'
  nw_uwnd.units = 'm/s'

  nw_vwnd = nw.createVariable('vwnd','d',('bulk_time', 'eta_v', 'xi_v'))
  nw_vwnd.long_name = 'v-wind'
  nw_vwnd.units = 'm/s'

  #result = endef(nw)

  #
  # Create global attributes
  #
  nw.title = title
  nw.date = datetime.today().strftime('%Y-%m-%d') 
  nw.grd_file = grdname
  nw.type = 'CROCO heat flux bulk forcing file'

  #
  # Write time variables
  #
  for tndx in range(len(bulkt)):
    if tndx%20 == 0:
      print(['Time Step Bulk: '+str(tndx)+' of '+str(len(bulkt))])
    
    nw['bulk_time'][tndx] = bulkt[tndx]
  
  nw.close()
  return