import numpy as np
import netCDF4 as netcdf
import sys

#######################################################
#Transfert a field at rho points to psi points
#######################################################
def rho2psi(var_rho):
    if np.ndim(var_rho)<3:
        var_psi = rho2psi_2d(var_rho)
    else:
        var_psi = rho2psi_3d(var_rho)
    return var_psi
###########
def rho2psi_2d(var_rho):
    var_psi = 0.25*(var_rho[1:,1:]+var_rho[1:,:-1]+var_rho[:-1,:-1]+var_rho[:-1,1:])
    return var_psi
###########
def rho2psi_3d(var_rho):
    var_psi = 0.25*(var_rho[1:,1:,:]+var_rho[1:,:-1,:]+var_rho[:-1,:-1,:]+var_rho[:-1,1:,:])
    return var_psi

#######################################################
#Transfert a field at rho points to u points
#######################################################
def rho2u(var_rho):
    if np.ndim(var_rho)==1:
        var_u = 0.5*(var_rho[1:]+var_rho[:-1])
    elif np.ndim(var_rho)==2:
        var_u = rho2u_2d(var_rho)
    elif np.ndim(var_rho)==3:
        var_u = rho2u_3d(var_rho)
    else:
        var_u = rho2u_4d(var_rho)
    return var_u
###########
def rho2u_2d(var_rho):
    var_u = 0.5*(var_rho[:,1:]+var_rho[:,:-1])
    return var_u
###########
def rho2u_3d(var_rho):
    var_u = 0.5*(var_rho[:,:,1:]+var_rho[:,:,:-1])
    return var_u
###########
def rho2u_4d(var_rho):
    var_u = 0.5*(var_rho[:,:,:,1:]+var_rho[:,:,:,:-1])
    return var_u

#######################################################
#Transfert a field at rho points to v points
#######################################################
def rho2v(var_rho):
    if np.ndim(var_rho)==1:
        var_v = 0.5*(var_rho[1:]+var_rho[:-1])
    elif np.ndim(var_rho)==2:
        var_v = rho2v_2d(var_rho)
    elif np.ndim(var_rho)==3:
        var_v = rho2v_3d(var_rho)
    else:
        var_v = rho2v_4d(var_rho)
    return var_v
###########
def rho2v_2d(var_rho):
    var_v = 0.5*(var_rho[1:,:]+var_rho[:-1,:])
    return var_v
###########
def rho2v_3d(var_rho):
    var_v = 0.5*(var_rho[:,1:,:]+var_rho[:,:-1,:])
    return var_v
###########
def rho2v_4d(var_rho):
    var_v = 0.5*(var_rho[:,:,1:,:]+var_rho[:,:,:-1,:])
    return var_v

#######################################################
#Transfert a field at u points to the rho points
#######################################################
def v2rho(var_v):
    if np.ndim(var_v)<3:
        var_rho = v2rho_2d(var_v)
    elif np.ndim(var_v)==3:
        var_rho = v2rho_3d(var_v)
    else:
        var_rho = v2rho_4d(var_v)
    return var_rho
###########
def v2rho_2d(var_v):
    [L,Mp]=var_v.shape
    Lp=L+1
    Lm=LM-1
    var_rho=np.zeros((Lp,Mp))
    var_rho[1:L,:]=0.5*(var_v[0:Lm:,:]+var_v[1:L,:])
    var_rho[0,:]=var_rho[1,:]
    var_rho[Lp-1,:]=var_rho[L-1,:]
    return var_rho
###########
def v2rho_3d(var_v):
    [N,L,Mp]=var_v.shape
    Lp=L+1
    Lm=L-1
    var_rho=np.zeros((N,Lp,Mp))
    var_rho[:,1:L,:]=0.5*(var_v[:,0:Lm,:]+var_v[:,1:L,:])
    var_rho[:,0,:]=var_rho[:,1,:]
    var_rho[:,Lp-1,:]=var_rho[:,L-1,:]
    return var_rho
##########
def v2rho_4d(var_v):
    [T,N,L,Mp]=var_v.shape
    Lp=L+1
    Lm=L-1
    var_rho=np.zeros((T,N,Lp,Mp))
    var_rho[:,:,1:L,:]=0.5*(var_v[:,:,0:Lm,:]+var_v[:,:,1:L,:])
    var_rho[:,:,0,:]=var_rho[:,:,1,:]
    var_rho[:,:,Lp-1,:]=var_rho[:,:,L-1,:]
    return var_rho

#######################################################
#Transfert a 2 or 2-D field at u points to the rho points
#######################################################
def u2rho(var_u):
    if np.ndim(var_u)<3:
        var_rho = u2rho_2d(var_u)
    elif np.ndim(var_u)==3:
        var_rho = u2rho_3d(var_u)
    else:
        var_rho = u2rho_4d(var_u)
    return var_rho
##########
def u2rho_2d(var_u):
    [Lp,M]=var_u.shape
    Mp=M+1
    Mm=M-1
    var_rho=np.zeros((Lp,Mp))
    var_rho[:,1:M]=0.5*(var_u[:,0:Mm]+var_u[:,1:M])
    var_rho[:,0]=var_rho[:,1]
    var_rho[:,Mp-1]=var_rho[:,M-1]
    return var_rho
##########
def u2rho_3d(var_u):
    [N,Lp,M]=var_u.shape
    Mp=M+1
    Mm=M-1
    var_rho=np.zeros((N,Lp,Mp))
    var_rho[:,:,1:M]=0.5*(var_u[:,:,0:Mm]+var_u[:,:,1:M])
    var_rho[:,:,0]=var_rho[:,:,1]
    var_rho[:,:,Mp-1]=var_rho[:,:,M-1]
    return var_rho
##########
def u2rho_4d(var_u):
    [T,N,Lp,M]=var_u.shape
    Mp=M+1
    Mm=M-1
    var_rho=np.zeros((T,N,Lp,Mp))
    var_rho[:,:,:,1:M]=0.5*(var_u[:,:,:,0:Mm]+var_u[:,:,:,1:M])
    var_rho[:,:,:,0]=var_rho[:,:,:,1]
    var_rho[:,:,:,Mp-1]=var_rho[:,:,:,M-1]
    return var_rho

#######################################################
#Transfert a field at psi points to rho points
#######################################################
def psi2rho(var_psi):
    if np.ndim(var_psi)<3:
        var_rho = psi2rho_2d(var_psi)
    else:
        var_rho = psi2rho_3d(var_psi)
    return var_rho
###########
def psi2rho_2d(var_psi):
    [L,M]=var_psi.shape
    Mp=M+1
    Lp=L+1
    Mm=M-1
    Lm=L-1
    var_rho=np.zeros((Lp,Mp))
    var_rho[1:L,1:M]=0.25*(var_psi[0:Lm,0:Mm]+var_psi[0:Lm,1:M]+var_psi[1:L,0:Mm]+var_psi[1:L,1:M])
    var_rho[:,0]=var_rho[:,1]
    var_rho[:,Mp-1]=var_rho[:,M-1]
    var_rho[0,:]=var_rho[1,:]
    var_rho[Lp-1,:]=var_rho[L-1,:]
    return var_rho
###########
def psi2rho_3d(var_psi):
    [Nz,Lz,Mz]=var_psi.shape
    var_rho=np.zeros((Nz,Lz+1,Mz+1))
    for iz in range(0, Nz, 1):
        var_rho[iz,:,:]=psi2rho_2d(var_psi[iz:,:])
    return var_rho


#######################################################
#Transfert a 3-D field from verical w points to vertical rho-points
#######################################################
def w2rho(var_w):
    [N,L,M]=var_w.shape
    print( '[N,L,M]',[N,L,M])
    var_rho = np.zeros((N-1,L,M))
    for iz in range(1,N-2):
        var_rho[iz,:,:]  = 0.5625*(var_w[iz+1,:,:] + var_w[iz,:,:]) -0.0625*(var_w[iz+2:,:] + var_w[iz-1,:,:])
    var_rho[0,:,:]  = -0.125*var_w[1,:,:] + 0.75*var_w[1,:,:] +0.375*var_w[0,:,:]
    var_rho[N-2,:,:]  = -0.125*var_w[N-3,:,:] + 0.75*var_w[N-2,:,:] +0.375*var_w[N-1,:,:]
    return var_rho


### FUNCTION ZLEVS #####################################################
#
#

def zlevs(h,zeta,theta_s,theta_b,hc,N,type,vtransform):
#
# function z=zlevs(h,zeta,theta_s,theta_b,hc,N,type,vtransform)
#
#  this function compute the depth of rho or w points for CROCO
#
#  On Input:
#
#    type    'r': rho point 'w': w point 
#    vtransform  1=> old v transform (Song, 1994); 
#                2=> new v transform (Shcheptekin, 2006)
#  On Output:
#
#    z       Depths (m) of RHO- or W-points (3D matrix).
#
########################################################################
#
# Test the number of dimension for h
#
  Ndim=np.size(np.shape(h))
  if Ndim==2:  
    M, L = np.squeeze(np.shape(h))
  elif Ndim==1:
    L = np.squeeze(np.shape(h))
  else:
    print('zlevs: error - incorrect dimension for h')
    return
    
  hmin = h.min()
  hmax = h.max()
#
# Set S-Curves in domain [-1 < sc < 0] at vertical W- and RHO-points.
#
  ds = 1.0 / float(N)      
  if type=='w': 
    sc = ds * (np.arange(0,N+1) - N)
    Nmax=N+1   
  elif type=='r':  
    sc = ds * (np.arange(1,N+1) - N - 0.5)
    Nmax=N
  else:
    print('Problem with type = ',type)
    sys.exit()

  if vtransform==1:  
    print('OLD_S_COORD')
    cff1 = 1.0 / np.sinh(theta_s)
    cff2 = 0.5 / np.tanh(0.5*theta_s)
    Cs = (1. - theta_b) * cff1 * np.sinh(theta_s * sc) + \
           theta_b * (cff2 * np.tanh(theta_s * (sc + 0.5)) - 0.5)
  elif vtransform==2:  
    print('NEW_S_COORD')
    Cs = get_csf(sc,theta_s,theta_b)         
  else:
    print('Problem with vtransform = ',vtransform)
    sys.exit()
#
#  Set vertical grid
#
  if Ndim==2:  
    z = np.zeros([Nmax, M, L])
  elif Ndim==1:
    z = np.zeros([Nmax, L])

  for k in np.arange(0,Nmax):
    if vtransform==1: 
      cff = hc*(sc[k]-Cs[k])
      cff1 = Cs[k]
      z0 = cff + cff1*h
      hinv=1./h  
      if Ndim==2:  
        z[k,:,:]=z0+zeta*(1.+z0*hinv)
      elif Ndim==1:
        z[k,:]=z0+zeta*(1.+z0*hinv)
    elif vtransform==2: 
      cff = hc*sc[k]
      cff1 = Cs[k]
      z0 = cff + cff1*h 
      hinv = 1./(h+hc)
      if Ndim==2:  
       z[k,:,:]=z0*h*hinv+zeta*(1.+z0*hinv)
      elif Ndim==1:
       z[k,:]=z0*h*hinv+zeta*(1.+z0*hinv)
  return z
#
#
### END FUNCTION ZLEVS #################################################




##############################################################
##############################################################
##############################################################

def vintegr(var,zw,zr,z01,z02):

#function [V,h0]=vintegr2(var,zw,zr,z01,z02)
#
# Vertically integrate a CROCO variable (var) from a constant depth 
# z01 (ex z01=-4000 m) to a constant depth z02 (ex z02=-2000m).
#
# If z01 = NaN : perform the integration from the bottom to z02.
# If z02 = NaN : perform the integration from z01 to the surface.
# If they are both NaNs perform the integration from the bottom
# to the surface.
#
# Input :
#
# var : CROCO variable at RHO-points (3D matrix)
# zw  : Depth of the W-points (3D matrix)
# zr  : Depth of the RHO-points (3D matrix)
# z01 : lower limit of integration (scalar)
# z02 : upper limit of integration (scalar)
#
# Output :
#
# V   : intgrated value (2D matrix)
# h0  : layer thickness (2D matrix)
#
# Pierrick Penven 2005
#
    if z02 <= z01:
        print('vintegr2:  z02 <= z01')
        sys.exit()

    [Np,M,L]=np.shape(zw)
    N=Np-1
    mask=np.zeros([M,L]) + 1.

    i1=np.arange(0,L)
    j1=np.arange(0,M)
    [i2,j2]=np.meshgrid(i1,j1)

    za=np.reshape(zw,Np*M*L)
    vara=np.reshape(var,N*M*L)
#
# Get the matrix of the variables above of z01 and below z02
# 
    if (np.isfinite(z01) & np.isfinite(z02)):
        isgood=np.int_(zw[1:,:,:]>z01) * np.int_(zw[1:,:,:]<z02)
    elif np.isfinite(z01):
        isgood=np.int_(zw[1:,:,:]>z01)
    elif np.isfinite(z02):
        isgood=np.int_(zw[1:,:,:]<z02)
    else:
        isgood=np.int_(var==var)
#
    if np.isfinite(z01):
#
# Find the bottom limit of the corresponding grid cells
#
        a=np.int_(zw<z01)
        levs=np.sum(a,axis=0)-1
        mask=np.zeros([M,L]) + 1.
        mask[np.where(levs<0)]=np.nan
        mask[np.where( (levs<0) | (levs==N) )]=np.nan
        levs[np.where(levs==N)]=1

        pos = L*M*levs + L*j2 + i2
        z1=mask*za[pos]
#
# Compute the vertical size of the partial step to add at the bottom
#
        dzbot=z1-z01
        dzbot[np.where(np.isnan(mask))]=0.
#
# Get the value of the variable in the partial step to add at the bottom
#
        Vbot=vara[pos]
    else:
        dzbot=0
        Vbot=0

    if np.isfinite(z02):
#
# Find the top positions
#
        a=np.int_(zw<z02)
        levs=np.sum(a,axis=0)-1
        mask=np.zeros([M,L]) + 1.
        mask[np.where( (levs<0) | (levs==N) )]=np.nan
        levs[np.where(levs==N)]=1

        pos = L*M*levs + L*j2 + i2
        z2=mask*za[pos]
#
# Compute the vertical size of the partial step to add at the top
#
        dztop=z02-z2
        dztop[np.where(np.isnan(mask))]=0.
#
# Get the value of the variable in the partial step to add at the bottom
#
        Vtop=vara[pos]
    else:
        dztop=0
        Vtop=0
#
# Perform the vertical integration
    dz=zw[1:,:,:]-zw[:-1,:,:]
    V=np.sum(dz*isgood*var,axis=0) + dzbot*Vbot + dztop*Vtop
#
# Get the depth
#
    h0=np.sum(dz*isgood,axis=0) + dzbot + dztop

    V[np.where(h0==0)]=0
    h0[np.where(h0==0)]=0

    return V,h0

########
def vintegr4D(var,zw,zr,z01,z02):

#function [V,h0]=vintegr2(var,zw,zr,z01,z02)
#
# Vertically integrate a CROCO variable (var) from a constant depth
# z01 (ex z01=-4000 m) to a constant depth z02 (ex z02=-2000m).
#
# If z01 = NaN : perform the integration from the bottom to z02.
# If z02 = NaN : perform the integration from z01 to the surface.
# If they are both NaNs perform the integration from the bottom
# to the surface.
#
# Input :
#
# var : CROCO variable at RHO-points (3D matrix)
# zw  : Depth of the W-points (3D matrix)
# zr  : Depth of the RHO-points (3D matrix)
# z01 : lower limit of integration (scalar)
# z02 : upper limit of integration (scalar)
#
# Output :
#
# V   : intgrated value (2D matrix)
# h0  : layer thickness (2D matrix)
#
# Pierrick Penven 2005
#
    if z02 <= z01:
        print('vintegr2:  z02 <= z01')
        sys.exit()

    [T,Np,M,L]=np.shape(zw)
    N=Np-1
    mask=np.zeros([T,M,L]) + 1.

    i1=np.arange(0,L)
    j1=np.arange(0,M)
    t1=np.arange(0,T)
    [j2,t2,i2]=np.meshgrid(j1,t1,i1)
#    [i2,j2]=np.meshgrid(i1,j1)

    za=np.reshape(zw,T*Np*M*L)
    vara=np.reshape(var,T*N*M*L)
#
# Get the matrix of the variables above of z01 and below z02
#
    if (np.isfinite(z01) & np.isfinite(z02)):
        isgood=np.int_(zw[:,1:,:,:]>z01) * np.int_(zw[:,1:,:,:]<z02)
    elif np.isfinite(z01):
        isgood=np.int_(zw[:,1:,:,:]>z01)
    elif np.isfinite(z02):
        isgood=np.int_(zw[:,1:,:,:]<z02)
    else:
        isgood=np.int_(var==var)
#
    if np.isfinite(z01):
#
# Find the bottom limit of the corresponding grid cells
#
        a=np.int_(zw<z01)
        levs=np.sum(a,axis=1)-1
        mask=np.zeros([T,M,L]) + 1.
        mask[np.where(levs<0)]=np.nan
        mask[np.where( (levs<0) | (levs==N) )]=np.nan
        levs[np.where(levs==N)]=1

#        pos = L*M*levs + L*j2 + i2
        pos= L*M*T*levs + L*M*j2 + L*j2 + i2
        z1=mask*za[pos]
#
# Compute the vertical size of the partial step to add at the bottom
#
        dzbot=z1-z01
        dzbot[np.where(np.isnan(mask))]=0.
#
# Get the value of the variable in the partial step to add at the bottom
#
        Vbot=vara[pos]
    else:
        dzbot=0
        Vbot=0

    if np.isfinite(z02):
#
# Find the top positions
#
        a=np.int_(zw<z02)
        levs=np.sum(a,axis=1)-1
        mask=np.zeros([T,M,L]) + 1.
        mask[np.where( (levs<0) | (levs==N) )]=np.nan
        levs[np.where(levs==N)]=1

#        pos = L*M*levs + L*j2 + i2
        pos= L*M*T*levs + L*M*j2 + L*j2 + i2
        z2=mask*za[pos]
#
# Compute the vertical size of the partial step to add at the top
#
        dztop=z02-z2
        dztop[np.where(np.isnan(mask))]=0.
#
# Get the value of the variable in the partial step to add at the bottom
#
        Vtop=vara[pos]
    else:
        dztop=0
        Vtop=0
#
# Perform the vertical integration
    dz=zw[:,1:,:,:]-zw[:,:-1,:,:]
    V=np.sum(dz*isgood*var,axis=0) + dzbot*Vbot + dztop*Vtop
#
# Get the depth
#
    h0=np.sum(dz*isgood,axis=1) + dzbot + dztop

    V[np.where(h0==0)]=0
    h0[np.where(h0==0)]=0

    return V,h0



##############################################################
##############################################################
##############################################################
def progressbar(it, prefix="", size=60, file=sys.stdout):
    count = len(it)
    def show(j):
        x = int(size*j/count)
        file.write("%s[%s%s] %i/%i\r" % (prefix, "#"*x, "."*(size-x), j, count))
        file.flush()
    show(0)
    for i, item in enumerate(it):
        yield item
        show(i+1)
    file.write("\n")
    file.flush()

##############################################################
##############################################################
##############################################################

def geo_idx(dd, dd_array):

    """
     - dd - the decimal degree (latitude or longitude)
     - dd_array - the list of decimal degrees to search.
     search for nearest decimal degree in an array of decimal degrees and return the index.
     np.argmin returns the indices of minium value along an axis.
     so subtract dd from all values in dd_array, take absolute value and find index of minium.
   """
    geo_idx = (np.abs(dd_array - dd)).argmin()
    return geo_idx

def indx_bound(x, x0):
    """
    Conversion of fortran tools indx_bound
    """
    n=x.shape[0]
    if x0 < x[0] :
        i=0                      # if x0 is outside the full range
    elif x0 > x[-1] :            # of x(1) ... x(n), then return
        i=n                      # i=0 or i=n.
    else:
        i=int( ( x[-1]-x0 +n*(x0-x[0]) )/(x[-1]-x[0]) )
        if x[i+1]<x0 :
            while x[i+1] <x0 :  # This algorithm computes "i" as
                i=i+1           # linear interpolation between x(1)
                                # and x(n) which should yield the
        elif x[i] > x0 :        # correct value for "i" right a way
            while x[i] > x0 :    # because array elements x(i) are
                i=i-1           # equidistantly spaced.  The while
                                # loops are here merely to address
                                # possible roundoff errors.

        if x[i+1]-x0 < 0 or x0-x[i] < 0 :
            print('### ERROR: indx_bound :: ',x[i], x0, x[i+1], x0-x[i], x[i+1]-x0)
            exit()
    indx_bound=i
    return indx_bound

##############################################################
##############################################################
##############################################################

def read_nc(filename,varname, indices="[:]"):
    try:
        with netcdf.Dataset(filename,'r') as nc:
            var = eval(''.join(("nc.variables[varname]", indices)))
    except Exception:
        raise
   #
    if 'float32' in var.dtype.name:
        return var.astype(np.float64)
    else:
        return var


def read_nc_mf(filename,varname,indices="[:]",time_dim='time'):
    # read multiple netcdf. filename need to be glob.glob(...)
    try:
        try:
            # Load over unlimited dimension
            with netcdf.MFDataset(filename) as nc:
                    var =  eval(''.join(("nc.variables[varname]", indices)))
        except Exception:
            try:
                # Load over time dimension
                with netcdf.MFDataset(filename, aggdim=time_dim) as nc:
                    var =  eval(''.join(("nc.variables[varname]", indices)))
            except Exception:
                print("Vars can not be loaded along %s. Please specify time_dim in read_nc_mf." % time_dim)
    except Exception:
        raise
    if 'float32' in var.dtype.name:
        return var.astype(np.float64)
    else:
        return var



