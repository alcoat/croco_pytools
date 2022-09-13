import numpy as np
from netCDF4 import Dataset  as netcdf
from scipy.interpolate import griddata
import os
import sys
import time
from scipy.spatial import Delaunay
import tools
from tools import progressbar

###############
# Get interpolation weight
###############
def get_tri_coef(X, Y, newX, newY, verbose=0):

    """
    Inputs:
        origin lon and lat 2d arrays (X,Y)
        child lon and lat 2d arrays (newX,newY)
    Ouputs:
        elem - pointers to 2d gridded data (at lonp,latp locations) from
            which the interpolation is computed (3 for each child point)
        coef - linear interpolation coefficients
    Use:
        To subsequently interpolate data from Fp to Fc, the following
        will work:      Fc  = sum(coef.*Fp(elem),3);  This line  should come in place of all
        griddata calls. Since it avoids repeated triangulations and tsearches (that are done
        with every call to griddata) it should be much faster.
    """

    Xp = np.array([X.ravel(), Y.ravel()]).T
    Xc = np.array([newX.ravel(), newY.ravel()]).T


    #Compute Delaunay triangulation
    if verbose==1: tstart = tm.time()
    tri = Delaunay(Xp)
    if verbose==1: print('Delaunay Triangulation', tm.time()-tstart)

    #Compute enclosing simplex and barycentric coordinate (similar to tsearchn in MATLAB)
    npts = Xc.shape[0]
    p = np.zeros((npts,3))

    points = tri.points[tri.vertices[tri.find_simplex(Xc)]]
    if verbose==1: tstart = tm.time()
    for i in progressbar(range(npts),'  Get_tri_coef: ', 40):

        if verbose==1: print(np.float(i)/npts)

        if tri.find_simplex(Xc[i])==-1:  #Point outside triangulation
             p[i,:] = p[i,:] * np.nan

        else:

            if verbose==1: tstart = tm.time()
            A = np.append(np.ones((3,1)),points[i] ,axis=1)
            if verbose==1: print('append A', tm.time()-tstart)

            if verbose==1: tstart = tm.time()
            B = np.append(1., Xc[i])
            if verbose==1: print('append B', tm.time()-tstart)

            if verbose==1: tstart = tm.time()
            p[i,:] = np.linalg.lstsq(A.T,B.T)[0]
            if verbose==1: print('solve', tm.time()-tstart)


    if verbose==1: print('Coef. computation 1', tm.time()-tstart)

    if verbose==1: tstart = tm.time()
    elem = np.reshape(tri.vertices[tri.find_simplex(Xc)],(newX.shape[0],newY.shape[1],3))
    coef = np.reshape(p,(newX.shape[0],newY.shape[1],3))
    if verbose==1: print('Coef. computation 2', tm.time()-tstart)

    return(elem,coef)

#######################################

def get_delaunay_bry(lon_bry,lat_bry,inputfile,bdy):
#   This function computes the delaunay matrices for the interpolations 
#   at the boundaies
#
#   Input:
#     lon_bry      Longitudes of the boundary (vector).
#     lat_bry      Latitudes of the boundary (vector).
#     inputfile    netcdf structure poiting to the input file
#     bdy          which boundary is done
#
#   Output:
#     LonT_bry,LatT_bry,iminT_bry,imaxT_bry,jminT_bry,jmaxT_bry,elemT_bry,coefT_bry,
#     LonU_bry,LatU_bry,iminU_bry,imaxU_bry,jminU_bry,jmaxU_bry,elemU_bry,coefU_bry,
#     LonV_bry,LatV_bry,iminV_bry,imaxV_bry,jminV_bry,jmaxV_bry,elemV_bry,coefV_bry
#
    comp_delaunay=1

#
# get grid positions
#
    LonU_bry,LatU_bry   = eval(''.join(("inputfile.lonU"+bdy))),eval(''.join(("inputfile.latU"+bdy)))
    LonV_bry,LatV_bry   = eval(''.join(("inputfile.lonV"+bdy))),eval(''.join(("inputfile.latV"+bdy)))
    LonT_bry,LatT_bry   = eval(''.join(("inputfile.lonT"+bdy))),eval(''.join(("inputfile.latT"+bdy))) 

#
# Get the 2D interpolation coefficients
#

    if comp_delaunay==1:

      print('\nCompute Delaunay triangulation from T-points to CROCO rho_points...')
      print('-------------------------------------------------------------------')
      [elemT_bry,coefT_bry] = get_tri_coef(LonT_bry,LatT_bry,lon_bry,lat_bry)
      coefnorm=np.sum(coefT_bry,axis=2)
      coefT_bry=coefT_bry/coefnorm[:,:,np.newaxis]

      print('\nCompute Delaunay triangulation from U-points to CROCO rho_points...')
      print('-------------------------------------------------------------------')
      [elemU_bry,coefU_bry] = get_tri_coef(LonU_bry,LatU_bry,lon_bry,lat_bry)
      coefnorm=np.sum(coefU_bry,axis=2)
      coefU_bry=coefU_bry/coefnorm[:,:,np.newaxis]

      print('\nCompute Delaunay triangulation from V-points to CROCO rho_points...')
      print('-------------------------------------------------------------------')
      [elemV_bry,coefV_bry] = get_tri_coef(LonV_bry,LatV_bry,lon_bry,lat_bry)
      coefnorm=np.sum(coefV_bry,axis=2)
      coefV_bry=coefV_bry/coefnorm[:,:,np.newaxis]

      # Save the Delaunay triangulation matrices
      np.savez('coeffs_bry'+bdy+'.npz',\
                coefT=coefT_bry,elemT=elemT_bry,\
                coefU=coefU_bry,elemU=elemU_bry,\
                coefV=coefV_bry,elemV=elemV_bry)

    else:
#
# Load the Delaunay triangulation matrices
#
      print('Load Delaunay triangulation...')
      data=np.load('coeffs_bry'+bdy+'.npz')
      coefT_bry = data['coefT_bry']
      elemT_bry = data['elemT_bry']
      coefU_bry = data['coefU_bry']
      elemU_bry = data['elemU_bry']
      coefV_bry = data['coefV_bry']
      elemV_bry = data['elemV_bry']

    print('Delaunay triangulation done')

    return(elemT_bry,coefT_bry,elemU_bry,coefU_bry,elemV_bry,coefV_bry)




################################################################
#

##################################################################

def ztosigma(vin,Z,zcroco):
# Do a vertical interpolation from z levels to sigma CROCO levels
    if len(zcroco.shape)>3:
        [T,N,M,L]=np.shape(zcroco)
        four_dim=True
    else:
        [N,M,L]=np.shape(zcroco)
        four_dim=False

    [Nz]=np.shape(Z)
#
# Find the grid position of the nearest vertical levels
#
    i1=np.arange(0,L)
    j1=np.arange(0,M)
    if four_dim:
        t1=np.arange(0,T)
        [jmat,tmat,imat]=np.meshgrid(j1,t1,i1)
        VAR=np.reshape(vin,T*Nz*M*L)
        vout=np.zeros((T,N,M,L))
    else:
        [imat,jmat]=np.meshgrid(i1,j1)
        VAR=np.reshape(vin,Nz*M*L)
        vout=np.zeros((N,M,L))
    for ks in progressbar(range(N),' Sigma layer : ', 40):
        if four_dim:
            sigmalev=zcroco[:,ks,:,:]
            thezlevs=np.zeros((T,M,L),dtype=int)-1
        else:
            sigmalev=zcroco[ks,:,:]
            thezlevs=np.zeros((M,L),dtype=int)-1

        for kz in range(Nz):
            thezlevs[np.where(sigmalev>Z[kz])]=thezlevs[np.where(sigmalev>Z[kz])]+1

        if four_dim:
            pos= L*M*Nz*tmat+ L*M*thezlevs + L*jmat + imat
        else:
            pos= L*M*thezlevs + L*jmat + imat
  
        z1=Z[thezlevs]
        z2=Z[thezlevs+1]
        v1=VAR[pos]
        v2=VAR[pos+L*M]
       
        if four_dim:
            vout[:,ks,:,:]=(((v1-v2)*sigmalev+v2*z1-v1*z2)/(z1-z2))
        else:
            vout[ks,:,:]=(((v1-v2)*sigmalev+v2*z1-v1*z2)/(z1-z2))
    return vout

####################################################################

def add2layers(vin):
# Add a layer below the bottom and above the surface to avoid
# vertical extrapolations when doing a vertical interpolation
    if len(np.shape(vin))==3:
        [Nz,M,L]=np.shape(vin)
        vout=np.zeros((Nz+2,M,L))

        vout[1:-1,:,:]=vin
        vout[0,:,:]=vin[0,:,:]
        vout[-1,:,:]=vin[-1,:,:]

    elif len(np.shape(vin))==4:
        [T,Nz,M,L]=np.shape(vin)
        vout=np.zeros((T,Nz+2,M,L))

        vout[:,1:-1,:,:]=vin
        vout[:,0,:,:]=vin[:,0,:,:]
        vout[:,-1,:,:]=vin[:,-1,:,:]
    return vout

######################################################################
# Contrary to what the name says, this function interpolates the 2-d fields (not tracers temp,salt,...)
def interp_tracers(inputfile,vname,l,k,coef,elem):
#  Remove the missing values from a gridded 2D field
#  and do an horizontal interpolation using Delaunay matrices (coef and elem)

# 0: Read input data informations
    nc        = inputfile.ncglo
    varinp    = inputfile.var
    if vname in [ 'u','ubar' ]:
        Lon,Lat   = inputfile.lonU,inputfile.latU
    elif  vname in [ 'v','vbar' ]:
        Lon,Lat   = inputfile.lonV,inputfile.latV
    else:
        Lon,Lat   = inputfile.lonT,inputfile.latT
   
# 1: Read data
    Vin = inputfile.var_periodicity(vname,l,k) #np.array(nc[varinp[vname]][l,jmin:jmax,imin:imax]) 
    try:
        fvalue=nc[vname][varinp[vname]]._FillValue
    except:
        try:
            fvalue=nc[vname][varinp[vname]].missing_value
        except:
            fvalue=None
            pass

# 2: Remove bad values (using nearest values)
    if fvalue == None:
        igood = np.where(np.isnan(Vin)==False)
        ibad  = np.where(np.isnan(Vin))
    else:
        igood = np.where(Vin!=fvalue)
        ibad  = np.where(Vin==fvalue)

    NzGood=np.size(igood) 
    Nbad=np.size(ibad)
 
    if NzGood==0:
        print('\nWarning: no good data')
        Vin[:]=np.nan
    elif NzGood<10:
        print('\nWarning: less than 10 good values')
        Vin[:]=np.mean(Vin[igood])
    else:
        Vin[ibad] = griddata((Lon[igood],Lat[igood]),Vin[igood],(Lon[ibad],Lat[ibad]),method='nearest')

# 3: 2D interpolation
    Vout = np.sum(coef*Vin.ravel()[elem],2)
    return Vout,NzGood

#################
def interp_tracers3D(inputfile,vname,k,coef,elem,bdy,dtmin,dtmax,prev,nxt): # interp tracers 3D for bdy
#  Remove the missing values from a gridded 2D field
#  and do an horizontal interpolation using Delaunay matrices (coef and elem)

# 0: Read input data informations
    nc        = inputfile.ncglo
    varinp    = inputfile.var
    if vname in [ 'u','ubar' ]:
        Lon,Lat   = eval(''.join(("inputfile.lonU"+bdy))),eval(''.join(("inputfile.latU"+bdy)))
    elif  vname in [ 'v','vbar' ]:
        Lon,Lat   = eval(''.join(("inputfile.lonV"+bdy))),eval(''.join(("inputfile.latV"+bdy)))
    else:
        Lon,Lat   = eval(''.join(("inputfile.lonT"+bdy))),eval(''.join(("inputfile.latT"+bdy)))

# 1: Read data
    if dtmin != dtmax:
        l=np.arange(dtmin,dtmax+1)
    else:
        l=dtmin

    Vin = inputfile.var_periodicity(vname,l,k,bdy=bdy)

    if dtmin==dtmax:
        Vin=Vin[np.newaxis,:]

    if prev == 1: # create overlap before when no data avail
        Vtmp=np.zeros([Vin.shape[0]+1,Vin.shape[1],Vin.shape[2]])
        Vtmp[1:,:]=np.copy(Vin)
        Vtmp[0,:]=np.copy(Vin[0,:])
        Vin=Vtmp
        del Vtmp
    if nxt==1: # create overlap after when no data avail
        Vtmp=np.zeros([Vin.shape[0]+1,Vin.shape[1],Vin.shape[2]])
        Vtmp[:-1,:]=np.copy(Vin)
        Vtmp[-1,:]=np.copy(Vin[-1,:])
        Vin=Vtmp
        del Vtmp

#                if nxt == 1:

    try:
        fvalue=nc[vname][varinp[vname]]._FillValue
    except:
        try:
            fvalue=nc[vname][varinp[vname]].missing_value
        except:
            fvalue=None
            pass

# 2: Remove bad values (using nearest values)
    if fvalue == None:
        igood = np.where(np.isnan(Vin)==False)
        ibad  = np.where(np.isnan(Vin))
    else:
        igood = np.where(Vin!=fvalue)
        ibad  = np.where(Vin==fvalue)

    NzGood=np.size(igood) 
    Nbad=np.size(ibad)  
    for tt in range(Vin.shape[0]):
        if NzGood==0:
#            print('\nWarning: no good data')
            Vin[:]=np.nan
        elif NzGood<10:
#            print('\nWarning: less than 10 good values')
            Vin[igood]=np.mean(Vin[igood])
        else:
            Vin[tt,ibad[1],ibad[2]] = griddata((Lon[igood[1:]],Lat[igood[1:]])\
                                               ,Vin[tt,igood[1],igood[2]],\
                                               (Lon[ibad[1:]],Lat[ibad[1:]]),method='nearest')

# 3: 2D interpolation
    Vout=np.zeros([Vin.shape[0],coef.shape[0],coef.shape[1]]) 
    for k in range(Vout.shape[0]):
        var=Vin[k,:]
        Vout[k,:] = np.sum(coef*var[:].ravel()[elem],2)

    return Vout,NzGood


#######################################################


def interp3d(inputfile,vname,tndx_glo,Nzgoodmin,z_rho,coef,elem):
#  Do a full interpolation of a 4d variable from GLORYS to a CROCO sigma grid
#  1 - Horizontal interpolation on each GLORYS levels
#  2 - Vertical Interpolation from z to CROCO sigma levels
    
    [N,M,L]=np.shape(z_rho)
    depth= inputfile.depth
    [Nz]=np.shape(depth)

    comp_horzinterp = 1 # if 1 compute horizontal interpolations - 0 use saved matrices (for debugging)
    if comp_horzinterp==1:
        print('Horizontal interpolation over z levels')
        t3d=np.zeros((Nz,M,L))
        kgood=-1
        for k in progressbar(range(Nz),vname+': ', 40):#range(Nz)
            (t2d,Nzgood) = interp_tracers(inputfile,vname,tndx_glo,k,coef,elem)
            if Nzgood>Nzgoodmin:
                kgood=kgood+1
                t3d[kgood,:,:]=t2d
                
        t3d=t3d[0:kgood,:,:]
        depth=depth[0:kgood]
        np.savez('t3d.npz',t3d=t3d,depth=depth)

    else:
        print('Load matrix...')
        data=np.load('t3d.npz')
        t3d = data['t3d']
        depth = data['depth']

    [Nz]=np.shape(depth)
    Z=-depth

#  Vertical interpolation
    print('Vertical interpolation')

# Add a layer below the bottom and above the surface to avoid vertical extrapolations
# and flip the matrices upside down (Z[Nz]=surface)

    Z=np.flipud(np.concatenate(([100.],Z,[-10000.])))
    [Nz]=np.shape(Z)
    t3d=np.flipud(add2layers(t3d))

# Do the vertical interpolations
    vout=ztosigma(t3d,Z,z_rho)
    return vout

####################################################################


def interp3d_uv(inputfile,tndx_glo,Nzgoodmin,z_rho,cosa,sina,\
        coefU,elemU,coefV,elemV):
#  1 - Horizontal interpolation on each GLORYS levels
#  2 - Vertical Interpolation from z to CROCO sigma levels

    [N,M,L]=np.shape(z_rho)
    depth=inputfile.depth
    dz=np.gradient(depth)
    [Nz]=np.shape(depth)
    comp_horzinterp = 1 # if 1 compute horizontal interpolations - 0 use saved matrices (for debugging)
    if comp_horzinterp==1:
        print('Horizontal interpolation of u and v over z levels')
        u3d=np.zeros((Nz,M,L-1))
        v3d=np.zeros((Nz,M-1,L))
        ubar=np.zeros((M,L-1))
        vbar=np.zeros((M-1,L))
        zu  =ubar
        zv  =vbar
        kgood=-1
        for k in progressbar(range(Nz),' uv : ', 40):
            (u2d,Nzgood_u) = interp_tracers(inputfile,'u',tndx_glo,k,coefU,elemU)
            (v2d,Nzgood_v) = interp_tracers(inputfile,'v',tndx_glo,k,coefV,elemV)
            Nzgood=np.min((Nzgood_u,Nzgood_v))
            if Nzgood>Nzgoodmin:
                kgood=kgood+1
#
# Rotation and put to u-points and v-points 
#
                u3d[kgood,:,:]=tools.rho2u(u2d*cosa+v2d*sina)
                v3d[kgood,:,:]=tools.rho2v(v2d*cosa-u2d*sina)

                ubar = ubar + tools.rho2u((u2d*dz[kgood])*cosa+(v2d*dz[kgood])*sina)
                zu   = zu   + dz[kgood]*np.ones(ubar.shape)
                vbar = vbar + tools.rho2v((v2d*dz[kgood])*cosa-(u2d*dz[kgood])*sina)
                zv   = zv   + dz[kgood]*np.ones(vbar.shape)


        u3d=u3d[0:kgood,:,:]
        v3d=v3d[0:kgood,:,:]
        ubar=ubar/zu
        vbar=vbar/zv
        depth=depth[0:kgood]
        np.savez('u3d.npz',u3d=u3d,v3d=v3d,depth=depth,ubar=ubar,vbar=vbar)

    else:
        print('Load matrices...')
        data=np.load('u3d.npz')
        u3d = data['u3d']
        v3d = data['v3d']
        depth = data['depth']
        ubar = data['ubar']
        vbar = data['vbar']

    [Nz]=np.shape(depth)
    Z=-depth
#--------------------------------------------------
#  Vertical interpolation
#----------------------------------------------------
#
    print('Vertical interpolation')
#
# Add a layer below the bottom and above the surface to avoid vertical extrapolations
# and flip the matrices upside down (Z[Nz]=surface)
#
    Z=np.flipud(np.concatenate(([100.],Z,[-10000.])))
    [Nz]=np.shape(Z)
    u3d=np.flipud(add2layers(u3d))
    v3d=np.flipud(add2layers(v3d))
#
# Do the vertical interpolations 
#
    uout=ztosigma(u3d,Z,tools.rho2u(z_rho))
    vout=ztosigma(v3d,Z,tools.rho2v(z_rho))

    return uout,vout,ubar,vbar



############################################################
###################Interp 4D ###############################
############################################################

def interp4d(inputfile,vname,Nzgoodmin,z_rho,coef,elem,bdy,dtmin,dtmax,prev,nxt):
#  Do a full interpolation of a 3d variable from GLORYS to a CROCO sigma grid
#  1 - Horizontal interpolation on each GLORYS levels
#  2 - Vertical Interpolation from z to CROCO sigma levels
    if np.ndim(z_rho)==3:
        z_rho=z_rho[np.newaxis,:]
    [T,N,M,L]=np.shape(z_rho)
    depth= inputfile.depth
    [Nz]=np.shape(depth)

    comp_horzinterp = 1 # if 1 compute horizontal interpolations - 0 use saved matrices (for debugging)
    if comp_horzinterp==1:
        print('Horizontal interpolation over z levels')
        t4d=np.zeros((T,Nz,M,L))
        kgood=-1
        for k in progressbar(range(Nz),vname+': ', 40):#range(Nz)
            (t3d,Nzgood) = interp_tracers3D(inputfile,vname,k,coef,elem,bdy,dtmin,dtmax,prev,nxt)
            if Nzgood>Nzgoodmin:
                kgood=kgood+1
                t4d[:,kgood,:,:]=t3d
                
        t4d=t4d[:,0:kgood,:,:]
        depth=depth[0:kgood]
        np.savez('t4d.npz',t4d=t4d,depth=depth)

    else:
        print('Load matrix...')
        data=np.load('t4d.npz')
        t3d = data['t4d']
        depth = data['depth']

    [Nz]=np.shape(depth)
    Z=-depth

#  Vertical interpolation
    print('Vertical interpolation')

# Add a layer below the bottom and above the surface to avoid vertical extrapolations
# and flip the matrices upside down (Z[Nz]=surface)

    Z=np.flipud(np.concatenate(([100.],Z,[-10000.])))
    [Nz]=np.shape(Z)
    t4d=np.flip(add2layers(t4d),axis=1)

# Do the vertical interpolations
    vout=ztosigma(t4d,Z,z_rho)
    return vout

####################################################################


def interp4d_uv(inputfile,Nzgoodmin,z_rho,cosa,sina,\
        coefU,elemU,coefV,elemV,bdy,dtmin,dtmax,prev,nxt):
#  1 - Horizontal interpolation on each GLORYS levels
#  2 - Vertical Interpolation from z to CROCO sigma levels
    if np.ndim(z_rho)==3:
        z_rho=z_rho[np.newaxis,:]
    [T,N,M,L]=np.shape(z_rho)

    cosa3d=np.tile(cosa, (T, 1, 1))
    sina3d=np.tile(sina, (T, 1, 1))

    depth=inputfile.depth
    dz=np.gradient(depth)
    [Nz]=np.shape(depth)
    comp_horzinterp = 1 # if 1 compute horizontal interpolations - 0 use saved matrices (for debugging)
    if comp_horzinterp==1:
        print('Horizontal interpolation of u and v over z levels')
        u4d=np.zeros((T,Nz,M,L-1))
        v4d=np.zeros((T,Nz,M-1,L))
        ubar=np.zeros((T,M,L-1))
        vbar=np.zeros((T,M-1,L))
        zu  =ubar
        zv  =vbar
        kgood=-1
        for k in progressbar(range(Nz),' uv : ', 40):
            (u3d,Nzgood_u) = interp_tracers3D(inputfile,'u',k,coefU,elemU,bdy,dtmin,dtmax,prev,nxt)
            (v3d,Nzgood_v) = interp_tracers3D(inputfile,'v',k,coefV,elemV,bdy,dtmin,dtmax,prev,nxt)
            Nzgood=np.min((Nzgood_u,Nzgood_v))
            if Nzgood>Nzgoodmin:
                kgood=kgood+1
#
# Rotation and put to u-points and v-points 
#
                u4d[:,kgood,:,:]=tools.rho2u(u3d*cosa3d+v3d*sina3d)
                v4d[:,kgood,:,:]=tools.rho2v(v3d*cosa3d-u3d*sina3d)

                ubar = ubar + tools.rho2u((u3d*dz[kgood])*cosa3d+(v3d*dz[kgood])*sina3d)
                zu   = zu   + dz[kgood]*np.ones(ubar.shape)
                vbar = vbar + tools.rho2v((v3d*dz[kgood])*cosa3d-(u3d*dz[kgood])*sina3d)
                zv   = zv   + dz[kgood]*np.ones(vbar.shape)


        u4d=u4d[:,0:kgood,:,:]
        v4d=v4d[:,0:kgood,:,:]
        ubar=ubar/zu
        vbar=vbar/zv
        depth=depth[0:kgood]
        np.savez('u4d.npz',u4d=u4d,v4d=v4d,depth=depth,ubar=ubar,vbar=vbar)

    else:
        print('Load matrices...')
        data=np.load('u4d.npz')
        u4d = data['u4d']
        v4d = data['v4d']
        depth = data['depth']
        ubar = data['ubar']
        vbar = data['vbar']

    [Nz]=np.shape(depth)
    Z=-depth
#--------------------------------------------------
#  Vertical interpolation
#----------------------------------------------------
#
    print('Vertical interpolation')
#
# Add a layer below the bottom and above the surface to avoid vertical extrapolations
# and flip the matrices upside down (Z[Nz]=surface)
#
    Z=np.flipud(np.concatenate(([100.],Z,[-10000.])))
    [Nz]=np.shape(Z)
    u4d=np.flip(add2layers(u4d),axis=1)
    v4d=np.flip(add2layers(v4d),axis=1)
#
# Do the vertical interpolations 
#
    uout=ztosigma(u4d,Z,tools.rho2u(z_rho))
    vout=ztosigma(v4d,Z,tools.rho2v(z_rho))

    return uout,vout,ubar,vbar
