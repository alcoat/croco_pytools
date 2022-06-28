import numpy as np
import xarray as xr

import gridop as gop
from scipy.sparse import lil_matrix


# relative vorticity
def relative_vorticity(model, ds=None, xgrid=None, u=None, v=None, f=None):
    
    """
    Compute relative vorticity normalized by f
    input: 
        model : instance of the Model class defined in the model.py module
        ds    : xarray DataSet
        xgrid : xgcm grid
        u     : xarray DataArray x-current component
        v     : xarray DataArray y-current component
    """
    
    xgrid = model.xgrid if xgrid is None else xgrid
    ds = model.ds if ds is None else ds    
    
    u = ds.u if u is None else u
    hgrid = gop.get_grid_point(u)
    if hgrid != 'u': u = gop.x2u(ds, u, xgrid)
        
    v = ds.v if v is None else v
    hgrid = gop.get_grid_point(v)
    if hgrid != 'v': v = gop.x2v(ds, v, xgrid)
        
    f = ds.f if f is None else f
    hgrid = gop.get_grid_point(f)
    if hgrid != 'f': f = gop.x2f(ds, f, xgrid)
    
    xi = ((-xgrid.derivative(u, 'y') 
           + xgrid.derivative(v, 'x')
          )/f).rename('vorticity')
    return xi.rename('relvort')


def ertel_pv(model, ds=None, xgrid=None, u=None, v=None, w=None, z=None, typ='ijk'):
    """
    #
    #   epv    - The ertel potential vorticity with respect to property 'lambda'
    #
    #                                       [ curl(u) + f ]
    #   -  epv is given by:           EPV = --------------- . del(lambda)
    #                                            rho
    #
    #   -  pvi,pvj,pvk - the x, y, and z components of the potential vorticity.
    #
    #   -  Ertel PV is calculated on horizontal rho-points, vertical w-points.
    #
    #
    #   tindex   - The time index at which to calculate the potential vorticity.
    #   depth    - depth
    #
    # Adapted from rob hetland.
    #
    """

    # Grid parameters
    ds = model.ds if ds is None else ds
    xgrid = model.xgrid if xgrid is None else xgrid

    for var in ['f', 'rho']:
        if var in ds:
            locals()[var] = ds[var]  
        else:
            print(var + ' not found in the dataset')
            return None
    rho0 = 1027 if 'rho0' not in ds else ds.rho0

    # 3D variables
    if z is  None : z = gop.get_z(model, ds)
    dz = xgrid.diff(z,'z')
    u = ds.xcur if u is None else u
    v = ds.ycur if v is None else v
    w = ds.zcur if w is None else w

    if 'k' in typ:

        # Ertel potential vorticity, term 1: [f + (dv/dx - du/dy)]*drho/dz
        # Compute d(v)/d(xi) at PSI-points.
        dvdxi = xgrid.derivative(v,'x')

        # Compute d(u)/d(eta) at PSI-points.
        dudeta = xgrid.derivative(u,'y')

        # Compute d(rho)/d(z) at horizontal RHO-points and vertical W-points
        drhodz = xgrid.diff(rho,'z') / dz

        #  Compute Ertel potential vorticity <k hat> at horizontal RHO-points and
        #  vertical W-points. 
        omega = dvdxi - dudeta
        omega = f + gop.x2rho(ds, omega, xgrid)
        pvk = xgrid.interp(omega,'z') * drhodz
        del dvdxi, dudeta, drhodz, omega
    else:
        pvk = 0.

    if 'i' in typ:

        #  Ertel potential vorticity, term 2: (dw/dy - dv/dz)*(drho/dx)
        #  Compute d(w)/d(y) at horizontal V-points and vertical RHO-points
        dwdy = xgrid.derivative(w,'y')

        #  Compute d(v)/d(z) at horizontal V-points and vertical W-points
        dz_v = xgrid.interp(dz,'y')
        dvdz = xgrid.diff(v,'z') / dz_v

        #  Compute d(rho)/d(xi) at horizontal U-points and vertical RHO-points
        drhodx = xgrid.derivative(rho,'x')

        #  Add in term 2 contribution to Ertel potential vorticity at horizontal RHO-points and
        #  vertical W-points.
        pvi = (gop.x2w(ds, dwdy, xgrid) - gop.x2w(ds, dvdz, xgrid)) * gop.x2w(ds, drhodx, xgrid)
        del dwdy, dz_v, dvdz, drhodx
    else:
        pvi = 0.

    if 'j' in typ:

        #  Ertel potential vorticity, term 3: (du/dz - dw/dx)*(drho/dy)
        #  Compute d(u)/d(z) at horizontal U-points and vertical W-points
        dz_u = xgrid.interp(dz, 'x')
        dudz = xgrid.diff(u,'z') / dz_u

        #  Compute d(w)/d(x) at horizontal U-points and vertical RHO-points
        dwdx = xgrid.derivative(w,'x')

        #  Compute d(rho)/d(eta) at horizontal V-points and vertical RHO-points
        drhodeta = xgrid.derivative(rho,'y')

        #  Add in term 3 contribution to Ertel potential vorticity at horizontal RHO-points and
        #  vertical W-points..
        pvj =  (gop.x2w(ds,dudz,xgrid)-gop.x2w(ds,dwdx,xgrid)) * gop.x2w(ds,drhodeta,xgrid)
        del dz_u, dudz, dwdx, drhodeta

    else:
        pvj = 0.

    #
    #
    # Sum potential vorticity components, and divide by rho0
    #
    pvi = pvi / rho0
    pvj = pvj / rho0
    pvk = pvk / rho0
    pv = pvi + pvj + pvk

    pv = pv.assign_coords(coords={"lon":rho.lon})
    pv = pv.assign_coords(coords={"lat":rho.lat})
    z_w = xgrid.interp(z, 'z') #.fillna(0.)
    pv = pv.assign_coords(coords={"z_w":z_w})

    return pv.squeeze().rename('pv')

def dtempdz(model, ds=None, xgrid=None, temp=None, z=None):
    """
    Compute dT/dz at horizontal rho point/vertical w point
    ds : dataset, containing T field
    z : xarray.DataArray, z in meters at rho points
    time : int, time index 
    """
    
    # Grid parameters
    if ds is None: ds = model.ds 
    if xgrid is None: xgrid = model.xgrid
    if z is  None: z = gop.get_z(model, ds=ds, z_sfc=ds.z_sfc)
    
    # Initialize temperature
    if temp is None: 
        try: 
            temp=ds.temp
        except:
            print("temp not not")
            return None

    # compute z coordinates at w level
    z_w = xgrid.interp(z,'z').squeeze()

    dtempdz = (xgrid.diff(temp,'z') / xgrid.diff(z,'z')).squeeze()
    if 'lon' in temp.coords: dtempdz = dtempdz.assign_coords(coords={"lon":temp.lon})
    if 'lat' in temp.coords: dtempdz = dtempdz.assign_coords(coords={"lat":temp.lat})
    dtempdz = dtempdz.assign_coords(coords={"z_w":z_w})
    return dtempdz.rename('dtdz')

def richardson(model, ds=None, u=None, v=None, rho=None, z=None, xgrid=None):
    """
    Ri is given by:      N²/((du/dz)² - (dv/dz)²)
         with N = sqrt(-g/rho0 * drho/dz)
    Ri is calculated at RHO-points and w level

    ds : dataset, which contains 3D u,v and rho fields
    z : xarray datarray, z depths at rho levels
    time : int, time index
    """

    # Grid parameters    
    if ds is None: ds=model.ds
    if xgrid is None: xgrid=model.xgrid
    try:
        rho0 = ds.rho0.values
    except:
        rho0 = 1027.
    try:
        g = ds.g.values 
    except:
        g = 9.81
    if u is None:
        try: 
            u=ds.xcur
        except:
            print("xcur not found in the dataset")
            return None
    if v is None:
        try: 
            v=ds.ycur
        except:
            print("ycur not found in the dataset")
            return None
    if rho is None:
        try: 
            rho=ds.rho
        except:
            print("rho not found in the dataset")
            return None

    # compute Z
    if z is None: z = gop.get_z(model, ds=ds, z_sfc=ds.z_sfc)
    z_w = xgrid.interp(z,'z').squeeze()

    N2 = get_N2(model, ds=ds, rho=rho, z=z, g=g)
    dudz = xgrid.diff(u,'z') / xgrid.diff(gop.x2u(ds,z,xgrid),'z')
    dvdz = xgrid.diff(v,'z') / xgrid.diff(gop.x2v(ds,z,xgrid),'z')

    Ri = xr.ufuncs.log10(N2 / (gop.x2w(ds,dudz,xgrid)**2 +  gop.x2w(ds,dvdz,xgrid)**2)).squeeze()
    if 'lon' in rho.coords: Ri = Ri.assign_coords(coords={"lon":rho.lon})
    if 'lat' in rho.coords: Ri = Ri.assign_coords(coords={"lat":rho.lat})
    Ri = Ri.assign_coords(coords={"z":z_w})
    return Ri.rename('Ri')

def get_N2(model, ds=None, rho=None, z=None, rho0=None, g=None, xgrid=None):
    """ Compute square buoyancy frequency N2 
    ... doc to be improved
    """
    if ds is None: ds=model.ds
    if xgrid is None: xgrid = model.xgrid
    try:
        rho0 = ds.rho0.values
    except:
        rho0 = 1027.
    try:
        g = ds.g.values
    except:
        g = 9.81
    if rho is None: 
        try: 
            rho=ds.rho
        except:
            print("rho not found in the dataset")
            return None
    print(g)
    N2 = -g/rho0 * xgrid.diff(rho, 'z', boundary='fill', fill_value=np.NaN) \
            / xgrid.diff(z, 'z', boundary='fill', fill_value=np.NaN)
    # cannot find a solution with xgcm, weird
    N2 = N2.fillna(N2.shift(s_w=-1))
    N2 = N2.fillna(N2.shift(s_w=1))
    return N2.rename('N2')

def get_streamfunction(model,pm,pn,pv,verbo=False):
    """
    Compute the stream function from the relative vorticity
    Invert the laplacian to solve the poisson equation Ax=b
    A is the horizontal laplacian, b is the vorticity
    Input:
        - pm : (DataArray) 1/dx metric
        - pn : (DataArray) 1/dy metric
        - pv : (DataArray) relative vorticity
        - verbo : (Boolean) verbose mode
    Output:
        (DataArray) the computed streamfunction 
    """

    if np.any(np.isnan(pv)): 
        print("Can't inverse the laplacian, non compact domain, pv contains nan values")
        return None

    #######################################################
    #Create matrix A
    #######################################################
    if verbo: print('creating matrix A')
    A = gop.poisson_matrix(pm.values,pn.values)

    #######################################################
    #Solve matrix A
    A = A.tocsr()
    #######################################################

    if verbo: print('creating matrix b')
    b = -1. * pv.values.flatten() # right hand side
    ml = ruge_stuben_solver(A)                # construct the multigrid hierarchy
    if verbo: print(ml)                             # print hierarchy information
    x = ml.solve(b, tol=1e-8)                       # solve Ax=b to a tolerance of 1e-8     
    #x = solve(A,b,verb=False,tol=1e-8)

    if verbo: print("residual: ", np.linalg.norm(b-A*x))          # compute norm of residual vector

    chi = xr.DataArray(
        data=x.reshape(pm.shape),
        dims=["y_rho", "x_rho"],
        coords={'nav_lon_rho':pv.nav_lon_rho, 'nav_lat_rho':pv.nav_lat_rho}
        )
    return chi.rename('streamfct')


def get_p(model, rho, z_w, z_r, ds=None, g=None, rho0=None, xgrid=None):
    """ 
    Compute (not reduced) pressure by integration from the surface, 
    taking rho at rho points and giving results on rho points (z grid).

    Parameters
    ----------
    grid : xgcm grid
    rho  : Density (DataArray)
    z_w  : depth at vertical w points (DataArray)
    z_r  : depth at vertical rho points (DataArray)
    g    : acceleration of gravity (float)
    rho0 : mean density (float)

    """
    # useful parameters
    eps = 1.0e-10
    if ds is None: ds=model.ds
    if xgrid is None: xgrid=model.xgrid
        
    if g is None: 
        try:
            g=ds.g.values
        except:
            g=9.81
    if rho0 is None: 
        try:
            rho0=ds.rho0.values
        except:
            rho0=1000.
                
    GRho=g/rho0
    HalfGRho=0.5*GRho
    N = rho.sizes['s']
    OneFifth = 1.0/5.0
    OneTwelfth = 1.0/12.0

    # dz and drho on w levels
    dR = xgrid.diff(rho,'z', boundary='extend').rename('dR')
    dZ = xgrid.diff(z_r,'z', boundary='extend').rename('dZ')

    # modified dz and dr on w levels
    dZi = 2. * ( dZ * dZ.shift(s_w=1,fill_value=0) 
            / (dZ + dZ.shift(s_w=1,fill_value=0)) )
    dZi = xr.concat([dZ.isel(s_w=0), dZi.isel(s_w=slice(1,None))], dim='s_w')
    dRi = 2. * ( dR * dR.shift(s_w=1,fill_value=0) 
            / (dR + dR.shift(s_w=1,fill_value=0)) )
    dRi = xr.concat([dR.isel(s_w=0), dRi.isel(s_w=slice(1,None))], dim='s_w')

    # Pressure at the surface on rho level
    Pn = (g*z_w.isel(s_w=-1) + GRho*( rho.isel(s=-1) 
                                   + 0.5*(rho.isel(s=-1)-rho.isel(s=-2)) 
                                   * (z_w.isel(s_w=-1)-z_r.isel(s=-1)) 
                                   / (z_r.isel(s=-1)-z_r.isel(s=-2))
                                  ) * (z_w.isel(s_w=-1)-z_r.isel(s=-1))
         )

    # Pressure of each slice
    rhoz = (rho + rho.shift(s=-1))*(z_r.shift(s=-1) - z_r)
    dRz = ((xgrid.diff(dRi,'z'))*(z_r.shift(s=-1) - z_r
                                 -OneTwelfth*(xgrid.interp(dZi, 'z'))) )
    dZr = ((xgrid.diff(dZi,'z'))*(rho.shift(s=-1) - rho
                                 -OneTwelfth*(xgrid.interp(dRi,'z'))) )
    Pk = (HalfGRho*( rhoz - OneFifth*(dRz - dZr)))

    # replace pressure at the surface
    Pk = xr.concat([Pk.isel(s=slice(0,-1)), Pn], dim='s')

    # integrate from the surface to bottom
    P = (Pk
            .sortby(Pk.s, ascending=False)                
            .cumsum("s")
            .sortby(Pk.s, ascending=True)
            .assign_coords(z_r=z_r)
        )

    return P.rename('P')