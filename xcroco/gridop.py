import numpy as np
import xarray as xr
from functools import partial
import pandas as pd
from xgcm import Grid
import intake

import time
import os.path as path
import itertools
from collections import OrderedDict

from math import radians, cos, sin, asin, sqrt

def open_files(model, gridname, filenames, 
               grid_metrics=1,
               drop_variables=[], 
               time=None,
               chunks={'t':1},
               suffix='', 
               remove_ghost_pts=True,
               xperiodic=False,
               yperiodic=False,
               verbose=False
              ):
    """
    open Netcdf files or zarr archive
    input: 
        model : instance of the Model class defined in the model.py module
        gridname : path to the grid file
        filenames : path to the Netcdf files or to the zarr archive
        grid_metrics : type of xgcm grid 
                       0: no metrics
                       1: horizontal metrics
                       2: horizontal and vertical metrics
        drop_variables : list of variables to drop
        chunks: chunks to override those in the yaml file
        suffix: to remove suffixes in dimensions or variables
        verbose: verbose mode
    outputs:
        ds: an xarray dataset
        grid: the associated xgcm grid
                       
    """
    if time is not None:
        def _preprocess(ds,time):
            return ds.sel(time_counter=time, method='nearest')
        partial_func = partial(_preprocess, time=time)
    else:
        partial_func = None
              
    # convert filenames to list of strings if string
    filenames = filenames.tolist() if isinstance(filenames,str) else filenames
    # find time dimension name
    concat_dim = [k for k,v in model.rename_vars.items() if v == 't'][0]
    open_kwargs = {'concat_dim': concat_dim,
                   'combine': 'nested',
                   'coords': 'minimal',
                   'parallel': True,
                   'compat': 'override'
                  }
    try : 
        # zarr archive
        ds = xr.open_zarr(filenames[0], drop_variables=drop_variables)
    except:
        try : 
            # list of explicit filenames
            ds = xr.open_mfdataset(filenames, drop_variables=drop_variables,
                                   preprocess=partial_func, **open_kwargs)  
        except :
            try :
                # list of files with wildcards
                ds = xr.open_mfdataset(filenames[0], drop_variables=drop_variables,
                                       preprocess=partial_func, **open_kwargs)
            except:
                print('open_files: unknown format: only Netcdf or Zarr')
                print('or filenames do not exist')
                return
        
    # delete suffix from dimension or variable names
    if suffix != '' : ds = del_name_suffix(ds,suffix)
    # change the names in the dataset according to the model instance
    model.ds = adjust_grid(model, ds)
    
    # add the grid and the xgcm grid to the dataset
    ds, grid = add_grid(model, gridname, grid_metrics=grid_metrics, suffix=suffix,
                        remove_ghost_pts=remove_ghost_pts,xperiodic=xperiodic, yperiodic=yperiodic)
    model.ds = ds.chunk(chunks=chunks).squeeze()
    return model.ds, grid

def open_catalog(model, gridname, catalog, source=None,
                 grid_metrics=1,
                 chunks={},
                 suffix='', 
                 remove_ghost_pts=True,
                 xperiodic=False,
                 yperiodic=False,
                 verbose=False,
                ):
    """
    open files through an intake catalog
    input: 
        model : instance of the Model class defined in the model.py module
        gridname : path to the grid file
        catalog : path to the intake yaml catalog
        source : source to open in the catalog (if None, the first)
        grid_metrics : type of xgcm grid 
                       0: no metrics
                       1: horizontal metrics
                       2: horizontal and vertical metrics
        chunks: chunks to override those in the yaml file
        suffix: to remove suffixes in dimensions or variables
        verbose: verbose mode
    outputs:
        ds: an xarray dataset
        grid: the associated xgcm grid
                       
    """
    
    try : 
        # open intake catalog
        cat = intake.open_catalog(catalog)
        if verbose: print('Available sources are: ',list(cat))
    except:
        print('open_catalog: yaml catalog not found')
        return
    
    # find the first source of the catalog if source is None
    source = list(cat)[0] if source is None else source
    if verbose: 
        print('Source :', source)
        print('   ', cat[source])
    # open the source as a dataset
    try:
        ds = cat[source].to_dask() # chunks={'time_counter': 50})
    except:
        print('May be the source is not found in the yaml catalog')
        return
    
    # delete suffix from dimension or variable names
    if suffix != '' : ds = del_name_suffix(ds,suffix)
    # change the names in the dataset according to the model instance
    model.ds = adjust_grid(model, ds)
    
    # add the grid and the xgcm grid to the dataset
    ds, grid = add_grid(model, gridname, grid_metrics=grid_metrics, suffix=suffix,
                        remove_ghost_pts=remove_ghost_pts,xperiodic=xperiodic, yperiodic=yperiodic)
    model.ds = ds.chunk(chunks=chunks)
    return model.ds.squeeze(), grid

def del_name_suffix(ds,suffix):
    """ remove a suffix from dimensions, coordinates and variables in the dataset """
    var_suffix = [v for v in itertools.chain(ds.dims,
                                          ds.data_vars.keys(),
                                          ds.coords.keys(),
                                          ) if suffix in v]
    if var_suffix:
        # remove duplicates
        var_suffix = list(set(var_suffix))
        var = [v.replace('_'+suffix,'').replace(suffix,'') for v in var_suffix]
        for (old,new)in zip(var_suffix,var):
            ds = ds.rename({old: new})
    return ds

def find_var(model,varname,ds,gd):
    
    def good_type(var):
        if isinstance(var,xr.DataArray) or isinstance(var,np.ndarray) or \
        isinstance(var,np.float32) or isinstance(var,np.float64) or \
        isinstance(var,np.int32) or isinstance(var,np.float64):
            return True
        else:
            return False
        
    def to_dataarray(model,varname, var):
        if  isinstance(var,np.ndarray):
            var = xr.DataArray(data=var,
                            dims=model.dims_var[varname]
            )
                         
        return var
        
    if varname in gd and good_type(gd[varname]):
        return to_dataarray(model,varname,gd[varname])
    elif varname in ds and good_type(ds[varname]):
        return to_dataarray(model,varname,ds[vvarnamear])
    elif varname in gd.attrs and good_type(gd.attrs[varname]):
        return to_dataarray(model,varname,gd.attrs[varname])
    elif varname in ds.attrs and good_type(ds.attrs[varname]):
        return to_dataarray(model,varname,ds.attrs[varname])
    else:
        return None

def get_cs(model, ds, gd, vgrid):
    """ get croco vertical grid  stretching 
    https://www.myroms.org/wiki/Vertical_S-coordinate
    """
    # search vtransform 
    if 'vtransform' in gd:
        vtransform = gd.vtransform.values
    elif 'vtransform' in ds:
        vtransform = ds.vtransform.values
    elif 'vtransform' in gd.attrs:
        vtransform = gd.attrs['vtransform']
    elif 'vtransform' in ds.attrs:
        vtransform = ds.attrs['vtransform']
    else:
        print("Can't find vtransform neither in filename nor in gridname  ")
        return None
    
    # search theta_s 
    if find_var(model,'theta_s',ds,gd) is not None: 
        theta_s = find_var(model,'theta_s',ds,gd).values
    else:
        print("Can't find theta_s neither in filename nor in gridname  ")
        return None
    # search theta_b 
    if find_var(model,'theta_b',ds,gd) is not None: 
        theta_b = find_var(model,'theta_b',ds,gd).values
    else:
        print("Can't find theta_b neither in filename nor in gridname  ")
        return None
    
    sc = ds.sc_r.values if vgrid=='r' else ds.sc_w.values
    
    # New coordinates
    if vtransform == 2 or vtransform.lower()=='new':
        if theta_s>0:
            csf = (1-np.cosh(theta_s*sc)) / (np.cosh(theta_s)-1.)
        else:
            csf = sc**2
        if theta_b>0:
            cs = (np.exp(theta_s*csf)-1.) / (1-np.exp(-theta_b))
        else:
            cs = csf
    # Old coordinates
    else:
        cs = (1-theta_b)*np.sinh(theta_s*sc)/np.sinh(theta_s) \
             + theta_b*0.5 \
               *(np.tanh((sc+0.5)*theta_s)/np.tanh(0.5*theta_s)-1.)
    return cs

def add_grid(model, gridname, grid_metrics=1, suffix='', remove_ghost_pts=True,
             xperiodic=False, yperiodic=False):
        
    # open grid file
    try : 
        gd = xr.open_zarr(gridname).squeeze()
    except:
        try : 
            gd = xr.open_dataset(gridname).squeeze() 
        except :
            print('add_grid: unknown format for grid : only Netcdf or Zarr')
            
    # Rename variable according model
    if suffix != '' : gd = del_name_suffix(gd,suffix)
    gd = adjust_grid(model, gd)
    ds = model.ds

    if find_var(model,'hc',ds,gd) is not None: ds['hc'] = find_var(model,'hc',ds,gd)
    if find_var(model,'h',ds,gd) is not None: ds['h'] = find_var(model,'h',ds,gd)
    if find_var(model,'pm',ds,gd) is not None: ds['pm']   = find_var(model,'pm',ds,gd)
    if find_var(model,'pn',ds,gd) is not None: ds['pn']   = find_var(model,'pn',ds,gd)
    try:
        N = ds.dims['s']
        if 'sc_r' not in ds:
            if find_var(model,'sc_r',ds,gd) is not None: 
                ds['sc_r'] = find_var(model,'sc_r',ds,gd)
            else:
                ds['sc_r'] = xr.DataArray(np.arange(-1.+1./(N+1),0., 1/(N+1)), dims='s')  
        if 'sc_w' not in ds:      
            if find_var(model,'sc_w',ds,gd) is not None: 
                ds['sc_w'] = find_var(model,'sc_w',ds,gd)
            else:
                ds['sc_w'] = xr.DataArray(np.arange(-1.,0., 1/(N+1)), dims='s_w')
        if 'Cs_r' not in ds:
            if  find_var(model,'Cs_r',ds,gd) is not None: 
                ds['Cs_r'] = find_var(model,'Cs_r',ds,gd)
            else:
                ds['Cs_r'] = get_cs(model,ds, gd, 'r')
        if 'Cs_w' not in ds:
            if  find_var(model,'Cs_w',ds,gd) is not None: 
                ds['Cs_w'] = find_var(model,'Cs_w',ds,gd)
            else:
                ds['Cs_w'] = get_cs(model,ds, gd, 'w')
    except:
        pass        
        
    if find_var(model,'angle',ds,gd) is not None: ds['angle'] = find_var(model,'angle',ds,gd)
    if find_var(model,'mask',ds,gd) is not None: ds['mask'] = find_var(model,'mask',ds,gd)
    if find_var(model,'lon',ds,gd) is not None: ds['lon'] = find_var(model,'lon',ds,gd)
    if find_var(model,'lat',ds,gd) is not None: ds['lat'] = find_var(model,'lat',ds,gd)
    if find_var(model,'f',ds,gd) is not None: ds['f'] = find_var(model,'f',ds,gd)
    if find_var(model,'rho0',ds,gd) is not None: ds['rho0'] = find_var(model,'rho0',ds,gd)
    if find_var(model,'g',ds,gd) is not None: ds['g'] = find_var(model,'g',ds,gd)
    
    
    # coords = [c for c in ds.coords if c not in ['t','s','s_w']]
    coords = [c for c in ds.coords if c in ['t','s','s_w','lat','lon']]
    ds = ds.reset_coords()
    ds = ds.set_coords(coords)
        
    # remove ghost points
    if remove_ghost_pts:
        ds = remove_ghost_points(model, ds, xperiodic=xperiodic, yperiodic=yperiodic)
    model.ds = ds
    
    # On crée la grille xgcm
    ds, grid = xgcm_grid(model, grid_metrics=grid_metrics, 
                         xperiodic=xperiodic, yperiodic=yperiodic)
    
    return ds, grid

def remove_ghost_points(model, ds, xperiodic=False, yperiodic=False):
    """
    remove ghost points
    """
    ds = ds.isel(x=slice(1,-1),y=slice(1,-1))
    if xperiodic:
        ds = ds.isel(x_u=slice(0,-1))
    if yperiodic:
        ds = ds.isel(y_v=slice(0,-1))
    return ds

def xgcm_grid(model, grid_metrics=1, xperiodic=False, yperiodic=False):
        
        # Create xgcm grid without metrics
        coords={}
        if all(d in model.ds.dims for d in ['x','x_u']):
            if xperiodic:
                coords.update({'x': {'center':'x', 'left':'x_u'}})
            else:
                coords.update({'x': {'center':'x', 'outer':'x_u'}})
        if all(d in model.ds.dims for d in ['y','y_v']):
            if yperiodic:
                coords.update({'y': {'center':'y', 'left':'y_v'}} )
            else:
                coords.update({'y': {'center':'y', 'outer':'y_v'}} )
        if 's' in model.ds.dims:
            coords.update({'z': {'center':'s', 'outer':'s_w'}})
            
        grid = Grid(model.ds, 
                  coords=coords,
                  periodic=False,
                  boundary='extend')
        
        if grid_metrics==0:           
            model.xgrid = grid
            return model.ds, grid
        
        # compute horizontal coordinates

        ds = model.ds
        if 'lon_u' not in ds and 'x_u' in ds.dims: ds['lon_u'] = grid.interp(ds.lon,'x')
        if 'lat_u' not in ds and 'y'   in ds.dims: ds['lat_u'] = grid.interp(ds.lat,'x')
        if 'lon_v' not in ds and 'x'   in ds.dims: ds['lon_v'] = grid.interp(ds.lon,'y')
        if 'lat_v' not in ds and 'y_v' in ds.dims: ds['lat_v'] = grid.interp(ds.lat,'y')
        if 'lon_f' not in ds and 'x_u' in ds.dims: ds['lon_f'] = grid.interp(ds.lon_v,'x')
        if 'lat_f' not in ds and 'y_v' in ds.dims: ds['lat_f'] = grid.interp(ds.lat_u,'y')
        _coords = [d for d in ds.data_vars.keys() if d.startswith(tuple(['lon','lat']))]
        ds = ds.set_coords(_coords)
        
        
        # add horizontal metrics for u, v and psi point
        if 'pm' in ds and 'pn' in ds:
            ds['dx'] = 1/ds['pm']
            ds['dy'] = 1/ds['pn']
        else: # backward compatibility, hack
            dlon = grid.interp(grid.diff(ds.lon,'x'),'x')
            dlat =  grid.interp(grid.diff(ds.lat,'y'),'y')
            ds['dx'], ds['dy'] = dll_dist(dlon, dlat, ds.lon, ds.lat)
        # dlon = grid.interp(grid.diff(ds.lon_u,'x'),'x')
        # dlat = grid.interp(grid.diff(ds.lat_u,'y'),'y')
        # ds['dx_u'], ds['dy_u'] = dll_dist(dlon, dlat, ds.lon_u, ds.lat_u)
        # dlon = grid.interp(grid.diff(ds.lon_v,'x'),'x')
        # dlat = grid.interp(grid.diff(ds.lat_v,'y'),'y')
        # ds['dx_v'], ds['dy_v'] = dll_dist(dlon, dlat, ds.lon_v, ds.lat_v)
        # dlon = grid.interp(grid.diff(ds.lon_f,'x'),'x')
        # dlat = grid.interp(grid.diff(ds.lat_f,'y'),'y')
        # ds['dx_psi'], ds['dy_psi'] = dll_dist(dlon, dlat, ds.lon_f, ds.lat_f)
        ds["dx_u"] = grid.interp(ds["dx"], "x")
        # ds["dy_u"] = grid.interp(ds["dy"], "x")
        # ds["dx_v"] = grid.interp(ds["dx"], "y")
        ds["dy_v"] = grid.interp(ds["dy"], "y")
        # ds["dx_psi"] = grid.interp(ds["dx_u"], "y")
        # ds["dy_psi"] = grid.interp(ds["dy_v"], "x")

        # add areas metrics for rho,u,v and psi points
        # ds['rAr'] = ds.dx_psi * ds.dy_psi
        # ds['rAu'] = ds.dx_v * ds.dy_v
        # ds['rAv'] = ds.dx_u * ds.dy_u
        # ds['rAf'] = ds.dx * ds.dy
        
        metrics={}    
        # if all(d in model.ds.dims for d in ['x','x_u']):
            # metrics.update({('x',): ['dx', 'dx_u', 'dx_v', 'dx_psi']})
        metrics.update({('x',): ['dx', 'dx_u']})
        # if all(d in model.ds.dims for d in ['y','y_v']):
            # metrics.update({('y',): ['dy', 'dy_u', 'dy_v', 'dy_psi']})
        metrics.update({('y',): ['dy', 'dy_v']})
        # if all(d in model.ds.dims for d in ['x','x_u','y','y_v']):
            # metrics.update({('x', 'y'): ['rAr', 'rAu', 'rAv', 'rAf']})
         
        if grid_metrics==1:
            # generate xgcm grid
            grid = Grid(ds,
                        coords=coords,
                        periodic=False,
                        metrics=metrics,
                        boundary='extend')
            model.xgrid = grid
            model.ds = ds
            return ds, grid
        
        # compute z coordinate at rho/w points
        if 'z_sfc' in [v for v in ds.data_vars] and \
           's' in [d for d in ds.dims.keys()] and \
            ds['s'].size>1:
            ds['is3D'] = True
            z = get_z(model, z_sfc=ds.z_sfc, xgrid=grid).fillna(0.)
            z_w = get_z(model, z_sfc=ds.z_sfc, xgrid=grid, vgrid='w').fillna(0.)
            ds['z'] = z
            ds['z_w'] = z_w
            ds['z_u'] = grid.interp(z,'x')
            ds['z_v'] = grid.interp(z,'y')
            ds['z_f'] = grid.interp(ds.z_u,'y')
            # set as coordinates in the dataset
            _coords = ['z','z_w','z_u','z_v','z_f']
            ds = ds.set_coords(_coords)
        else:
            ds['is3D'] = False

        # add vertical metrics for u, v, rho and psi points
        # if 'z' in [v for v in ds.coords]:
        if ds['is3D']:
            ds['dz'] = grid.interp(grid.diff(z,'z'),'z')
            ds['dz_w'] = grid.interp(grid.diff(z_w,'z'),'z')
            # ds['dz'] = grid.interp(grid.diff(ds.z,'z'),'z')
            # ds['dz_w'] = grid.interp(grid.diff(ds.z_w,'z'),'z')
            # ds['dz_u'] = grid.interp(grid.diff(ds.z_u,'z'),'z')
            # ds['dz_v'] = grid.interp(grid.diff(ds.z_v,'z'),'z')
            # ds['dz_f'] = grid.interp(grid.diff(ds.z_f,'z'),'z')
            
        # add coords and metrics for xgcm for the vertical direction
        # if 'z' in ds:
        if ds['is3D']:
# #             coords.update({'z': {'center':'s', 'outer':'s_w'}})
#             metrics.update({('z',): ['dz', 'dz_u', 'dz_v', 'dz_f', 'dz_w']}), # Z distances
            metrics.update({('z',): ['dz', 'dz_w']}), # Z distances
        # generate xgcm grid
        grid = Grid(ds,
                    coords=coords,
                    periodic=False,
                    metrics=metrics,
                    boundary='extend')

        model.xgrid = grid
        model.ds = ds

        return ds, grid

def fast_xgcm_grid(ds, grid_metrics=1, xperiodic=False, yperiodic=False):
    
    """
    Create the xgcm grid without computing any metrics. Just use those which are already 
    in the dataset.
    Input:
        ds: (Xarray Dataset) the dataset to create the xgcm grid
        grid_metrics: (integer) 0:no metrics, 1:horizontal metrics, 2:hor + vert metrics
        xperiodic: True if ds periodic in x
        yperiodic: True if ds periodic in y
        
    Return:
        grid: the xgcm grid
        
    """
    
    # Create xgcm grid without metrics
    coords={}
    if all(d in ds.dims for d in ['x','x_u']):
        if xperiodic:
            coords.update({'x': {'center':'x', 'left':'x_u'}})
        else:
            coords.update({'x': {'center':'x', 'outer':'x_u'}})
    if all(d in ds.dims for d in ['y','y_v']):
        if yperiodic:
            coords.update({'y': {'center':'y', 'left':'y_v'}} )
        else:
            coords.update({'y': {'center':'y', 'outer':'y_v'}} )
    if all(d in ds.dims for d in ['s','s_w']):
        coords.update({'z': {'center':'s', 'outer':'s_w'}})
        
    grid = Grid(ds, 
              coords=coords,
              periodic=False,
              boundary='extend')

    if grid_metrics==0: return grid

    # set all lon/lat variables as coordinates
    _coords = [d for d in ds.data_vars.keys() if d.startswith(tuple(['lon','lat']))]
    ds = ds.set_coords(_coords)
             
    # set horizontal metrics
    # move horizontal metrics from global attributes to variables
    attrs = [k for k in ds.attrs.keys() if k.startswith(('dx','dy','rA'))]
    if attrs is not None:
        for k in attrs: ds[k] = ds.attrs[k]
    metrics={}   
    # add dx metrics
    if all(d in ds.dims for d in ['x','x_u']):
        dx = [v for v in ds.data_vars if v in ['dx','dx_u','dx_v','dx_psi']]
        metrics.update({('x',): dx})
    # add dy metrics
    if all(d in ds.dims for d in ['y','y_v']):
        dy = [v for v in ds.data_vars if v in ['dy','dy_u','dy_v','dy_psi']]
        metrics.update({('y',): dy})
    # add area metrics
    if all(d in ds.dims for d in ['x','x_u','y','y_v']):        
        rA = [v for v in ds.data_vars if v in ['rAr','rAu','rAv','rAf']]
        metrics.update({('x', 'y'): rA})
    
    if grid_metrics==1:
        # generate xgcm grid
        grid = Grid(ds,
                    coords=coords,
                    periodic=False,
                    metrics=metrics,
                    boundary='extend')
        return grid

    # Set z variables as coordinates
    _coords = [d for d in ds.data_vars.keys() if d in ['z','z_w','z_u','z_v','z_f']]
    ds = ds.set_coords(_coords)

    # add vertical metrics
    # move vertical metrics from global attributes to variables
    attrs = [k for k in ds.attrs.keys() if k.startswith(('dz'))]
    if attrs is not None:
        for k in attrs: ds[k] = ds.attrs[k]
    # add dz metrics
    if all(d in ds.dims for d in ['s','s_w']):        
        dz = [v for v in ds.data_vars if v in ['dz', 'dz_u', 'dz_v', 'dz_f', 'dz_w']]
        metrics.update({('z',): dz})
    
    # generate xgcm grid
    grid = Grid(ds,
                coords=coords,
                periodic=False,
                metrics=metrics,
                boundary='extend')

    return grid
    
def dll_dist(dlon, dlat, lon, lat):
    """
    Converts lat/lon differentials into distances in meters
    PARAMETERS
    ----------
    dlon : xarray.DataArray longitude differentials 
    dlat : xarray.DataArray latitude differentials 
    lon : xarray.DataArray longitude values
    lat : xarray.DataArray latitude values
    RETURNS
    -------
    dx : xarray.DataArray distance inferred from dlon 
    dy : xarray.DataArray distance inferred from dlat 
    """
    distance_1deg_equator = 111000.0
    # dx = dlon * xr.ufuncs.cos(xr.ufuncs.deg2rad(lat)) * distance_1deg_equator 
    dx = dlon * np.cos(np.deg2rad(lat)) * distance_1deg_equator 
    dy = ((lon * 0) + 1) * dlat * distance_1deg_equator
    return dx, dy



def adjust_grid(model, ds):
    ''' 
    Change the names in the dataset according to the model
    Input : model: Instance of the model class
            ds : dataset to change
    Output : changed dataset
    '''
    # change names in dims, coordinates and variables
    for k,v in model.rename_vars.items():
        if (k in ds and v not in ds) or \
            k in ds.dims.keys():
            ds = ds.rename({k: v})
    # change names in attributes
    for k,v in model.rename_vars.items():
        if (k in ds.attrs and v not in ds.attrs):
            ds.attrs[v] = ds.attrs.pop(k)
    return ds
    

def get_spatial_dims(v):
    """ Return an ordered dict of spatial dimensions in the s/z, y, x order
    """
    if isinstance(v, xr.DataArray):
        dims = OrderedDict( (d, next((x for x in v.dims if x[0]==d), None))
                        for d in ['s','y','x'] )
    elif isinstance(v, xr.Dataset):
        # dims = OrderedDict( (d, next((x for x in v.dims if x==d), None))
        dims = OrderedDict( (d, [x for x in v.dims if x[0]==d])
                        for d in ['s','y','x'] ) 
        # convert empty list to None
        dims = {k: None if not d else d for k, d in dims.items() }
    else:
        print('get_spatial_dims: ERROR!!! the argument must be a DataArray or a Dataset')
    return dims


def get_spatial_coords(v):
    """ Return an ordered dict of spatial dimensions in the s/z, y, x order
    """
    coords = OrderedDict( (d, next((x for x in v.coords if x.startswith(d)), None))
                       for d in ['z','lat','lon'] )
    for k,c in coords.items():
        if c is not None and v.coords[c].size==1: coords[k]= None
    return coords

def reorder_dims(da):    
    # reorder spatial dimensions and place them last
    sdims = list(get_spatial_dims(da).values())
    sdims = tuple(filter(None,sdims)) # delete None values
    reordered_dims = tuple(d for d in da.dims if d not in sdims) + sdims
    return da.transpose(*reordered_dims, transpose_coords=True)

def x2rho(ds,v, grid):
    """ Interpolate from any grid to rho grid
    """
    dims = get_spatial_dims(v)
    coords = get_spatial_coords(v)
    vout = v.copy()
    if coords['lon']: lonout = v[coords['lon']].copy()
    if coords['lat']: latout = v[coords['lat']].copy()
    if coords['z']: zout = v[coords['z']].copy()
    if dims['x'] == 'x_u' and v.x_u.size>1:
        vout = grid.interp(vout, 'x')
        if coords['lon'] and 'x_u' in lonout.dims: lonout = grid.interp(lonout, 'x')
        if coords['lat'] and 'x_u' in latout.dims: latout = grid.interp(latout, 'x')
        if coords['z'] and 'x_u' in zout.dims: zout = grid.interp(zout, 'x')
    if dims['y'] == 'y_v' and v.y_v.size>1:
        vout = grid.interp(vout, 'y')
        if coords['lon'] and 'y_v' in lonout.dims: lonout = grid.interp(lonout, 'y')
        if coords['lat'] and 'y_v' in latout.dims: latout = grid.interp(latout, 'y')
        if coords['z'] and 'y_v' in zout.dims: zout = grid.interp(zout, 'y')
    if dims['s'] == 's_w' and v.s_w.size>1:
        vout = grid.interp(vout, 'z')
        if coords['z'] and 's_w' in zout.dims: zout = grid.interp(zout, 'z')
    # assign coordinates
    if coords['lon']: vout = vout.assign_coords(coords={'lon':lonout})
    if coords['lat']: vout = vout.assign_coords(coords={'lat':latout})
    if coords['z']: vout = vout.assign_coords(coords={'z':zout})

    return vout

def x2u(ds,v, grid):
    """ Interpolate from any grid to u grid
    """
    dims = get_spatial_dims(v)
    coords = get_spatial_coords(v)
    vout = v.copy()
    if coords['lon']: lonout = v[coords['lon']].copy()
    if coords['lat']: latout = v[coords['lat']].copy()
    if coords['z']: zout = v[coords['z']].copy()
    if dims['x'] == 'x' and v.x.size>1:
        vout = grid.interp(vout, 'x')
        if coords['lon'] and 'x' in lonout.dims: lonout = grid.interp(lonout, 'x')
        if coords['lat'] and 'x' in latout.dims: latout = grid.interp(latout, 'x')
        if coords['z'] and 'x' in zout.dims: zout = grid.interp(zout, 'x')
    if dims['y'] == 'y_v' and v.y_v.size>1:
        vout = grid.interp(vout, 'y')
        if coords['lon'] and 'y_v' in lonout.dims: lonout = grid.interp(lonout, 'y')
        if coords['lat'] and 'y_v' in latout.dims: latout = grid.interp(latout, 'y')
        if coords['z'] and 'y_v' in zout.dims: zout = grid.interp(zout, 'y')
    if dims['s'] == 's_w' and v.s_w.size>1:
        vout = grid.interp(vout, 'z')
        if coords['z'] and 's_w' in zout.dims: zout = grid.interp(zout, 'z')
    # assign coordinates
    if coords['lon']: vout = vout.assign_coords(coords={'lon_u':lonout})
    if coords['lat']: vout = vout.assign_coords(coords={'lat_u':latout})
    if coords['z']:vout = vout.assign_coords(coords={'z_u':zout})
    return vout

def x2v(ds,v, grid):
    """ Interpolate from any grid to v grid
    """
    dims = get_spatial_dims(v)
    coords = get_spatial_coords(v)
    vout = v.copy()
    if coords['lon']: lonout = v[coords['lon']].copy()
    if coords['lat']: latout = v[coords['lat']].copy()
    if coords['z']: zout = v[coords['z']].copy()
    if dims['x'] == 'x_u' and v.x_u.size>1:
        vout = grid.interp(vout, 'x')
        if coords['lon'] and 'x_u' in lonout.dims: lonout = grid.interp(lonout, 'x')
        if coords['lat'] and 'x_u' in latout.dims: latout = grid.interp(latout, 'x')
        if coords['z'] and 'x_u' in zout.dims: zout = grid.interp(zout, 'x')
    if dims['y'] == 'y' and v.y.size>1:
        vout = grid.interp(vout, 'y')
        if coords['lon'] and 'y' in lonout.dims: lonout = grid.interp(lonout, 'y')
        if coords['lat'] and 'y' in latout.dims: latout = grid.interp(latout, 'y')
        if coords['z'] and 'y' in zout.dims: zout = grid.interp(zout, 'y')
    if dims['s'] == 's_w' and v.s_w.size>1:
        vout = grid.interp(vout, 'z')
        if coords['z'] and 's_w' in zout.dims: zout = grid.interp(zout, 'z')
    # assign coordinates
    if coords['lon']: vout = vout.assign_coords(coords={'lon_v':lonout})
    if coords['lat']: vout = vout.assign_coords(coords={'lat_v':latout})
    if coords['z']:vout = vout.assign_coords(coords={'z_v':zout})
    return vout

def x2w(ds,v, grid):
    """ Interpolate from any grid to w grid
    """
    dims = get_spatial_dims(v)
    coords = get_spatial_coords(v)
    vout = v.copy()
    if coords['lon']: lonout = v[coords['lon']].copy()
    if coords['lat']: latout = v[coords['lat']].copy()
    if coords['z']: zout = v[coords['z']].copy()
    if dims['x'] == 'x_u' and v.x_u.size>1:
        vout = grid.interp(vout, 'x')
        if coords['lon'] and 'x_u' in lonout.dims: lonout = grid.interp(lonout, 'x')
        if coords['lat'] and 'x_u' in latout.dims: latout = grid.interp(latout, 'x')
        if coords['z'] and 'x_u' in zout.dims: zout = grid.interp(zout, 'x')
    if dims['y'] == 'y_v' and v.y_v.size>1:
        vout = grid.interp(vout, 'y')
        if coords['lon'] and 'y_v' in lonout.dims: lonout = grid.interp(lonout, 'y')
        if coords['lat'] and 'y_v' in latout.dims: latout = grid.interp(latout, 'y')
        if coords['z'] and 'y_v' in zout.dims: zout = grid.interp(zout, 'y')
    if dims['s'] == 's' and v.s.size>1:
        vout = grid.interp(vout, 'z')
        if coords['z'] and 's' in zout.dims: zout = grid.interp(zout, 'z')
    # assign coordinates
    if coords['lon']: vout = vout.assign_coords(coords={'lon':lonout})
    if coords['lat']: vout = vout.assign_coords(coords={'lat':latout})
    if coords['z']:vout = vout.assign_coords(coords={'z_w':zout})
    return vout


def x2f(ds,v, grid):
    """ Interpolate from any grid to psi grid
    """
    dims = get_spatial_dims(v)
    coords = get_spatial_coords(v)
    vout = v.copy()
    if coords['lon']: lonout = v[coords['lon']].copy()
    if coords['lat']: latout = v[coords['lat']].copy()
    if coords['z']: zout = v[coords['z']].copy()
    if dims['x'] == 'x' and v.x.size>1:
        vout = grid.interp(vout, 'x')
        if coords['lon'] and 'x' in lonout.dims: lonout = grid.interp(lonout, 'x')
        if coords['lat'] and 'x' in latout.dims: latout = grid.interp(latout, 'x')
        if coords['z'] and 'x' in zout.dims: zout = grid.interp(zout, 'x')
    if dims['y'] == 'y' and v.y.size>1:
        vout = grid.interp(vout, 'y')
        if coords['lon'] and 'y' in lonout.dims: lonout = grid.interp(lonout, 'y')
        if coords['lat'] and 'y' in latout.dims: latout = grid.interp(latout, 'y')
        if coords['z'] and 'y' in zout.dims: zout = grid.interp(zout, 'y')
    if dims['s'] == 's_w' and v.s_w.size>1:
        vout = grid.interp(vout, 'z')
        if coords['z'] and 's_w' in zout.dims: zout = grid.interp(zout, 'z')
    # assign coordinates
    if coords['lon']: vout = vout.assign_coords(coords={'lon_f':lonout})
    if coords['lat']: vout = vout.assign_coords(coords={'lat_f':latout})
    if coords['z']:vout = vout.assign_coords(coords={'z_f':zout})
    return vout

def x2x(ds,v, grid, target):
    if target in ['rho', 'r']:
        return x2rho(ds,v, grid)
    elif target == 'u':
        return x2u(ds,v, grid)
    elif target == 'v':
        return x2v(ds,v, grid)
    elif target == 'w':
        return x2w(ds,v, grid)
    elif target in ['psi', 'p', 'f']:
        return x2f(ds,v, grid)




def get_z(model, ds=None, z_sfc=None, h=None, xgrid=None, vgrid='r',
          hgrid='r', vtransform=2):
    ''' Compute vertical coordinates
        Spatial dimensions are placed last, in the order: s/s_w, y, x

        Parameters
        ----------
        ds: xarray dataset
        z_sfc: xarray.DataArray, optional
            Sea level data, default to 0 if not provided
            If you use slices, make sure singleton dimensions are kept, i.e do:
                z_sfc.isel(x=[i])
            and not :
                z_sfc.isel(x=i)
        h: xarray.DataArray, optional
            Water depth, searche depth in grid if not provided
        vgrid: str, optional
            Vertical grid, 'r'/'rho' or 'w'. Default is 'r'
        hgrid: str, optional
            Any horizontal grid: 'r'/'rho', 'u', 'v', 'f'. Default is 'r'
        vtransform: int, str, optional
            croco vertical transform employed in the simulation.
            1="old": z = z0 + (1+z0/_h) * _z_sfc  with  z0 = hc*sc + (_h-hc)*cs
            2="new": z = z0 * (_z_sfc + _h) + _z_sfc  with  z0 = (hc * sc + _h * cs) / (hc + _h)
    '''

    xgrid = model.xgrid if xgrid is None else xgrid
    ds = model.ds if ds is None else ds

    h = ds.h if h is None else h
    z_sfc = 0*ds.h if z_sfc is None else z_sfc

    # switch horizontal grid if needed
    if hgrid in ['u','v','f']:
        h = x2x(ds, h, xgrid, hgrid)
        z_sfc = x2x(ds, z_sfc, xgrid, hgrid)

    # align datasets (z_sfc may contain a slice along one dimension for example)
    h, z_sfc  = xr.align(h, z_sfc, join='inner')

    if vgrid in ['r', 'rho']:
        vgrid = 'r'
        sc = ds['sc_r']
        cs = ds['Cs_r']
    else:
        sc = ds['sc_'+vgrid]
        cs = ds['Cs_'+vgrid]

    hc = ds['hc']

    if vtransform == 1:
        z0 = hc*sc + (h-hc)*cs
        z = z0 + (1+z0/h) * z_sfc
    elif vtransform == 2:
        z0 = (hc * sc + h * cs) / (hc + h)
        z = z0 * (z_sfc + h) + z_sfc

    # reorder spatial dimensions and place them last
    sdims = list(get_spatial_dims(z).values())
    sdims = tuple(filter(None,sdims)) # delete None values
    reordered_dims = tuple(d for d in z.dims if d not in sdims) + sdims
    z = z.transpose(*reordered_dims, transpose_coords=True).rename('z_'+hgrid)
    z.name = z.name.replace('z_r','z_'+vgrid)
    
    return z.fillna(0.) #.rename('z_'+hgrid).replace('z_r','z_'+vgrid)



def rotuv(model, ds=None, xgrid=None, u=None, v=None, angle=None):
    '''
    Rotate winds or u,v to lat,lon coord -> result on rho grid by default
    '''

    import timeit
    
    xgrid = model.xgrid if xgrid is None else xgrid
    ds = model.ds if ds is None else ds
        
    u = ds.u if u is None else u
    hgrid,vgrid = get_grid_point(u)
    if hgrid != 'r': u = x2rho(ds, u, xgrid)
    #u = ds_hor_chunk(u, wanted_chunk=100)
        
    v = ds.v if v is None else v
    hgrid,vgrid = get_grid_point(v)
    if hgrid != 'r': v = x2rho(ds, v, xgrid)
    #v = ds_hor_chunk(v, wanted_chunk=100)
        
    angle = ds.angle if angle is None else angle
    hgrid,vgrid = get_grid_point(angle)
    if hgrid != 'r': angle = x2rho(ds, angle, xgrid)
    
    cosang = np.cos(angle)
    sinang = np.sin(angle)

    # All the program statements
    urot = (u*cosang - v*sinang)
    
    #start = timeit.default_timer()
    vrot = (u*sinang + v*cosang)
    #stop = timeit.default_timer()
    #print("time vrot: "+str(stop - start))
    
    # assign coordinates to urot/vrot
    dims = get_spatial_dims(u)
    coords = get_spatial_coords(u)
    for k,c in coords.items(): 
        if c is not None: urot = urot.assign_coords(coords={c:u[c]})
    dims = get_spatial_dims(v)
    coords = get_spatial_coords(v)
    for k,c in coords.items(): 
        if c is not None: vrot = vrot.assign_coords(coords={c:v[c]})

    return [urot,vrot]



def hgrad(
    q,
    xgrid,
    which="both",
    z=None,
    hcoord=None,
    scoord=None,
    hboundary="extend",
    hfill_value=None,
    sboundary="extend",
    sfill_value=None,
    attrs=None,
):
    """Return gradients of property q accounting for s coordinates.

    Note that you need the 3D metrics for horizontal derivatives for ROMS, so ``include_3D_metrics=True`` in ``xroms.roms_dataset()``.

    Parameters
    ----------

    q: DataArray
        Property to take gradients of.
    xgrid: xgcm.grid
        Grid object associated with q.
    which: string, optional
        Which components of gradient to return.
        * 'both': return both components of hgrad.
        * 'xi': return only xi-direction.
        * 'eta': return only eta-direction.
    z: DataArray, ndarray, optional
        Depth [m]. If None, use z coordinate attached to q.
    hcoord: string, optional
        Name of horizontal grid to interpolate output to.
        Options are 'rho', 'psi', 'u', 'v'.
    scoord: string, optional
        Name of vertical grid to interpolate output to.
        Options are 's_rho', 's_w', 'rho', 'w'.
    hboundary: string, optional
        Passed to `grid` method calls; horizontal boundary selection
        for calculating horizontal derivatives of q. This same value
        will be used for all horizontal grid changes too.
        From xgcm documentation:
        A flag indicating how to handle boundaries:
        * None:  Do not apply any boundary conditions. Raise an error if
          boundary conditions are required for the operation.
        * 'fill':  Set values outside the array boundary to fill_value
          (i.e. a Neumann boundary condition.)
        * 'extend': Set values outside the array to the nearest array
          value. (i.e. a limited form of Dirichlet boundary condition.
    hfill_value: float, optional
        Passed to `grid` method calls; horizontal boundary selection
        fill value.
        From xgcm documentation:
        The value to use in the boundary condition with `boundary='fill'`.
    sboundary: string, optional
        Passed to `grid` method calls; vertical boundary selection
        for calculating horizontal derivatives of q. This same value will
        be used for all vertical grid changes too.
        From xgcm documentation:
        A flag indicating how to handle boundaries:
        * None:  Do not apply any boundary conditions. Raise an error if
          boundary conditions are required for the operation.
        * 'fill':  Set values outside the array boundary to fill_value
          (i.e. a Neumann boundary condition.)
        * 'extend': Set values outside the array to the nearest array
          value. (i.e. a limited form of Dirichlet boundary condition.
    sfill_value: float, optional
        Passed to `grid` method calls; vertical boundary selection
        fill value.
        From xgcm documentation:
        The value to use in the boundary condition with `boundary='fill'`.
    attrs: dict, optional
        Dictionary of attributes to add to resultant arrays. Requires that
        q is DataArray.

    Returns
    -------
    DataArray(s) of dqdxi and/or dqdeta, the gradients of q
    in the xi- and eta-directions with attributes altered to reflect calculation.

    Notes
    -----
    dqdxi = dqdx*dzdz - dqdz*dzdx

    dqdeta = dqdy*dzdz - dqdz*dzdy

    Derivatives are taken in the ROMS curvilinear grid native xi- and eta- directions.

    These derivatives properly account for the fact that ROMS vertical coordinates are
    s coordinates and therefore can vary in time and space.

    The xi derivative will alter the number of points in the xi and s dimensions.
    The eta derivative will alter the number of points in the eta and s dimensions.

    Examples
    --------
    >>> dtempdxi, dtempdeta = xroms.hgrad(ds.temp, xgrid)
    """

    assert isinstance(q, xr.DataArray), "var must be DataArray"

    if not [dim for dim in q.dims if dim.startswith('s')]:
        is3D = False
    else:
        is3D = True
        
    if is3D and z is None:
        try:
            coords = list(q.coords)
            z_coord_name = coords[[coord[:2] == "z_" for coord in coords].index(True)]
            z = q[z_coord_name]
            is3D = True
        except:
            # if we get here that means that q doesn't have z coords (like zeta)
            print("!!! Missing z coordinate, only horizontal gradient")
            is3D = False

    if which in ["both", "x"]:

        if is3D:
            dqdx = xgrid.interp(
                xgrid.derivative(q, "x", boundary=hboundary, fill_value=hfill_value),
                "z",
                boundary=sboundary,
                fill_value=sfill_value,
            )
            dqdz = xgrid.interp(
                xgrid.derivative(q, "z", boundary=sboundary, fill_value=sfill_value),
                "x",
                boundary=hboundary,
                fill_value=hfill_value,
            )
            dzdx = xgrid.interp(
                xgrid.derivative(z, "x", boundary=hboundary, fill_value=hfill_value),
                "z",
                boundary=sboundary,
                fill_value=sfill_value,
            )
            dzdz = xgrid.interp(
                xgrid.derivative(z, "z", boundary=sboundary, fill_value=sfill_value),
                "x",
                boundary=hboundary,
                fill_value=hfill_value,
            )

            dqdx = dqdx * dzdz - dqdz * dzdx

        else:  # 2D variables
            dqdx = xgrid.derivative(q, "x", boundary=hboundary, fill_value=hfill_value)

        if attrs is None and isinstance(q, xr.DataArray):
            attrs = q.attrs.copy()
            attrs["name"] = "d" + q.name + "dx"
            attrs["units"] = "1/m * " + attrs.setdefault("units", "units")
            attrs["long_name"] = "horizontal xi derivative of " + attrs.setdefault(
                "long_name", "var"
            )

    if which in ["both", "y"]:

        if is3D:
            dqdy = xgrid.interp(
                xgrid.derivative(q, "y", boundary=hboundary, fill_value=hfill_value),
                "z",
                boundary=sboundary,
                fill_value=sfill_value,
            )
            dqdz = xgrid.interp(
                xgrid.derivative(q, "z", boundary=sboundary, fill_value=sfill_value),
                "y",
                boundary=hboundary,
                fill_value=hfill_value,
            )
            dzdy = xgrid.interp(
                xgrid.derivative(z, "y", boundary=hboundary, fill_value=hfill_value),
                "z",
                boundary=sboundary,
                fill_value=sfill_value,
            )
            dzdz = xgrid.interp(
                xgrid.derivative(z, "z", boundary=sboundary, fill_value=sfill_value),
                "y",
                boundary=hboundary,
                fill_value=hfill_value,
            )

            dqdy = dqdy * dzdz - dqdz * dzdy

        else:  # 2D variables
            dqdy = xgrid.derivative(
                q, "y", boundary=hboundary, fill_value=hfill_value
            )

        if attrs is None and isinstance(q, xr.DataArray):
            attrs = q.attrs.copy()
            attrs["name"] = "d" + q.name + "dy"
            attrs["units"] = "1/m * " + attrs.setdefault("units", "units")
            attrs["long_name"] = "horizontal eta derivative of " + attrs.setdefault(
                "long_name", "var"
            )

    if which == "both":
        return dqdx, dqdy
    elif which == "x":
        return dqdx
    elif which == "y":
        return dqdy
    else:
        print("nothing being returned from hgrad")
        

def get_grid_point(var):
    dims = var.dims
    # horizontal point
    hpoint='r'
    if "x_u" in dims:
        if "y" in dims:
            hpoint='u'
        else:
            hpoint='f'
    elif "y_v" in dims:
        hpoint='v'
    if 's' in dims:
        vpoint='r'
    else:
        vpoint='w'
    return hpoint,vpoint
        
    
def slices(model, var, z, ds=None, xgrid=None, longitude=None, latitude=None, depth=None):
    """
    #
    #
    # This function interpolate a 3D variable on slices at constant depths/longitude/latitude
    # This function use xcgm transform method and needs xgcm.Grid to be defined over the 3 axes.
    # !!! For now, it works only with curvilinear coordinates !!!
    #
    # On Input:
    #
    #    ds      dataset to find the grid
    #    var     (dataArray) Variable to process (3D matrix).
    #    z       (dataArray) Depths at the same point than var (3D matrix).
    #    longitude   (scalar,list or ndarray) longitude of the slice (scalar meters, negative).
    #    latitude    (scalar,list or ndarray) latitude of the slice (scalar meters, negative).
    #    depth       (scalar,list or ndarray) depth of the slice (scalar meters, negative).
    #
    # On Output:
    #
    #    vnew    (dataArray) Horizontal slice
    #
    #
    """
    from matplotlib.cbook import flatten
  
    xgrid = model.xgrid if xgrid is None else xgrid
    ds = model.ds if ds is None else ds
    
    if longitude is None and latitude is None and depth is None:
        "Longitude or latitude or depth must be defined"
        return None

    # check typ of longitude/latitude/depth
    # longitude = longitude.tolist() if isinstance(longitude,np.ndarray) else longitude
    # longitude = [longitude] if (isinstance(longitude,int) or isinstance(longitude,float)) else longitude
    longitude = np.asarray(longitude) if isinstance(longitude,list) else longitude
    longitude = np.asarray([longitude]) if (isinstance(longitude,int) or isinstance(longitude,float)) else longitude

    # latitude = latitude.tolist() if isinstance(latitude,np.ndarray) else latitude
    # latitude = [latitude] if (isinstance(latitude,int) or isinstance(latitude,float)) else latitude
    latitude = np.asarray(latitude) if isinstance(latitude,list) else latitude
    latitude = np.asarray([latitude]) if (isinstance(latitude,int) or isinstance(latitude,float)) else latitude

    # depth = depth.tolist() if isinstance(depth,np.ndarray) else depth
    # depth = [depth] if (isinstance(depth,int) or isinstance(depth,float)) else depth
    depth = np.asarray(depth) if isinstance(depth,list) else depth
    depth = np.asarray([depth]) if (isinstance(depth,int) or isinstance(depth,float)) else depth

     # Find dimensions and coordinates of the variable
    dims = get_spatial_dims(var)
    coords = get_spatial_coords(var)
    if dims['s'] is not None and coords['z'] is None: 
        var = var.assign_coords(coords={'z':z})
        coords = get_spatial_coords(var)
    # hgrid,vgrid = get_grid_point(var)

    if longitude is not None:
        axe = 'x'
        coord_ref = coords['lon']
        coord_x = coords['lat']
        coord_y = coords['z']
        slices_values = longitude
    elif latitude is not None:
        axe = 'y'
        coord_ref = coords['lat']
        coord_x = coords['lon']
        coord_y = coords['z']
        slices_values = latitude
    else:
        axe = 'z'
        coord_ref = coords['z']
        coord_x = coords['lon']
        coord_y = coords['lat']
        slices_values = depth

    # Recursively loop over time if needed
    if len(var.squeeze().dims) == 4:
        vnew = [slices(model, var.isel(t=t), z.isel(t=t), ds=ds, xgrid=xgrid,
                      longitude=longitude, latitude=latitude, depth=depth)
                      for t in range(len(var.t))]
        vnew = xr.concat(vnew, dim='t')
    else:
        vnew = xgrid.transform(var, axe, slices_values,
                               target_data=var[coord_ref]).squeeze()
    # Do the linear interpolation
    if not depth:
        x = xgrid.transform(var[coord_x], axe, slices_values,
                                   target_data=var[coord_ref]).squeeze() #\
                     #.expand_dims({dims['s']: len(var[dims['s']])})            
        vnew = vnew.assign_coords(coords={coord_x:x})

        #y = xgrid.transform(var[coord_y], axe, slices_values,
        if dims['s'] is not None:
            y = xgrid.transform(z, axe, slices_values,
                                       target_data=var[coord_ref]).squeeze()
            # Add the coordinates to dataArray
            vnew = vnew.assign_coords(coords={coord_y:y})
    else:
        # Add the coordinates to dataArray
        vnew = vnew.assign_coords(coords={coord_x:var[coord_x]})
        vnew = vnew.assign_coords(coords={coord_y:var[coord_y]})

#     return vnew.squeeze().unify_chunks().fillna(0.)  #unify_chunks() 
    return vnew.squeeze().fillna(0.)  #unify_chunks()   

def cross_section(grid, da, lon1, lat1, lon2, lat2, dlon=None):
    
    # check input parameters
    if not isinstance(grid,Grid): print('grid must be a xgcm grid'); return None
    if not isinstance(da,xr.DataArray): print('da must be a xarray DataArray'); return None
    dims = get_spatial_dims(da)
    coords = get_spatial_coords(da)
    if coords['lon'] is None or coords['lat'] is None:
        print('da must have longitude AND latitude coordinates')
        return None
    if not isinstance(lon1,(int,float)): print('lon1 must be a float'); return None
    if not isinstance(lat1,(int,float)): print('lat1 must be a float'); return None
    if not isinstance(lon2,(int,float)): print('lon2 must be a float'); return None
    if not isinstance(lat2,(int,float)): print('lat2 must be a float'); return None
    if dlon is not None and not isinstance(dlon,(int,float)): print('dlon must be a number'); return None
           
    # compute the linear function from the two points
    a = (lat2 - lat1)/(lon2 - lon1)
    b = lat1 - a*lon1
    
    # get the longitude interval of the new grid , compute the new longitude grid
    if dlon is None:
        dlon = ((da[coords['lon']].max().values - da[coords['lon']].min().values) /
                da[coords['lon']].sizes[dims['x']])
    longrd = np.arange(lon1,lon2,dlon)
    
    # compute the latitude coordinates of the new grid with the linear function
    latgrd = a * longrd + b

    # interpolate on the regular longitude grid
    newda = auto_chunk(da, keep_complete_dim='x', wanted_chunk=200)
    newda = grid.transform(newda,'x', longrd, target_data=da[coords['lon']])
    newda = newda.rename({coords['lon']:dims['x']})
    newlat = auto_chunk(da[coords['lat']], keep_complete_dim='x', wanted_chunk=200)
    newlat = grid.transform(newlat, 'x', longrd, target_data=da[coords['lon']])              
    newlat = newlat.rename({coords['lon']:dims['x']})
    if coords['z'] is not None:
        newz = auto_chunk(da[coords['z']], keep_complete_dim='x', wanted_chunk=200)
        newz = grid.transform(newz,'x', longrd, target_data=da[coords['lon']]).fillna(0.)
        newz = newz.rename({coords['lon']:dims['x']})
        
    # interpolate on a new latitude grid
    newda = auto_chunk(newda, keep_complete_dim='y', wanted_chunk=200)
    newda = grid.transform(newda,'y',latgrd,target_data=newlat)
    newda = newda.rename({coords['lat']:dims['y']})
    if coords['z'] is not None:
        newz = auto_chunk(newz, keep_complete_dim='y', wanted_chunk=200)
        newz = grid.transform(newz,'y',latgrd,target_data=newlat).fillna(0.)
        newz = newz.rename({coords['lat']:dims['y']}).fillna(0.)
    
    # extract the cross section
    crossda = []; crossz=[]
    for lon,lat in zip(longrd,latgrd):
        crossda.append(newda.loc[{dims['x']:lon}].loc[{dims['y']:lat}])
        if coords['z'] is not None:
            crossz.append(newz.loc[{dims['x']:lon}].loc[{dims['y']:lat}])

    cross = xr.concat(crossda,dim=dims['x'])
    if coords['z'] is not None:
        crossz = xr.concat(crossz,dim=dims['x'])
    
    # assign the coordinates lon/lat/z to the section
    cross = cross.assign_coords({coords['lon']:cross[dims['x']],
                                 coords['lat']:cross[dims['y']]})
    if coords['z'] is not None:
        cross = cross.assign_coords({coords['z']:crossz})    
        
    return reorder_dims(cross.fillna(0.))

def interp_regular(da,grid,axis,tgrid,rgrid=None):
    """
    interpolate on a regular grid
    imputs:
        - da (DataArray) : variable to interpolate
        - grid (xgcm grid): xgcm grid 
        - axis (str): axis of the xgcm grid for the interpolation ('x', 'y' or 'z')
        - tgrid (numpy vector): target relular grid space
        - rgrid (numpy array or DataArray): reference grid of da
    return:
        - (DataArray): regurlarly interpolated variable
    """
    
    # check axis
    if axis not in ['x','y','z']: 
        print('axis must be x, y or z')
        return None
    
    # corresponding keys between spatial coords/dims and axes of the xgcm grid
    refc = {'x':'lon', 'y':'lat', 'z':'z'}
    refd = {'x':'x', 'y':'y', 'z':'s'}
    
    # find spatial coordinates/dims of da
    coords = get_spatial_coords(da)
    coord = coords[refc[axis]]
    dims = get_spatial_dims(da)
    dim = dims[refd[axis]]
    
    # initialize the reference coordinate
    if rgrid is None and coord is not None:
        rgrid = da[coord]
    else:
        print('the reference grid is missing along the axis of interpolation')
        return None

    # interpolate da on the regular grid
    newvar = grid.transform(da,axis,tgrid,target_data=rgrid).rename({coord:dim})
    
    return reorder_dims(newvar )


def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    # Radius of earth in kilometers is 6371
    km = 6371* c
    return km

    
# ----------------------------- grid rechunk -----------------------------

def auto_chunk(ds, keep_complete_dim=None, wanted_chunk=150):
    """
    Rechunk Dataset or DataArray such as each partition size is about 150Mb
    Input:
        - ds : (Dataset or DataArray) object to rechunk
        - keep_complete_dim : (character) Horizontal axe to keep with no chunk ('x','y','s')
        - wanted_chunk : (integer) size of each partition in Mb
    Output:
        - object rechunked
    """

    #check input parameters
    if not isinstance(ds, (xr.Dataset,xr.DataArray)):
        print('argument must be a xarray.DataArray or xarray.Dataset')
        return
    if keep_complete_dim is not None:
        keep_complete_dim = list(keep_complete_dim) if not isinstance(keep_complete_dim,list) else keep_complete_dim
        if not all(item in ['s','y','x'] for item in keep_complete_dim):
            print('keep_complete_dim must equal x or y or s')
            return

    # get horizontal dimensions names of the Dataset/DataArray
    dname = get_spatial_dims(ds)
    # remove None values
    dname = {k: v for k, v in dname.items() if v is not None}
    chunks_name = dname.copy()
    
    # get max dimensions sizes of the Dataset/DataArray
    chunks_size={}
    for k,v in chunks_name.items():
        if isinstance(v,list):       # for a dataset
            chunks_size[k] = max([ds.sizes[d] for d in v])
        else:
            chunks_size[k] = ds.sizes[v]

    # always chunk in time
    if 't' in ds.dims: chunks_size['t'] = 1
        
    if keep_complete_dim:
        # remove keep_complete_dim from the dimensions of the Dataset/DatAarray
        for d in keep_complete_dim: del chunks_name[d]
        
    # reduce chunks size  beginning by 's' then 'y' then 'x' if necessary
    for k in chunks_name.keys():
        for d in range(chunks_size[k],0,-1):
            # chunk_size = (chunks_size['x']*chunks_size['y']*chunks_size['s']*4 / 1.e6)
            chunk_size = 4 / 1.e6
            for chunk in chunks_size.values():
                chunk_size = chunk_size*chunk
            if chunk_size > wanted_chunk:
                chunks_size[k] = d
            else:
                break
        if chunk_size > wanted_chunk:
            break

    if isinstance(ds,xr.Dataset):
        # set chunk for all the dimensions of the dataset (ie : x and x_u)
        for c in list(itertools.product(dname.keys(), ds.dims.keys())):
            if c[1].startswith(c[0]):
                chunks_size[c[1]] = chunks_size[c[0]]
    else:
        # rename the dimension name by the right values (ie: x_u instead of x)
        for key in dname.keys():
            chunks_size[dname[key]] = chunks_size.pop(key)

        
    return ds.chunk(chunks_size)

