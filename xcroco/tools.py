import os

import dask
from dask.distributed import progress
import xarray as xr


import numpy as np 

# try to force flushing of memory
import gc, ctypes

#---------- cluster tools

def wait_cluster_ready(cluster, nworkers):
    """
    Wait for the client to be ready (all workers started)
    """
    import time, sys
    Nw = 0.
    while(len(cluster.scheduler.workers)<nworkers):
        p = len(cluster.scheduler.workers)/nworkers
        time.sleep(2)
        # print('{} % of the workers started'.format(int(p*100)))
        if p>=Nw:
            Nw = int(p*10)/10
            Nw+=0.1
    print('100 % of the workers started: {} workers'.format(nworkers))
    
def trim_memory() -> int:
    libc = ctypes.CDLL("libc.so.6")
    return libc.malloc_trim(0)

def dask_compute_batch(computations, client, batch_size=None):
    """ breaks down a list of computations into batches
    """
    from dask.diagnostics import ProgressBar
    ProgressBar().register()

    # compute batch size according to number of workers
    if batch_size is None:
        # batch_size = len(client.scheduler_info()["workers"])
        batch_size = sum(list(client.nthreads().values()))
    # find batch indices
    total_range = range(len(computations))
    splits = max(1, np.ceil(len(total_range)/batch_size))
    batches = np.array_split(total_range, splits)
    # launch computations
    outputs = []
    for b in batches:
        out = dask.compute(*computations[slice(b[0], b[-1]+1)])
        progress(out)

        # try to manually clean up memory
        # https://coiled.io/blog/tackling-unmanaged-memory-with-dask/
        # client.run(gc.collect)
        client.run(trim_memory)  # should not be done systematically
        outputs.append(out)
    return outputs
    
#--------- Storage tools

def is_zarr_archive(zarr_dir, key):
    """ Utils, test existence of a zarr archive
    """
    zarr_archive = path.join(zarr_dir, key+'.zarr')
    return path.isdir(zarr_archive)

def store_zarr(ds, 
               zarr_archive,
               chunks={}, 
               auto_rechunk=True,
               mode='w',
               append_dim='time_counter',
               **kwargs):
    """ writes data in zarr archives

    Parameters
    ----------
    chunks: dict, optional
        Dictionary with output keys as keys and dimension chunk sizes 
        as values
    auto_rechunk: boolean, optional
        Activate automatic rechunking which will ensure chunks are larger than
        _chunk_size_threshold (see postp.py). Default is True.
    **kwargs:
        Passed to xr.to_zarr method
    """

#       change singleton variables and coords to attrs
    ds = move_singletons_as_attrs(ds)

    # fix encoding inplace for nonchunked data
    fix_nochunk_encoding(ds)
    #
    if mode=='a':
        ds.to_zarr(zarr_archive, mode=mode, append_dim=append_dim, **kwargs)
    else:
        ds.to_zarr(zarr_archive, **kwargs)
    # print('- {} stored'.format(ds.attrs['name']))
    print_zarr_archive_info(zarr_archive)
        
def store_netcdf(ds, filename, chunks={},  auto_rechunk=False, mode='w',
               append_dim='time_counter', **kwargs):
    """ 
    writes data in a unique netcdf file

    Parameters
    ----------
    chunks: dict, optional
        Dictionary with output keys as keys and dimension chunk sizes 
        as values
    auto_rechunk: boolean, optional
        Activate automatic rechunking which will ensure chunks are larger than
        _chunk_size_threshold (see postp.py). Default is True.
    **kwargs:
        Passed to ds.to_netcdf method
    """
    #       change singleton variables and coords to attrs
    ds = move_singletons_as_attrs(ds)

    # fix encoding inplace for nonchunked data
    fix_nochunk_encoding(ds)
    #
    ds.to_netcdf(filename, **kwargs)
    print('- {} stored'.format(filename))
    print_netcdf_file_info(filename)
        
def move_singletons_as_attrs(ds, ignore=[]):
    """ change singleton variables and coords to attrs
    This seems to be required for zarr archiving
    """
    for c,co in ds.coords.items():
        if len(co.dims)==0:
            ds = ds.drop_vars(c).assign_attrs({c: co.values})
    for k,v in ds.data_vars.items():
        if len(v.dims)==0:
            ds = ds.drop_vars(k).assign_attrs({k: v.values})
    return ds

def fix_nochunk_encoding(da):
    ''' Fix in place the encoding for nonchunked arrays such that zarr 
    writing does not automatically chunked arrays.
    
    Parameters
    ----------
    da: xr.DataArray, xr.Dataset
        variable or dataset to fix
    '''
    if isinstance(da, xr.Dataset):
        for v in da:
            fix_nochunk_encoding(da[v])
    if isinstance(da, xr.DataArray):
        if not da.chunks:
            da.encoding['chunks'] = -1
            
def print_zarr_archive_info(zarr_archive):
    """ Print basic information about a zarr archive
    """
    print('  location: {} '.format(zarr_archive))
    # get archive size
    arch_size = get_dir_size(zarr_archive)
    print('  size:     {:.1f} GB'.format(arch_size/1e9))
    #
    ds = xr.open_zarr(zarr_archive)
    # get largest item typical chunks
    n_dim_max = 0
    for v in ds:
        if ds[v].chunks and ds[v].ndim>n_dim_max:
            _size = list(ds[v].sizes.values())
            _chunks = [np.max(d) for d in ds[v].chunks]
            n_dim_max = ds[v].ndim
    if n_dim_max>0:
        print('  typical chunks: ('
              +','.join('{}'.format(c) for c in _chunks)
              +') for size ('
              +','.join('{}'.format(c) for c in _size)
              +')'
             )
    else:
        print('  data is not chunked')

        
def print_netcdf_file_info(netcdf_file):
    """ Print basic information about a netcdf file
    """
    print('  location: {} '.format(netcdf_file))
    # get file size
    file_size = get_dir_size(netcdf_file)
    print('  size:     {:.1f} GB'.format(file_size/1e9))
    #
    ds = xr.open_dataset(netcdf_file)
    n_dim_max = 0
    for v in ds:
        if ds[v].chunks and ds[v].ndim>n_dim_max:
            _size = list(ds[v].sizes.values())
            _chunks = [np.max(d) for d in ds[v].chunks]
            n_dim_max = ds[v].ndim
    if n_dim_max>0:
        print('  typical chunks: ('
              +','.join('{}'.format(c) for c in _chunks)
              +') for size ('
              +','.join('{}'.format(c) for c in _size)
              +')'
             )
    else:
        print('  data is not chunked')
    # # get largest item typical chunks
    # chks = dsg.chunks
    # print('  chunks:')
    # for k,v in chks.items():
    #     print('     ',k,': ',v[0])
    # print('  size of a 2D chunk: {:.1f} MB'.format(chks['time_counter'][0]*
    #                                                chks['x_rho'][0]*
    #                                                chks['y_rho'][0]* 4 / 1.e6))
    # if 's_rho' in ds.dims.keys():
    #     print('  size of a 3D chunk: {:.1f} MB '.format(chks['time_counter'][0]*
    #                                                     chks['x_rho'][0]*
    #                                                     chks['y_rho'][0]*
    #                                                     chks['s_rho'][0]* 4 / 1.e6))
        

def get_dir_size(dir_path):
    ''' Returns the size of a directory in bytes
    '''
    process = os.popen('du -s '+dir_path)
    size = int(process.read().split()[0]) # du returns kb
    process.close()
    return size*1e3