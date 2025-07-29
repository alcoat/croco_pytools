"""xcroco in/out module"""

import os
from functools import partial

import intake
import numpy as np
import xarray as xr

import gridop


def open_files(
    model,
    gridname,
    filenames,
    grid_metrics=2,
    drop_variables=None,
    chunks=None,
    remove_ghost_pts=True,
    xperiodic=False,
    yperiodic=False,
    verbose=False,
    **prepro_kwargs,
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
        verbose: verbose mode
    outputs:
        ds: an xarray dataset
        grid: the associated xgcm grid

    """

    def _preprocess(ds, **prepro_kwargs):
        for key, value in prepro_kwargs.items():
            try:
                ds = ds.sel({key: value}, method="nearest")
            except:
                ds = ds.isel({key: value})
        return ds

    if chunks is None:
        chunks = {"t": 1}

    partial_func = partial(_preprocess, **prepro_kwargs)

    # convert filenames to list of strings if string
    filenames = filenames.tolist() if isinstance(filenames, str) else filenames
    # find time dimension name
    concat_dim = [k for k, v in model.rename_vars.items() if v == "t"][0]
    open_kwargs = {
        "concat_dim": concat_dim,
        "combine": "nested",
        "coords": "minimal",
        "parallel": True,
        "compat": "override",
    }
    try:
        # zarr archive
        ds = xr.open_zarr(filenames[0], drop_variables=drop_variables)
    except:
        try:
            # list of explicit filenames
            ds = xr.open_mfdataset(
                filenames,
                drop_variables=drop_variables,
                preprocess=partial_func,
                **open_kwargs,
            )
        except Exception as exc:
            print(exc)
            try:
                # list of files with wildcards
                ds = xr.open_mfdataset(
                    filenames[0],
                    drop_variables=drop_variables,
                    preprocess=partial_func,
                    **open_kwargs,
                )
            except Exception as exc:
                print(exc)
                print("open_files: unknown format: only Netcdf or Zarr")
                print("or filenames do not exist")
                return None

    # change the names in the dataset according to the model instance
    model.ds = gridop.adjust_grid(model, ds)

    # add the grid and the xgcm grid to the dataset
    ds, grid = gridop.add_grid(
        model,
        gridname,
        grid_metrics=grid_metrics,
        remove_ghost_pts=remove_ghost_pts,
        xperiodic=xperiodic,
        yperiodic=yperiodic,
    )
    ds = force_cf_convention(ds)
    model.ds = ds.chunk(chunks=chunks).squeeze()
    return model.ds, grid


def open_catalog(
    model,
    gridname,
    catalog,
    source=None,
    grid_metrics=1,
    chunks=None,
    remove_ghost_pts=True,
    xperiodic=False,
    yperiodic=False,
    verbose=False,
    **prepro_kwargs,
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
        verbose: verbose mode
    outputs:
        ds: an xarray dataset
        grid: the associated xgcm grid

    """

    def _preprocess(ds, **prepro_kwargs):
        for key, value in prepro_kwargs.items():
            try:
                ds = ds.sel({key: value}, method="nearest")
            except:
                ds = ds.isel({key: value})
        return ds

    if chunks is None:
        chunks = {}

    try:
        # open intake catalog
        cat = intake.open_catalog(catalog)
        if verbose:
            print("Available sources are: ", list(cat))
    except:
        print("open_catalog: yaml catalog not found")
        return None

    # find the first source of the catalog if source is None
    source = list(cat)[0] if source is None else source
    if verbose:
        print("Source :", source)
        print("   ", cat[source])
    # open the source as a dataset
    try:
        ds = cat[source].to_dask().pipe(lambda x: _preprocess(x, **prepro_kwargs))
    except:
        print("May be the source is not found in the yaml catalog")
        return None

    # change the names in the dataset according to the model instance
    model.ds = gridop.adjust_grid(model, ds)

    # add the grid and the xgcm grid to the dataset
    ds, grid = gridop.add_grid(
        model,
        gridname,
        grid_metrics=grid_metrics,
        remove_ghost_pts=remove_ghost_pts,
        xperiodic=xperiodic,
        yperiodic=yperiodic,
    )
    ds = force_cf_convention(ds)
    model.ds = ds.chunk(chunks=chunks)
    return model.ds.squeeze(), grid


def force_cf_convention(ds):
    """Force CF convention attributes of dimensions and coordinates for using cf_xarray

    Args:
    ----
        ds (dataset): input xarray dataset

    Returns:
    -------
        ds (dataset): xarray dataset with CF convention for dimensions and coordinates

    Examples:
    --------
    >>> ds = force_cf_convention(ds)
    """
    # make sure dimensions have axis attribute
    tdims = [dim for dim in ds.dims if dim.startswith("x")]
    for dim in tdims:
        ds[dim] = (dim, np.arange(ds.sizes[dim]), {"axis": "X"})
    tdims = [dim for dim in ds.dims if dim.startswith("y")]
    for dim in tdims:
        ds[dim] = (dim, np.arange(ds.sizes[dim]), {"axis": "Y"})
    tdims = [dim for dim in ds.dims if dim.startswith("s")]
    for dim in tdims:
        if dim in ds.coords:
            ds[dim].attrs["axis"] = "Z"
        else:
            ds[dim] = (dim, np.arange(ds.sizes[dim]), {"axis": "Z"})

    tdims = [dim for dim in ds.dims if dim.startswith("t")]
    for dim in tdims:
        if dim in ds.coords:
            ds[dim].attrs["axis"] = "T"
        else:
            ds[dim] = (dim, np.arange(ds.sizes[dim]), {"axis": "T"})

    # make sure lon/lat/depth/time have standard names
    tcoords = [coord for coord in ds.coords if coord.startswith("lon")]
    for coord in tcoords:
        ds[coord].attrs["standard_name"] = "longitude"
    tcoords = [coord for coord in ds.coords if coord.startswith("lat")]
    for coord in tcoords:
        ds[coord].attrs["standard_name"] = "latitude"

    tcoords = [coord for coord in ds.coords if coord.startswith("z")]
    for coord in tcoords:
        ds[coord].attrs["standard_name"] = "vertical"
    tcoords = [coord for coord in ds.coords if coord.startswith("t")]
    for coord in tcoords:
        ds[coord].attrs["standard_name"] = "time"

    return ds


def find_var(model, varname, ds, gd):
    """Find a variable in the gridname or history files variables or attributes

    Args:
        model (string): model class
        varname (string): variable name to find
        ds (dataset):  dataset of history file
        gd (dataset): dataset of the grid
        value (optional): value of the variable if not found. Defaults to None.
    Returns:
        (DataArray): the DataArray corresponding to varname
    """

    def good_type(var):
        return isinstance(
            var,
            (xr.DataArray, np.ndarray, np.float32, np.float64, np.int32, np.float64),
        )

    def to_dataarray(model, varname, var):
        if isinstance(var, np.ndarray):
            var = xr.DataArray(data=var, dims=model.dims_var[varname])

        return var

    if varname in gd and good_type(gd[varname]):
        return to_dataarray(model, varname, gd[varname])
    if varname in gd.attrs and good_type(gd.attrs[varname]):
        return to_dataarray(model, varname, gd.attrs[varname])
    if varname in ds and good_type(ds[varname]):
        return to_dataarray(model, varname, ds[varname])
    if varname in ds.attrs and good_type(ds.attrs[varname]):
        return to_dataarray(model, varname, ds.attrs[varname])

    return None


# --------- Storage tools


def is_zarr_archive(zarr_dir, key):
    """Utils, test existence of a zarr archive"""
    zarr_archive = os.path.join(zarr_dir, key + ".zarr")
    return os.path.isdir(zarr_archive)


def store_zarr(
    ds,
    zarr_archive,
    chunks=None,
    auto_rechunk=True,
    mode="w",
    append_dim="time_counter",
    **kwargs,
):
    """writes data in zarr archives

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

    if chunks is None:
        chunks = {}
    #       change singleton variables and coords to attrs
    ds = move_singletons_as_attrs(ds)

    # fix encoding inplace for nonchunked data
    fix_nochunk_encoding(ds)
    #
    if mode == "a":
        ds.to_zarr(zarr_archive, mode=mode, append_dim=append_dim, **kwargs)
    else:
        ds.to_zarr(zarr_archive, **kwargs)
    # print('- {} stored'.format(ds.attrs['name']))
    print_zarr_archive_info(zarr_archive)


def store_netcdf(
    ds,
    filename,
    chunks=None,
    auto_rechunk=False,
    mode="w",
    append_dim="time_counter",
    **kwargs,
):
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
    if chunks is None:
        chunks = {}
    #       change singleton variables and coords to attrs
    ds = move_singletons_as_attrs(ds)

    # fix encoding inplace for nonchunked data
    fix_nochunk_encoding(ds)
    #
    ds.to_netcdf(filename, **kwargs)
    print(f"- {filename} stored")
    print_netcdf_file_info(filename)


def move_singletons_as_attrs(ds, ignore=None):
    """change singleton variables and coords to attrs
    This seems to be required for zarr archiving
    """
    for c, co in ds.coords.items():
        if len(co.dims) == 0:
            if isinstance(co.values, np.datetime64):
                ds = ds.drop_vars(c).assign_attrs({c: str(co.values)})
            else:
                ds = ds.drop_vars(c).assign_attrs({c: co.values})
    for k, v in ds.data_vars.items():
        if len(v.dims) == 0:
            if isinstance(v.values, np.datetime64):
                ds = ds.drop_vars(k).assign_attrs({k: str(v.values)})
            else:
                ds = ds.drop_vars(k).assign_attrs({k: v.values})
    return ds


def fix_nochunk_encoding(da):
    """Fix in place the encoding for nonchunked arrays such that zarr
    writing does not automatically chunked arrays.

    Parameters
    ----------
    da: xr.DataArray, xr.Dataset
        variable or dataset to fix
    """
    if isinstance(da, xr.Dataset):
        for v in da:
            fix_nochunk_encoding(da[v])
    if isinstance(da, xr.DataArray):
        if not da.chunks:
            da.encoding["chunks"] = -1


def print_zarr_archive_info(zarr_archive):
    """Print basic information about a zarr archive"""
    print(f"  location: {zarr_archive} ")
    # get archive size
    arch_size = get_dir_size(zarr_archive)
    print(f"  size:     {arch_size / 1e9:.1f} GB")
    #
    ds = xr.open_zarr(zarr_archive)
    # get largest item typical chunks
    n_dim_max = 0
    for v in ds:
        if ds[v].chunks and ds[v].ndim > n_dim_max:
            _size = list(ds[v].sizes.values())
            _chunks = [np.max(d) for d in ds[v].chunks]
            n_dim_max = ds[v].ndim
    if n_dim_max > 0:
        print(
            "  typical chunks: ("
            + ",".join(f"{c}" for c in _chunks)
            + ") for size ("
            + ",".join(f"{c}" for c in _size)
            + ")"
        )
    else:
        print("  data is not chunked")


def print_netcdf_file_info(netcdf_file):
    """Print basic information about a netcdf file"""
    print(f"  location: {netcdf_file} ")
    # get file size
    file_size = get_dir_size(netcdf_file)
    print(f"  size:     {file_size / 1e9:.1f} GB")
    #
    ds = xr.open_dataset(netcdf_file)
    n_dim_max = 0
    for v in ds:
        if ds[v].chunks and ds[v].ndim > n_dim_max:
            _size = list(ds[v].sizes.values())
            _chunks = [np.max(d) for d in ds[v].chunks]
            n_dim_max = ds[v].ndim
    if n_dim_max > 0:
        print(
            "  typical chunks: ("
            + ",".join(f"{c}" for c in _chunks)
            + ") for size ("
            + ",".join(f"{c}" for c in _size)
            + ")"
        )
    else:
        print("  data is not chunked")
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
    """Returns the size of a directory in bytes"""
    process = os.popen("du -s " + dir_path)
    size = int(process.read().split()[0])  # du returns kb
    process.close()
    return size * 1e3
