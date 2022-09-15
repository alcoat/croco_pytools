import os

import numpy as np

from dask import delayed

from matplotlib import pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import cartopy.crs as ccrs

# check if python run in batch mode or from jupyter (for plot in batch mode)
import psutil
def running_in_jupyter():
    """Return True if any of our parent processes is jupyter"""
    parent_names = [parent.name() for parent in psutil.Process().parents()]
    return any('jupyter' in string for string in parent_names)
if not running_in_jupyter():
    import matplotlib
    matplotlib.use('AGG')

import scipy.io
from collections import OrderedDict


import gridop as gop
from tools import dask_compute_batch

# -------------------------------- Colormap -------------------------------

def get_cmap_colors(Nc, cmap='plasma'):
    """ load colors from a colormap to plot lines
    
    Parameters
    ----------
    Nc: int
        Number of colors to select
    cmap: str, optional
        Colormap to pick color from (default: 'plasma')
    """
    scalarMap = cmx.ScalarMappable(norm=colors.Normalize(vmin=0, vmax=Nc),
                                   cmap=cmap)
    return [scalarMap.to_rgba(i) for i in range(Nc)]

def DefCmap():
    """ construct a python colormap from matlab mat file, stored in the same directory as this module"""
    matfile = scipy.io.loadmat(os.path.join(os.path.dirname(__file__),'map_64_wc.mat'))
    return array2cmap(np.array(matfile['cm']))


def array2cmap(X):
    N = X.shape[0]

    r = np.linspace(0., 1., N + 1)
    r = np.sort(np.concatenate((r, r)))[1:-1]

    rd = np.concatenate([[X[i, 0], X[i, 0]] for i in range(N)])
    gr = np.concatenate([[X[i, 1], X[i, 1]] for i in range(N)])
    bl = np.concatenate([[X[i, 2], X[i, 2]] for i in range(N)])

    rd = tuple([(r[i], rd[i], rd[i]) for i in range(2 * N)])
    gr = tuple([(r[i], gr[i], gr[i]) for i in range(2 * N)])
    bl = tuple([(r[i], bl[i], bl[i]) for i in range(2 * N)])

    cdict = {'red': rd, 'green': gr, 'blue': bl}
    return colors.LinearSegmentedColormap('my_colormap', cdict, N)#

# -------------------------------- Images --------------------------------------



def plotfig(da, numimage=0, fig_dir=None, fig_suffix=None, date=' ', save=False, 
                cmap=None, figsize=(10,8), dpi=150, **kwargs):
    '''
    Plot an 2d xarray DataArray
    '''
    # Init, create the directory where to generate the figure
    if fig_dir is None:
        try:
            fig_dir = os.environ['SCRATCH']+'/figs/'
        except:
            fig_dir = os.environ['PWD']+'/figs/'
    if not os.path.isdir(fig_dir):
        os.mkdir(fig_dir)
    
    # Init the figure name
    if fig_suffix is None:
        if hasattr(da,'name'): fig_suffix = ''.join(da.name) 
    figname = fig_dir+fig_suffix+'_t%05d' %(numimage)+'.png'

    # Define the colormap
    if cmap is None: cmap = DefCmap()

    # retrieve dimensions and coordinates of the variable
    dims = gop.get_spatial_dims(da.squeeze())
    # delete None dimension
    dims = OrderedDict([(k,v) for k,v in dims.items() if v is not None])
    coords = gop.get_spatial_coords(da.squeeze())
    
    if 'x' in dims.keys() and 'y' in dims.keys():
        # horizontal section
        coordx = coords['lon']
        coordy = coords['lat']
        # prepare the plot
        # First the Map Projection
        projection = ccrs.Mercator()
        # Specify the CRS (coordinate reference system)
        crs = ccrs.PlateCarree()
        fig = plt.figure(figsize=figsize)
        ax = plt.axes(projection = projection)
        gl = ax.gridlines(crs=crs, draw_labels=True, linewidth=.6,
                          color='gray', alpha=0.5, linestyle='dashed')
        gl.top_labels = False
        gl.right_labels = False 
        ax.coastlines() 
        da.plot(x=coordx, y=coordy, ax=ax, cmap=cmap, transform=crs, **kwargs) 
    else:
        # vertical section
        if 'x' in dims.keys(): coordx = coords['lon']
        if 'y' in dims.keys(): coordx = coords['lat']
        coordy = coords['z']
        fig, ax = plt.subplots(figsize=figsize)
        da.plot(x=coordx, y=coordy, ax=ax, cmap=cmap, **kwargs) 
        ax.grid(color='gray', alpha=0.5, linestyle='dashed')
    
    # put the title
    if 't' in da.coords: date = np.datetime_as_string(da.t, unit='m')
    title = fig_suffix+', date = %s'%(date)
    ax.set_title(title)
    
    # save in a file
    if save: 
        fig.savefig(figname, dpi=dpi)
        plt.close()

    return None

# -------------------------------- movies --------------------------------------

def movie_wrapper(da, client, fig_dir=None, fig_suffix=None, figsize=(10,8),  
                      dpi=150, fps=5, **kwargs):

    if fig_dir is None:
        try:
            fig_dir = os.environ['SCRATCH']+'/figs/'
        except:
            fig_dir = os.environ['PWD']+'/figs/'
    if not os.path.isdir(fig_dir):
        os.mkdir(fig_dir)
    
    if fig_suffix is None:
        if hasattr(da,'name'): fig_suffix = ''.join(da.name) 
        
    # plotfig_delayed = delayed(plotfig)
    tasks = [
        delayed(plotfig)(da.isel(t=i), i, fig_dir=fig_dir, 
                         fig_suffix=fig_suffix, 
                         save=True, dpi=dpi, figsize=figsize, **kwargs,
                         ) for i in range(da.t.size)
    ]

    dask_compute_batch(tasks, client)
    
    fig_name = fig_dir+fig_suffix+'_t%05d.png'
    movie_name = fig_dir+fig_suffix+'.mp4'
    commande = "ffmpeg -framerate "+str(fps)+"  -i "+fig_name+" -vcodec mpeg4 -y "+movie_name
    # print(commande)
    os.system(commande)

    figure_name = fig_dir+fig_suffix+'_*.png'
    commande = "rm -rf "+figure_name
    os.system(commande)