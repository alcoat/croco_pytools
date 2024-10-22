
import numpy as np
import matplotlib
import matplotlib.pyplot
#import make_grid_param as param

from cartopy import crs as ccrs, feature as cfeature
import cartopy.io.shapereader as shpreader


def plot_outline(outputs, figure):

    def outline(lon, lat):
        '''
        Return lon, lat of perimeter around the grid
        '''
        def func(var):
            return np.hstack([var[:, 0], var[-1, 1:-1],
                              var[::-1, -1], var[0, ::-1][1:]])
        return func(lon), func(lat)

    x_rho = outputs.lon_rho

    ox_psi, oy_psi = outline(
        outputs.lon_psi[1:-1, 1:-1], outputs.lat_psi[1:-1, 1:-1])

    if min(x_rho.ravel()) < 180 and max(x_rho.ravel()) > 180:
        # check if grid cross lon=180 and change central lon for cartopy
        mid = 180
        x_rho = x_rho-180
        ox_psi = ox_psi-180
    else:
        mid = 0

    ax = figure.add_axes(rect=[0.1, 0.04, 0.80, 0.9],
                         projection=ccrs.PlateCarree(central_longitude=mid))
    ax.plot(ox_psi, oy_psi, 'r', zorder=5)  # , ax=figure.gca(), zorder=3)

    return ax


#def plot_grid(outputs, figure, ax, zview='grid outline', shp_file=param.shp_file, plot_shape=True):
def plot_grid(outputs, figure, ax, shp_file, zview='grid outline', plot_shape=True):
    x_rho, y_rho = outputs.lon_rho, outputs.lat_rho
    x_psi, y_psi = outputs.lon_psi, outputs.lat_psi

    # Colorbar position
    axpos = ax.get_position()
    pos_x = axpos.x0  # + 0.25*axpos.width
    pos_y = axpos.y0-0.08  # -axpos.height - 0.02
    cax_width = axpos.width
    cax_height = 0.04
    # pos_cax = figure.add_axes([pos_x,pos_y,cax_width,cax_height])

    if 'grid points' in zview:
        figure.suptitle('Grid points: rho=green, psi=blue')
        ax.scatter(x_rho, y_rho, s=5, c='g', edgecolor='g', zorder=3)
        ax.scatter(x_psi, y_psi, s=5, c='b', edgecolor='b', zorder=3)
    elif 'topo' in zview:
        cb = ax.pcolormesh(x_rho, y_rho, outputs.hraw,
                           shading='auto', zorder=2)
        pos_cax = figure.add_axes([pos_x, pos_y, cax_width, cax_height])
        figure.colorbar(cb, cax=pos_cax, orientation='horizontal')
        #M.scatter(x_rho, y_rho, s=2, c='w', edgecolor='w', ax=figure.gca(), zorder=3)
    elif '1/pm' in zview:
        cb = ax.pcolormesh(x_rho, y_rho, 1 / outputs.pm, zorder=2)
        #M.scatter(x_rho, y_rho, s=2, c='w', edgecolor='w', ax=figure.gca(), zorder=3)
        pos_cax = figure.add_axes([pos_x, pos_y, cax_width, cax_height])
        figure.colorbar(cb, cax=pos_cax, orientation='horizontal')

    elif '1/pn' in zview:
        cb = ax.pcolormesh(x_rho, y_rho, 1 / outputs.pn, zorder=2)
        #M.scatter(x_rho, y_rho, s=2, c='w', edgecolor='w', ax=figure.gca(), zorder=3)
        pos_cax = figure.add_axes([pos_x, pos_y, cax_width, cax_height])
        figure.colorbar(cb, cax=pos_cax, orientation='horizontal')

    elif 'angle' in zview:
        cb = ax.pcolormesh(x_rho, y_rho, outputs.angle, zorder=2)
        #M.scatter(x_rho, y_rho, s=2, c='w', edgecolor='w', ax=figure.gca(), zorder=3)
        pos_cax = figure.add_axes([pos_x, pos_y, cax_width, cax_height])
        figure.colorbar(cb, cax=pos_cax, orientation='horizontal')

    elif 'mask' in zview:
        cb = ax.pcolormesh(x_rho, y_rho, outputs.mask_rho,
                           cmap=matplotlib.cm.Spectral, zorder=2)
    #M.scatter(x_rho, y_rho, s=2, c='k', edgecolor='k', ax=figure.gca(), zorder=3)

    if not 'mask' in zview:
        if plot_shape:
            tch = shpreader.Reader(shp_file)
            # here ccrs define the projection of the data (coastline)
            shape_feature = cfeature.ShapelyFeature(
                tch.geometries(), ccrs.PlateCarree())
            ax.add_feature(shape_feature, facecolor='#929591',
                           edgecolor='Gainsboro', zorder=4)
        gl = ax.gridlines(draw_labels=True, linewidth=2,
                          color='gray', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False

    return ax


def plot_topo(outputs, figure, ax):
    x_rho, y_rho = outputs.lon_rho, outputs.lat_rho

    # Colorbar position
    axpos = ax.get_position()
    pos_x = axpos.x0  # + 0.25*axpos.width
    pos_y = axpos.y0-0.08  # -axpos.height - 0.02
    cax_width = axpos.width
    cax_height = 0.04

    ax.pcolormesh(x_rho, y_rho, np.ma.masked_where(outputs.mask_rho > 0, outputs.mask_rho),
                  vmin=0, vmax=1, zorder=2, cmap=matplotlib.pyplot.cm.copper_r)

    plotTopo = ax.pcolormesh(x_rho, y_rho, np.ma.masked_where(
        outputs.mask_rho < 1, outputs.h), zorder=2)
    
    
    pos_cax = figure.add_axes([pos_x, pos_y, cax_width, cax_height])
    
    cb = figure.colorbar(plotTopo, cax=pos_cax, orientation='horizontal')
    
    #gl = ax.gridlines(draw_labels=True, linewidth=2,
    #                  color='gray', alpha=0.5, linestyle='--')
    #gl.top_labels = False
    #gl.right_labels = False
    
    # Add gridlines
    ax.grid(True, which='both', linewidth=2, color='gray', alpha=0.5, linestyle='--')

    # Optionally adjust grid labels
    ax.xaxis.set_major_formatter(matplotlib.pyplot.FuncFormatter(lambda x, _: f'{x:.1f}'))
    ax.yaxis.set_major_formatter(matplotlib.pyplot.FuncFormatter(lambda y, _: f'{y:.1f}'))
    

    return ax, cb
