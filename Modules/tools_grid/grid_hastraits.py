import numpy as np
import matplotlib
# We want matplotlib to use a wxPython backend
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_wx import NavigationToolbar2Wx

import scipy.interpolate as itp

from cartopy import crs as ccrs, feature as cfeature
import cartopy.io as cio
import cartopy.io.shapereader as shpreader

from threading import Thread
from time import sleep
import wx

from traits.api import *
from traitsui.api import View, Item, Group, HSplit, Handler, EnumEditor, FileEditor,DirectoryEditor
from traitsui.menu import NoButtons
from traitsui.wx.editor import Editor
#from traitsui.wx.basic_editor_factory import BasicEditorFactory
from traitsui.basic_editor_factory import BasicEditorFactory

import matplotlib.pyplot as plt
import tools
import tools_topo
import toolsf
import netCDF4 as netcdf
from datetime import datetime

class Inputs(HasTraits):
    """
    Inputs object
    """
    zview = Enum('grid outline', 'grid points', 'topo', '1/pm', '1/pn', 'angle', 'mask',
        desc="the data to view",
        label="View", )

    tra_lon = CFloat(15,
        desc="a central longitude",
        label="longitude", )

    tra_lat = CFloat(-32,
        desc="a central latitude",
        label="latitude", )

    size_x = CFloat(1556,
        desc="the mean distance along xi",
        label="x size [km]", )

    size_y = CFloat(1334,
        desc="the mean distance along eta",
        label="y size [km]", )

    rot = CFloat(0,
        desc="rotation about longitude, latitude",
        label="rotation [deg]", )

    nx = CInt(62,
        desc="the number of points along xi",
        label="nx", )

    ny = CInt(53,
        desc="the number of points along eta",
        label="ny", )

class Inputs_smth(HasTraits):
    """
    Inputs object for smoothing
    """
    depthmin=CFloat( 50,
            desc="minimum depth",
            label="Minimum depth [m]",)

    depthmax=CFloat( 6000,
            desc="maximum depth",
            label="Maximum depth [m]",)

    smthr=CFloat( 2,
            desc="smoothing radius",
            label="Smth radius [nb points]",)

    rfact=CFloat( 0.2,
            desc="maximum r-factor",
            label="r-factor",)

    smooth= Enum('smooth', 'lsmooth', 'lsmooth_legacy', 'lsmooth2', 'lsmooth1', 'cond_rx0_topo',
        desc="smoothing method",
        label="Smoothing method", )

class Inputs_smth_c2c(HasTraits):
    """
    Inputs object for smoothing
    """
    smthr=CFloat( 2,
            desc="smoothing radius",
            label="Smth radius [nb points]",)

    rfact=CFloat( 0.2,
            desc="maximum r-factor",
            label="r-factor",)

    smooth= Enum('smooth', 'lsmooth', 'lsmooth_legacy', 'lsmooth2', 'lsmooth1', 'cond_rx0_topo',
        desc="smoothing method",
        label="Smoothing method", )



class Inputs_zm(HasTraits):
    """
    Inputs object
    """
    tra_lon = CFloat(18,
        desc="a central longitude",
        label="longitude", )

    tra_lat = CFloat(-33,
        desc="a central latitude",
        label="latitude", )

    size_x = CFloat(550,
        desc="the mean distance along xi",
        label="x size [km]", )

    size_y = CFloat(550,
        desc="the mean distance along eta",
        label="y size [km]", )

    rot = CFloat(0,
        desc="rotation about longitude, latitude",
        label="rotation [deg]", )

    nx = CInt(55,
        desc="the number of points along xi",
        label="nx", )

    ny = CInt(55,
        desc="the number of points along eta",
        label="ny", )


class Inputs_c2c(HasTraits):
    """
    Inputs object
    """
    coef = CInt(3,
        desc="Refinement coefficient",
        label="Refinement coef")

    imin = CInt(35,
        desc="Parent imin",
        label="imin", )

    imax = CInt(55,
        desc="Parent imax",
        label="imax", )

    jmin = CInt(8,
        desc="Parent jmin",
        label="jmin", )

    jmax = CInt(28,
        desc="Parent jmax",
        label="jmax", )

    
class Outputs(HasTraits):
    """
    Outputs object
    """
    lon_rho = CArray()
    lat_rho = CArray()
    lon_u = CArray()
    lat_u = CArray()
    lon_v = CArray()
    lat_v = CArray()
    h = CArray()
    hraw = CArray()
    pm = CArray()
    pn = CArray()
    angle = CArray()
    f = CArray()
    mask_rho = CArray()

class BmapOptions(HasTraits):
    """
    Coastline options object
    """
    bmap_res = Enum('Crude', 'Low', 'Intermediate', 'High', 'Full',
        desc="the coastline resolution",
        label="Coastline resolution", )


class GetTopo(HasTraits):
    """
    GetTopo object.  At present this class will identify and read nc files from:
      ETOPO
      GEBCO
      Romstools (etopo2.nc; http://www.romsagrif.org)
    
    The default topography is an Etopo5 file (etopo5.nc).
    A file selector is available to point to choose alternative topo products.
    """
    def nearest1d(self, point, array):
        '''
        Return index of nearest point in array
        '''
        return np.argmin(np.abs(array - point))

    def topo(self,outputs, topo_file, smooth=None,hmin=None,hmax=None):

        if smooth is not None:
            rd=smooth.smthr
        else:
            rd = 1 

        topo_type=topo_reader.topo_file_id(topo_file)
        if topo_type['srtm']:
            srtm_file='/'.join(topo_file.split('/')[:-1])
            topo=toolsf.srtopo(srtm_file,outputs.lon_rho,outputs.lat_rho,outputs.pm,outputs.pn,rd)
        else:          
            lonmin,lonmax,latmin,latmax=toolsf.roms_grid_geo_bounds(outputs.lon_rho,outputs.lat_rho,rd)
            topo_lon,topo_lat,topo=tools_topo.topo_periodicity(topo_file,[lonmin,lonmax,latmin,latmax])
            print('Interpolating topography to ROMS grid')
            topo=toolsf.compute_hraw(topo_lon,topo_lat,topo.T,outputs.lon_rho,outputs.lat_rho,outputs.pm,outputs.pn,rd)
            print('Finished interpolating')

        if smooth is not None:
            outputs.hraw = topo
            if hmin is not None and hmax is not None: #Means you are in zoom AGRIF
                topo=eval(''.join(("toolsf.",smooth.smooth,"(topo,hmin,hmax, \
                                    smooth.rfact,outputs.mask_rho)")))
            else:
                topo=eval(''.join(("toolsf.",smooth.smooth,"(topo,smooth.depthmin,smooth.depthmax, \
                                    smooth.rfact,outputs.mask_rho)")))
            outputs.h=topo
            return outputs
        else:
            outputs.hraw = topo
            return outputs

    def match_topo(self,prt,outputs,WESN):
        WEST=True if 'West' in WESN[0] else False
        EAST=True if 'East' in WESN[0] else False
        SOUTH=True if 'South' in WESN[0] else False
        NORTH=True if 'North' in WESN[0] else False
        mask_tmp=np.copy(outputs.mask_rho)
        topo=toolsf.r2r_match_topo(WEST,EAST,SOUTH,NORTH,WESN[1],outputs.lon_rho.T,outputs.lat_rho.T,outputs.h.T,mask_tmp.T,\
                                  prt.lon_rho.T,prt.lat_rho.T,prt.h.T)
        outputs.h=topo.T
        return outputs



class GetMask(HasTraits):
     def outline(lon, lat):
        '''
        Return lon, lat of perimeter around the grid
        '''
        def func(var):
            return np.hstack([var[:, 0], var[-1, 1:-1],
                              var[::-1, -1], var[0, ::-1][1:]])
        return func(lon), func(lat)

     def mask(self, outputs,bmapoptions,gfile,sgl_connect=None):

         llcrnrlon = outputs.lon_rho[1:-1, 1:-1].min()
         urcrnrlon = outputs.lon_rho[1:-1, 1:-1].max()
         llcrnrlat = outputs.lat_rho[1:-1, 1:-1].min()
         urcrnrlat = outputs.lat_rho[1:-1, 1:-1].max()
    
         bmapoptions_dic = {'Crude':'c', 'Low':'l', 'Intermediate':'l',
                           'High':'h', 'Full':'f'}
         resolution = bmapoptions_dic[bmapoptions.bmap_res]
         
         rmask=toolsf.gshhs_to_roms_mask(outputs.lon_rho,outputs.lat_rho,gfile+'/gshhs_'+resolution+'.b')
         outputs.mask_rho=np.zeros(rmask.shape)
         outputs.mask_rho[rmask==0]=1
         if sgl_connect is not None:
             outputs.mask_rho=toolsf.single_connect(sgl_connect[1],sgl_connect[2],outputs.mask_rho.T).T
         return outputs

class EasyGrid(HasTraits):
    """
    EasyGrid object. Implements both the easygrid computation, and
    the picture acquisition.
    """
    '''def easygrid(self, inputs, outputs):
    
        lon_rho = np.linspace(inputs.tra_lon-size_x/111.,
                              inputs.tra_lon+size_x/111.,inputs.nx)
        lat_rho = np.linspace(inputs.tra_lat-size_y/111.,
                              inputs.tra_lat+size_y/111.,inputs.ny)
                            
        lon_rho, lat_rho = np.meshgrid(lon_rho,lat_rho)
        
        outputs.lon_rho = lon_rho
        outputs.lat_rho = lat_rho
        
        return outputs'''


    
    def easygrid(self, inputs, outputs):
        """
        Easy grid makes a rectangular, orthogonal grid with minimal gridsize variation
        It uses a Mercator projection around the equator and then rotates the sphere around its three axes to position the grid wherever it is desired.
   
        Inputs:
          nx:      Number of grid point in the x direction
          ny:      Number of grid point in the y direction
          size_x:  Domain size in x-direction
          size_y:  Domain size in y-direction
          tra_lon: Desired longitude of grid center
          tra_lat: Desired latitude of grid center
          rot:     Rotation of grid direction (0: x direction is west-east)

        Example:  > lon_rho,lat_rho,pm,pn,ang,lon_psi,lat_psi = easy_grid(30,20,4e6,3e6,-180,59,0)

        Translated to Python from EGRID (c) (Matlab) by Jeroen Molemaker, UCLA, 2008
        """
        def tra_sphere(lon1, lat1, tra):
            """
            Translate sphere about its y-axis
            Part of easy grid
            (c) 2008, Jeroen Molemaker, UCLA
            """

            n, m = lon1.shape
            tra = np.deg2rad(tra) # translation in latitude direction

            # translate into x,y,z
            # conventions:  (lon,lat) = (0,0)  corresponds to (x,y,z) = ( 0,-r, 0)
            #               (lon,lat) = (0,90) corresponds to (x,y,z) = ( 0, 0, r)
            x1 = np.sin(lon1) * np.cos(lat1)
            y1 = np.cos(lon1) * np.cos(lat1)
            z1 = np.sin(lat1)

            """
            We will rotate these points around the small circle defined by
            the intersection of the sphere and the plane that
            is orthogonal to the line through (lon,lat) (90,0) and (-90,0)

            The rotation is in that plane around its intersection with
            aforementioned line.

            Since the plane is orthogonal to the x-axis (in my definition at least),
            rotations in the plane of the small circle maintain constant x and are around
            (x,y,z) = (x1,0,0)
            """

            rp1 = np.sqrt(y1 ** 2 + z1 ** 2)

            ap1 = 0.5 * np.pi * np.ones((n, m))
            tmpi = np.abs(y1) > 1.e-5 
            ap1[tmpi] = np.arctan(np.abs(z1[tmpi] /
                                         y1[tmpi]))
            tmpi = y1 < 0.
            ap1[tmpi] = np.pi - ap1[tmpi]
            tmpi = z1 < 0.
            ap1[tmpi] = -ap1[tmpi]

            ap2 = ap1 + tra
            x2  = x1.copy()
            y2  = rp1 * np.cos(ap2)
            z2  = rp1 * np.sin(ap2)

            # transformation from (x,y,z) to (lat,lon)
            lon2 = 0.5 * np.pi * np.ones((n, m))
            tmpi = np.abs(y2) > 1.e-5
            lon2[tmpi] = np.arctan(np.abs(x2[tmpi] /
                                          y2[tmpi]))
            tmpi = y2 < 0.
            lon2[tmpi] = np.pi - lon2[tmpi]
            tmpi = x2 < 0.
            lon2[tmpi] = -lon2[tmpi]

            pr2 = np.sqrt(x2 ** 2 + y2 ** 2)
            lat2 = 0.5 * np.pi * np.ones((n, m))
            tmpi = np.abs(pr2) > 1.e-5
            lat2[tmpi] = np.arctan(np.abs(z2[tmpi] /
                                         pr2[tmpi]))
            tmpi = z2 < 0.
            lat2[tmpi] = -lat2[tmpi]

            return lon2, lat2
        
        
        def rot_sphere(lon1, lat1, rot):
            """
            Rotate sphere around its y-axis
            Part of Easy Grid
            (c) 2008, Jeroen Molemaker, UCLA
            """
            n, m = lon1.shape
            rot = np.deg2rad(rot)
            # translate into x,y,z
            # conventions:  (lon,lat) = (0,0)  corresponds to (x,y,z) = ( 0,-r, 0)
            #                   (lon,lat) = (0,90) corresponds to (x,y,z) = ( 0, 0, r)
            x1 = np.sin(lon1) * np.cos(lat1)
            y1 = np.cos(lon1) * np.cos(lat1)
            z1 = np.sin(lat1)
            """
            We will rotate these points around the small circle defined by
            the intersection of the sphere and the plane that
            is orthogonal to the line through (lon,lat) (0,0) and (180,0)

            The rotation is in that plane around its intersection with
            aforementioned line.

            Since the plane is orthogonal to the y-axis (in my definition at least),
            rotations in the plane of the small circle maintain constant y and are around
            (x,y,z) = (0,y1,0)
            """
            rp1 = np.sqrt(x1 ** 2 + z1 ** 2)

            ap1 = 0.5 * np.pi * np.ones((n, m))
            ap1[np.abs(x1) > 1.e-5] = np.arctan(np.abs(z1[np.abs(x1) > 1.e-5] /
                                                       x1[np.abs(x1) > 1.e-5]))
            ap1[x1 < 0.] = np.pi - ap1[x1 < 0.]
            ap1[z1 < 0.] = -ap1[z1 < 0.]

            ap2 = ap1 + rot
            x2 = rp1 * np.cos(ap2)
            y2 = y1.copy()
            z2 = rp1 * np.sin(ap2)

            lon2 = 0.5 * np.pi * np.ones((n, m))
            lon2[np.abs(y2) > 1.e-5] = np.arctan(np.abs(x2[np.abs(y2) > 1.e-5] /
                                                        y2[np.abs(y2) > 1.e-5]))
            lon2[y2 < 0.] = np.pi - lon2[y2 < 0.]
            lon2[x2 < 0.] = -lon2[x2 < 0.]

            pr2 = np.hypot(x2, y2)
            lat2 = 0.5 * np.pi * np.ones((n, m))
            lat2[np.abs(pr2) > 1.e-5] = np.arctan(np.abs(z2[np.abs(pr2) > 1.e-5] /
                                                        pr2[np.abs(pr2) > 1.e-5]))
            lat2[z2 < 0] = -lat2[z2 < 0]

            return lon2, lat2
        
        
        def gc_dist(lon1, lat1, lon2, lat2):
            '''
            Use Haversine formula to calculate distance
            between one point and another
            '''
            dlat = lat2 - lat1
            dlon = lon2 - lon1
            dang = 2. * np.arcsin(np.sqrt(np.power(np.sin(0.5 * dlat),2) + \
                        np.cos(lat2) * np.cos(lat1) * np.power(np.sin(0.5 * dlon),2)))
            return r_earth * dang # distance


        r_earth = 6371315. # Mean earth radius in metres (from scalars.h)
        
        size_x = inputs.size_x * 1000. # convert to meters
        size_y = inputs.size_y * 1000. # convert to meters

        # Mercator projection around the equator
        if size_y > size_x:
            length = np.float64(size_y)
            nl = np.float64(inputs.ny)
            width = np.float64(size_x)
            nw = np.float64(inputs.nx)
        else:
            length = np.float64(size_x)
            nl = np.float64(inputs.nx)
            width = np.float64(size_y)
            nw = np.float64(inputs.ny)
    
        dlon = length / r_earth
        lon1d = (dlon * np.arange(-0.5, nl + 1.5, 1.) / nl ) - ( 0.5 * dlon)

        mul = 1.
        dlat = width / r_earth

        for it in range(100):
            y1 = np.log(np.tan((0.25 * np.pi) - (0.25 * dlat)))
            y2 = np.log(np.tan((0.25 * np.pi) + (0.25 * dlat)))
            y = ((y2 - y1) * np.arange(-0.5, nw + 1.5, 1.) / nw ) + y1
            
            lat1d = np.arctan(np.sinh(y))
            dlat_cen = 0.5 * (lat1d[np.int32(np.round(0.5 * nw))] -  
                              lat1d[np.int32(np.round(0.5 * nw) - 2)])
            dlon_cen = dlon / nl
            mul = (dlat_cen / dlon_cen) * (length/width) * (nw / nl)
            dlat /= mul

        lon1de = (dlon * np.arange(-1., nl + 2., 1.) / nl ) - (0.5 * dlon)
        ye = ((y2-y1) * np.arange(-1., nw + 2., 1.) / nw ) + y1
        lat1de = np.arctan(np.sinh(ye))
        lat1de /= mul
        
        lon1, lat1 = np.meshgrid(lon1d, lat1d)
        lone, late = np.meshgrid(lon1de, lat1de)
        lonu = 0.5 * (lon1[:, :-1] + lon1[:, 1:])
        latu = 0.5 * (lat1[:, :-1] + lat1[:, 1:])
        lonv = 0.5 * (lon1[:-1] + lon1[1:])
        latv = 0.5 * (lat1[:-1] + lat1[1:])
        
        if size_y > size_x:
            lon1, lat1 = rot_sphere(lon1, lat1, 90.)
            lonu, latu = rot_sphere(lonu, latu, 90.)
            lonv, latv = rot_sphere(lonv, latv, 90.)
            lone, late = rot_sphere(lone, late, 90.)
            
            lon1 = lon1[::-1].T
            lat1 = lat1[::-1].T
            lone = lone[::-1].T
            late = late[::-1].T

            lonu_tmp = lonv[::-1].T
            latu_tmp = latv[::-1].T
            lonv = lonu[::-1].T
            latv = latu[::-1].T
            lonu = lonu_tmp
            latu = latu_tmp

        lon2, lat2 = rot_sphere(lon1, lat1, inputs.rot)
        lonu, latu = rot_sphere(lonu, latu, inputs.rot)
        lonv, latv = rot_sphere(lonv, latv, inputs.rot)
        lone, late = rot_sphere(lone, late, inputs.rot)

        lon3, lat3 = tra_sphere(lon2, lat2, inputs.tra_lat)
        lonu, latu = tra_sphere(lonu, latu, inputs.tra_lat)
        lonv, latv = tra_sphere(lonv, latv, inputs.tra_lat)
        lone, late = tra_sphere(lone, late, inputs.tra_lat)

        lon4 = lon3 + np.deg2rad(inputs.tra_lon)
        lonu += np.deg2rad(inputs.tra_lon)
        lonv += np.deg2rad(inputs.tra_lon)
        lone += np.deg2rad(inputs.tra_lon)
        
        pi2 = 2 * np.pi 
        lon4[lon4 < -np.pi] += pi2
        lonu[lonu < -np.pi] += pi2
        lonv[lonv < -np.pi] += pi2
        lone[lone < -np.pi] += pi2
        lat4 = lat3.copy()

        # Compute pm and pn
        pmu = gc_dist(lonu[:, :-1], latu[:, :-1],
                      lonu[:, 1:], latu[:, 1:])
        pm = np.zeros_like(lon4)
        pm[:, 1:-1] = pmu
        pm[:, 0] = pm[:, 1]
        pm[:, -1] = pm[:, -2]
        pm = 1. / pm
        
        pnv = gc_dist(lonv[:-1], latv[:-1],
                      lonv[1:], latv[1:])
        pn = np.zeros_like(lon4)
        pn[1:-1] = pnv
        pn[0] = pn[1]
        pn[-1] = pn[-2]
        pn = 1. / pn

        # Compute angles of local grid positive x-axis relative to east
        dellat = latu[:, 1:] - latu[:, :-1]
        dellon = lonu[:, 1:] - lonu[:, :-1]
        dellon[dellon > np.pi] -= pi2
        dellon[dellon < -np.pi] += pi2
        dellon = dellon * np.cos(0.5 * (latu[:, 1:] + latu[:, :-1]) )

        ang = np.zeros_like(lon4)
        ang_s = np.arctan(dellat / (dellon + 1.e-16))
        deli = np.logical_and(dellon < 0., dellat < 0.)
        ang_s[deli] -= np.pi
        
        deli = np.logical_and(dellon < 0., dellat >= 0.)
        ang_s[deli] += np.pi
        ang_s[ang_s > np.pi] -= np.pi
        ang_s[ang_s < -np.pi] += np.pi

        ang[:, 1:-1] = ang_s
        ang[:, 0] = ang[:, 1]
        ang[:, -1] = ang[:, -2]

        outputs.lon_rho = np.rad2deg(lon4)
        outputs.lat_rho = np.rad2deg(lat4)
        outputs.lon_u = np.rad2deg(lonu)
        outputs.lat_u = np.rad2deg(latu)
        outputs.lon_v = np.rad2deg(lonv)
        outputs.lat_v = np.rad2deg(latv)
        outputs.pm = pm
        outputs.pn = pn
        outputs.angle = ang
        outputs.lon_psi = np.rad2deg(lone)
        outputs.lat_psi = np.rad2deg(late)
        
        return outputs

    def AGRIFgrid(self,prt_grd,inputs,outputs):
        def gc_dist(lon1, lat1, lon2, lat2):
            r_earth=6371315. # Mean earth radius in metres (from scalars.h
            '''
            Use Haversine formula to calculate distance
            between one point and another
            '''
            dlat = lat2 - lat1
            dlon = lon2 - lon1
            dang = 2. * np.arcsin(np.sqrt(np.power(np.sin(0.5 * dlat),2) + \
                        np.cos(lat2) * np.cos(lat1) * np.power(np.sin(0.5 * dlon),2)))
            return r_earth * dang # distance

        [Mp,Lp]=prt_grd.lon_rho.shape
        
        igrdp=np.arange(0,Lp-1);jgrdp=np.arange(0,Mp-1)
        igrdr=np.arange(0,Lp);jgrdr=np.arange(0,Mp)
        igrdu=np.arange(0,Lp-1);jgrdu=np.arange(0,Mp)
        igrdv=np.arange(0,Lp);jgrdv=np.arange(0,Mp-1)


        ipch=np.arange(inputs.imin,inputs.imax+0.5/inputs.coef,1/inputs.coef)
        jpch=np.arange(inputs.jmin,inputs.jmax+0.5/inputs.coef,1/inputs.coef)
        ################
        irch=np.arange(inputs.imin+0.5-0.5/inputs.coef,inputs.imax+0.5+0.75/inputs.coef,1/inputs.coef)
        jrch=np.arange(inputs.jmin+0.5-0.5/inputs.coef,inputs.jmax+0.5+0.75/inputs.coef,1/inputs.coef)

        [ichildgrd_p,jchildgrd_p]=np.meshgrid(ipch,jpch)
        [ichildgrd_r,jchildgrd_r]=np.meshgrid(irch,jrch)
        [ichildgrd_u,jchildgrd_u]=np.meshgrid(ipch,jrch)
        [ichildgrd_v,jchildgrd_v]=np.meshgrid(irch,jpch)
         
        spline_lonp    = itp.RectBivariateSpline(jgrdp,igrdp,prt_grd.lon_psi)
        spline_latp    = itp.RectBivariateSpline(jgrdp,igrdp,prt_grd.lat_psi)
        outputs.lon_psi = spline_lonp(jchildgrd_p,ichildgrd_p,grid=False)
        outputs.lat_psi = spline_latp(jchildgrd_p,ichildgrd_p,grid=False)
        ###############
        spline_lonr    = itp.RectBivariateSpline(jgrdr,igrdr,prt_grd.lon_rho)
        spline_latr    = itp.RectBivariateSpline(jgrdr,igrdr,prt_grd.lat_rho)
        outputs.lon_rho = spline_lonr(jchildgrd_r,ichildgrd_r,grid=False)
        outputs.lat_rho = spline_latr(jchildgrd_r,ichildgrd_r,grid=False)
        ###############
        spline_lonu    = itp.RectBivariateSpline(jgrdu,igrdu,prt_grd.lon_u)
        spline_latu    = itp.RectBivariateSpline(jgrdu,igrdu,prt_grd.lat_u)
        outputs.lon_u = spline_lonr(jchildgrd_u,ichildgrd_u,grid=False)
        outputs.lat_u = spline_latr(jchildgrd_u,ichildgrd_u,grid=False)
        ###############
        spline_lonv    = itp.RectBivariateSpline(jgrdv,igrdv,prt_grd.lon_v)
        spline_latv    = itp.RectBivariateSpline(jgrdv,igrdv,prt_grd.lat_v)
        outputs.lon_v = spline_lonr(jchildgrd_v,ichildgrd_v,grid=False)
        outputs.lat_v = spline_latr(jchildgrd_v,ichildgrd_v,grid=False)

        #################

        # Compute pm and pn
        pmu = gc_dist(np.deg2rad(outputs.lon_u[:, :-1]), np.deg2rad(outputs.lat_u[:, :-1]),
                      np.deg2rad(outputs.lon_u[:, 1:]), np.deg2rad(outputs.lat_u[:, 1:]))
        pm = np.zeros_like(outputs.lon_rho)
        pm[:, 1:-1] = pmu
        pm[:, 0] = pm[:, 1]
        pm[:, -1] = pm[:, -2]
        pm = 1. / pm

        pnv = gc_dist(np.deg2rad(outputs.lon_v[:-1]),np.deg2rad(outputs.lat_v[:-1]),
                      np.deg2rad(outputs.lon_v[1:]), np.deg2rad(outputs.lat_v[1:]))
        pn = np.zeros_like(outputs.lon_rho)
        pn[1:-1] = pnv
        pn[0] = pn[1]
        pn[-1] = pn[-2]
        pn = 1. / pn

        # Compute angles of local grid positive x-axis relative to east
        pi2 = 2 * np.pi
        dellat = np.deg2rad(outputs.lat_u[:, 1:]) - np.deg2rad(outputs.lat_u[:, :-1])
        dellon = np.deg2rad(outputs.lon_u[:, 1:]) - np.deg2rad(outputs.lon_u[:, :-1])
        dellon[dellon > np.pi] -= pi2
        dellon[dellon < -np.pi] += pi2
        dellon = dellon * np.cos(0.5 * (np.deg2rad(outputs.lat_u[:, 1:]) +np.deg2rad(outputs.lat_u[:, :-1])) )

        ang = np.zeros_like(outputs.lon_rho)
        ang_s = np.arctan(dellat / (dellon + 1.e-16))
        deli = np.logical_and(dellon < 0., dellat < 0.)
        ang_s[deli] -= np.pi

        deli = np.logical_and(dellon < 0., dellat >= 0.)
        ang_s[deli] += np.pi
        ang_s[ang_s > np.pi] -= np.pi
        ang_s[ang_s < -np.pi] += np.pi

        ang[:, 1:-1] = ang_s
        ang[:, 0] = ang[:, 1]
        ang[:, -1] = ang[:, -2]

        outputs.pm = pm
        outputs.pn = pn
        outputs.angle = ang

        return outputs


class Save2Netcdf(HasTraits):
    """
    Save2Netcdf object
    """
    def save2netcdf(self,output_dir, inputs, outputs,prt_grd=None):
    
        """
        Create and save a new CROCO grid file
        """
#        prt_grd=[AGRIF,prt_file,coef,imi,imax,jmin,jmax]
        if prt_grd is not None:
            if prt_grd[0]==True: # Means we are in AGRIF
                lev=prt_grd[1][-1]
                if not lev.isnumeric():
                    grid_name='croco_grd.nc.1'
                else:
                    print(prt_grd[1][0:-1],lev)
                    grid_name=prt_grd[1][0:-1]+str(int(lev)+1)

            else :  # We create offline zoom
                grid_name='croco_chd_grd.nc'
        else:
            grid_name='croco_grd.nc'
         
        nc = netcdf.Dataset(output_dir+grid_name, 'w', format='NETCDF4')

        # create global variables
        nc.created = datetime.now().isoformat()
        nc.type = 'ROMS grid file produced by easygrid_python.py'
        nc.VertCoordType = 'NEW';

        nc.createDimension('one', 1)
        if prt_grd is not None and prt_grd[0]: #AGRIF case
            nc.nx=outputs.h.shape[1]-2
            nc.ny=outputs.h.shape[0]-2

            # create dimensions
            nc.createDimension('xi_rho', outputs.h.shape[1])
            nc.createDimension('eta_rho', outputs.h.shape[0])
            nc.createDimension('xi_u', outputs.h.shape[1] - 1)
            nc.createDimension('eta_v', outputs.h.shape[0] - 1)
            nc.createDimension('xi_psi', outputs.h.shape[1] - 1)
            nc.createDimension('eta_psi', outputs.h.shape[0] - 1)
            nc.createDimension('four', 4)

            # Some empty variables in AGRIF
            nc.createVariable('xl', 'f8', ('one'))
            nc.variables['xl'].long_name = 'domain length in the XI-direction'
            nc.variables['xl'].units = 'meters'

            nc.createVariable('el', 'f8', ('one'))
            nc.variables['el'].long_name = 'domain length in the ETA-direction'
            nc.variables['el'].units = 'meters'
            nc.variables['el'][:] = inputs.ny

        else: # Usual case

            nc.nx = np.int32(inputs.nx)
            nc.ny = np.int32(inputs.ny)
            nc.size_x = inputs.size_x
            nc.size_y = inputs.size_y
            nc.tra_lon = inputs.tra_lon
            nc.tra_lat = inputs.tra_lat
            nc.rotation = inputs.rot
            
            # create dimensions
            nc.createDimension('xi_rho', inputs.nx + 2)
            nc.createDimension('eta_rho', inputs.ny + 2)
            nc.createDimension('xi_u', inputs.nx + 1)
            nc.createDimension('eta_v', inputs.ny + 1)
            nc.createDimension('xi_psi', inputs.nx + 1)
            nc.createDimension('eta_psi', inputs.ny + 1)

            nc.createVariable('xl', 'f8', ('one'))
            nc.variables['xl'].long_name = 'domain length in the XI-direction'
            nc.variables['xl'].units = 'meters'
            nc.variables['xl'][:] = inputs.nx

            nc.createVariable('el', 'f8', ('one'))
            nc.variables['el'].long_name = 'domain length in the ETA-direction'
            nc.variables['el'].units = 'meters'
            nc.variables['el'][:] = inputs.ny

        
        # create variables and attributes
        nc.createVariable('spherical', 'S1', ('one'))
        nc.variables['spherical'].long_name = 'Grid type logical switch'
        nc.variables['spherical'].option_T = 'spherical'
        nc.variables['spherical'][:] = 'T'
        
        nc.createVariable('angle', 'f8', ('eta_rho', 'xi_rho'))
        nc.variables['angle'].long_name = 'angle between xi axis and east'
        nc.variables['angle'].units = 'radians' 
        nc.variables['angle'][:] = outputs.angle

        nc.createVariable('h', 'f8', ('eta_rho', 'xi_rho'))
        nc.variables['h'].long_name = 'Final bathymetry at RHO-points'
        nc.variables['h'].units = 'meter'
        nc.variables['h'][:] = outputs.h

        nc.createVariable('hraw', 'f8', ('eta_rho', 'xi_rho'))
        nc.variables['hraw'].long_name = 'Working bathymetry at RHO-points'
        nc.variables['hraw'].units = 'meter'
        nc.variables['hraw'][:] = outputs.hraw

        nc.createVariable('f', 'f8', ('eta_rho', 'xi_rho'))
        nc.variables['f'].long_name = 'Coriolis parameter at RHO-points'
        nc.variables['f'].units = 'second-1'
        nc.variables['f'][:] = (4 * np.pi * np.sin(np.deg2rad(outputs.lat_rho)) /
                                (23.9344699 * 3600))

        nc.createVariable('pm', 'f8', ('eta_rho', 'xi_rho'))
        nc.variables['pm'].long_name = 'curvilinear coordinate metric in XI'
        nc.variables['pm'].units = 'meter-1'
        nc.variables['pm'][:] = outputs.pm

        nc.createVariable('pn', 'f8', ('eta_rho', 'xi_rho'))
        nc.variables['pn'].long_name = 'curvilinear coordinate metric in ETA'
        nc.variables['pn'].units = 'meter-1'
        nc.variables['pn'][:] = outputs.pn

        nc.createVariable('lon_rho', 'f8', ('eta_rho', 'xi_rho'))
        nc.variables['lon_rho'].long_name = 'longitude of RHO-points'
        nc.variables['lon_rho'].units = 'degree_east'
        nc.variables['lon_rho'][:] = outputs.lon_rho

        nc.createVariable('lat_rho', 'f8', ('eta_rho', 'xi_rho'))
        nc.variables['lat_rho'].long_name = 'latitude of RHO-points'
        nc.variables['lat_rho'].units = 'degree_north'
        nc.variables['lat_rho'][:] = outputs.lat_rho

        nc.createVariable('mask_rho', 'f8', ('eta_rho', 'xi_rho'))
        nc.variables['mask_rho'].long_name = 'mask on RHO-points'
        nc.variables['mask_rho'].option_0 = 'land'
        nc.variables['mask_rho'].option_1 = 'water'
        nc.variables['mask_rho'][:] = outputs.mask_rho

        # Extraneous variables should be placed at the end (ensures no
        # later problems with e.g., partit
        nc.createVariable('lon_psi', 'f8', ('eta_psi', 'xi_psi'))
        nc.variables['lon_psi'].long_name = 'longitude of PSI-points'
        nc.variables['lon_psi'].units = 'degree_east'
        if prt_grd is not None and prt_grd[0]:
            nc.variables['lon_psi'][:] = outputs.lon_psi
        else:
            nc.variables['lon_psi'][:] = outputs.lon_psi[1:-1, 1:-1]

        nc.createVariable('lat_psi', 'f8', ('eta_psi', 'xi_psi'))
        nc.variables['lat_psi'].long_name = 'latitude of PSI-points'
        nc.variables['lat_psi'].units = 'degree_north'
        if prt_grd is not None and prt_grd[0]:
            nc.variables['lat_psi'][:] = outputs.lat_psi
        else:
            nc.variables['lat_psi'][:] = outputs.lat_psi[1:-1, 1:-1]

        nc.createVariable('lon_u', 'f8', ('eta_rho', 'xi_u'))
        nc.variables['lon_u'].long_name = 'longitude of U-points'
        nc.variables['lon_u'].units = 'degree_east'
        nc.variables['lon_u'][:] = outputs.lon_u

        nc.createVariable('lat_u', 'f8', ('eta_rho', 'xi_u'))
        nc.variables['lat_u'].long_name = 'latitude of U-points'
        nc.variables['lat_u'].units = 'degree_north'
        nc.variables['lat_u'][:] = outputs.lat_u

        nc.createVariable('lon_v', 'f8', ('eta_v', 'xi_rho'))
        nc.variables['lon_v'].long_name = 'longitude of V-points'
        nc.variables['lon_v'].units = 'degree_east'
        nc.variables['lon_v'][:] = outputs.lon_v

        nc.createVariable('lat_v', 'f8', ('eta_v', 'xi_rho'))
        nc.variables['lat_v'].long_name = 'latitude of RHO-points'
        nc.variables['lat_v'].units = 'degree_north'
        nc.variables['lat_v'][:] = outputs.lat_v

        if prt_grd is not None and prt_grd[0]: 
            nc.createVariable('refine_coef', 'i', ('one'))
            nc.variables['refine_coef'].long_name ='Grid refinement coefficient'
            nc.variables['refine_coef'][:]=prt_grd[2]

            nc.createVariable('grd_pos','i',('four'))
            nc.variables['grd_pos'].long_name='Subgrid location in the parent grid: psi corner points (imin imax jmin jmax)'
            nc.variables['grd_pos'][:]=prt_grd[3:]


        nc.close()
        print('Writting '+grid_name+' done')

        if prt_grd is not None and prt_grd[0]:
            print('Create an AGRIF_FixedGrids.in file')
            fname='AGRIF_FixedGrids.in'
            fid=open(fname,'w')
            fid.write('    1\n')#'%s\n','    1');
            fid.write('    '+str(prt_grd[3])+ \
                               '    '+str(prt_grd[4])+ \
                               '    '+str(prt_grd[5])+\
                               '    '+str(prt_grd[6])+\
                               '    '+str(prt_grd[2])+\
                               '    '+str(prt_grd[2])+\
                               '    '+str(prt_grd[2])+\
                               '    '+str(prt_grd[2]))
            fid.write('\n    0')
            fid.write('\n# number of children per parent')
            fid.write('\n# imin imax jmin jmax spacerefx spacerefy timerefx timerefy')
            fid.write('\n# [all coordinates are relative to each parent grid!]')
            fid.write('\n~')
            fid.close()

