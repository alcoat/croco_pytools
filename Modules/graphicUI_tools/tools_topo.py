import numpy as np
import netCDF4 as netcdf
import topo_reader
from datetime import datetime

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
            while x[i] > x0 :   # because array elements x(i) are
                i=i-1           # equidistantly spaced.  The while
                                # loops are here merely to address
                                # possible roundoff errors.

        if x[i+1]-x0 < 0 or x0-x[i] < 0 :
            print('### ERROR: indx_bound :: ',x[i], x0, x[i+1], x0-x[i], x[i+1]-x0)
            exit()
    indx_bound=i
    return indx_bound


def topo_periodicity(topo_file, geolim):
    '''
    topo_periodicity checks whether domain is inside the topo file.
    If so, check if there is a need for create a periodicity between
    the last and first longitude points ( for global data).
    It is returning lon/lat/topo adapted to the desired domain
    geolim = [lonmin,lonmax,latmin,latmax]
    '''

    topo_type = topo_reader.topo_file_id(topo_file)

    print('Reading topography file:', topo_file)
    nc = netcdf.Dataset(topo_file)
    topo_lon = nc.variables[topo_type['lon']][:]
    topo_lat = nc.variables[topo_type['lat']][:]
    if topo_lon.size==2: # gebco is a bit different
        topo_lon = np.linspace(topo_lon[0],
                               topo_lon[1], num=nc.variables['dimension'][0])
        topo_lat = np.linspace(topo_lat[0],
                               topo_lat[1], num=nc.variables['dimension'][1])[::-1]
        gebco = True
    else:
        gebco = False

    for i in range(1,topo_lon.shape[0]): # Fix etopo5 discontinuity
        if topo_lon[i]<topo_lon[i-1]:    # between 180/-180 in the
            topo_lon[i]=topo_lon[i]+360  # middle

####
    jmin=indx_bound(topo_lat, geolim[2])
    jmax=indx_bound(topo_lat, geolim[-1])
    if 0 < jmin and jmin < topo_lat.shape[0] and 0 < jmax and jmax < topo_lat.shape[0] :
        if jmin > 1 :
            jmin=jmin-1
        jmax=jmax+2
    else:
        print('North-south extents of the dataset ',topo_lat[0],topo_lat[-1],' are not sufficient to cover the entire model grid.')
        exit()

    imin=indx_bound(topo_lon, geolim[0])
    imax=indx_bound(topo_lon, geolim[1])
    
    if 0 < imin and imin < topo_lon.shape[0] and 0 < imax and imax < topo_lon.shape[0] :
        if imax > 1:
            imin=imin-1
        imax=imax+2
        shft_west=0 ; shft_east=0
        print('Single region dataset imin/imax=',imin,imax, )
    else:
        ######
        ptest=topo_lon[-1]-topo_lon[0]-360
        dx=(topo_lon[-1]-topo_lon[0])/(topo_lon.shape[0]-1)
        epsil=0.01*abs(dx)
        if abs(ptest) < epsil :
            period=topo_lon.shape[0]-1
        elif abs(ptest+dx) < epsil :
            period=topo_lon.shape[0]
        else:
            period=0

        if period>0:
            print('Identified periodicity domain in data of ', period,' points out of', topo_lon.shape[0])
        else :
            print('ERROR: The data does not cover the entire grid. Change your grid definition')
            exit()
        ##
        shft_west=0
        if imin==0 :
            shft_west=-1
            imin=indx_bound(topo_lon, geolim[0]+360)
        elif imin==topo_lon.shape[0] :
            shft_west=+1
            imin=indx_bound(topo_lon, geolim[0]-360)
        ##
        shft_east=0
        if imax == 0:
            shft_east=-1
            imax=indx_bound(topo_lon, geolim[1]+360)
        elif imax == topo_lon.shape[0]:
            shft_east=+1
            imax=indx_bound(topo_lon, geolim[1]-360)

        if 0<imin and imin <topo_lon.shape[0] and 0<imax and imax<topo_lon.shape[0] :
            if imin>1:
                imin=imin-1
            imax=imax+1
        else:
            print('ERROR: Data longitude covers 360 degrees, but still cannot find  starting and ending indices.')
            exit()
    
    print('Bounding indices of the relevant part to be extracted from the entire dataset:\n', \
          'imin,imax =', imin,imax,'out of', topo_lon.shape[0],'jmin,jmax =',jmin,jmax, 'out of',topo_lat.shape[0])
    ny_lat=jmax-jmin+1
    start2=jmin ; end2=start2+ny_lat; count2=ny_lat
    lat_tmp=np.zeros([ny_lat])
    for j in range(0,ny_lat):
        lat_tmp[j]=topo_lat[j+jmin-1]
 
    #####

    if imin < imax :
        nx_lon=imax-imin+1
        start1=imin ; end1=start1+nx_lon ; count1=nx_lon
        if gebco:
            topo = nc.variables[topo_type['topo']][:]
            topo = np.reshape(topo, (topo_lat.size, topo_lon.size))
            topo = topo[start2:end2, start1:end1]
        else:
            topo = nc.variables[topo_type['topo']][start2:end2, start1:end1]
        nc.close()

        ishft=imin-1
        lon_tmp=np.zeros([topo.shape[1]])
        if shft_west>0 and shft_east>0:
            for i in range(0,nx_lon):
                lon_tmp[i]=topo_lon[i+ishft] +360
        elif shft_west<0 and shft_east<0:
            for i in range(0,nx_lon):
                 lon_tmp[i]=topo_lon[i+ishft]-360
        elif shft_west== 0 and shft_east==0:
            for i in range(0,nx_lon) :
                lon_tmp[i]=topo_lon[i+ishft]
        else:
            print('Error in shifting algoritm')
            exit()

    elif imin>imax:
        print('Reading topography in two separate parts adjacent through 360-degree periodicity, first...' )

        nx_lon=imax+period-imin+1
        htopo = np.zeros([ny_lat,nx_lon])
        xtmp  = np.zeros([nx_lon])
        start1=0 ; end1=start1+nx_lon; count1=imax
        if gebco:
            topo = nc.variables[topo_type['topo']][:]
            topo = np.reshape(topo, (topo_lat.size, topo_lon.size))
            topo = topo[start2:end2, start1:end1]
        else:
            topo = nc.variables[topo_type['topo']][start2:end2, start1:end1]
        for j in range(0,count2):
            for i in range(0,count1):
                htopo[j,nx_lon-imax+i-1]=topo[j,i]
        del topo

        ishft=nx_lon-count1
        if shft_east>0:
            for i in range(0,count1):
                xtmp[i+ishft]=topo_lon[i] +360
        elif shft_east<0:
            for i in range(0,count1):
                xtmp[i+ishft]=topo_lon[i] -360
        else:
            for i in range(0,count1):
                xtmp[i+ishft]=topo_lon[i]

        print('second...')
        start1=imin ; count1=period-imin; end1=start1+count1
        if gebco:
            topo = nc.variables[topo_type['topo']][:]
            topo = np.reshape(topo, (topo_lat.size, topo_lon.size))
            topo = topo[start2:end2, start1:end1]
        else:
            topo = nc.variables[topo_type['topo']][start2:end2, start1:end1]
        nc.close()

        for j in range(0,count2):
            for i in range(0,count1):
                htopo[j,i]=topo[j,i]
        del topo
        ishft=imin-1
        if shft_west>0:
            for i in range(0,count1):
                xtmp[i]=topo_lon[i+ishft] +360
        elif shft_west<0 :
            for i in range(0,count1):
                xtmp[i]=topo_lon[i+ishft] -360
        else:
            for i in range(0,count1):
                xtmp[i]=topo_lon[i+ishft]
        lon_tmp=np.zeros([xtmp.shape[0]])
        for i in range(0,nx_lon):
            lon_tmp[i]=xtmp[i]

        topo=np.copy(htopo)

    del topo_lon,topo_lat
    topo_lon=np.copy(lon_tmp)
    topo_lat=np.copy(lat_tmp)

    return topo_lon,topo_lat,topo


##################################
def savegrid(inputs,outputs):
    """
    Create and save a new CROCO grid file
    """
    nc = netcdf.Dataset('./croco_grd.nc', 'w', format='NETCDF4')

    # create global variables
    nc.created = datetime.now().isoformat()
    nc.type = 'ROMS grid file produced by easygrid_python.py'
    nc.VertCoordType = 'NEW';
    nc.nx = np.int32(inputs.nx)
    nc.ny = np.int32(inputs.ny)
    nc.size_x = inputs.size_x
    nc.size_y = inputs.size_y
    nc.lonmin = inputs.lonmin
    nc.lonmax = inputs.lonmax
    nc.latmin = inputs.latmin
    nc.latmax = inputs.latmax
    nc.rotation = inputs.rot

    # create dimensions
    nc.createDimension('xi_rho', inputs.nx)
    nc.createDimension('eta_rho', inputs.ny)
    nc.createDimension('xi_u', inputs.nx - 1)
    nc.createDimension('eta_v', inputs.ny - 1)
    nc.createDimension('xi_psi', inputs.nx - 1)
    nc.createDimension('eta_psi', inputs.ny - 1)
    nc.createDimension('one', 1)

    # create variables and attributes
    nc.createVariable('spherical', 'S1', ('one'))
    nc.variables['spherical'].long_name = 'Grid type logical switch'
    nc.variables['spherical'].option_T = 'spherical'
    nc.variables['spherical'][:] = 'T'

    nc.createVariable('xl', 'f8', ('one'))
    nc.variables['xl'].long_name = 'domain length in the XI-direction'
    nc.variables['xl'].units = 'meters'
    nc.variables['xl'][:] = inputs.nx

    nc.createVariable('el', 'f8', ('one'))
    nc.variables['el'].long_name = 'domain length in the ETA-direction'
    nc.variables['el'].units = 'meters'
    nc.variables['el'][:] = inputs.ny

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
    nc.variables['lon_psi'][:] = outputs.lon_psi

    nc.createVariable('lat_psi', 'f8', ('eta_psi', 'xi_psi'))
    nc.variables['lat_psi'].long_name = 'latitude of PSI-points'
    nc.variables['lat_psi'].units = 'degree_north'
    nc.variables['lat_psi'][:] = outputs.lat_psi

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

    nc.close()

#############################
### Class for normal mode ###
class inputs():
    '''
    Inputs to locate grid
    '''
    def __init__(self,tra_lon,tra_lat,size_x,size_y,nx,ny,rot):

        self.tra_lon = tra_lon
        self.tra_lat = tra_lat
        self.size_x  = size_x
        self.size_y  = size_y
        self.rot     = rot
        self.nx      = nx
        self.ny      = ny

class inputs_smth():
    '''
    Inputs for smoothing
    '''
    def __init__(self,depthmin,depthmax,smthr,rfact,smooth):
        self.depthmin  = depthmin
        self.depthmax  = depthmax
        self.smthr     = smthr
        self.rfact     = rfact
        self.smooth    = smooth
#######################
### Read parent grid ##

class topo_prt():

    def __init__(self,prt_file):
       nc=netcdf.Dataset(prt_file,'r')
       self.lon_rho  = nc.variables['lon_rho'][:]
       self.lat_rho  = nc.variables['lat_rho'][:]
       self.lon_psi  = nc.variables['lon_psi'][:]
       self.lat_psi  = nc.variables['lat_psi'][:]
       self.lon_u    = nc.variables['lon_u'][:]
       self.lat_u    = nc.variables['lat_u'][:]
       self.lon_v    = nc.variables['lon_v'][:]
       self.lat_v    = nc.variables['lat_v'][:]
       self.mask_rho = nc.variables['mask_rho'][:]
       self.h        = nc.variables['h'][:]
       nc.close()
