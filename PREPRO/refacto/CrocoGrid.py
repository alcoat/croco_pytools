import numpy as np
import sigmagrid_tools as sig_tools
import netcdf_tools as nc_tools
import netCDF4 as netcdf

class CrocoGrid:

    def __init__(self, filename, sigma_params=None):
        print('Reading CROCO grid: %s' %filename )
        self.grid_file = filename
        nc=netcdf.Dataset(filename,'r')
        self.lon = nc_tools.read_nc(filename,'lon_rho')
        self.lat = nc_tools.read_nc(filename,'lat_rho')
        self.lonu =nc_tools.read_nc(filename,'lon_u')
        self.latu = nc_tools.read_nc(filename,'lat_u')
        self.lonv = nc_tools.read_nc(filename,'lon_v')
        self.latv = nc_tools.read_nc(filename,'lat_v')
        self.pm  = nc_tools.read_nc(filename,'pm')
        self.pn  = nc_tools.read_nc(filename,'pn')
        self.maskr = nc_tools.read_nc(filename,'mask_rho')
        self.angle = nc_tools.read_nc(filename,'angle')
        self.h = nc_tools.read_nc(filename,'h')
        self.hraw = nc_tools.read_nc(filename,'hraw')
        self.f = nc_tools.read_nc(filename,'f')
        self.umask= self.maskr[:,0:-1]*self.maskr[:,1:]
        self.vmask= self.maskr[0:-1,:]*self.maskr[1:,:]
        self.pmask= self.umask[0:-1,:]*self.umask[1:,:]
        if sigma_params is not None:
            self.theta_s = np.double(sigma_params['theta_s'])
            self.theta_b = np.double(sigma_params['theta_b'])
            self.hc = np.double(sigma_params['hc'])
            self.N = np.double(sigma_params['N'])
        self.sc_r = None
        nc.close
    def mask3d(self):
        return np.tile(self.maskr, (np.int(self.N), 1, 1))

    def umask3d(self):
        return np.tile(self.umask, (np.int(self.N), 1, 1))

    def vmask3d(self):
        return np.tile(self.vmask, (np.int(self.N), 1, 1))

    def lonmin(self):
        return np.min(self.lon)

    def lonmax(self):
        return np.max(self.lon)

    def latmin(self):
        return np.min(self.lat)

    def latmax(self):
        return np.max(self.lat)

    def scoord2z_r(self, zeta=0., bdy="", scoord='new2008'):
        '''
        Depths at vertical rho points
        '''
        return sig_tools.scoord2z('r', zeta=zeta, topo=eval(''.join(("self.h",bdy))), theta_s=self.theta_s, theta_b=self.theta_b,\
                N=self.N,hc=self.hc,scoord=scoord)[0]


    def Cs_r(self, zeta=0., bdy="", scoord='new2008'):
        '''
        S-coordinate stretching curves at rho points
        '''
        return sig_tools.scoord2z('r', zeta=zeta, topo=eval(''.join(("self.h",bdy))), theta_s=self.theta_s, theta_b=self.theta_b,\
                N=self.N,hc=self.hc,scoord=scoord)[1]


    def scoord2z_w(self, zeta=0., bdy="", scoord='new2008'):
        '''
        Depths at vertical w points
        '''
        return sig_tools.scoord2z('w', zeta=zeta, topo=eval(''.join(("self.h",bdy))), theta_s=self.theta_s, theta_b=self.theta_b,\
                N=self.N,hc=self.hc,scoord=scoord)[0]

    def Cs_w(self, zeta=0., bdy="", scoord='new2008'):
        '''
        S-coordinate stretching curves at w points
        '''
        return sig_tools.scoord2z('w', zeta=zeta, topo=eval(''.join(("self.h",bdy))), theta_s=self.theta_s, theta_b=self.theta_b,\
                N=self.N,hc=self.hc,scoord=scoord)[1]

    def WEST_grid(self,indices="[:,0:2]"):
        '''
        Defines some vars for Western grid
        '''
        self.h_west=eval(''.join(('self.h',indices)))
        self.lon_west=eval(''.join(('self.lon',indices)))
        self.lat_west=eval(''.join(('self.lat',indices)))
        self.lonu_west=eval(''.join(('self.lonu',indices)))
        self.latu_west=eval(''.join(('self.latu',indices)))
        self.lonv_west=eval(''.join(('self.lonv',indices)))
        self.latv_west=eval(''.join(('self.latv',indices)))
        self.maskr_west=eval(''.join(('self.maskr',indices)))
        self.umask_west=eval(''.join(('self.umask',indices)))
        self.vmask_west=eval(''.join(('self.vmask',indices)))
        self.angle_west=eval(''.join(('self.angle',indices)))

        return self

    def EAST_grid(self,indices="[:,-2:]"):
        '''
        Defines some vars for Western grid
        '''
        self.h_east=eval(''.join(('self.h',indices)))
        self.lon_east=eval(''.join(('self.lon',indices)))
        self.lat_east=eval(''.join(('self.lat',indices)))
        self.lonu_east=eval(''.join(('self.lonu',indices)))
        self.latu_east=eval(''.join(('self.latu',indices)))
        self.lonv_east=eval(''.join(('self.lonv',indices)))
        self.latv_east=eval(''.join(('self.latv',indices)))
        self.maskr_east=eval(''.join(('self.maskr',indices)))
        self.umask_east=eval(''.join(('self.umask',indices)))
        self.vmask_east=eval(''.join(('self.vmask',indices)))
        self.angle_east=eval(''.join(('self.angle',indices)))

        return self

    def SOUTH_grid(self,indices="[0:2,:]"):
        '''
        Defines some vars for Western grid
        '''
        self.h_south=eval(''.join(('self.h',indices)))
        self.lon_south=eval(''.join(('self.lon',indices)))
        self.lat_south=eval(''.join(('self.lat',indices)))
        self.lonu_south=eval(''.join(('self.lonu',indices)))
        self.latu_south=eval(''.join(('self.latu',indices)))
        self.lonv_south=eval(''.join(('self.lonv',indices)))
        self.latv_south=eval(''.join(('self.latv',indices)))
        self.maskr_south=eval(''.join(('self.maskr',indices)))
        self.umask_south=eval(''.join(('self.umask',indices)))
        self.vmask_south=eval(''.join(('self.vmask',indices)))
        self.angle_south=eval(''.join(('self.angle',indices)))

        return self

    def NORTH_grid(self,indices="[-2:,:]"):
        '''
        Defines some vars for Western grid
        '''
        self.h_north=eval(''.join(('self.h',indices)))
        self.lon_north=eval(''.join(('self.lon',indices)))
        self.lat_north=eval(''.join(('self.lat',indices)))
        self.lonu_north=eval(''.join(('self.lonu',indices)))
        self.latu_north=eval(''.join(('self.latu',indices)))
        self.lonv_north=eval(''.join(('self.lonv',indices)))
        self.latv_north=eval(''.join(('self.latv',indices)))
        self.maskr_north=eval(''.join(('self.maskr',indices)))
        self.umask_north=eval(''.join(('self.umask',indices)))
        self.vmask_north=eval(''.join(('self.vmask',indices)))
        self.angle_north=eval(''.join(('self.angle',indices)))

        return self


