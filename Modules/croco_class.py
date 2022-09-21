import numpy as np
import Cgrid_transformation_tools as grd_tools
import sigmagrid_tools as sig_tools
import netcdf_tools as nc_tools
import netCDF4 as netcdf
from datetime import datetime
from collections import OrderedDict

class CROCO_grd(object):

    def __init__(self, filename, sigma_params):
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
        self.umask= grd_tools.rho2u(self.maskr)
        self.vmask= grd_tools.rho2v(self.maskr)
        self.pmask= grd_tools.rho2psi(self.maskr)
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

class CROCO(CROCO_grd):

    def create_ini_nc(self, filename, grdobj, created_by='make_ini.py'):#fillval
        # Global attributes
        nc = netcdf.Dataset(filename, 'w', format='NETCDF4')
        nc.created = datetime.now().isoformat()
        nc.type  = 'CROCO initial file produced by %s' % created_by
        nc.grd_file = grdobj.grid_file
        nc.hc = grdobj.hc
        nc.theta_s = grdobj.theta_s
        nc.theta_b = grdobj.theta_b
        nc.Tcline = grdobj.hc
        nc.Cs_r = grdobj.Cs_r()
        nc.Cs_w = grdobj.Cs_w()
        nc.VertCoordType = grdobj.scoord

        # Dimensions
        nc.createDimension('xi_rho', grdobj.lon.shape[1])
        nc.createDimension('xi_u', grdobj.lon.shape[1] - 1)
        nc.createDimension('eta_rho', grdobj.lon.shape[0])
        nc.createDimension('eta_v', grdobj.lon.shape[0] - 1)
        nc.createDimension('s_rho', grdobj.N)
        nc.createDimension('s_w', grdobj.N + 1)
        nc.createDimension('time', None)
        nc.createDimension('one', 1)

        # Create the variables and write...
        nc.createVariable('theta_s', 'f', ('one'), zlib=True)
        nc.variables['theta_s'].long_name = 'S-coordinate surface control parameter'
        nc.variables['theta_s'].units = 'nondimensional'
        nc.variables['theta_s'][:] = grdobj.theta_s

        nc.createVariable('theta_b', 'f', ('one'), zlib=True)
        nc.variables['theta_b'].long_name = 'S-coordinate bottom control parameter'
        nc.variables['theta_b'].units = 'nondimensional'
        nc.variables['theta_b'][:] = grdobj.theta_b

        nc.createVariable('Tcline', 'f', ('one'), zlib=True)
        nc.variables['Tcline'].long_name = 'S-coordinate surface/bottom layer width'
        nc.variables['Tcline'].units = 'meters'
        nc.variables['Tcline'][:] = grdobj.hc

        nc.createVariable('hc', 'f', ('one'), zlib=True)
        nc.variables['hc'].long_name = 'S-coordinate parameter, critical depth'
        nc.variables['hc'].units = 'meters'
        nc.variables['hc'][:] = grdobj.hc

        nc.createVariable('sc_r', 'f8', ('s_rho'))
        nc.variables['sc_r'].long_name = 'S-coordinate at RHO-points'
        nc.variables['sc_r'].units = 'nondimensional'
        nc.variables['sc_r'].valid_min = -1.
        nc.variables['sc_r'].valid_max = 0.

        nc.createVariable('Cs_r', 'f8', ('s_rho'), zlib=True)
        nc.variables['Cs_r'].long_name = 'S-coordinate stretching curves at RHO-points'
        nc.variables['Cs_r'].units = 'nondimensional'
        nc.variables['Cs_r'].valid_min = -1.
        nc.variables['Cs_r'].valid_max = 0.
        nc.variables['Cs_r'][:] = grdobj.Cs_r()

        nc.createVariable('Cs_w', 'f8', ('s_w'), zlib=True)
        nc.variables['Cs_w'].long_name = 'S-coordinate stretching curves at w-points'
        nc.variables['Cs_w'].units = 'nondimensional'
        nc.variables['Cs_w'].valid_min = -1.
        nc.variables['Cs_w'].valid_max = 0.
        nc.variables['Cs_w'][:] = grdobj.Cs_w()

        nc.createVariable('ocean_time', 'f8', ('time'), zlib=True)
        nc.variables['ocean_time'].long_name = 'time since initialization'
        nc.variables['ocean_time'].units     = 'seconds'

        nc.createVariable('scrum_time', 'f8', ('time'), zlib=True)
        nc.variables['scrum_time'].long_name = 'time since initialization'
        nc.variables['scrum_time'].units     = 'seconds'

        nc.createVariable('tstart', 'f8', ('one'), zlib=True)
        nc.variables['tstart'].long_name = 'start processing day'
        nc.variables['tstart'].units     = 'days'

        nc.createVariable('tend', 'f8', ('one'), zlib=True)
        nc.variables['tend'].long_name = 'end processing day'
        nc.variables['tend'].units     = 'days'

        # dictionary for the prognostic variables
        prog_vars = OrderedDict()
        prog_vars['temp'] = ['rho3d',
                             'initial potential temperature',
                             'Celsius']
        prog_vars['salt'] = ['rho3d',
                             'initial salinity',
                             'psu']
        prog_vars['u']    = ['u3d',
                             'initial u-momentum component',
                             'meters second-1']
        prog_vars['v']    = ['v3d',
                             'initial v-momentum component',
                             'meters second-1']
        prog_vars['ubar'] = ['u2d',
                             'initial vertically integrated u-momentum component',
                             'meters second-1']
        prog_vars['vbar'] = ['v2d',
                             'initial vertically integrated v-momentum component',
                             'meters second-1']
        prog_vars['zeta'] = ['rho2d',
                             'initial sea surface height',
                             'meters']

        for varname, value in zip(prog_vars.keys(), prog_vars.values()):

            if 'rho3d' in value[0]:
                dims = ('time', 's_rho', 'eta_rho', 'xi_rho')

            elif 'u3d' in value[0]:
                dims = ('time', 's_rho', 'eta_rho', 'xi_u')

            elif 'v3d' in value[0]:
                dims = ('time', 's_rho', 'eta_v', 'xi_rho')

            elif 'u2d' in value[0]:
                dims = ('time', 'eta_rho', 'xi_u')

            elif 'v2d' in value[0]:
                dims = ('time', 'eta_v', 'xi_rho')

            elif 'rho2d' in value[0]:
                dims = ('time', 'eta_rho', 'xi_rho')

            else: error

            nc.createVariable(varname, 'f8', dims, zlib=True)#fill_value=fillval
            nc.variables[varname].long_name = value[1]
            nc.variables[varname].units     = value[2]

        nc.close()

    def create_bry_nc(self,filename, grdobj, obc_dict, cycle, created_by='make_bdy.py'):#fillval
        # Global attributes
        nc = netcdf.Dataset(filename, 'w', format='NETCDF4')
        nc.created = datetime.now().isoformat()
        nc.type = 'ROMS boundary file produced by %s' %created_by
        nc.grd_file = grdobj.grid_file
        nc.hc = grdobj.hc
        nc.theta_s = grdobj.theta_s
        nc.theta_b = grdobj.theta_b
        nc.Tcline = grdobj.hc
        nc.Cs_r = grdobj.Cs_r()
        nc.Cs_w = grdobj.Cs_w()
        nc.VertCoordType = 'NEW'

        # Dimensions
        nc.createDimension('xi_rho', grdobj.lon.shape[1])
        nc.createDimension('xi_u', grdobj.lon.shape[1]-1)
        nc.createDimension('eta_rho', grdobj.lon.shape[0])
        nc.createDimension('eta_v', grdobj.lon.shape[0]-1)
        nc.createDimension('s_rho', grdobj.N)
        nc.createDimension('s_w', grdobj.N+1)
        nc.createDimension('bry_time', None)
        nc.createDimension('one', 1)

        # Create the variables and write...
        nc.createVariable('theta_s', 'f', ('one'))
        nc.variables['theta_s'].long_name = 'S-coordinate surface control parameter'
        nc.variables['theta_s'].units = 'nondimensional'
        nc.variables['theta_s'][:] = grdobj.theta_s

        nc.createVariable('theta_b', 'f', ('one'))
        nc.variables['theta_b'].long_name = 'S-coordinate bottom control parameter'
        nc.variables['theta_b'].units = 'nondimensional'
        nc.variables['theta_b'][:] = grdobj.theta_b

        nc.createVariable('Tcline', 'f', ('one'))
        nc.variables['Tcline'].long_name = 'S-coordinate surface/bottom layer width'
        nc.variables['Tcline'].units = 'meters'
        nc.variables['Tcline'][:] = grdobj.hc

        nc.createVariable('hc', 'f', ('one'))
        nc.variables['hc'].long_name = 'S-coordinate parameter, critical depth'
        nc.variables['hc'].units = 'meters'
        nc.variables['hc'][:] = grdobj.hc

        nc.createVariable('sc_r', 'f8', ('s_rho'))
        nc.variables['sc_r'].long_name = 'S-coordinate at RHO-points'
        nc.variables['sc_r'].units = 'nondimensional'
        nc.variables['sc_r'].valid_min = -1.
        nc.variables['sc_r'].valid_max = 0.
        nc.variables['sc_r'][:] = grdobj.sc_r

        nc.createVariable('Cs_r', 'f8', ('s_rho'))
        nc.variables['Cs_r'].long_name = 'S-coordinate stretching curves at RHO-points'
        nc.variables['Cs_r'].units = 'nondimensional'
        nc.variables['Cs_r'].valid_min = -1.
        nc.variables['Cs_r'].valid_max = 0.
        nc.variables['Cs_r'][:] = grdobj.Cs_r()

        nc.createVariable('Cs_w', 'f8', ('s_w'))
        nc.variables['Cs_w'].long_name = 'S-coordinate stretching curves at w-points'
        nc.variables['Cs_w'].units = 'nondimensional'
        nc.variables['Cs_w'].valid_min = -1.
        nc.variables['Cs_w'].valid_max = 0.
        nc.variables['Cs_w'][:] = grdobj.Cs_w()

        nc.createVariable('bry_time', 'f8', ('bry_time'), zlib=True)
        nc.variables['bry_time'].long_name = 'time for boundary data'
        nc.variables['bry_time'].units = 'days'

        # dictionary for the prognostic variables
        prog_vars = OrderedDict()
        prog_vars['temp_'] = ['rho2',
                              ' boundary potential temperature',
                              'Celsius']
        prog_vars['salt_'] = ['rho2',
                              ' boundary salinity',
                              'psu']
        prog_vars['u_']    = ['u2',
                              ' boundary u-momentum component',
                              'meters second-1']
        prog_vars['v_']    = ['v2',
                              ' boundary v-momentum component',
                              'meters second-1']
        prog_vars['ubar_'] = ['u1',
                              ' boundary vertically integrated u-momentum component',
                              'meters second-1']
        prog_vars['vbar_'] = ['v1',
                              ' boundary vertically integrated v-momentum component',
                              'meters second-1']
        prog_vars['zeta_'] = ['rho1',
                              ' boundary sea surface height',
                              'meters']
        # Loop over boundary
        for boundary, flag in zip(obc_dict.keys(), obc_dict.values()):
            if flag:
                varlabel = '%sern'   % boundary
                for key, value in zip(prog_vars.keys(), prog_vars.values()):
                    if 'rho2' in value[0]:
                        if boundary=='east' or boundary=='west':
                            dims = ('bry_time', 's_rho', 'eta_rho')
                        elif boundary=='south' or boundary=='north':
                            dims = ('bry_time', 's_rho', 'xi_rho')
                    elif 'u2' in value[0]:
                        if boundary=='south' or boundary=='north':
                            dims = ('bry_time', 's_rho', 'xi_u')
                        else:
                            dims = ('bry_time', 's_rho', 'eta_rho')
                    elif 'v2' in value[0]:
                        if boundary=='east' or boundary=='west':
                            dims = ('bry_time', 's_rho', 'eta_v')
                        else:
                            dims = ('bry_time', 's_rho', 'xi_rho')
                    elif 'u1' in value[0]:
                        if boundary=='south' or boundary=='north':
                            dims = ('bry_time', 'xi_u')
                        else:
                            dims = ('bry_time', 'eta_rho')
                    elif 'v1' in value[0]:
                        if boundary=='east' or boundary=='west':
                            dims = ('bry_time', 'eta_v')
                        else:
                            dims = ('bry_time', 'xi_rho')
                    elif 'rho1' in value[0]:
                        if boundary=='east' or boundary=='west':
                            dims = ('bry_time', 'eta_rho')
                        elif boundary=='south' or boundary=='north':
                            dims = ('bry_time', 'xi_rho')
                    else: error

                    varname  = ''.join((key, '%s' % boundary))
                    nc.createVariable(varname, 'f8', dims, zlib=True)
                    nc.variables[varname].long_name = varlabel + value[1]
                    nc.variables[varname].units     = value[2]

        nc.close()







