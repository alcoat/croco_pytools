__author__ = 'Mathieu Le Corre'
__email__  = 'mathieu.le.corre@shom.fr'
__date__   = '2022-09'
__license__='GPL3'
'''
===========================================================================
Further Information:  
  http://www.croco-ocean.org
  
LThis file is part of CROCO_PYTOOLS

Create a CROCO initial file
In the current state the script can handle:
    - mercator (glorys)
    - soda
    - eccov4  (V4.3)

To add a new dataset you just have to go in Readers/inputs_readers.py and
create a dico with correspondance between croco_var(ssh,u,v,temp,salt)
and dataset variable names.

The periodicity of a dataset (if it is -180/180 or 0/360) is handle by the
script. Just download the dataset and it should be fine

The script works as follow:
    - reads croco_grd
    - reads input data and restricts the area of spatial coverage
    - creates croco_ini.nc
    - Loop on var with:
        * horizontal interpolation
        * vertical interpolation
    - Writes data in netcdf
===========================================================================
'''
#--- Dependencies ---------------------------------------------------------
import namelist_ini as params
from CrocoGrid import CrocoGrid
import inputs_reader
from InputData import InputData
import interp_tools
import sigmagrid_tools as sig_tools
import Cgrid_transformation_tools as grd_tools
from croco_files_def import create_ini_nc
import pylab as plt
import numpy as np
from datetime import datetime
import netCDF4 as netcdf
#--------------------------------------------------------------------------

if __name__ == '__main__':

    # --- Read croco grid --------------------------------------------------

    print('Read croco grid')
    crocogrd = CrocoGrid(filename=''.join((params.croco_dir, params.croco_grd)), sigma_params=params.sigma_params)

    # --- Get input data --------------------------------------------------

    # edit ini_filename to add starting date
    ini_filename = params.ini_filename.replace('.nc', '_%s_Y%sM%s.nc' %(params.datatype, params.Yini, params.Mini))

    # read variables names through reader dictionnary
    print('Read data from ' + params.datatype)
    dicovar = inputs_reader.lookvar(params.datatype)

    tracers = ['temp', 'salt']
    vectors = [('u', 'v')]

    # First get SSH
    for var in ['ssh']:
        print('Start with ssh: ')
        if params.multi_files == False: input_file = params.input_file
        else: input_file = params.input_file[var]

        myvar = InputData(var=dicovar[var], input_file=input_file, grid='r')
        print('----- Get ' + var + ' -----')
        myvar.get_data()

        print('----- Interpolate ' + var + ' -----')
        ssh, Ngood = interp_tools.interp_horiz(myvar.data, crocogrd)

        # compute z of sigma coords for vertical interpolations
        print('----- Compute depth of sigma coordinates (will be used for vertical interpolations) -----')
        z_rho = crocogrd.scoord2z_r(zeta = ssh)
        z_w = crocogrd.scoord2z_w(zeta = ssh)

    # Then tracers
    for var in tracers:
        # TRACERS: 
        #    - horiz interp for each z level
        #    - concatenate all levels
        #    - add 2 levels (above and under) to avoid vertical extrap.
        #    - z to sigma interp of the 3D field

        print('Then get tracers: ')
        if params.multi_files == False: input_file = params.input_file
        else: input_file = params.input_file[var]

        myvar = InputData(var=dicovar[var], input_file=input_file, grid='r')
        print('----- Get ' + var + ' -----')
        myvar.get_data()

        print('----- Interpolate ' + var + ' on each z level -----')
        vertvar, Z = interp_tools.interp_horiz_var3D(myvar, crocogrd)
        print('----- Interpolate ' + var + ' on sigma levels -----')
        myinterpvar = sig_tools.ztosigma(vertvar, Z, z_rho)

     # Then vectors
     for var in vectors:
        # VECTORS:
        #    - horiz interp for each z level
        #    - concatenate all levels
        #    - add 2 levels (above and under) to avoid vertical extrap.
        #    - move to u, v grids
        #    - z to sigma interp on u, v points
        #    - compute ubar, vbar (on u, v points)

        print('Then get vectors: ')
        if params.multi_files == False: input_file = params.input_file
        else: input_file = params.input_file[var]

        print('----- Get ' + var + ' -----')
        uvec = InputData(var=dicovar[var[0]], input_file=input_file, grid='u')
        vvec = InputData(var=dicovar[var[1]], input_file=input_file, grid='v')
        uvec.get_data()
        vvec.get_data()

        print('----- Interpolate ' + var + ' on each z level -----')
        uinterp, Z = interp_tools.interp_horiz_var3D(uvec, crocogrd)
        vinterp, Z = interp_tools.interp_horiz_var3D(vvec, crocogrd)


        # Rotation and put to u-points and v-points 
        print('----- Eventually rotate according to grid angle, and put on u-v-points  -----')
        cosa3d = np.tile(np.cos(crocogrd.angle), (np.shape(uinterp)[0], 1, 1))
        sina3d = np.tile(np.sin(crocogrd.angle), (np.shape(uinterp)[0], 1, 1))
        uz = grd_tools.rho2u(uinterp*cosa3d + vinterp*sina3d)
        vv = grd_tools.rho2v(vinterp*cosa3d - uinterp*sina3d)

        # z to sigma interp
        print('----- Interpolate ' + var + ' on sigma levels -----')
        usigma = sig_tools.ztosigma(uz, Z, grd_tools.rho2u(z_rho))
        vsigma = sig_tools.ztosigma(vz, Z, grd_tools.rho2v(z_rho))

        # compute ubar, vbar
        print('----- Compute ' + var + ' integrated value over the water column -----')
        dz = np.gradient(zw)
        ubar = np.nansum(usigma*grd_tools.rho2u(dz)) / np.nansum(grd_tools.rho2u(dz))
        vbar = np.nansum(vsigma*grd_tools.rho2v(dz)) / np.nansum(grd_tools.rho2v(dz))

        # Correct the horizontal transport i.e. remove the intergrated tranport and add the OGCM transport
        u = usigma + np.tile(ubar, (z_rho.shape[0],1,1))
        v = vsigma + np.tile(vbar, (z_rho.shape[0],1,1))


    # --- Create croco_ini file -----------------------------------------

    print(' ')
    print(' Making initial file: ' + ini_filename)
    print(' ')

    create_ini_nc(None, ''.join((params.croco_dir + ini_filename)), crocogrd)

    # Handle initial time
    ini_date_num = datetime(int(params.Yini), int(params.Mini), int(params.Dini))
    ini_date_num = plt.date2num(ini_date_num) + 0.5
    day_zero_num = datetime(int(params.Yzer), int(params.Mzer), int(params.Dzer))
    day_zero_num = plt.date2num(day_zero_num)
    tstart=0
    if ini_date_num != day_zero_num:
        tstart = ini_date_num - day_zero_num # days
    scrumt = tstart*3600*24 # convert in second
    oceant = tstart*3600*24
    tend=0.

    nc=netcdf.Dataset(params.croco_dir+ini_filename, 'a')
    nc.variables['ocean_time'][:] = oceant
    nc.variables['scrum_time'][:] = scrumt
    nc.variables['scrum_time'].units='seconds since %s-%s-%s 00:00:00' %(params.Yzer,params.Mzer,params.Dzer)
    nc.variables['tstart'][:] = tstart
    nc.variables['tend'][:] = tend
    nc.Input_data_type = params.datatype # ?

    for vars in ['ssh','tracers','velocity']:
        print('\nWriting *%s*' %vars)
        if vars == 'ssh' :
            nc.variables['zeta'][0,:,:] = zeta * crocogrd.maskr

        elif vars == 'tracers':
            print('\nIn tracers writing Temp')
            nc.variables['temp'][0,:,:,:] = temp * crocogrd.mask3d()
            print('\nIn tracers writing Salt')
            nc.variables['salt'][0,:,:,:] = salt * crocogrd.mask3d()

        elif vars == 'velocity':
            nc.variables['u'][0,:,:,:] = u * crocogrd.umask3d()
            nc.variables['v'][0,:,:,:] = v * crocogrd.vmask3d()
            nc.variables['ubar'][0,:,:] = ubar * crocogrd.umask
            nc.variables['vbar'][0,:,:] = vbar * crocogrd.vmask

    nc.close()

