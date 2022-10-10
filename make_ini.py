__author__ = 'Mathieu Le Corre'
__email__  = 'mathieu.le.corre@shom.fr'
__date__   = '2022-09'
__license__='GPL3'
'''
===========================================================================
Further Information:  
  http://www.croco-ocean.org
  
This file is part of CROCOTOOLS

Create a CROCO initial file
In the current state the script can handle:
    - mercator (glorys)
    - soda
    - eccov4  (V4.3)

To add a new dataset you just have to go in Modules/inputs_readers.py and
create a dico with correspondance between croco_var(ssh,u,v,temp,salt)
and dataset variable names.

The periodicity of a dataset (if it is -180/180 or 0/360) is handle by the
script. Just download the dataset and it should be fine

The script works as follow:
    - reads croco_grd
    - reads input data and restricts the area of spatial coverage
    - creates croco_ini.nc
    - computes coefficients for horizontal interpolation on each grid (rho,u,v)
    - Loop on var with:
        * horizontal interpolation
        * vertical interpolation
    - Writes data in netcdf
===========================================================================
'''
#--- Dependencies ---------------------------------------------------------
import netCDF4 as netcdf
import pylab as plt
import numpy as np
import glob as glob
from datetime import datetime
import sys
sys.path.append("./Modules/")
import Cgrid_transformation_tools as grd_tools
import interp_tools
import sigmagrid_tools as sig_tools

import croco_class as Croco
import input_class as Inp
#--------------------------------------------------------------------------


#--- USER CHANGES ---------------------------------------------------------

# input informations
inputdata='mercator'   # At hte current time can handle mercator,soda,eccov4
input_dir = './'
input_prefix='raw_motu_mercator_'
multi_files=False # If variables are in different netcdf

# input files
Yini,Mini,Dini  = '2005','01','01' # Month and days need to be 2-digits format
date_str = (Yini, Mini)
input_file  = input_dir + input_prefix + 'Y%sM%s.nc' % date_str

if multi_files: # Mutiple files
    input_file = { 'ssh'  : input_dir + input_prefix + 'ETAN.%s.nc' % date_str, \
                   'temp' : input_dir + input_prefix + 'THETA.%s.nc' % date_str, \
                   'salt' : input_dir + input_prefix + 'SALT.%s.nc' % date_str, \
                   'u'    : input_dir + input_prefix + 'EVEL.%s.nc' % date_str, \
                   'v'    : input_dir + input_prefix + 'NVEL.%s.nc' % date_str\
                }
tndx=0 # time index in the file

# CROCO path and filename informations
croco_dir = './'
croco_grd = 'croco_grd.nc'
sigma_params = dict(theta_s=7, theta_b=2, N=32, hc=75) # Vertical streching, sig_surf/sig_bot/ nb level/critical depth

Yzer,Mzer,Dzer = Yini, Mini, Dini # reference time (default = ini time).Month and days need to be 2-digits format

# inifile informations
ini_filename    = 'croco_ini.nc'

# create delaunay weight
comp_delaunay=1

Nzgoodmin = 4  # default value to consider a z-level fine to be used

#--- END USER CHANGES -----------------------------------------------------

if __name__ == '__main__':
    # edit ini_filename to add starting date
    ini_filename = ini_filename.replace('.nc', '_%s_Y%sM%s.nc' %(inputdata,Yini, Mini))

    # Load croco_grd
    crocogrd = Croco.CROCO_grd(''.join((croco_dir, croco_grd)), sigma_params)

    # --- Load input (restricted to croco_grd) ----------------------------

    inpdat=Inp.getdata(inputdata,input_file,crocogrd,multi_files)

    print(' ')
    print(' Making initial file: '+ini_filename)
    print(' ')

    # --- Create the initial file -----------------------------------------

    Croco.CROCO.create_ini_nc(None,''.join((croco_dir + ini_filename)),crocogrd)

    # --- Handle initial time ---------------------------------------------
    ini_date_num = datetime(int(Yini), int(Mini), int(Dini))
    ini_date_num = plt.date2num(ini_date_num) + 0.5

    day_zero_num = datetime(int(Yzer), int(Mzer), int(Dzer))
    day_zero_num = plt.date2num(day_zero_num)
    
    tstart=0
    
    if ini_date_num != day_zero_num:
        tstart = ini_date_num - day_zero_num # days

    scrumt = tstart*3600*24 # convert in second
    oceant = tstart*3600*24
    tend=0.


   # --- Get the 2D interpolation coefficients ----------------------------

    if comp_delaunay==1:
    # Compute the Delaunay triangulation matrices (long but only done once)
    # (u and v are interpolated on croco rho_points because we may need to rotate them)

        print('\nCompute Delaunay triangulation from GLORYS T-points to CROCO rho_points...')
        print('--------------------------------------------------------------------------')

        [elemT,coefT] = interp_tools.get_tri_coef(inpdat.lonT,inpdat.latT,crocogrd.lon,crocogrd.lat)
        coefnorm=np.sum(coefT,axis=2)
        coefT=coefT/coefnorm[:,:,np.newaxis]

        print('\nCompute Delaunay triangulation from GLORYS U-points to CROCO rho_points...')
        print('--------------------------------------------------------------------------')

        [elemU,coefU] = interp_tools.get_tri_coef(inpdat.lonU,inpdat.latU,crocogrd.lon,crocogrd.lat)
        coefnorm=np.sum(coefU,axis=2)
        coefU=coefU/coefnorm[:,:,np.newaxis]

        print('\nCompute Delaunay triangulation from GLORYS V-points to CROCO rho_points...')
        print('--------------------------------------------------------------------------')

        [elemV,coefV] = interp_tools.get_tri_coef(inpdat.lonV,inpdat.latV,crocogrd.lon,crocogrd.lat)
        coefnorm=np.sum(coefV,axis=2)
        coefV=coefV/coefnorm[:,:,np.newaxis]

    # Save the Delaunay triangulation matrices
        np.savez('coeffs.npz',coefT=coefT,elemT=elemT,\
             coefU=coefU,elemU=elemU,coefV=coefV,elemV=elemV)
    else:
    # Load the Delaunay triangulation matrices
        print('Load Delaunay triangulation...')
        data=np.load('coeffs.npz')
        coefT = data['coefT']
        elemT = data['elemT']
        coefU = data['coefU']
        elemU = data['elemU']
        coefV = data['coefV']
        elemV = data['elemV']

        print('Delaunay triangulation done')

   #  --- Compute and save variables on CROCO grid ---------------

    for vars in ['ssh','tracers','velocity']:
        print('\nProcessing *%s*' %vars)
        nc=netcdf.Dataset(croco_dir+ini_filename, 'a')
        if vars == 'ssh' :
            (zeta,NzGood) = interp_tools.interp_tracers(inpdat,vars,tndx,-1,coefT,elemT)
            nc.variables['zeta'][0,:,:] = zeta*crocogrd.maskr

            nc.variables['ocean_time'][:] = oceant
            nc.variables['scrum_time'][:] = scrumt
            nc.variables['tstart'][:] = tstart
            nc.variables['tend'][:] = tend

            z_rho = crocogrd.scoord2z_r(zeta=zeta)
            z_w   = crocogrd.scoord2z_w(zeta=zeta)
            
        elif vars == 'tracers':
            print('\nIn tracers processing Temp')
            temp= interp_tools.interp3d(inpdat,'temp',tndx,Nzgoodmin,z_rho,coefT,elemT)
            nc.variables['temp'][0,:,:,:] = temp*crocogrd.mask3d()
            print('\nIn tracers processing Salt')

            salt= interp_tools.interp3d(inpdat,'salt',tndx,Nzgoodmin,z_rho,coefT,elemT)
            nc.variables['salt'][0,:,:,:] = salt*crocogrd.mask3d()

        elif vars == 'velocity':

            cosa=np.cos(crocogrd.angle)
            sina=np.sin(crocogrd.angle)

            [u,v,ubar,vbar]=interp_tools.interp3d_uv(inpdat,tndx,Nzgoodmin,z_rho,cosa,sina,\
                                   coefU,elemU,coefV,elemV)
              
            conserv=1  # Correct the horizontal transport i.e. remove the intergrated tranport and add the OGCM transport          
            if conserv == 1:
                (ubar_croco,h0)=sig_tools.vintegr(u,grd_tools.rho2u(z_w),grd_tools.rho2u(z_rho),np.nan,np.nan)/grd_tools.rho2u(crocogrd.h)
                (vbar_croco,h0)=sig_tools.vintegr(v,grd_tools.rho2v(z_w),grd_tools.rho2v(z_rho),np.nan,np.nan)/grd_tools.rho2v(crocogrd.h)

                u = u - ubar_croco ; u = u + np.tile(ubar,(z_rho.shape[0],1,1))
                v = v - vbar_croco ; v = v + np.tile(vbar,(z_rho.shape[0],1,1))
           
            nc.variables['u'][0,:,:,:] = u *crocogrd.umask3d()
            nc.variables['v'][0,:,:,:] = v * crocogrd.vmask3d()
            nc.variables['ubar'][0,:,:] = ubar *crocogrd.umask
            nc.variables['vbar'][0,:,:] = vbar * crocogrd.vmask

   
    nc.close()



