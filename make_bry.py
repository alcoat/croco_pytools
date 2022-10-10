__author__ = 'Mathieu Le Corre'
__email__  = 'mathieu.le.corre@shom.fr'
__date__   = '2022-09'
__license__='GPL3'

'''
===========================================================================
Further Information:  
  http://www.croco-ocean.org
  
This file is part of CROCOTOOLS

Create a CROCO bounday files
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
    - creates croco_bry.nc
    - computes coefficients for horizontal interpolation on each grid (rho,u,v)
      and for each open boudary
    - Loop on open boundaries with:
          check if data for before or after for continuity, if not duplicate first or last
          loop on var with:
              * horizontal interpolation
              * vertical interpolation
    - Writes data in netcdf

===========================================================================
'''

#--- Dependencies ---------------------------------------------------------
import netCDF4 as netcdf
import xarray as xr
import pylab as plt
import numpy as np
import glob as glob
from dateutil.relativedelta import relativedelta
import sys
sys.path.append("./Modules/")
import interp_tools
import sigmagrid_tools as sig_tools
import Cgrid_transformation_tools as grd_tools
import croco_class as Croco
import input_class as Inp
#--------------------------------------------------------------------------


#--- USER CHANGES ---------------------------------------------------------

# input informations
inputdata='mercator'   # At hte current time can handle mercator,soda,eccov4
input_dir = '/local/tmp/3/'
input_prefix='raw_motu_mercator_*'# Please use * to include all files

multi_files=False
if multi_files: # Multiple data files. Time is read in ssh file
    input_file = { 'ssh'  : sorted(glob.glob(input_dir + input_prefix + 'ETAN.%s.nc' % date_str)), \
                   'temp' : sorted(glob.glob(input_dir + input_prefix + 'THETA.%s.nc' % date_str)), \
                   'salt' : sorted(glob.glob(input_dir + input_prefix + 'SALT.%s.nc' % date_str)), \
                   'u'    : sorted(glob.glob(input_dir + input_prefix + 'EVEL.%s.nc' % date_str)), \
                   'v'    : sorted(glob.glob(input_dir + input_prefix + 'NVEL.%s.nc' % date_str))\
                }
else:  # glob all files
    input_file  = sorted(glob.glob(input_dir + input_prefix))


# CROCO path and filename informations
croco_dir = './' 
croco_grd = 'croco_grd.nc'
sigma_params = dict(theta_s=7, theta_b=2, N=32, hc=75) # Vertical streching, sig_surf/sig_bot/ nb level/critical depth

# bryfile informations
bry_filename    = 'croco_bry.nc'
obc_dict = dict(south=1, west=1, east=1, north=1) # open boundaries (1=open , [S W E N])
output_file_format="MONTHLY" # How outputs are spit (MONTHLY,YEARLY,FULL)
cycle_bry=0.

Yorig=2005 # year origin of time : days since Yorig-01-01
Ystart,Mstart = '2005', '01'   # Starting month
Yend,Mend  = '2005','02'       # Ending month 

#  create delaunay weight 
comp_delaunay=1

Nzgoodmin=4 # default value to consider a z-level fine to be used
#--- END USER CHANGES -----------------------------------------------------

#--- START MAIN SCRIPT ----------------------------------------------------


if __name__ == '__main__':

    # Put origin date to the right format
    day_zero   = str(Yorig)+'0101'    
    day_zero_num = plt.datetime.datetime(int(day_zero[:4]),
                                         int(day_zero[4:6]),
                                         int(day_zero[6:8]))
    day_zero_num = plt.date2num(day_zero_num)

    # Put start and end date to the right format
    start_date = Ystart+Mstart+'01'+'12'  # defaut start day is 1st
   

    dtstrdt = plt.datetime.datetime(int(start_date[:4]),
                                    int(start_date[4:6]),
                                    int(start_date[6:8]),
                                    int(start_date[8:]))

    dtenddt = plt.datetime.datetime(int(Yend),int(Mend),1,12) \
            + relativedelta(months=1,days=-1) # Last day of the ending month


    dtstr, dtend = plt.date2num(dtstrdt), plt.date2num(dtenddt)

    # --- Load croco_grd --------------------------------------------------

    crocogrd = Croco.CROCO_grd(''.join((croco_dir, croco_grd)), sigma_params)

    # --- Initialize boundary vars ----------------------------------------

    crocogrd.WEST_grid()
    crocogrd.EAST_grid()
    crocogrd.SOUTH_grid()
    crocogrd.NORTH_grid()
    
    # --- Initialize input data class -------------------------------------

    inpdat = Inp.getdata(inputdata,input_file,crocogrd,multi_files,bdy=[obc_dict,cycle_bry])

    # --- Get the 2D interpolation coefficients ---------------------------

    if comp_delaunay==1:
        for boundary, is_open in zip(obc_dict.keys(), obc_dict.values()):
            if is_open:
                print('\n--- Processing %sern boundary' % boundary)
                print('---------------------------------------')

            if 'west' in boundary and is_open:
                (elemT_west,coefT_west,elemU_west,coefU_west,elemV_west,coefV_west)\
                =interp_tools.get_delaunay_bry(crocogrd.lon_west,crocogrd.lat_west,inpdat,'W')

            elif 'east' in boundary and is_open:
                (elemT_east,coefT_east,elemU_east,coefU_east,elemV_east,coefV_east)\
                =interp_tools.get_delaunay_bry(crocogrd.lon_east,crocogrd.lat_east,inpdat,'E')

            elif 'south' in boundary and is_open:
                (elemT_south,coefT_south,elemU_south,coefU_south,elemV_south,coefV_south)\
                =interp_tools.get_delaunay_bry(crocogrd.lon_south,crocogrd.lat_south,inpdat,'S')

            elif 'north' in boundary and is_open:
                (elemT_north,coefT_north,elemU_north,coefU_north,elemV_north,coefV_north)\
                =interp_tools.get_delaunay_bry(crocogrd.lon_north,crocogrd.lat_north,inpdat,'N')
    else:
    # Load the Delaunay triangulation matrices
        print('Load Delaunay triangulation...')
        for boundary, is_open in zip(obc_dict.keys(), obc_dict.values()):
            if 'west' in boundary and is_open:
                data=np.load('coeffs_bry'+boundary[0].upper()+'.npz')
                coefT_west = data['coefT']; elemT_west = data['elemT']
                coefU_west = data['coefU']; elemU_west = data['elemU']
                coefV_west = data['coefV']; elemV_west = data['elemV']
            if 'east' in boundary and is_open:
                data=np.load('coeffs_bry'+boundary[0].upper()+'.npz')
                coefT_east = data['coefT']; elemT_east = data['elemT']
                coefU_east = data['coefU']; elemU_east = data['elemU']
                coefV_east = data['coefV']; elemV_east = data['elemV']
            if 'south' in boundary and is_open:
                data=np.load('coeffs_bry'+boundary[0].upper()+'.npz')
                coefT_south = data['coefT']; elemT_south = data['elemT']
                coefU_south = data['coefU']; elemU_south = data['elemU']
                coefV_south = data['coefV']; elemV_south = data['elemV']
            if 'north' in boundary and is_open:
                data=np.load('coeffs_bry'+boundary[0].upper()+'.npz')
                coefT_north = data['coefT']; elemT_north = data['elemT']
                coefU_north = data['coefU']; elemU_north = data['elemU']
                coefV_north = data['coefV']; elemV_north = data['elemV']


    # --- Work on date format for the loop in time ------------------------

    startloc=plt.datetime.datetime(int(start_date[:4]),
                                   int(start_date[4:6]),
                                   1)
    if output_file_format.upper() == "MONTHLY":
        endloc= startloc+relativedelta(months=1,days=-1,hours=12)
    elif output_file_format.upper() == "YEARLY":
        if plt.datetime.datetime.strptime(str(startloc), "%Y-%m-%d %H:%M:%S").year == int(Yend) :
            endloc=plt.num2date(dtend).replace(tzinfo=None)
        else:
            endloc= plt.datetime.datetime(int(start_date[:4]), 12,31,12)

    elif output_file_format.upper() == "FULL":
        endloc=plt.num2date(dtend).replace(tzinfo=None)
    else:
        print("\n Output file format \"%s\" is not setup. Pease change it to MONTHLY, YEARLY or FULL")
        sys.exit()
 
    # --- Start time loop loop in time ------------------------------------

    while plt.date2num(endloc) <= dtend:

        # Load full time dataset
        time = plt.date2num(inpdat.ncglo['time'].values)
        # find index for the time range 
        ind= np.where((time>plt.date2num(startloc)) & (time<=plt.date2num(endloc)))
 
        if len(ind[0])==0 :
            print('\nData is missing for range %s to %s' % (startloc ,endloc))
            sys.exit()
            
        [dtmin,dtmax]=np.min(ind),np.max(ind)
        # create monthly file
        tmp_date = plt.datetime.datetime.strptime(str(startloc), "%Y-%m-%d %H:%M:%S")
            # file name depending on format chosen
        if output_file_format.upper() == "MONTHLY":
            bdy_filename = bry_filename.replace('.nc', '_%s_Y%sM%02i.nc' %(inputdata,tmp_date.year,tmp_date.month))
        elif output_file_format.upper() == "YEARLY":
            bdy_filename = bry_filename.replace('.nc', '_%s_Y%s.nc' %(inputdata,tmp_date.year))
        elif output_file_format.upper() == "FULL":
            bdy_filename = bry_filename.replace('.nc', '_%s.nc' %(inputdata))

        Croco.CROCO.create_bry_nc(None,bdy_filename,crocogrd,obc_dict,cycle_bry)
        #
        print('\n-----------------------------------')
        if output_file_format.upper() == "MONTHLY":
            print(' Processing Year %s - Month %02i' %(tmp_date.year,tmp_date.month))
        elif output_file_format.upper() == "YEARLY":
            print(' Processing Year %s' %(tmp_date.year))
        elif output_file_format.upper() == "FULL": 
            tmp_end_date = plt.datetime.datetime.strptime(str(endloc), "%Y-%m-%d %H:%M:%S")
            print(' Processing from Year %s - Month %02i  to Year %s - Month %02i' %(tmp_date.year,tmp_date.month,tmp_end_date.year,tmp_end_date.month))
            del tmp_end_date

        # Check if there is at least 1 point by month when doing Yearly or full
        if output_file_format.upper() != "MONTHLY":
            tmp_str_date=startloc
            tmp_end_date=tmp_str_date+relativedelta(months=1,days=-1,hours=12)
            while(plt.date2num(tmp_end_date) <= plt.date2num(endloc)):
                ind_tmp= np.where((time>plt.date2num(tmp_str_date)) & (time<=plt.date2num(tmp_end_date)))
                tmp_date=plt.datetime.datetime.strptime(str(tmp_str_date), "%Y-%m-%d %H:%M:%S")
                if len(ind_tmp[0])==0:
                    print('\nLacking %s data for  Y%s - M%02i' %(inputdata,tmp_date.year,tmp_date.month))
                    sys.exit()            
                tmp_str_date=tmp_end_date+relativedelta(days=1,hours=-12)
                tmp_end_date=tmp_str_date+relativedelta(months=1,days=-1,hours=12)

            del tmp_str_date,tmp_end_date

        # --- Check if data availabla for the surrounded months -----------
        print('-----------------------------------')
        tmp_str_date=startloc ; tmp_end_date = endloc
        prev_month_str = tmp_str_date+relativedelta(months=-1)
        prev_month_end = tmp_str_date+relativedelta(days=-1,hours=12)
        ind_prev       = np.where((time>plt.date2num(prev_month_str)) & (time<plt.date2num(prev_month_end)))
        if len(ind_prev[0])==0:
            print('   No data for the previous month: using current month')
            prev=1
        else:
            prev=0
            dtmin=dtmin-1 # create overlap before (in this case it takes the previous month)
        del prev_month_str,prev_month_end,ind_prev
        #
        next_month_str = tmp_end_date+relativedelta(days=1)
        next_month_end = next_month_str+relativedelta(months=1,days=-1,hours=12)
        ind_next       = np.where((time>plt.date2num(next_month_str)) & (time<plt.date2num(next_month_end)))
        if len(ind_next[0])==0:
            print('   No data for the next month: using current month')
            nxt=1
        else:
            nxt=0
            dtmax=dtmax+1 # create overlap after (in this case it takes the next month)

        del next_month_str,next_month_end,ind_next,tmp_str_date,tmp_end_date
         
        if prev ^ nxt and np.std(np.gradient(time[dtmin:dtmax+1])) >=5: # Abnormal distribution of days
            Question = input( "Abnormal distribution of days (standart deviation > 5 days) \
                    \nThis can be due to the use of different time resolution dataset.\
                    \n Do you want to proceed?: y,n ")
            if Question.lower() == ("y") or Question.lower() == ("yes"):
                continue
            else:
                print('Aborting')
                sys.exit()
         
        ## --- Handle bry_time --------------------------------------------

        bry_time= time[dtmin:dtmax+1] - day_zero_num 
        if prev == 1 and len(bry_time)==1:
            prev_time = plt.date2num(plt.num2date(bry_time[0]) + relativedelta(days=-30))
            bry_time=np.append(prev_time,bry_time)
            del prev_time
        elif prev == 1 and len(bry_time)>1: 
            date_dt=np.gradient(bry_time)[0]
            prev_time = plt.date2num(plt.num2date(bry_time[0]) + relativedelta(days=-date_dt))
            bry_time=np.append(prev_time,bry_time)
            del prev_time

        if nxt == 1 and len(bry_time)==1:
            nxt_time = plt.date2num(plt.num2date(bry_time[-1]) + relativedelta(days=30))
            bry_time=np.append(bry_time,nxt_time)
            del nxt_time
        elif nxt == 1 and len(bry_time)>1:
            date_dt=np.gradient(bry_time)[-1]
            nxt_time = plt.date2num(plt.num2date(bry_time[-1]) + relativedelta(days=date_dt))
            bry_time=np.append(bry_time,nxt_time)
            del nxt_time
        
        nc=netcdf.Dataset(bdy_filename, 'a')

        nc.variables['bry_time'].cycle=cycle_bry
        nc.variables['bry_time'][:]=bry_time

        # --- Loop on boundaries ------------------------------------------

        for boundary, is_open in zip(obc_dict.keys(), obc_dict.values()):
            if is_open:
                for vars in ['ssh','tracers','velocity']:
                    print('\n     Processing *%s* for %sern boundary' %(vars, boundary))
                    print('     ------------------------------------------')
                    if vars == 'ssh': 
                        (zeta,NzGood) = interp_tools.interp_tracers3D(inpdat,vars,-1,eval(''.join(("coefT_"+boundary))),eval(''.join(("elemT_"+boundary))),dtmin,dtmax,prev,nxt,boundary[0].upper()) 
                        z_rho = crocogrd.scoord2z_r(zeta=zeta,bdy="_"+boundary)
                        z_w   = crocogrd.scoord2z_w(zeta=zeta,bdy="_"+boundary)
    
                    elif vars == 'tracers':
                        print('\nIn tracers processing Temp')
                        temp= interp_tools.interp4d(inpdat,'temp',Nzgoodmin,z_rho,\
                                                   eval(''.join(("coefT_"+boundary))),eval(''.join(("elemT_"+boundary))),dtmin,dtmax,prev,nxt,bdy=boundary[0].upper())
            
                        print('\nIn tracers processing Salt')
                        salt= interp_tools.interp4d(inpdat,'salt',Nzgoodmin,z_rho,\
                                                eval(''.join(("coefT_"+boundary))),eval(''.join(("elemT_"+boundary))),dtmin,dtmax,prev,nxt,boundary[0].upper())
        
                    elif vars == 'velocity':

                        cosa=np.cos(eval(''.join(('crocogrd.angle_',boundary))) )
                        sina=np.sin(eval(''.join(('crocogrd.angle_',boundary))) )

                        [u,v,ubar,vbar]=interp_tools.interp4d_uv(inpdat,Nzgoodmin,z_rho,cosa,sina,\
                                           eval(''.join(("coefU_"+boundary))),eval(''.join(("elemU_"+boundary))),\
                                           eval(''.join(("coefV_"+boundary))),eval(''.join(("elemV_"+boundary))),\
                                           dtmin,dtmax,prev,nxt,bdy=boundary[0].upper())

                        conserv=1  # Correct the horizontal transport i.e. remove the intergrated tranport and add the OGCM transport          
                        if conserv == 1:
                            ubar_croco=sig_tools.vintegr4D(u,grd_tools.rho2u(z_w),grd_tools.rho2u(z_rho),np.nan,np.nan)[0]/grd_tools.rho2u(eval(''.join(('crocogrd.h_'+boundary))))
                            vbar_croco=sig_tools.vintegr4D(v,grd_tools.rho2v(z_w),grd_tools.rho2v(z_rho),np.nan,np.nan)[0]/grd_tools.rho2v(eval(''.join(('crocogrd.h_'+boundary))))

                            u = u - ubar_croco ; u = u + np.tile(ubar[:,np.newaxis,:,:],(1,z_rho.shape[1],1,1))
                            v = v - vbar_croco ; v = v + np.tile(vbar[:,np.newaxis,:,:],(1,z_rho.shape[1],1,1))


    # --- Saving in netcdf ------------------------------------------------

                print('\nSaving %sern boundary in Netcdf' % boundary)
                print('----------------------------------')
  
                 # handle indices (as 2 points where taken next to bdy)
                if str(boundary) == 'west' and is_open:
                    indices3D="[:,:,:,0]" # T,N,J,i=0
                    indices2D="[:,:,0]"   # T,J,i=0
                elif str(boundary) == 'east' and is_open:
                    indices3D="[:,:,:,-1]" # T,N,J,i=last
                    indices2D="[:,:,-1]"   # T,J,i=last
                elif str(boundary) == 'south' and is_open:
                    indices3D="[:,:,0,:]" # T,N,j=0,I
                    indices2D="[:,0,:]"   # T,j=0,I
                elif str(boundary) == 'north' and is_open:
                    indices3D="[:,:,-1,:]" # T,N,j=last,I
                    indices2D="[:,-1,:]"   # T,j=last,I

                mask_zet = np.tile(eval(''.join(('crocogrd.maskr_',boundary))),[zeta.shape[0],1,1])
                mask_tra = np.tile(eval(''.join(('crocogrd.maskr_',boundary))),[temp.shape[0],temp.shape[1],1,1])
                mask_u   = np.tile(eval(''.join(('crocogrd.umask_',boundary))),[u.shape[0],u.shape[1],1,1]) 
                mask_v   = np.tile(eval(''.join(('crocogrd.vmask_',boundary))),[u.shape[0],u.shape[1],1,1])
                mask_ubar   = np.tile(eval(''.join(('crocogrd.umask_',boundary))),[u.shape[0],1,1])
                mask_vbar   = np.tile(eval(''.join(('crocogrd.vmask_',boundary))),[v.shape[0],1,1])
                
                nc.variables['zeta_'+str(boundary)][:]=eval(''.join(('zeta',indices2D)))*eval(''.join(('mask_zet',indices2D)))
                nc.variables['temp_'+str(boundary)][:]=eval(''.join(('temp',indices3D)))*eval(''.join(('mask_tra',indices3D)))
                nc.variables['salt_'+str(boundary)][:]=eval(''.join(('salt',indices3D)))*eval(''.join(('mask_tra',indices3D)))
                nc.variables['u_'+str(boundary)][:]   =eval(''.join(('u',indices3D)))*eval(''.join(('mask_u',indices3D)))
                nc.variables['v_'+str(boundary)][:]   =eval(''.join(('v',indices3D)))*eval(''.join(('mask_v',indices3D)))
                nc.variables['ubar_'+str(boundary)][:]=eval(''.join(('ubar',indices2D)))*eval(''.join(('mask_ubar',indices2D)))
                nc.variables['vbar_'+str(boundary)][:]=eval(''.join(('vbar',indices2D)))*eval(''.join(('mask_vbar',indices2D)))

                # handle prev and nxt + save

        nc.close()
    # --- END writting netcdf ---------------------------------------------
    
    # --- Preparing time for next loop ------------------------------------
        startloc=endloc+relativedelta(days=1,hours=-12)
        if output_file_format.upper() == "MONTHLY":
            endloc= startloc+relativedelta(months=1,days=-1,hour=12)
        elif output_file_format.upper() == "YEARLY":
            yearloc=plt.datetime.datetime.strptime(str(startloc), "%Y-%m-%d %H:%M:%S")
            if yearloc.year == int(Yend) :
                endloc=plt.num2date(dtend).replace(tzinfo=None) 
            else:
                endloc= plt.datetime.datetime(int(yearloc.year), 12,31,12)
        elif output_file_format.upper() == "FULL":
            endloc=startloc

