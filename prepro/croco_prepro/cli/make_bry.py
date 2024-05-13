__author__ = 'Mathieu Le Corre'
__email__  = 'mathieu.le.corre@shom.fr'
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

To add a new dataset you just have to go in Readers/inputs_readers.py and
create a dico with correspondance between croco_var(ssh,u,v,temp,salt)
and dataset variable names.

The periodicity of a dataset (if it is -180/180 or 0/360) is handle by the
script. Just download the dataset and it should be fine

The script works as follow:
    - reads croco_grd
    - reads input data and restricts the area of spatial coverage
    - creates croco_bry.nc
    - Loop on open boundaries with:
          check if there is data before or after for continuity, if not duplicate first or last
          loop on var with:
              * horizontal interpolation
              * vertical interpolation
    - Writes data in netcdf

===========================================================================
'''

#--- Dependencies ---------------------------------------------------------
import os
import glob as glob
import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

import netCDF4 as netcdf
import pylab as plt
import numpy as np
from dateutil.relativedelta import relativedelta
import pandas as pd

from .. import croco as Croco
from .. import ibc as Inp
from .. import interp
from .. import sigma
from .. import cgrid


# Dates
Yorig = 2013                    # year defining the origin of time as: days since Yorig-01-01
Ystart, Mstart = '2013', '01'   # Starting month
Yend, Mend  = '2013','03'       # Ending month


def get_parser():
    parser = ArgumentParser(
        description="Make CROCO bry files",
        formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "date_start", type=pd.Timestamp, help="date start rounded at month")
    parser.add_argument(
        "date_end", type=pd.Timestamp, help="date end rounded at month")


    parser.add_argument(
        "inputdata", choices=["mercator_croco", "occov4", "soda", "mercator"],
        help="input reader name", default="mercator_croco")

    parser.add_argument(
        "--input-dir", default=".",
        help="directory that contains input files")
    parser.add_argument(
        "--input-prefix", default="mercator_*",
        help="file name prefix with wildcard")
    parser.add_argument(
        "--multi-files", action="store_true" ,
        help="multiple data files with time from a ssh file")

    parser.add_argument(
        "--nzgoodmin", default=4, type=int,
        help="default value to consider a z-level fine to be used")

    parser.add_argument(
        "--croco-dir", default=".",
        help="the CROCO_FILES directory")
    parser.add_argument("--croco-grd", default="croco_grd.nc")

    parser.add_argument("--sigma-theta-n", default=32, type=int)
    parser.add_argument("--sigma-theta-s", default=7, type=float)
    parser.add_argument("--sigma-theta-b", default=2, type=float)
    parser.add_argument("--sigma-theta-hc", default=75, type=float)

    parser.add_argument(
        "--bry-file-name", default="croco_bry.nc",
        help="output file name saved in the croco_dir")
    parser.add_argument(
        "--bry-file-format", default="monthly",
        choices=["monthy", "yearly", "full"],
        help="how outputs are split")
    parser.add_argument("--bry-cycle", default=0., type="float")

    parser.add_argument(
        "--no-obc-south", action="store_false", dest="obc_south",
        help="deactivate south open boundary generation")
    parser.add_argument(
        "--no-obc-north", action="store_false", dest="obc_north",
        help="deactivate north open boundary generation")
    parser.add_argument(
        "--no-obc-west", action="store_false", dest="obc_west",
        help="deactivate west open boundary generation")
    parser.add_argument(
        "--no-obc-east", action="store_false", dest="obc_east",
        help="deactivate east open boundary generation")

    parser.add_argument(
        "--tracer", action="append", dest="tracers",
        help="add tracers", default=["temp", "salt"])


    parser.add_argument(
        "--year-orig", default=2023, type=int)

    return parser

def main(parser):

    parser = get_parser()
    args = parser.parse_args()

    # Post-process args
    sigma_params = dict(theta_s=args.sigma_theta_s, theta_b=args.sigma_theta_b, N=args.sigma_n, hc=args.sigma_hcS)
    obc_dict = dict(south=args.obc_south, west=args.obc_west, east=args.obc_east, north=args.obc_north)
    if args.multi_files: # Multiple data files. Time is read in ssh file
        input_file = {'ssh':sorted(glob.glob(args.input_dir+args.input_prefix+'ETAN.*.nc')),
                      'temp':sorted(glob.glob(args.input_dir+args.input_prefix+'THETA.*.nc')),
                      'salt':sorted(glob.glob(args.input_dir+args.input_prefix+'SALT.*.nc')),
                      'u':sorted(glob.glob(args.input_dir+args.input_prefix+'EVEL.*.nc')),
                      'v':sorted(glob.glob(args.input_dir+args.input_prefix+'NVEL.*.nc'))
                    }
    else:
        input_file = sorted(glob.glob(os.path.join(args.input_dir, args.input_prefix)))


    # Put origin date to the right format
    # day_zero   = str(Yorig)+'0101'
    # day_zero_num = plt.datetime.datetime(int(day_zero[:4]),
    #                                      int(day_zero[4:6]),
    #                                      int(day_zero[6:8]))
    # day_zero_num = plt.date2num(day_zero_num)

    # Put start and end date to the right format
    # start_date = Ystart+Mstart+'01'+'12'  # defaut start day is 1st


    # dtstrdt = plt.datetime.datetime(int(start_date[:4]),
    #                                 int(start_date[4:6]),
    #                                 int(start_date[6:8]),
    #                                 int(start_date[8:]))

    # dtenddt = plt.datetime.datetime(int(Yend),int(Mend),1,12) \
    #         + relativedelta(months=1,days=-1) # Last day of the ending month


    # dtstr, dtend = plt.date2num(dtstrdt), plt.date2num(dtenddt)


    # --- Load croco_grd --------------------------------------------------

    crocogrd = Croco.CROCO_grd(os.path.join(args.croco_dir, args.croco_grd), sigma_params)

    # --- Initialize boundary vars ----------------------------------------

    crocogrd.WEST_grid()
    crocogrd.EAST_grid()
    crocogrd.SOUTH_grid()
    crocogrd.NORTH_grid()

    # --- Initialize input data class -------------------------------------

    inpdat = Inp.getdata(args.inputdata,input_file,crocogrd,args.multi_files,bdy=[obc_dict,args.bry_cycle])

    # --- Work on date format for the loop in time ------------------------

    # startloc=plt.datetime.datetime(int(start_date[:4]),
    #                                int(start_date[4:6]),
    #                                1)
    dtend = pd.Timestamp(*args.date_start.timetuple()[:2]+(1, 12))
    dtend += pd.
    startloc = pd.Timestamp(args.date_start).floor("m")

    if args.bry_file_format.upper() == "MONTHLY":
        endloc= startloc+relativedelta(months=1,days=-1,hours=12)
    elif args.bry_file_format.upper() == "YEARLY":
        if plt.datetime.datetime.strptime(str(startloc), "%Y-%m-%d %H:%M:%S").year == int(Yend) :
            endloc=plt.num2date(dtend).replace(tzinfo=None)
        else:
            endloc= plt.datetime.datetime(int(start_date[:4]), 12,31,12)

    elif args.bry_file_format.upper() == "FULL":
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
        if args.bry_file_format.upper() == "MONTHLY":
            bdy_filename = args.croco_dir+args.bry_file_name.replace('.nc', '_%s_Y%sM%02i.nc' %(inputdata,tmp_date.year,tmp_date.month))
        elif args.bry_file_format.upper() == "YEARLY":
            bdy_filename = args.croco_dir+args.bry_file_name.replace('.nc', '_%s_Y%s.nc' %(inputdata,tmp_date.year))
        elif args.bry_file_format.upper() == "FULL":
            bdy_filename = args.croco_dir+args.bry_file_name.replace('.nc', '_%s.nc' %(inputdata))

        Croco.CROCO.create_bry_nc(None,bdy_filename,crocogrd,obc_dict,cycle_bry,tracers=args.tracers)
        #
        print('\n-----------------------------------')
        if args.bry_file_format.upper() == "MONTHLY":
            print(' Processing Year %s - Month %02i' %(tmp_date.year,tmp_date.month))
        elif args.bry_file_format.upper() == "YEARLY":
            print(' Processing Year %s' %(tmp_date.year))
        elif args.bry_file_format.upper() == "FULL":
            tmp_end_date = plt.datetime.datetime.strptime(str(endloc), "%Y-%m-%d %H:%M:%S")
            print(' Processing from Year %s - Month %02i  to Year %s - Month %02i' %(tmp_date.year,tmp_date.month,tmp_end_date.year,tmp_end_date.month))
            del tmp_end_date

        # Check if there is at least 1 point by month when doing Yearly or full
        if args.bry_file_format.upper() != "MONTHLY":
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

        if np.nanvar(np.gradient(time[dtmin:dtmax+1])) >=5: # Abnormal distribution of days
            Question = input( "Abnormal distribution of days (variance to high) \
                    \nThis may be due to the use of different temproral resolution dataset.\
                    \n Do you want to proceed?: y,[n] ") or 'no'
            if Question.lower() == ("n") or Question.lower() == ("no"):
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

        nc.Input_data_type=inputdata
        nc.variables['bry_time'].cycle=args.bry_cycle
        nc.variables['bry_time'][:]=bry_time
        if args.bry_cycle==0:
            nc.variables['bry_time'].units='days since %s-01-01 00:00:00' %(Yorig)
        # --- Loop on boundaries ------------------------------------------

        if len(tracers) == 0:
            var_loop = ['ssh','velocity']
        else:
            var_loop = ['ssh','tracers','velocity']

        for boundary, is_open in zip(obc_dict.keys(), obc_dict.values()):
            if is_open:
                for vars in var_loop:
                    print('\n     Processing *%s* for %sern boundary' %(vars, boundary))
                    print('     ------------------------------------------')
                    if vars == 'ssh':
                        (zeta,NzGood) = interp.interp_tracers(inpdat,vars,-1,crocogrd,dtmin,dtmax,prev,nxt,boundary[0].upper())
                        z_rho = crocogrd.scoord2z_r(zeta=zeta,bdy="_"+boundary)
                        z_w   = crocogrd.scoord2z_w(zeta=zeta,bdy="_"+boundary)

                    elif vars == 'tracers':
                        trac_dict = dict()
                        for trc in tracers:
                            print(f'\nIn tracers processing {trc}')
                            trac_dict[trc] = interp_tools.interp(inpdat,trc,Nzgoodmin,z_rho,crocogrd,dtmin,dtmax,prev,nxt,bdy=boundary[0].upper())

                    elif vars == 'velocity':

                        cosa=np.cos(eval(''.join(('crocogrd.angle_',boundary))) )
                        sina=np.sin(eval(''.join(('crocogrd.angle_',boundary))) )

                        [u,v,ubar,vbar]=interp.interp_uv(inpdat,args.nzgoodmin,z_rho,cosa,sina,crocogrd,dtmin,dtmax,prev,nxt,bdy=boundary[0].upper())

                        conserv=1  # Correct the horizontal transport i.e. remove the intergrated tranport and add the OGCM transport
                        if conserv == 1:
                            ubar_croco=sigma.vintegr4D(u,cgrid.rho2u(z_w),cgrid.rho2u(z_rho),np.nan,np.nan)[0]/cgrid.rho2u(eval(''.join(('crocogrd.h_'+boundary))))
                            vbar_croco=sigma.vintegr4D(v,cgrid.rho2v(z_w),cgrid.rho2v(z_rho),np.nan,np.nan)[0]/cgrid.rho2v(eval(''.join(('crocogrd.h_'+boundary))))

                            u = u - np.tile(ubar_croco[:,np.newaxis,:,:],(1,z_rho.shape[1],1,1))
                            u = u + np.tile(ubar[:,np.newaxis,:,:],(1,z_rho.shape[1],1,1))
                            v = v - np.tile(vbar_croco[:,np.newaxis,:,:],(1,z_rho.shape[1],1,1))
                            v = v + np.tile(vbar[:,np.newaxis,:,:],(1,z_rho.shape[1],1,1))


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
<<<<<<< HEAD:prepro/make_bry.py
                if "velocity" in var_loop:
                    mask_u   = np.tile(eval(''.join(('crocogrd.umask_',boundary))),[u.shape[0],u.shape[1],1,1])
                    mask_v   = np.tile(eval(''.join(('crocogrd.vmask_',boundary))),[u.shape[0],u.shape[1],1,1])
                    mask_ubar   = np.tile(eval(''.join(('crocogrd.umask_',boundary))),[u.shape[0],1,1])
                    mask_vbar   = np.tile(eval(''.join(('crocogrd.vmask_',boundary))),[v.shape[0],1,1])

=======
                mask_u   = np.tile(eval(''.join(('crocogrd.umask_',boundary))),[u.shape[0],u.shape[1],1,1])
                mask_v   = np.tile(eval(''.join(('crocogrd.vmask_',boundary))),[u.shape[0],u.shape[1],1,1])
                mask_ubar   = np.tile(eval(''.join(('crocogrd.umask_',boundary))),[u.shape[0],1,1])
                mask_vbar   = np.tile(eval(''.join(('crocogrd.vmask_',boundary))),[v.shape[0],1,1])

>>>>>>> 6ac27aa (mv prepro files to standard locations):prepro/croco_prepro/cli/make_bry.py
                nc.variables['zeta_'+str(boundary)][:]=eval(''.join(('zeta',indices2D)))*eval(''.join(('mask_zet',indices2D)))
                nc.variables['u_'+str(boundary)][:]   =eval(''.join(('u',indices3D)))*eval(''.join(('mask_u',indices3D)))
                nc.variables['v_'+str(boundary)][:]   =eval(''.join(('v',indices3D)))*eval(''.join(('mask_v',indices3D)))
                nc.variables['ubar_'+str(boundary)][:]=eval(''.join(('ubar',indices2D)))*eval(''.join(('mask_ubar',indices2D)))
                nc.variables['vbar_'+str(boundary)][:]=eval(''.join(('vbar',indices2D)))*eval(''.join(('mask_vbar',indices2D)))

                if 'tracers' in var_loop:
                    for varname, value in zip(trac_dict.keys(), trac_dict.values()):
                        mask_tra = np.tile(eval(''.join(('crocogrd.maskr_',boundary))),[value.shape[0],value.shape[1],1,1])
                        nc.variables[f"{varname}_{boundary}"][:] = eval(f'value{indices3D}')*eval(''.join(('mask_tra',indices3D)))

                # handle prev and nxt + save

        nc.close()
    # --- END writting netcdf ---------------------------------------------

    # --- Preparing time for next loop ------------------------------------
        startloc=endloc+relativedelta(days=1,hours=-12)
        if args.bry_file_format.upper() == "MONTHLY":
            endloc= startloc+relativedelta(months=1,days=-1,hour=12)
        elif args.bry_file_format.upper() == "YEARLY":
            yearloc=plt.datetime.datetime.strptime(str(startloc), "%Y-%m-%d %H:%M:%S")
            if yearloc.year == int(Yend) :
                endloc=plt.num2date(dtend).replace(tzinfo=None)
            else:
                endloc= plt.datetime.datetime(int(yearloc.year), 12,31,12)
        elif args.bry_file_format.upper() == "FULL":
            endloc=startloc

