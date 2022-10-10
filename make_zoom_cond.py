__author__ = 'Mathieu Le Corre'
__email__  = 'mathieu.le.corre@shom.fr'
__date__   = '2022-09'
__license__='GPL3'
'''
===========================================================================
Further Information:  
  http://www.croco-ocean.org
  
This file is part of CROCOTOOLS

Create a CROCO initial and bry file for zoom

===========================================================================
'''

#--- Dependencies ---------------------------------------------------------
import glob
import numpy as np
import sys
sys.path.append("./Modules/")
import toolsf
#--------------------------------------------------------------------------

#--- USER CHANGES ---------------------------------------------------------
prt_grd='croco_grd.nc'     # Parent grid file
chd_grd='croco_chd_grd.nc' # Child grid file
ts=6;tb=4;hc=75;n=50       # Child vertical coordinates parameters

### Ini file ###
make_ini=False    # Do you build ini file
prt_his_ini='croco_his_20050201_20050205.nc' # History file to start child simulation
rec=1    # record index in the ini file

### Bry file ###

# Create a list of all the files you desired to build you bry file
# It can handle all bash sign such as ?,*,[]
# The only constraint is that each string may be less than 160 characters
# Be aware that duplicated files are only used once. Use "sorted" intead of "set" 
#  if you really want to use mutiple imes one file

make_bry=True # Do you build bry file
prt_his_bry=['croco_his_2005021*','croco_his_2005020?*']  
obc_cond='SWEN' #SWEN First letters of the boundaries that are opened.

#--- END USER CHANGES -----------------------------------------------------

# --- Make ini ------------------------------------------------------------
if make_ini:
    toolsf.r2r_init(chd_grd,ts,tb,hc,n,prt_grd,prt_his_ini,rec)

# --- Make bry ------------------------------------------------------------
if make_bry:

    # --- Build input file list from prt_his_bry --------------------------

    inputfiles=[]
    for i in range(len(prt_his_bry)):
        if len(prt_his_bry[i])>160:
            print('The size of his path need to be less than 160. It pos %i it is %i' % (i,len(prt_his_bry[i])))
            sys.exit()
        inputfiles+=glob.glob(prt_his_bry[i])
    
    inputfiles=list(set(inputfiles))

    all_files = np.zeros((160,len(inputfiles) ), dtype='c')

    for i in range(len(inputfiles)):
        ll=len(inputfiles[i])
        all_files[0:ll,i]=inputfiles[i]

    # --- Create child bry ------------------------------------------------
    toolsf.r2r_bry(chd_grd,ts,tb,hc,n,obc_cond,prt_grd,all_files)
