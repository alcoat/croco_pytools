###################
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys
sys.path.append("./Modules/")
import tools
import tools_topo
sys.path.append("./Modules/tools_grid/")
from grid_main import *
#####################

tra_lon =  15
tra_lat = -32

#
# Grid resolution [km]
#
size_x = 1556
size_y = 1334
nx     = 62
ny     = 53
#
# Rotation [degree]
#
rot=0
#
# Smoothing params
hmin        = 20    # Minimum depth [m]
hmax        = 6000  # Maximum depth [m]
smth_rad    = 4     # Smoothing radius [nb points] (usually between 2 and 8)
rfact       = 0.2   # Maximum r-fact to reach (the lower it is, the smoother it will be)    
smooth_meth = 'smooth' # Smoothing method ('smooth', 'lsmooth', 'lsmooth_legacy', 'lsmooth2', 'lsmooth1', 'cond_rx0_topo') 
                    # See README.topo for more intel
#
# Topo file
#
topofile='./etopo5.nc'
#
# Mask resolution ( From basemap which is related to gshhs)
#
gshhs_file='./Modules/gshhs'
res="Crude" # 'Crude', 'Low', 'Intermediate', 'High', 'Full'
#
# Single Connect [Mask water not connected to the main domain]
#
sgl_connect=[True,20,50] # Precise True or false and a point index inside the main domain
#
# Output dir
#
output_dir="./"

###############################################
####### END USER CHANGES ######################
###############################################

r_earth=6371315 # Mean earth radius in meter
#dl = np.rad2deg((dm*1000)/r_earth)   # Grid resolution [degree]

####

###########################################################
###################     MAIN    ###########################
###########################################################
if __name__ == "__main__":


    Question = input( "Do you want to use interactive grid maker ? \
                      \n (e.g., for grid rotation or parameter adjustments) : y,[n] ")

    if Question.lower() == ("y") or Question.lower() == ("yes"):
        print("In interactive mode")
        MainWindow().configure_traits() # Go into Modules/tools_grid/grid_main.py
    elif Question.lower() == ("n") or Question.lower() == ("no"):
        print("In normal mode")
 
        from grid_hastraits import EasyGrid,Outputs,GetMask,GetTopo,BmapOptions,Save2Netcdf

        inputs=tools_topo.inputs(tra_lon,tra_lat,size_x,size_y,nx,ny,rot)
        inputs_smth=tools_topo.inputs_smth(hmin,hmax,smth_rad,rfact,smooth_meth)
        outputs=Outputs()

        EasyGrid.easygrid(None,inputs,outputs)
        ################
        bmap=BmapOptions()
        bmap.bmap_res=res 
        ###############################
        GetMask.mask(None,outputs,bmap,gshhs_file,sgl_connect=sgl_connect)
        GetTopo.topo(None,outputs,topofile,smooth=inputs_smth)

        ##############################
        
        print('Writing Topography')
        Save2Netcdf.save2netcdf(None,output_dir,inputs,outputs)
















