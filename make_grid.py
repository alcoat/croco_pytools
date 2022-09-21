__author__ = 'Mathieu Le Corre'
__email__  = 'mathieu.lecorre@univ-brest.fr'
__date__   = '2022-09'
__license__='GPL3'
'''
===========================================================================
Further Information:  
  http://www.croco-ocean.org
  
This file is part of CROCOTOOLS

Create a CROCO grid file
In the current state the script can handle:
    - etopo (5,2,1)
    - srtm30
    - gebco

To add a new dataset you just have to go in Modules/topo_readers.py and
create a dico with name of lon,lat,topo in the dataset.
At this time it only handle 1D lon/lat 

The script make a grid with mercator projection centred at the equator
before rotating the sphere to put the middle of the grid in tra_lon/tra_lat
position.

Then, it reads topo's dataset and apply the desired smoothing
The mask is generated using the coastline dataset GSHHS

The smoothing and mask generation are fortran functions to be faster
'''

#--- Dependencies ---------------------------------------------------------
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys
sys.path.append("./Modules/")
sys.path.append("./Modules/graphicUI_tools/")
import tools_topo
from grid_main import *
#--------------------------------------------------------------------------


#--- USER CHANGES ---------------------------------------------------------

tra_lon =  15 # Longitude of the grid center 
tra_lat = -32 # Latitude of the grid center

# Grid resolution [km]
size_x = 1556
size_y = 1334
nx     = 62
ny     = 53

# Rotation [degree]
rot=0

# Smoothing params
hmin        = 20    # Minimum depth [m]
hmax        = 6000  # Maximum depth [m]
smth_rad    = 4     # Smoothing radius [nb points] (usually between 2 and 8)
rfact       = 0.2   # Maximum r-fact to reach (the lower it is, the smoother it will be)    
smooth_meth = 'smooth' # Smoothing method ('smooth', 'lsmooth', 'lsmooth_legacy', 'lsmooth2', 'lsmooth1', 'cond_rx0_topo') 
                    # See README.topo for more intel

# Topo file
topofile='./etopo5.nc'

# Mask resolution ( From basemap which is related to gshhs)
gshhs_file='./Modules/gshhs'
res="Crude" # 'Crude', 'Low', 'Intermediate', 'High', 'Full'

# Single Connect [Mask water not connected to the main domain]
sgl_connect=[True,20,50] # Precise True or false and a point index inside the main domain

# Output dir
output_dir="./"

#--- END USER CHANGES -----------------------------------------------------


# CST
r_earth=6371315 # Mean earth radius in meter


if __name__ == "__main__":


    Question = input( "Do you want to use interactive grid maker ? \
                      \n (e.g., for grid rotation or parameter adjustments) : y,[n] ")

    if Question.lower() == ("y") or Question.lower() == ("yes"):

   # --- Building grid with graphicUI -------------------------------------

        print("In interactive mode")
        MainWindow().configure_traits()
    elif Question.lower() == ("n") or Question.lower() == ("no"):
   
   # --- Building grid without graphicUI ----------------------------------

        print("In normal mode")
 
        from grid_hastraits import EasyGrid,Outputs,GetMask,GetTopo,BmapOptions,Save2Netcdf

        # --- Create inputs and outputs class -----------------------------

        inputs=tools_topo.inputs(tra_lon,tra_lat,size_x,size_y,nx,ny,rot)
        inputs_smth=tools_topo.inputs_smth(hmin,hmax,smth_rad,rfact,smooth_meth)
        outputs=Outputs()

        # --- Create lon/lat grid -----------------------------------------

        EasyGrid.easygrid(None,inputs,outputs)
        
        # --- Build mask and topo -----------------------------------------

        bmap=BmapOptions()
        bmap.bmap_res=res 
        
        GetMask.mask(None,outputs,bmap,gshhs_file,sgl_connect=sgl_connect)
        GetTopo.topo(None,outputs,topofile,smooth=inputs_smth)

        # --- Save netcdf -------------------------------------------------
       
        print('Writing Topography')
        Save2Netcdf.save2netcdf(None,output_dir,inputs,outputs)
















