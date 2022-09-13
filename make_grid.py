# -*- coding: utf-8 -*-
# %run easygrid_python.py
ETS_TOOLKIT="wx"
import numpy as np
import matplotlib
# We want matplotlib to use a wxPython backend
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
#from mpl_toolkits.basemap    import Basemap, maskoceans

from threading import Thread
from time import sleep
import wx
#from scipy import *

import scipy.interpolate as si

from traits.api import *
from traitsui.api import View, Item, Group, HSplit, Handler, EnumEditor, FileEditor,DirectoryEditor
from traitsui.menu import NoButtons
from traitsui.wx.editor import Editor
#from traitsui.wx.basic_editor_factory import BasicEditorFactory
from traitsui.basic_editor_factory import BasicEditorFactory

import matplotlib.pyplot as plt
import sys
sys.path.append("./Modules/")
import tools
import tools_topo
#import toolsf

sys.path.append("./Modules/tools_grid/")

from grid_main import *

import netCDF4 as netcdf
from datetime import datetime
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
        Save2Netcdf.save2netcdf(None,inputs,outputs)
















