#--- Dependencies ---------------------------------------------------------
import numpy as np
import os,sys
os.environ['ETS_TOOLKIT'] = 'wx'
sys.path.append("./Modules/")
sys.path.append("./Readers/")
sys.path.append("./Modules/graphicUI_tools/")
#----------------------------------------------#

import tools_make_grid  # Importation des fonctions depuis tools_make_grid
from tools_make_grid import GetMask, GetTopo, topo_prt, Inputs_smth_c2c, Inputs_c2c
from croco_class import CROCO
from Modules.map_tools.map_tools import plot_grid, plot_outline, plot_topo

from scipy.spatial import distance
import matplotlib.pyplot as plt
import ipywidgets as widgets
from ipywidgets import interact_manual, FloatText, VBox
from IPython.display import display, clear_output

from matplotlib.widgets import RectangleSelector

from tools_create_grids import create_child_grid_AGRIF
from tools_grid_inputs import load_namelist_agrif

#---Code - body ----------------------------------------------------------

# ============================
# DATA LOADING
# ============================

#Load namelist in a dictionnary
config= load_namelist_agrif('config_grid_agrif.ini')

coef = config['coef']
imin = config['imin']
imax = config['imax']
jmin = config['jmin']
jmax = config['jmax']
hmin = config['hmin']
hmax = config['hmax']
interp_rad = config['interp_rad']
rfact = config['rfact']
smooth_meth = config['smooth_meth']
topofile = config['topofile']
shp_file = config['shp_file']
output_file = config['output_file']
sgl_connect = [config['sgl_connect'], config['sgl_index1'], config['sgl_index2']]
parent_grid= config['parent_grid']
open_boundary= {'NORTH': config['north'], 'SOUTH': config['south'],  'WEST': config['west'],  'EAST': config['east'],  'MERGING_AREA': config['merging_area'],  'dirs': config['dirs']}




# ============================
# GRID CREATION
# ============================

chd_grid, inputs, prt= create_child_grid_AGRIF(coef, imin, imax, jmin, jmax, hmin, hmax, interp_rad, rfact, smooth_meth, topofile, shp_file, sgl_connect[0], sgl_connect[1], sgl_connect[2], open_boundary,parent_grid)

# --- Save netcdf -------------------------------------------------

prt_grd=[True,parent_grid,coef,imin,imax,jmin,jmax]
CROCO.create_grid_nc(None,output_file,inputs,chd_grid, prt_grd)
