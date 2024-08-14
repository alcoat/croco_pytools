#--- Dependencies ---------------------------------------------------------
import numpy as np
import os,sys

os.environ['ETS_TOOLKIT'] = 'wx'
sys.path.append("./Modules/")
sys.path.append("./Readers/")
sys.path.append("./Modules/graphicUI_tools/")

from croco_class import CROCO
from tools_make_grid import *  # Importation des fonctions nécessaires de tools_make_grid
from tools_create_grids import create_child_grid_offline
from tools_grid_inputs import load_namelist


#---Code - body ----------------------------------------------------------

# ============================
# DATA LOADING
# ============================

#Load namelist in a dictionnary
config = load_namelist('config_grid_offline.ini', parent_grid= True)

#Data assiociated to variables from the '.ini' namelist
tra_lon = config['tra_lon']
tra_lat = config['tra_lat']
size_x = config['size_x']
size_y = config['size_y']
nx = config['nx']
ny = config['ny']
rot = config['rot']
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

chd_grid, inputs, prt= create_child_grid_offline(tra_lon, tra_lat, size_x, size_y, nx, ny, rot, hmin, hmax, interp_rad, rfact, smooth_meth, topofile, shp_file, sgl_connect[0], sgl_connect[1], sgl_connect[2], open_boundary, parent_grid)

# --- Save netcdf -------------------------------------------------
print('Writing Topography')
CROCO.create_grid_nc(None,output_file,inputs,chd_grid)