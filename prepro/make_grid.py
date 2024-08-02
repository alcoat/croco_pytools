#--- Dependencies ---------------------------------------------------------
import numpy as np
import os,sys
os.environ['ETS_TOOLKIT'] = 'wx'
sys.path.append("./Modules/")
sys.path.append("./Readers/")
sys.path.append("./Modules/graphicUI_tools/")
from main_window import *
from croco_class import CROCO

from tools_make_grid import *
from tools_create_grid import create_classic_grid

#---Code - body ----------------------------------------------------------

# ============================
# DATA LOADING
# ============================

#Load namelist in a dictionnary
config = load_namelist('config_grid.ini')

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
sgl_connect = [config['sgl_connect'], config['sgl_index1'], config['sgl_index2']

# ============================
# GRID CREATION
# ============================
               
#Merging topo/bathy,coastline mask and grid size and mesh data
grid= create_classic_grid(tra_lon, tra_lat, size_x, size_y, nx, ny, rot, hmin, hmax, interp_rad, rfact, smooth_meth, topofile, shp_file, sgl_connect, output_file)

# --- Save netcdf -------------------------------------------------
print('Writing Topography')
CROCO.create_grid_nc(None, output_file, inputs, grid)