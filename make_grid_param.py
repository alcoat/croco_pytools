import time

tra_lon =  15 # Longitude of the grid center
tra_lat = -32 # Latitude of the grid center

# Grid resolution [km]
size_x = 1556
size_y = 1334
nx     = 43
ny     = 44

# Rotation [degree]
rot = 0

# Smoothing params
hmin        = 20         # Minimum depth [m]
hmax        = 6000       # Maximum depth [m]
smth_rad    = 4          # Smoothing radius [nb points] (usually between 2 and 8)
rfact       = 0.2        # Maximum r-fact to reach (the lower it is, the smoother it will be)
smooth_meth = 'smooth'   # Smoothing method ('smooth', 'lsmooth', 'lsmooth_legacy', 'lsmooth2', 'lsmooth1', 'cond_rx0_topo')
# See README.topo for more intel

# Topo file
topofile = './data/input/etopo2.nc'

# Coastline file (for the mask)
shp_file = './data/input/gshhs/GSHHS_shp/f/GSHHS_f_L1.shp'

# Single Connect [Mask water not connected to the main domain]
# Precise True or false and a point index inside the main domain
sgl_connect = [False, 20, 20]

# Output dir
actual_time = int(time.time())

output_dir = "./data/output/"

# Remove / if has one
output_dir = output_dir.rstrip('/')

output_file = f"{output_dir}/{actual_time}_croco_grd.nc"

# FOR ZOOM 
# We do not know why it is used
# I's called here :
# Modules/graphicUI_tools/main_window.py:        self.compute_zm_thread.openb  = [self.checklist,self.merge]
# Modules/graphicUI_tools/main_window.py:        self.compute_c2c_thread.openb  = [self.checklist,self.merge]
zoom_output_file = f"{output_dir}/{actual_time}_croco_chd_grd.nc"
zoom_merge = 5
zoom_tra_lon = 18
zoom_tra_lat = -33
zoom_size_x = 550 
zoom_size_y =  550
zoom_rot = 0
zoom_nx = 55
zoom_ny = 55

# FOR AGRIF
c2c_coef = 3
c2c_imin = 35
c2c_imax = 55
c2c_jmin = 8
c2c_jmax = 28