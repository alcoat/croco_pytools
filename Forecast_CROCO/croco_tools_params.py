#######################################################################
#
# crocotools_param: common parameter file for the preprocessing
#                  of CROCO simulations using CROCOTOOLS
#
#                  This file is used by make_grid.m+ make_forcing.m+ 
#                  make_clim.m+ make_biol.m+ make_bry.m+ make_tides.m+
#                  make_NCEP.m+ make_OGCM.m+ make_...
# 
#
#######################################################################
#
# 1  - Configuration parameters
#      used by make_grid.m (and others..)
#
#######################################################################
#isoctave=exist('octave_config_info')
#
#  CROCO title names and directories
#
CROCO_title  = 'CARIBE_ANH'
CROCO_config = 'CARIBE_ANH'
#
# Grid dimensions:
#
lonmin = -98.5800    # Minimum longitude [degree east
lonmax = -58.4133    # Maximum longitude [degree east
latmin = 6.0500      # Minimum latitude  [degree north
latmax = 32.7899     # Maximum latitude  [degree north
#
# Grid resolution [degree
#
dl = 1/6
#
# Number of vertical Levels (! should be the same in param.h !)
#
N = 32
#
#  Vertical grid parameters (! should be the same in croco.in !)
#
theta_s    =  7.
theta_b    =  2.
hc         = 200.
vtransform =  2. # s-coordinate type (1: old-  2: new- coordinates)
                  # ! take care to define NEW_S_COORD cpp-key in cppdefs.h 
#
# Topography: choice of filter
#
topo_smooth =  1 # 1: old  2: new filter (better but slower)
#
# Minimum depth at the shore [m (depends on the resolution+
# rule of thumb: dl=1+ hmin=300+ dl=1/4+ hmin=150+ ...)
# This affect the filtering since it works on grad(h)/h.
#
hmin = 75
#
# Maximum depth at the shore [m (to prevent the generation
# of too big walls along the coast)
#
hmax_coast = 500
#
# Maximum depth [m (cut the topography to prevent
# extrapolations below WOA data)
#
hmax = 7000
#
# Slope parameter (r=grad(h)/h) maximum value for topography smoothing
#
rtarget = 0.25
#
# Number of pass of a selective filter to reduce the isolated
# seamounts on the deep ocean.
#
n_filter_deep_topo=4
#
# Number of pass of a single hanning filter at the end of the
# smooting procedure to ensure that there is no 2DX noise in the 
# topography.
#
n_filter_final=2
#
# Objective analysis decorrelation scale [m
# (if Roa=0: nearest extrapolation method crude but much cheaper)
#
#Roa=300e3
Roa=0
#
interp_method = 'cubic' # Interpolation method: 'linear', 'cubic' or 'quintic'
#
makeplot     = 0         # 1: create graphics after each preprocessing step
#
#######################################################################
#
# 2 - Generic file and directory names 
#
#######################################################################
#
#  CROCOTOOLS directory
#
CROCOTOOLS_dir = '/home/oscar/Documentos/ANH/croco_tools_py-main/'
#
#  Run directory
#
RUN_dir = CROCOTOOLS_dir
#
#  CROCO input netcdf files directory
#
CROCO_files_dir=RUN_dir+'CROCO_FILES/' 
#
#  Global data directory (etopo+ coads+ datasets download from ftp+ etc..)
#
DATADIR='/home/oscar/Documentos/ANH/DATASETS/' 
#
#  Forcing data directory (ncep+ quikscat+ datasets download with opendap+ etc..)
#
FORC_DATA_DIR = RUN_dir+'DATA/'
#
# CROCO file names (grid+ forcing+ bulk+ climatology+ initial)
#
grdname  = CROCO_files_dir+'croco_grd.nc'
frcname  = CROCO_files_dir+'croco_frc.nc'
blkname  = CROCO_files_dir+'croco_blk.nc'
clmname  = CROCO_files_dir+'croco_clm.nc'

bryname  = CROCO_files_dir+'croco_bry.nc'
ininame  = CROCO_files_dir+'croco_ini.nc'
#
# intermediate z-level data files (not used in simulations)
#
oaname   = CROCO_files_dir+'croco_oa.nc'    # for climatology data processing
Zbryname = CROCO_files_dir+'croco_bry_Z.nc' # for boundary data processing
#
# Generic forcing file root names for interannual simulations (NCEP/GFS)
#
frc_prefix=CROCO_files_dir+'croco_frc'     # forcing file name 
blk_prefix=CROCO_files_dir+'croco_blk'      # bulk file name
#
#######################################################################
#
#
#  Topography netcdf file name (ETOPO 2 or any other netcdf file
#  in the same format)
#
topofile = DATADIR+'Topo/etopo2.nc'
#
#######################################################################
#
# 3 - Surface forcing parameters
#     used by make_forcing.m and by make_bulk.m
#
#######################################################################
#
# COADS directory (for climatology runs)
#
coads_dir=DATADIR+'COADS05/'
#
# COADS time (for climatology runs)
#
#coads_time=(15:30:345) # days: middle of each month
coads_cycle=360        # repetition of a typical year of 360 days  
#
#coads_time=(15.2188:30.4375:350.0313) # year of 365.25 days in case
#coads_cycle=365.25                    # interannual QSCAT winds  
#                                       # are used with clim. heat flux
#
# Pathfinder SST data used by pathfinder_sst.m
#
pathfinder_sst_name=DATADIR+'SST_pathfinder/climato_pathfinder.nc'
#
#######################################################################
#
# 4 - Open boundaries and initial conditions parameters
#     used by make_clim.m+ make_biol.m+ make_bry.m
#             make_OGCM.m and make_OGCM_frcst.m
#
#######################################################################
#
#  Open boundaries switches (! should be consistent with cppdefs.h !)
#
obc = [0, 1, 1, 0] # open boundaries (1=open + [S E N W])
#
#  Level of reference for geostrophy calculation
#
zref = -7000
#
#  initial/boundary data options (1 = process)
#  (used in make_clim+ make_biol+ make_bry+
#   make_OGCM.m and make_OGCM_frcst.m)
#
makeini    = 1   # initial data
makeclim   = 1   # climatological data (for boundaries and nudging layers)
makebry    = 1   # lateral boundary data
#
makeoa     = 1   # oa data (intermediate file)
makeZbry   = 1   # boundary data in Z coordinate (intermediate file)
insitu2pot = 1   # transform in-situ temperature to potential temperature
#
#  Day of initialisation for climatology experiments (=0 : 1st january 0h)
#
tini=0  
#
# Select Climatology Atlas (temp+ salt and biological variables) from:
#    - World Ocean Atlas directory (WOA2009)  OR ...
#    - CARS2009 climatology directory (CARS2009)
#
woa_dir       = DATADIR+'WOA2009/'
cars2009_dir  = DATADIR+'CARS2009/'
climato_dir   = woa_dir
#
# Pisces biogeochemical seasonal climatology
#
woapisces_dir = DATADIR+'WOAPISCES/'  # only compatible with woa_dir
#
# Surface chlorophyll seasonal climatology (SeaWifs)
#
chla_dir=DATADIR+'SeaWifs/'
#
# Runoff monthly seasonal climatology (Dai and Trenberth)
#
global_clim_riverdir=DATADIR+'RUNOFF_DAI/'
global_clim_rivername=global_clim_riverdir+'Dai_Trenberth_runoff_global_clim.nc'
#
#  Set times and cycles for the boundary conditions: 
#   monthly climatology 
#
#woa_time=(15:30:345) # days: middle of each month
woa_cycle=360        # repetition of a typical year of 360 days  
#
#woa_time=(15.2188:30.4375:350.0313) # year of 365.25 days in case
#woa_cycle=365.25                    # interannual QSCAT winds are used 
#                                     # with clim. boundary conditions
#
#   For rivers setup : go in the routine Rivers/make_runoff.m to
#   setup your options
#        
#######################################################################
#
# 5 - Reference date and simulation times
#     (used for make_tides+ make_CFSR (or make_NCEP)+ make_OGCM)
#
#######################################################################
#
Yorig         = 2005          # reference time for vector time
                              # in croco initial and forcing files
#
Ymin          = 2005          # first forcing year
Ymax          = 2005          # last  forcing year
Mmin          = 1             # first forcing month
Mmax          = 1             # last  forcing month
#
Dmin          = 1             # Day of initialization
Hmin          = 0             # Hour of initialization
Min_min       = 0             # Minute of initialization
Smin          = 0             # Second of initialization
#
SPIN_Long     = 0             # SPIN-UP duration in Years
#
Mth_format    = '#02d'        # Number of digit for month on input files
#
#######################################################################
#
# 6 - Parameters for Interannual forcing (SODA+ ECCO+ CFSR+ NCEP+ ...)
#
#######################################################################
#
Download_data = 1   # Get data from OPENDAP sites  
level         = 0   # AGRIF level 0 = parent grid
#					  
#--------------------------------------------
# Options for make_OGCM or make_OGCM_mercator
#--------------------------------------------
#
OGCM        = 'mercator'        # Select OGCM: SODA+ ECCO+ mercator
#
OGCM_dir    = FORC_DATA_DIR+OGCM+'_'+CROCO_config+'/'  # OGCM data dir. croco format
#
bry_prefix  = CROCO_files_dir+'croco_bry_'+OGCM+'_'    # generic boundary file name
clm_prefix  = CROCO_files_dir+'croco_clm_'+OGCM+'_'    # generic climatology file name
ini_prefix  = CROCO_files_dir+'croco_ini_'+OGCM+'_'    # generic initial file name
OGCM_prefix = OGCM+'_'                                 # generic OGCM file name 

if OGCM == 'mercator':
    # For GLORYS 12 reanalysis extraction + download using python motuclient
    # ========================
    motu_url_reana='http://my.cmems-du.eu/motu-web/Motu'
    service_id_reana='GLOBAL_MULTIYEAR_PHY_001_030-TDS'
    product_id_reana='cmems_mod_glo_phy_my_0.083_P1D-m'
#
# Number of OGCM bottom levels to remove 
# (usefull if CROCO depth is shallower than OGCM depth)
#
rmdepth     = 2
#
# Overlap parameters : nb of records around each monthly sequence
#
itolap_a    = 1   # before
itolap_p    = 1   # after
#######################################################################
#
# 8 - Parameters for the forecast system
#
#     --> select OGCM name above (mercator ...)
#     --> don't forget to define in cppdefs.h:
#                    - ROBUST_DIAG
#                    - CLIMATOLOGY
#                    - BULK_FLUX
#                    - TIDES if you choose so+ but without TIDERAMP
#
#######################################################################
#
FRCST_dir = FORC_DATA_DIR+'Forecast/'  # path to local OGCM data directory
#
# Number of hindcast/forecast days
#
if OGCM=='ECCO':
    hdays=1
    fdays=6
elif OGCM=='mercator':
    hdays=1
    fdays=1
#
# Local time= UTC + timezone
#
timezone = 0
#
# Add tides
#
add_tides_fcst = 1       # 1: add tides
#
#  MERCATOR case: 
#  =============
#  To download data: set login/password (http://marine.copernicus.eu)
#  and path to croco's motuclient python package
#  or set pathMotu='' (empty) to use your own motuclient
#
#  Various sets of data are proposed in the 
#  Copernicus web site (Mercator+ UK Met Office ...)
#
user     = 'orivera1'
password = 'uvYHmM9VYpxPVh3b'

pathMotu =''

mercator_type=1   # 1 -->  1/12 deg Mercator forecast
                 # 2 -->  1/4  deg Met-Office forecast (GloSea5)
if mercator_type==1:
    motu_url_fcst='https://nrt.cmems-du.eu/motu-web/Motu'
    service_id_fcst='GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS'
    product_id_fcst='global-analysis-forecast-phy-001-024'
      
elif mercator_type==2:
    motu_url_fcst='http://nrt.cmems-du.eu/motu-web/Motu'
    service_id_fcst='GLOBAL_ANALYSISFORECAST_PHY_CPL_001_015-TDS'
    product_id_fcst='MetO-GLO-PHY-CPL-dm-TEM'