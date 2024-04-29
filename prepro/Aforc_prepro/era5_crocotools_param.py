#
# For ERA5 python crocotools parameters list
#
# CAUTION IT MUST BE CONSISTENT with your MATLAB CROCOTOOLS_PARAM.m file in Run directory
# *******************************************************************************
#                         U S E R  *  O P T I O N S
# *******************************************************************************
import sys
sys.path.append("./Readers/")
from aforc_reader import lookvar
from aforc_class import aforc_class, create_class
#
# General path
#
config_dir = 'path_to_my_run_dir/'            # must be the same than crocotools_param
config_name = 'my_config'

#
# Data origin : era_dataref (datarmor) / era_ecmwf (downwloaded by script) / arome
#
data_origin = 'era_ecmwf'

#
# Original ERA5 directory
#
#era5_dir_raw = config_dir + 'DATA/ERA5_native_' + config_name
if data_origin == 'era_dataref': era5_dir_raw = '/path_to_ERA5_dataref_data'
else: era5_dir_raw = '/path_to_downloaded_ERA5_data'
#

# Output ERA5 directory
#
#era5_dir_processed = config_dir + 'DATA/ERA5_' + config_name
#

#
# READ PATM
#
READ_PATM = False # must match with READ_PATM key (air_pressure_at_sea_level)


# extraction wave variables
#
wave_extract=False # True to extract wave variables
#
# Dates limits
#
year_start = 1980
month_start = 1
year_end = 1980
month_end =1
#
# Year origin of time
#
Yorig=1950
#
# Overlapping days (at the beginning/end of each month)
#
n_overlap = 0
#
# Request time (daily hours '00/01/.../23')
#
time = '00/01/02/03/04/05/06/07/08/09/10/11/12/13/14/15/16/17/18/19/20/21/22/23'
#
# Request variables (see available at ERA5_variables.json)
#variables = ['lsm','tp','strd','ssr','t2m','q','u10','v10'] #note lsm is land_sea_mask
#
# Request area ([north, west, south, east])
#
ownArea = 0 	# 0 if area from a crocotools_param.m file
                # 1 if own area

if ownArea == 0: 
    # To complete if ownArea==0
    paramFile = config_dir + 'crocotools_param.m' # path the crocotools_param file of the simulation
    
else:
    # To complete if ownArea==1
    lonmin=7
    lonmax=23
    latmin=-45
    latmax=-20

#
# Reader
#
variables = create_class(lookvar(data_origin))

if wave_extract:
    ## append waves variables
    wave_var=['swh', 'mwd', 'pp1d' ,'cdww'];variables.extend(wave_var)
    wave_conv_cff=[1.,  1., 1. , 1.];  conv_cff.extend(wave_conv_cff)
    wave_units=['m','Degrees true','s', 'dimensionless']; units.extend(wave_units)


# *******************************************************************************
#                         E N D     U S E R  *  O P T I O N S
# *******************************************************************************
