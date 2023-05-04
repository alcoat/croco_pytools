# CROCO path and filename informations
#-------------------------------------
croco_dir = 'CROCO_FILES/'
croco_grd = 'croco_grd.nc'
sigma_params = dict(theta_s=7, theta_b=2, N=32, hc=200) # Vertical streching, sig_surf/sig_bot/ nb level/critical depth
dl = 1. # number of degrees to consider around croco grid to extract data (for interpolations)

ini_filename    = 'croco_ini.nc'

# Dates
#-------------------------------------
#    initialization date chosen for CROCO
Yini, Mini, Dini  = '2019', '02', '01' # Month and days need to be 2-digits format
#    reference time chosen for date axis in croco files
Yzer, Mzer, Dzer = Yini, '01', Dini # reference time (default = ini time). Month and days need to be 2-digits format

# Input data informations
#-------------------------------------
datatype = 'mercator'   # At hte current time can handle mercator,soda,eccov4
#input_dir = 'DATA/'
input_dir = '/dataref/ref39/intranet/global-reanalysis-phy-001-030-daily/2019/02/'
#input_prefix = 'raw_motu_mercator_'
input_prefix='mercatorglorys12v1_gl12_mean_*'

multi_files = False # If variables are in different netcdf

#    Here check and eventually define you filename patterns
date_str = (Yini, Mini)
#input_file  = input_dir + input_prefix + 'Y%sM%s.nc' % date_str
input_file = '/dataref/ref39/intranet/global-reanalysis-phy-001-030-daily/2019/02/mercatorglorys12v1_gl12_mean_20190201_R20190206.nc'

if multi_files: # Mutiple files
    input_file = { 'ssh'  : input_dir + input_prefix + 'ETAN.%s.nc' % date_str, \
                   'temp' : input_dir + input_prefix + 'THETA.%s.nc' % date_str, \
                   'salt' : input_dir + input_prefix + 'SALT.%s.nc' % date_str, \
                   'u'    : input_dir + input_prefix + 'EVEL.%s.nc' % date_str, \
                   'v'    : input_dir + input_prefix + 'NVEL.%s.nc' % date_str\
                }

#     Here check the index of the time you want to use in the data file
tndx = 0 # time index in the file

# Default value to allow considering a z-level from the input data if enough number of values on this level are available (for interpolation issues)
#-------------------------------------
Nzgoodmin = 4  # default value to consider a z-level fine to be used

