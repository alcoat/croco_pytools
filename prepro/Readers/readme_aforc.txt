# Indications for filling the aforc_reader.py file

#--------------------------------------------------
# Construction of the file and information needed
#--------------------------------------------------
The var_info variable in the aforc_reader.py file is composed of 4 or 5 (if multi files) data :

     1                            2                    3                    4                            5
variable_name (see below) | variable_name (netcdf) | unit | conversion_to_access_unity_wanted | variable_file_name (if multi files)

Variable_name (1) and unity (3) are must be invariant in form (except for the relative humidity, unit accepted is % or (0-1). 
The list of the possible variables is given below.

Variable_name (netcdf) (3) is the name of the variable present in the user data.

Conversion_to_access_unity_wanted (4) is the transformation needed to go from the raw data unit to the expected unit. If the unit is already the right one, just write '1'.

Variable_file_name (5) is for multi files raw data cases and so optional. Put how the variable is written in the file name.

#--------------------------------------------------
# List of the possible variables
#--------------------------------------------------

Variable | Variable definition                  | Unit wanted    | Additional processing to metadata formatting  
--------------------------------------------------------------------------------------------------
lon      | longitude                            | degrees        | Used to properly read the file
lat      | latitude                             | degrees        | Used to properly read the file
tp       | total precipitation                  | cm/day         | Conversion into cm/day if necessary
ssr      | surface net solar radiation          | W m-2          | 
dswrf    | downward short wave rad flux surface | W m-2          | (if no ssr variable) calculate ssr with uswrf
uswrf    | upward short wave rad flux surface   | W m-2          | (if no ssr variable) calculate ssr with dswrf
t2m      | temperature at 2m                    | ÂC               Conversion K to Â° if necessary
u10m     | eastward wind at 10m                 | m s-1          |
v10m     | northward wind at 10m                | m s-1          |
strd     | surface thermal radiation downwards  | W m-2          |
str      | surface net thermal radiation        | W m-2          | (if no strd variable) calculate strd with sst
sst      | sea surface temperature              | K              | (if no strd variable) calculate strd with str
r        | relative humidity                    | (0-1)          | if % convert into (0-1) if necessary
q        | specfic humidity                     | kg kg-1        |
msl      | air pressure at sea level            | Pa             |

#--------------------------------------------------
# Specific operations to the data source
#--------------------------------------------------

era_ecmwf      |  - Inverse latitude
------------------------------------
cfsr           |  - Inverse latitude



