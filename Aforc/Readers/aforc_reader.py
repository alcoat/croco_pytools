import sys

#
# Variable names and conversion coefficients  
# TP: convert from accumlated m in a hour into   kg m-2 s-1
#
cff_tp=1000./3600. # m in 1 hour -> kg m-2 s-1
# Heat flux J m-2 in one hour into W m-2
#
cff_heat=1./3600.   # J m-2 in 1 hour -> W m-2
# Names, conversion coefficients and new units
#

def lookvar(data_origin):

    if data_origin == 'era_dataref':
        var_info = [ ['sst'  ,'sst'           ,'K'          ,1.           ],\
                     ['tp'   ,'tp'            ,'kg m-2 s-1' ,cff_tp       ],\
                     ['ssr'  ,'ssr'           ,'W m-2'      ,cff_heat     ],\
                     ['t2m'  ,'t2m'           ,'K'          ,1.           ],\
                     ['u10m' ,'u10'           ,'m s-1'      ,1.           ],\
                     ['v10m' ,'v10'           ,'m s-1'      ,1.           ],\
                     ['str'  ,'str'           ,'W m-2'      ,cff_heat     ],\
                     ['r'    ,'r'             ,'kg kg-1'    ,1.           ],\
                     ['msl'  ,'msl'           ,'Pa'         ,1.           ] \
                   ]

    elif data_origin == 'era_ecmwf':
        var_info = [ ['sst'  ,'sst'           ,'K'          ,1.           ],\
                     ['tp'   ,'tp'            ,'kg m-2 s-1' ,cff_tp       ],\
                     ['ssr'  ,'ssr'           ,'W m-2'      ,cff_heat     ],\
                     ['t2m'  ,'t2m'           ,'K'          ,1.           ],\
                     ['u10m' ,'u10'           ,'m s-1'      ,1.           ],\
                     ['v10m' ,'v10'           ,'m s-1'      ,1.           ],\
                     ['strd' ,'strd'          ,'W m-2'      ,cff_heat     ],\
                     ['q'    ,'q'             ,'kg kg-1'    ,1.           ],\
                     ['msl'  ,'msl'           ,'Pa'         ,1.           ] \
                   ]

    elif data_origin == 'formatted':
        var_info = [ ['sst'  ,'SST'           ,'K'          ,1.           ],\
                     ['tp'   ,'TP'            ,'kg m-2 s-1' ,1.           ],\
                     ['ssr'  ,'SSR'           ,'W m-2'      ,1.           ],\
                     ['t2m'  ,'T2M'           ,'K'          ,1.           ],\
                     ['u10m' ,'U10M'          ,'m s-1'      ,1.           ],\
                     ['v10m' ,'V10M'          ,'m s-1'      ,1.           ],\
                     ['strd' ,'STRD'          ,'W m-2'      ,1.           ],\
                     ['q'    ,'Q'             ,'kg kg-1'    ,1.           ],\
                     ['msl'  ,'MSL'           ,'Pa'         ,1.           ] \
                   ]




    else:
        sys.exit(print('No \'%s\' dico available in Modules/inputs_readers.py. Please add it an run the script' % data_origin))

    return var_info

