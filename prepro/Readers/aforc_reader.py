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


# T2M : Better to use temperature at 1000 hPa because humidity if given at 1000 hPa
# Except : CFSR : Humidity given at 2m -> better to use T2M

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

    elif data_origin == 'CFSR':
        var_info = [ ['tp'   ,'PRATE_L1_Avg_1','kg m-2 s-1' ,1.           ],\
                     ['dswrf','DSWRF_L1_Avg_1','W m-2'      ,1.           ],\
                     ['uswrf','USWRF_L1_Avg_1','W m-2'      ,1.           ],\
                     ['t2m'  ,'TMP_L103'      ,'K'          ,1.           ],\
                     ['u10m' ,'U_GRD_L103'    ,'m s-1'      ,1.           ],\
                     ['v10m' ,'V_GRD_L103'    ,'m s-1'      ,1.           ],\
                     ['strd' ,'DLWRF_L1_Avg_1','W m-2'      ,1.           ],\
                     ['q'    ,'SPF_H_L103'    ,'kg kg-1'    ,1.           ],\
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

    elif data_origin == 'arome':
        var_info = [ ['tp'   ,'rain'          ,'kg m-2 s-1' ,1.           ],\
                     ['ssr'  ,'swhf'          ,'W m-2'      ,1.           ],\
                     ['t2m'  ,'t2m'           ,'K'          ,1.           ],\
                     ['u10m' ,'u10m'          ,'m s-1'      ,1.           ],\
                     ['v10m' ,'v10m'          ,'m s-1'      ,1.           ],\
                     ['strd' ,'lwhf'          ,'W m-2'      ,1.           ],\
                     ['r'    ,'rh'            ,'kg kg-1'    ,1.           ],\
                     ['msl'  ,'pmer'          ,'Pa'         ,1.           ] \
                   ]



    else:
        sys.exit(print('No \'%s\' dico available in Modules/inputs_readers.py. Please add it an run the script' % data_origin))

    return var_info

