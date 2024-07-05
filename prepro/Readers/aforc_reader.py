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

# List of possible variables : see reader_list.txt


# T2M : Better to use temperature at 1000 hPa because humidity if given at 1000 hPa
# Except : CFSR : Humidity given at 2m -> better to use T2M


# Conversion pour arriver à l'unité inscrite à gauche

def lookvar(data_origin):

    if data_origin == 'era_dataref':
        var_info = [ ['lon'  ,'longitude025'  ,'degrees'    ,1.           ],\
                     ['lat'  ,'latitude025'   ,'degrees'    ,1.           ],\
                     ['sst'  ,'sst'           ,'K'          ,1.           ],\
                     ['tp'   ,'tp'            ,'kg m-2 s-1' ,cff_tp       ],\
                     ['ssr'  ,'ssr'           ,'W m-2'      ,cff_heat     ],\
                     ['t2m'  ,'t2m'           ,'K'          ,1.           ],\
                     ['u10m' ,'u10'           ,'m s-1'      ,1.           ],\
                     ['v10m' ,'v10'           ,'m s-1'      ,1.           ],\
                     ['str'  ,'str'           ,'W m-2'      ,cff_heat     ],\
                     ['r'    ,'r'             ,'%'          ,1.           ],\
                     ['msl'  ,'msl'           ,'Pa'         ,1.           ] \
                   ]

    elif data_origin == 'era_ecmwf':
        var_info = [ ['lon'  ,'longitude'     ,'degrees'    ,1.           ,' '       ],\
                     ['lat'  ,'latitude'      ,'degrees'    ,1.           ,' '       ],\
                     ['sst'  ,'sst'           ,'K'          ,1.           ,'SST'     ],\
                     ['tp'   ,'tp'            ,'kg m-2 s-1' ,cff_tp       ,'TP'      ],\
                     ['ssr'  ,'ssr'           ,'W m-2'      ,cff_heat     ,'SSR'     ],\
                     ['t2m'  ,'t2m'           ,'K'          ,1.           ,'T2M'     ],\
                     ['u10m' ,'u10'           ,'m s-1'      ,1.           ,'U10'     ],\
                     ['v10m' ,'v10'           ,'m s-1'      ,1.           ,'V10'     ],\
                     ['strd' ,'strd'          ,'W m-2'      ,cff_heat     ,'STRD'    ],\
                     ['q'    ,'q'             ,'kg kg-1'    ,1.           ,'Q'       ],\
                     ['msl'  ,'msl'           ,'Pa'         ,1.           ,'LSM'     ] \
                   ]

    elif data_origin == 'cfsr':
        var_info = [ ['lon'  ,'lon'           ,'degrees'    ,1.           ,' '       ],\
                     ['lat'  ,'lat'           ,'degrees'    ,1.           ,' '       ],\
                     ['tp'   ,'PRATE_L1_Avg_1','kg m-2 s-1' ,1.           ,'prate'   ],\
                     ['dswrf','DSWRF_L1_Avg_1','W m-2'      ,1.           ,'dswsfc'  ],\
                     ['uswrf','USWRF_L1_Avg_1','W m-2'      ,1.           ,'uswsfc'  ],\
                     ['t2m'  ,'TMP_L103'      ,'K'          ,1.           ,'tmp2m'   ],\
                     ['u10m' ,'U_GRD_L103'    ,'m s-1'      ,1.           ,'wnd10m'  ],\
                     ['v10m' ,'V_GRD_L103'    ,'m s-1'      ,1.           ,'wnd10m'  ],\
                     ['strd' ,'DLWRF_L1_Avg_1','W m-2'      ,1.           ,'dlwsfc'  ],\
                     ['q'    ,'SPF_H_L103'    ,'kg kg-1'    ,1.           ,'q2m'     ],\
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
                     ['r'    ,'rh'            ,'(0-1)'      ,1.           ],\
                     ['msl'  ,'pmer'          ,'Pa'         ,1.           ] \
                   ]
    elif data_origin == 'ncep':
        var_info = [ ['tp'   ,'pratesfc'      ,'kg m-2 s-1' ,1.           ],\
                     ['dswrf','dswrfsfc'      ,'W m-2'      ,1.           ],\
                     ['uswrf','uswrfsfc'      ,'W m-2'      ,1.           ],\
                     ['t2m'  ,'tmp2m'         ,'K'          ,1.           ],\
                     ['strd' ,'dlwrfsfc'      ,'W m-2'      ,1.           ],\
                     ['u10m' ,'ugrd10m'       ,'m s-1'      ,1.           ],\
                     ['v10m' ,'vgrd10m'       ,'m s-1'      ,1.           ],\
                     ['q'    ,'spfh2m'        ,'kg kg-1'    ,1.           ],\
                   ]
    elif data_origin == 'ncep_ftp':
        var_info = [ ['tp'   ,'prate'         ,'kg m-2 s-1' ,1.           ,'prate '  ],\
                     ['dswrf','dswrf'         ,'W m-2'      ,1.           ,'dswrf'   ],\
                     ['uswrf','uswrf'         ,'W m-2'      ,1.           ,'uswrf'   ],\
                     ['t2m'  ,'air'           ,'K'          ,1.           ,'air.2m'  ],\
                     ['strd' ,'dlwrf'         ,'W m-2'      ,1.           ,'dlwrf'   ],\
                     ['u10m' ,'uwnd'          ,'m s-1'      ,1.           ,'uwnd'    ],\
                     ['v10m' ,'vwnd'          ,'m s-1'      ,1.           ,'vwnd'    ],\
                     ['q'    ,'shum'          ,'kg kg-1'    ,1.           ,'shum'    ],\
                   ]


    else:
        sys.exit(print('No \'%s\' dico available in Modules/inputs_readers.py. Please add it an run the script' % data_origin))

    return var_info

