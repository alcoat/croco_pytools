#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 15:08:24 2025

@author: annelou
"""
import sys

#
# Variable names and conversion coefficients  
# Unity wanted is written in croco_variables.json

# TP: convert from accumlated m in x hour(s) into cm/day
cff_tp = 1000./3600. * 8640 # m in x hour(s) (-> kg m-2 s-1 ->) cm/day
# TP: convert from accumlated kg m-2 s-1 into cm/day
cff_tp2 = 8640
# TP: convert from accumlated kg m-2 in x hour(s) into cm/day
cff_tp3 = 1/3600 * 8640


# Heat flux J m-2 in x hour(s) = W m-2
cff_heat = 1/3600

# Convert temperature in Celsius
cff_temp = 273.15 # If cff_temp == 273.15 then it will be calculate data-cff_temp

# Names, conversion coefficients and new units
# List of possible variables : see reader_list.txt

# T2M : Better to use temperature at 1000 hPa because humidity if given at 1000 hPa
# Except : CFSR : Humidity given at 2m -> better to use T2M


def lookvar(data_origin):

    
#   ORDER         ['var dict', 'your_var', 'unit_conv', 'cumul or not', 'name multifile']
#                       1          2           3            4              5

    if data_origin == 'era_dataref':
        var_info = [ ['lon'  ,'longitude025'  ,1.         ,' '     ],\
                     ['lat'  ,'latitude025'   ,1.         ,' '     ],\
                     ['sst'  ,'sst'           ,1.         ,' '     ],\
                     ['tp'   ,'tp'            ,cff_tp     ,'cumul' ],\
                     ['ssr'  ,'ssr'           ,cff_heat   ,' '     ],\
                     ['t2m'  ,'t2m'           ,cff_temp   ,' '     ],\
                     ['u10m' ,'u10'           ,1.         ,' '     ],\
                     ['v10m' ,'v10'           ,1.         ,' '     ],\
                     ['str'  ,'str'           ,cff_heat   ,' '     ],\
                     ['r'    ,'r'             ,1./100.    ,' '     ],\
                     ['msl'  ,'msl'           ,1.         ,' '     ] \
                   ]

    elif data_origin == 'era_ecmwf':
        var_info = [ ['lon'  ,'longitude'     ,1.         ,' '     ,' '       ],\
                     ['lat'  ,'latitude'      ,1.         ,' '     ,' '       ],\
                     ['tp'   ,'tp'            ,cff_tp     ,'cumul' ,'TP'      ],\
                     ['ssr'  ,'ssr'           ,cff_heat   ,' '     ,'SSR'     ],\
                     ['t2m'  ,'t2m'           ,cff_temp   ,' '     ,'T2M'     ],\
                     ['u10m' ,'u10'           ,1.         ,' '     ,'U10'     ],\
                     ['v10m' ,'v10'           ,1.         ,' '     ,'V10'     ],\
                     ['strd' ,'strd'          ,cff_heat   ,' '     ,'STRD'    ],\
                     ['q'    ,'q'             ,1.         ,' '     ,'Q'       ],\
                     ['msl'  ,'msl'           ,1.         ,' '     ,'LSM'     ] \
                   ]

    elif data_origin == 'arpege':
        var_info = [ ['lon'  ,'lon'           ,1.         ,' '      ,' '       ],\
                     ['lat'  ,'lat'           ,1.         ,' '      ,' '       ],\
                     ['tp'   ,'tp'            ,cff_tp3    ,'cumul'  ,'tp'      ],\
                     ['ssr'  ,'ssr'           ,1.         ,'cumul'  ,'ssr'     ],\
                     ['dswrf','dswrf'         ,1.         ,'cumul'  ,'dswrf'   ],\
                     ['t2m'  ,'t2m'           ,1.         ,' '      ,'t2m'     ],\
                     ['u10m' ,'u10m'          ,1.         ,' '      ,'u10m'    ],\
                     ['v10m' ,'v10m'          ,1.         ,' '      ,'v10m'    ],\
                     ['str'  ,'str'           ,1.         ,'cumul'  ,'str'     ],\
                     ['sst'  ,'sst'           ,cff_temp   ,' '      ,'sst'     ],\
                     ['q'    ,'q'             ,1.         ,' '      ,'q'       ]\
                   ]
            
    elif data_origin == 'arpege_grib':
        var_info = [ ['lon'  ,'longitude'     ,1.         ,' '      ,' '        ],\
                     ['lat'  ,'latitude'      ,1.         ,' '      ,' '        ],\
                     ['tp'   ,'tp'            ,cff_tp3    ,'cumul'  ,'PRECIP'   ],\
                     ['ssr'  ,'p111.128'      ,1.         ,'cumul'  ,'FLSOLAIRE'],\
                     ['t2m'  ,'t2m'           ,1.         ,' '      ,'T2m'      ],\
                     ['u10m' ,'u10'           ,1.         ,' '      ,'UV10'     ],\
                     ['v10m' ,'v10'           ,1.         ,' '      ,'UV10'     ],\
                     ['str'  ,'p112.128'     ,1.          ,'cumul'  ,'FLTHERM'  ],\
                     ['sst'  ,'p112.128'      ,cff_temp   ,' '      ,'SST'      ],\
                     ['q'    ,'q'             ,1.         ,' '      ,'Q2m'      ]\
                   ]
        

    elif data_origin == 'cfsr':
        var_info = [ ['lon'  ,'lon'           ,1.          ,' '    ,' '       ],\
                     ['lat'  ,'lat'           ,1.          ,' '    ,' '       ],\
                     ['tp'   ,'PRATE_L1_Avg_1',1.          ,' '    ,'prate'   ],\
                     ['dswrf','DSWRF_L1_Avg_1',1.          ,' '    ,'dswsfc'  ],\
                     ['uswrf','USWRF_L1_Avg_1',1.          ,' '    ,'uswsfc'  ],\
                     ['t2m'  ,'TMP_L103'      ,cff_temp    ,' '    ,'tmp2m'   ],\
                     ['u10m' ,'U_GRD_L103'    ,1.          ,' '    ,'wnd10m'  ],\
                     ['v10m' ,'V_GRD_L103'    ,1.          ,' '    ,'wnd10m'  ],\
                     ['strd' ,'DLWRF_L1_Avg_1',1.          ,' '    ,'dlwsfc'  ],\
                     ['q'    ,'SPF_H_L103'    ,1.          ,' '    ,'q2m'     ],\
                   ]

    elif data_origin == 'formatted':
        var_info = [ ['sst'  ,'SST'           ,1.          ,' ' ],\
                     ['tp'   ,'TP'            ,1.          ,' ' ],\
                     ['ssr'  ,'SSR'           ,1.          ,' ' ],\
                     ['t2m'  ,'T2M'           ,1.          ,' ' ],\
                     ['u10m' ,'U10M'          ,1.          ,' ' ],\
                     ['v10m' ,'V10M'          ,1.          ,' ' ],\
                     ['strd' ,'STRD'          ,1.          ,' ' ],\
                     ['q'    ,'Q'             ,1.          ,' ' ],\
                     ['msl'  ,'MSL'           ,1.          ,' ' ] \
                   ]

    elif data_origin == 'output_MF':
        var_info = [ ['lon'  ,'lon'           ,1.          ,' '    ,' '      ],\
                     ['lat'  ,'lat'           ,1.          ,' '    ,' '      ],\
                     ['sst'  ,'ts'            ,1.          ,' '    ,'ts'     ],\
                     ['tp'   ,'pr'            ,cff_tp2     ,' '    ,'pr'     ],\
                     ['dswrf','rsds'          ,1.          ,' '    ,'rsds'   ],\
                     ['uswrf','rsus'          ,1.          ,' '    ,'rsus'   ],\
                     ['t2m'  ,'tas'           ,cff_temp    ,' '    ,'tas'    ],\
                     ['u10m' ,'uas'           ,1.          ,' '    ,'uas'    ],\
                     ['v10m' ,'vas'           ,1.          ,' '    ,'vas'    ],\
                     ['strd' ,'rlds'          ,1.          ,' '    ,'rlds'   ],\
                     ['q'    ,'hus950'        ,1.          ,' '    ,'hus950' ],\
                   ]
    elif data_origin == 'arome':
        var_info = [ ['tp'   ,'rain'          ,1.          ,' '   ],\
                     ['ssr'  ,'swhf'          ,1.          ,' '   ],\
                     ['t2m'  ,'t2m'           ,cff_temp    ,' '   ],\
                     ['u10m' ,'u10m'          ,1.          ,' '   ],\
                     ['v10m' ,'v10m'          ,1.          ,' '   ],\
                     ['strd' ,'lwhf'          ,1.          ,' '   ],\
                     ['r'    ,'rh'            ,1.          ,' '   ],\
                     ['msl'  ,'pmer'          ,1.          ,' '   ] \
                   ]
    elif data_origin == 'ncep':
        var_info = [ ['tp'   ,'pratesfc'      ,1.          ,' '   ],\
                     ['dswrf','dswrfsfc'      ,1.          ,' '   ],\
                     ['uswrf','uswrfsfc'      ,1.          ,' '   ],\
                     ['t2m'  ,'tmp2m'         ,cff_temp    ,' '   ],\
                     ['strd' ,'dlwrfsfc'      ,1.          ,' '   ],\
                     ['u10m' ,'ugrd10m'       ,1.          ,' '   ],\
                     ['v10m' ,'vgrd10m'       ,1.          ,' '   ],\
                     ['q'    ,'spfh2m'        ,1.          ,' '   ],\
                   ]
    elif data_origin == 'ncep_ftp':
        var_info = [ ['tp'   ,'prate'         ,1.           ,' '   ,'prate '  ],\
                     ['dswrf','dswrf'         ,1.           ,' '   ,'dswrf'   ],\
                     ['uswrf','uswrf'         ,1.           ,' '   ,'uswrf'   ],\
                     ['t2m'  ,'air'           ,cff_temp     ,' '   ,'air.2m'  ],\
                     ['strd' ,'dlwrf'         ,1.           ,' '   ,'dlwrf'   ],\
                     ['u10m' ,'uwnd'          ,1.           ,' '   ,'uwnd'    ],\
                     ['v10m' ,'vwnd'          ,1.           ,' '   ,'vwnd'    ],\
                     ['q'    ,'shum'          ,1.           ,' '   ,'shum'    ],\
                   ]


    else:
        sys.exit(print('No \'%s\' dico available in Modules/aforc_reader.py. Please add it an run the script' % data_origin))

    return var_info

