Build CROCO initial and boundary conditions
--------------------------------------------

This section presents the use of ``make_aforc.py``
to create atmospheric forcing for
realistic regional configurations with CROCO. 

The different steps to build atmospheric forcing are:

#. Select area of interest
#. Format metadata
#. Save data by variable on a daily or monthly files

Fill the Reader
^^^^^^^^^^^^^^^

In this tutorial, we will use as amtospheric data the 
`ERA5 Reanalysis (single levels) product <https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=overview>`_.

You can check the variables name from these input data by doing:
::

  ncdump -h filename.nc

You thus, first need to check/edit the atmospheric reader, ``Readers/aforc_reader.py``, for these data. 
Here we use the ``era_ecmwf`` style of formatting:

::

    elif data_origin == 'era_ecmwf':
        var_info = [ ['lon'  ,'longitude'     ,1.           ,' '      ],\
                     ['lat'  ,'latitude'      ,1.           ,' '      ],\
                     ['tp'   ,'tp'            ,cff_tp       ,'cumul'  ],\
                     ['ssr'  ,'ssr'           ,cff_heat     ,'cumul'  ],\
                     ['t2m'  ,'t2m'           ,cff_temp     ,' '      ],\
                     ['u10m' ,'u10'           ,1.           ,' '      ],\
                     ['v10m' ,'v10'           ,1.           ,' '      ],\
                     ['strd' ,'strd'          ,cff_heat     ,' '      ],\
                   ]


.. note:: 
    If the different variables are in separate files, you need to add how the variables are named to the file names. 

::

     elif data_origin == 'cfsr':
        var_info = [ ['lon'  ,'lon'           ,1.           ,' '      ,' '       ],\
                     ['lat'  ,'lat'           ,1.           ,' '      ,' '       ],\
                     ['tp'   ,'PRATE_L1_Avg_1',1.           ,'cumul'  ,'prate'   ],\
                     ['dswrf','DSWRF_L1_Avg_1',1.           ,'cumul'  ,'dswsfc'  ],\
                     ['uswrf','USWRF_L1_Avg_1',1.           ,'cumul'  ,'uswsfc'  ],\
                     ['t2m'  ,'TMP_L103'      ,cff_temp     ,' '      ,'tmp2m'   ],\
                     ['u10m' ,'U_GRD_L103'    ,1.           ,' '      ,'wnd10m'  ],\
                     ['v10m' ,'V_GRD_L103'    ,1.           ,' '      ,'wnd10m'  ],\
                     ['strd' ,'DLWRF_L1_Avg_1',1.           ,'cumul'  ,'dlwsfc'  ],\
                     ['q'    ,'SPF_H_L103'    ,1.           ,' '      ,'q2m'     ],\
                   ]

 
Using make_aforc
^^^^^^^^^^^^^^^^

Before running ``make_aforc.y``, **USER OPTIONS** section needs to be filled. 
There are several parts in it:


* **INPUT**:

::

  data_origin = 'era_ecmwf'
  input_dir = '/path/in/'
  input_prefix = 'era_5-copernicus__*'
  multi_files = False 

``data_origin`` input data dictionnary as defined in the Readers/aforc_reader.py.

``input_prefix`` for multifiles, if the name of the file begin with the variable name, write '*' before sufix.

``multi_files`` if one file per variable in input, set True.

* **OUTPUT**:

::

  output_dir = '/path/out/'
  output_file_format = "MONTHLY"

``output_file_format`` How output files are split (MONTHLY,DAILY)


* **Grid size**:

::

  ownArea = 0 # 0 if area from croco_grd.nc
              # 1 if own area
  if ownArea == 0:
      croco_grd = '/pathin/to/your/croco/grid/croco_grd.nc'
  else:
      lon_min,lon_max,lat_min,lat_max = 4,27,-40,-24


* **Dates limits**:

::

  Yorig = 1950                 # year defining the origin of time as: days since Yorig-01-01
  Ystart, Mstart = 1980,1   # Starting month
  Yend, Mend  = 1980,1  # Ending month

.. note:: 

    Origin time and initial time can be different.

* **Options**

::

  cumul_step = 1
  READ_PATM = False
  extrapolation_sst = True

``cumul_step`` data accumulation time in hours.

``READ_PATM`` to convert the atmospheric pressure, set True.

``extrapolation_sst`` is for the case there is no STRD variable in raw data. It will be calculated with STR and SST. SST may or may not be extrapolated to the coast. If it is not, this may result in temperature spike at the coast at certain points. However, note that extrapolation increases pre-processing time. If STRD is in raw data, extrapolation_sst will no be considered. If your variable is not the sea surface temperature but the surface temperature, no need to set True.


To use ``make_aforc.py``, do:
::

  python make_aforc.py
