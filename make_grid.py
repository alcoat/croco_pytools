__author__ = 'Mathieu Le Corre'
__email__ = 'mathieu.le.corre@shom.fr'
__date__ = '2022-09'
__license__ = 'GPL3'
'''
===========================================================================
Further Information:  
  http://www.croco-ocean.org
  
This file is part of CROCOTOOLS

Create a CROCO grid file
In the current state the script can handle:
    - etopo (5,2,1)
    - srtm30
    - gebco

To add a new dataset you just have to go in Modules/topo_readers.py and
create a dico with name of lon,lat,topo in the dataset.
At this time it only handle 1D lon/lat 

The script makes a grid with mercator projection centred at the equator
before rotating the sphere to put the middle of the grid in tra_lon/tra_lat
position.

Then, it reads topo's dataset and apply the desired smoothing
The mask is generated using a shapefile (.shp. default is gshhs dataset)

Smoothing use fortran routines for better performance
'''


if __name__ == "__main__":

    Question = input("Do you want to use interactive grid maker ? \
                      \n (e.g., for grid rotation or parameter adjustments) : y,[n] ")

    if Question.lower() == ("y") or Question.lower() == ("yes"):

        # --- Building grid with graphicUI -------------------------------------
        print("In interactive mode")

        # Auto set ETS_TOOLKIT environnement var
        import os
        os.environ["ETS_TOOLKIT"] = "wx"

        from Modules.graphicUI_tools.main_window import MainWindow
        MainWindow().configure_traits()

    elif Question.lower() == ("n") or Question.lower() == ("no") or Question.lower() == (""):

        # --- Building grid without graphicUI ----------------------------------
        print("In normal mode")

        import make_grid_param as param
        from Modules.tools_make_grid import inputs, inputs_smth, EasyGrid, GetTopo
        # --- Create inputs -----------------------------
        inputs = inputs(param.tra_lon, param.tra_lat, param.size_x,
                        param.size_y, param.nx, param.ny, param.rot)
        inputs_smth = inputs_smth(
            param.hmin, param.hmax, param.smth_rad, param.rfact, param.smooth_meth)

        # --- Create lon/lat grid -----------------------------------------
        outputs = EasyGrid.easygrid(None, inputs)

        # --- Build mask and topo -----------------------------------------
        GetTopo.topo(None, outputs, param.topofile, param.shp_file,
                     smooth=inputs_smth, sgl_connect=param.sgl_connect)

        # --- Save netcdf -------------------------------------------------
        print('Writing Topography')
        from Modules.croco_class import CROCO
        CROCO.create_grid_nc(None, param.output_file, inputs, outputs)
