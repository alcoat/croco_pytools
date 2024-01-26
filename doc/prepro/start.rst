=====
Start
=====

croco_pytools
-------------

croco_pytools are available on gitlab: https://gitlab.inria.fr/croco-ocean/croco_pytools
or can be downloaded from `CROCO website <https://www.croco-ocean.org/download/>`_.

Dependencies
------------

Python
======

croco_pytools is using the following package:

.. list-table::
   :widths: 10 90
  
   * - `matplotlib <https://matplotlib.org/>`_
     - Matplotlib is a comprehensive library for creating static, animated,
       and interactive visualizations in Python.
   * - `cartopy <https://scitools.org.uk/cartopy/docs/latest/>`_
     - Cartopy is a Python package designed for geospatial data processing 
       in order to produce maps and other geospatial data analyses.   
   * - `wxpython <https://wxpython.org/>`_
     - wxPython is a cross-platform GUI toolkit for the Python 
       programming language. 
   * - `geopandas <https://geopandas.org/en/stable/>`_
     - GeoPandas is an open source project to make working with geospatial 
       data in python easier.
   * - `regionmask <https://regionmask.readthedocs.io/en/stable/>`_
     - Create masks of geographical regions
   * - `pyinterp <https://pangeo-pyinterp.readthedocs.io/en/latest/>`_
     - Python library for optimized geo-referenced interpolation.
   * - `pandas <https://pandas.pydata.org/>`_
     - pandas is a fast, powerful, flexible and easy to use open source
       data analysis and manipulation tool, built on top of the
       Python programming language.
   * - `scipy <https://www.scipy.org/scipylib/index.html>`_
     - Scipy provides many user-friendly and efficient numerical routines,
       such as routines for numerical integration, interpolation,
       optimization, linear algebra, and statistics.
   * - `xarray <http://xarray.pydata.org/en/stable/>`_
     - xarray is an open source project and Python package that makes working
       with labelled multi-dimensional arrays simple, efficient, and fun

Fortran
=======

The scripts used to build the CROCO grid and the initial/boundary 
conditions for nests are using fortran routines which have been 
interfaced with python through `f2py <https://numpy.org/doc/stable/f2py/>`_
and must therefore be compiled.

Before compiling make sure that you have the following:

* Open MP-capable Fortran compiler, Ifort or Gfortran, may be others.
* netCDF library capable of handling netCDF-4/hdf5 format.
* working "ncview" (optional, but very useful).

