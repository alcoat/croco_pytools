CROCO_PYTOOLS preprocessing philosophy
######################################

CROCO_PYTOOLS/PREPRO remains in the same line as the matlab CROCO_TOOLS by separating all operation 
relative to creation of the simulation. To avoid having a too long namelist file, and each routines 
being independant, parameters need to be specified in the header of each routine.

Steps for creating a configuration (with these routines) are:

* build the grid
* build the lateral boundary conditions (3D currents, temperature and salinity, barotropic 
  currents, surface elevation)
* build the initial conditions

And eventually:

* build tides forcing
* build rivers field

Structure of directories
========================

* **doc**
  
  This is the Sphinx source code for this documentation

* **env**
  
  Contains **yml** files to install conda environment decidated to croco_pytools/PREPRO

* **Modules**

  Contains all python routines to run the croco_pytools/PREPRO

* **Readers**

  This contains **readers** to decode the input datasets


Detail of **PREPRO** routines
=============================

The main directory contains the following files:

.. list-table::
   
   * - __init__.py
     - File containing informations to build conda environment and compile fortran routines
   * - install.py
     - Routine to facilitate the installation of conda environnement and the compilation of fortran routines
   * - make_grid.py
     - Routines to build CROCO grid and associated nests
   * - make_bry.py
     - Build the latteral boundary conditions (surface elevation, 3D currents, barotropic currents, temperature and salinity, other tracers)
   * - make_ini.py
     - Build initial boundary condition
   * - make_zoom_cond.py
     - Build initial and boundary conditions for nested grid (offline nest or AGRIF nest)
   * - make_tides.py
     - Build tidal forcing (amplitude and phase) for elevation and barotropic current
   * - make_rivers.py
     - Create netcdf file containing runoff flows

   











