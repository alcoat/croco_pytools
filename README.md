# CROCO_PYTOOLS - DEVELOPMENT BRANCH (not stable)

This is a python redesign of the matlab pre- and post-processing tools for CROCO
(croco_tools) developed over the last 20 years.

It's a work in progress, which needs to be improved and further tested for some
of its features.

We provide separate DOIs (Digital Object Identifier) for the static versions of
the code at https://doi.org/10.5281/zenodo.11047402

You can also subscribe at CROCO forum at https://forum.croco-ocean.org

## Reminders

These scripts come with no warranty, and are provided on a free community
development basis.

## Contents:

The repository consists of several parts :

- **_prepro_**: a set of python routines, interface with fortran to process the
  grid and forcing files (initialization, open-boundary conditions, tides
  forcing, rivers forcing)
- **_croco_pyvisu_**: a vizualisation GUI
- **_xcroco_**: for analysis based on xarray and xgcm

#

### prepro

**croco_pytools/prepro** remains in the footsteps of the matlab **croco_tools**
with the following teps for creating a configuration:

- `make_grid`: build the grid, its mask and bathymetry (with chosen smoothing)
- `make_bry`: build the lateral boundary conditions (3D currents, temperature
  and salinity, barotropic currents, surface elevation)
- `make_ini`: build the initial conditions

And eventually:

- `make_tides`: build tidal forcing
- `make_rivers`: build river forcing

No routine is available for managing surface forcing at the moment, and the
users are suggested to use the `ONLINE` cppkey in CROCO, which performs the
surface forcing interpolation, directly during the run. See the usual
croco_tools tutorial, and tools for downloading and formating the data.

These pre-processing tools are taking advantage of several fortran routines
created by Alexander Shchepetkin, which are interfaced within python with f2py.

**Readers** are available for the various type of data, to ease the use and
addition of different type of input data.

**croco_pytools/prepro** installation required to set-up a python environment
(see prepro/env/environment_tools.yml), and compilation of fortran routines
interfaced with python routines. Scripts for easy installation are provided:

- general script: prepro/install.py
- script for conda environment: prepro/env/conda_install.py

See the full
[documentation](https://croco-ocean.gitlabpages.inria.fr/croco_pytools/prepro)
for more details.

#

### croco_pyvisu

Croco_pyvisu is a GUI tool written in python to visualize history files
genarated by the Croco model.

The croco_pyvisu directory is provided by the CROCO_PYTOOLS.  
See the
[documentation](https://croco-ocean.gitlabpages.inria.fr/croco_pytools/croco_pyvisu)
for the installation and the way to use it.

#

### xcroco

Xcroco is a library written in python to analyse history files genarated by the
Croco model.

The Xcroco directory is provided by the CROCO_PYTOOLS. See
[documentation](https://croco-ocean.gitlabpages.inria.fr/croco_pytools/xcroco)
for the installation and the way to use it.
