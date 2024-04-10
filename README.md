# croco_pytools - v1.0

This is a python redesign of the matlab pre- and post-processing tools for CROCO (croco_tools) developed over the last 20 years.

It's a work in progress, which needs to be improved and further tested for some of its features.

## Release notes: 


## GIT checkout:

`git clone --branch release https://gitlab.inria.fr/croco-ocean/croco_pytools.git`

## Reminders

These scripts come with no warranty, and are provided on a free community development basis.

## Contents:
The repository consists of several parts :
  - ***PREPO***: a set of python routines, interface with fortran to process the grid and forcing files (initialization, open-boundary conditions, atmospheric forcing, tides forcing, rivers forcing)
  - ***croco_pyvisu***: a vizualisation GUI
  - ***xcroco***: for analysis based on xarray and xgcm

#
### PREPRO
See README.INSTALL for installation of the required python environement (using conda) and compilation of python librairies interfaced with fortran

See README.topo for grid and coastline management

- ***to be continued*** 

#
### croco_pyvisu
Croco_pyvisu is a  GUI tool written in python to visualize history files genarated by the Croco model.

The croco_pyvisu directory is provided by the CROCO_PYTOOLS.  
See the [documentation](https://gitlab.inria.fr/croco-ocean/croco_pytools/-/blob/release/doc/croco_pyvisu) for the installation and the way to use it.


#
### xcroco
Xcroco is a library written in python to analyse history files genarated by the Croco model.

The Xcroco directory is provided by the CROCO_PYTOOLS.  
For the installation, see [Xcroco main page](https://gitlab.inria.fr/croco-ocean/croco_pytools/-/blob/release/xcroco)  
For the content of the library, see the [documentation](https://gitlab.inria.fr/croco-ocean/croco_pytools/-/blob/release/doc/xcroco)

