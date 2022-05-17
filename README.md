XCROCO
====== 




Installation
=============

Install miniconda:
Download Miniconda3 (i.e. for python3) from the [conda website](https://conda.io/miniconda.html) and run:
```
./Miniconda3-latest-Linux-x86_64.sh
```

Download the repository:
```
git clone https://gitlab.inria.fr/croco-ocean/croco_tools.git
```

Install an appropriate conda-environment:
```
conda update -y conda
conda create -n xcroco -c conda-forge -y python=3.10 
conda activate xcroco
conda install -y -c conda-forge dask dask-jobqueue \
            xarray zarr netcdf4 jupyterlab ipywidgets cartopy \
            geopandas nodejs intake-xarray xgcm numba
cd croco-tools/xcroco; pip install -e .

jupyter labextension install @jupyter-widgets/jupyterlab-manager \
                             @pyviz/jupyterlab_pyviz \
                             jupyter-leaflet
```

see also [conda doc](doc/conda.md)


Tutorials:
=========
You can find tutorials in the directory:
. tuto_xcroco.ipynb : several tools/diagnostics available in the library
. tuto_movie.ipynb : a way to make a movie