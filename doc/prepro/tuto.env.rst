Dedicated environment
=====================

Installation and compilation
----------------------------

* Get the source code as described here :doc:`start`

.. note:: 
  
  Be sure to have conda install on your computer. If not, 
  refer to `Conda website <https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_ to install it.

* If necessary edit `__init__.py` to change compilers or netcdf library to use.

.. note:: 
  
  On datarmor you can use gfortan and gcc for compilers and use 
  the library `NetCDF/4.4.1.1__gcc-6.3.0__nop`

* Execute the installation script ::

    python install.py

* Answer the question on the screen

  ::

    Which venv do you want to install (leave empty for no installation)? 
    [tools,doc,full]: full
    Do you want to compile fortran tools? 
    y,[n]: y


If no error is raised your environment and compilation were successful.
Now you can load your conda environment

::

  conda activate croco_pyenv
 
That's it, you are ready to use the preprocessing tools!
Remind to activate this envrironment each time you want to use the preprocessing tools.


Common errors
-------------

You might face some errors while trying to compile fortran tools. Here is a list 
of what have been already encountered by some users and th associated solution.

* ifort can raise wn error while compiling. In **__init__.py** try to add "--fcompiler=intelem" in **ENV_FFLAGS**.

* ImportError means you have missing librairies. In your terminal do nf-config --flibs 
  and check that you have -lnetcdff -lnetcdf. 

