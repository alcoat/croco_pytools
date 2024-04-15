

import os
import sys

#: Directory of the python library
ENV_CP = os.path.dirname(__file__)

#: Modules dir
ENV_MOD = os.path.join(ENV_CP, 'Modules')

#: Conda env dir
ENV_CONDA = os.path.join(ENV_CP, 'env')


ENV_FCOMP_TPL = os.path.join(ENV_MOD,'tools_fort_routines/Makedefs.tpl')
ENV_CXX = 'cpp'
ENV_F90 = 'gfortran'
ENV_CC = 'gcc'

ENV_FFLAGS = '--f77flags="-std=legacy" --f90flags="-g -check all -CA -CB -CS -fopenmp" -lgomp'

ENV_NETCDFINC = '$(shell nf-config --includedir)'
ENV_NETCDFLIB = '-L$(shell nf-config --includedir)/../lib -lnetcdff -lnetcdf'

