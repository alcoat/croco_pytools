import os
import sys
import subprocess


THIS_DIR = os.path.dirname(__file__)
ROOT_DIR = os.path.dirname(THIS_DIR)

sys.path.insert(0, ROOT_DIR)
from __init__ import (
    ENV_MOD,
    ENV_CONDA,
    ENV_FCOMP_TPL,
    ENV_CXX,
    ENV_F90,
    ENV_CC,
    ENV_NETCDFINC,
    ENV_NETCDFLIB,
    ENV_FFLAGS
)


f = open(ENV_FCOMP_TPL,'r')
filedata = f.read()
f.close()

old_var = ['env_CXX','env_F90','env_CC','env_NETCDFINC','env_NETCDFLIB','env_FFLAGS']
new_var = [ENV_CXX,ENV_F90,ENV_CC,ENV_NETCDFINC,ENV_NETCDFLIB,ENV_FFLAGS]

for i,var in enumerate(old_var):
    filedata = filedata.replace(var,str(new_var[i]))

f = open(ENV_FCOMP_TPL[:-4],'w')
f.write(filedata)
f.close()

print('Compiling fortran tools')

os.chdir(ENV_MOD+'/tools_fort_routines')
subprocess.run(["make", "clean"], check=True)
boo = subprocess.run(["make"], check=True)
print(boo)
#except:
#    print('Compilation failled.... You won\'t be able to use make_grid and make_zoom_cond')

