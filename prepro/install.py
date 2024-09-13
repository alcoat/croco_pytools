#!/usr/bin python3
# -*- coding: utf-8 -*-
"""
Main croco_pytools install execuable script
"""


import os
import subprocess
import shutil
from __init__ import (
    ENV_MOD,
    ENV_CONDA,
)


def install_conda_environment(venv_conda):
    # Check if mamba is available
    mamba_available = shutil.which("mamba") is not None
    
    # Determine the installation command
    if mamba_available:
        print("Mamba is available. Using mamba to create the environment.")
        install_cmd = f"mamba env create -q -f {venv_conda}"
    else:
        print("Mamba is not available. Using conda to create the environment.")
        install_cmd = f"conda env create -q -f {venv_conda}"
    
    # Attempt to install the environment
    try:
        subprocess.run(install_cmd.split(), check=True)
        print(f"Environment installation succeeded using {'mamba' if mamba_available else 'conda'}.")
    except subprocess.CalledProcessError as e:
        print(f"Failed to install environment from {venv_conda}.")
        print(e)

# Function to check if a conda environment exists
def conda_env_exists(env_name):
    try:
        result = subprocess.run(['conda', 'env', 'list'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return env_name in result.stdout
    except subprocess.CalledProcessError:
        return False


# Ask if the user wants to install a conda environment
install_venv_conda = input("Do you want to install conda environment? [y,[n]]: ").lower()
if install_venv_conda in ['y', 'yes']:
    venv_conda = ''.join((ENV_CONDA,'/environment_tools.yml'))
    install_conda_environment(venv_conda)
else:
    print("Skipping environment installation.")

# Ask if the user wants to compile Fortran tools
fortran_compilation = input("Do you want to compile Fortran tools? [y,[n]]: ").lower()
if fortran_compilation in ['y', 'yes']:
    print("Compiling Fortran tools...")
    tools_path = os.path.join(ENV_MOD, 'tools_fort_routines')
    env_fine = 'croco_pyenv'
    comp_cmd = None
    if conda_env_exists(env_fine) :
        print(f"Using conda environment {env_fine} to compile Fortran tools.")
        comp_cmd = f"conda run -n {env_fine} make"
    else:
        print(f"Conda environment {env_fine} not found.")
        try:
            fc = subprocess.run(["nc-config --fc"],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            shell=True,
                            check=True)
            if "gfortran" in str(fc.stdout):
                print('Found a NETCDF library compiled with gfortran. Trying to compile tools')
                comp_cmd = "make"
            else:
                print("Found NETCDF library but not compiled with gfortran. Fortran tools can not be compiled")
            
        except subprocess.CalledProcessError:
            print('No NETCDF library found, fortran tools can not be compiled')
              
    if comp_cmd is not None:
        # Change to the tools directory
        os.chdir(tools_path)
 
        # Clean previous build
        subprocess.run(["make", "clean"], check=True)
        # Compile Fortran tools
        try:
            log = subprocess.run(comp_cmd, 
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 text=True,
                                 shell=True,
                                 check=True)
            print(log.stdout)
            print(f"Compilation was a success")
        except subprocess.CalledProcessError:
            print(f"Compilation failed")
else:
    print("Skipping Fortran tools compilation.")


