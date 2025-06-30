#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main croco_pytools install executable script
"""

import os
import sys
import subprocess
import shutil
import argparse
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
        print(
            f"Environment installation succeeded using {'mamba' if mamba_available else 'conda'}."
        )
    except subprocess.CalledProcessError as e:
        print(f"Failed to install environment from {venv_conda}.")
        print(e)
        sys.exit(1)


# Function to check if a conda environment exists
def conda_env_exists(env_name):
    try:
        result = subprocess.run(
            ["conda", "env", "list"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        return env_name in result.stdout
    except subprocess.CalledProcessError:
        return False


def compile_fortran_tools():
    """Compile Fortran tools using conda environment or system compiler."""
    print("Compiling Fortran tools...")
    tools_path = os.path.join(ENV_MOD, "tools_fort_routines")
    env_fine = "croco_pyenv"
    comp_cmd = None
    if conda_env_exists(env_fine):
        print(f"Using conda environment {env_fine} to compile Fortran tools.")
        comp_cmd = f"conda run -n {env_fine} make"
    else:
        print(f"Conda environment {env_fine} not found.")
        try:
            fc = subprocess.run(
                ["nc-config --fc"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                shell=True,
                check=True,
            )
            if "gfortran" in str(fc.stdout):
                print(
                    "Found a NETCDF library compiled with gfortran. Trying to compile tools"
                )
                comp_cmd = "make"
            else:
                print(
                    "Found NETCDF library but not compiled with gfortran. Fortran tools can not be compiled"
                )

        except subprocess.CalledProcessError:
            print("No NETCDF library found, fortran tools can not be compiled")

    if comp_cmd is not None:
        # Change to the tools directory
        os.chdir(tools_path)

        # Clean previous build
        subprocess.run(["make", "clean"], check=True)
        # Compile Fortran tools
        try:
            log = subprocess.run(
                comp_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                shell=True,
                check=True,
            )
            print(log.stdout)
            print("Compilation was a success")
        except subprocess.CalledProcessError as e:
            print("Compilation failed")
            print(e)
            sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Install croco_pytools and compile Fortran tools."
    )
    parser.add_argument(
        "--install-conda", action="store_true", help="Install conda environment"
    )
    parser.add_argument(
        "--compile-fortran", action="store_true", help="Compile Fortran tools"
    )
    args = parser.parse_args()

    if not args.install_conda and not args.compile_fortran:
        install_conda = (
            input("Do you want to install the conda environment? (yes/no): ")
            .strip()
            .lower()
            == "yes"
        )
        compile_fortran = (
            input("Do you want to compile Fortran tools? (yes/no): ").strip().lower()
            == "yes"
        )
    else:
        install_conda = args.install_conda
        compile_fortran = args.compile_fortran

    if install_conda:
        venv_conda = "".join((ENV_CONDA, "/environment_tools.yml"))
        install_conda_environment(venv_conda)
    else:
        print("Skipping environment installation.")

    if compile_fortran:
        compile_fortran_tools()
    else:
        print("Skipping Fortran tools compilation.")


if __name__ == "__main__":
    main()
