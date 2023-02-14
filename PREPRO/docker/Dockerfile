# syntax=docker/dockerfile:1.2
# this syntax is understand only by the new 
# BuildKit the new build engine shipped with Docker
# this engine is handle caching.
# For exemple you do not need to dl apt dep or conda dep
# at all buildtime

# Container for building the environment
# It condaforge/mambaforge use minial Ubuntu install
# We are using this env to build the code with f2py
FROM condaforge/mambaforge:22.9.0-3 as conda

# Here define witch conda environment file you want to use
# exemple : conda-lock.yml or environment.yml
ARG CONDA_ENV_FILE=environment.yml
ENV CONDA_ENV_FILE=environment.yml

ARG HOME=/root
ENV HOME=$HOME

ENV PYTHON_TOOLS_DIR=$HOME/Python_tools
ENV SRC_DIR=$PYTHON_TOOLS_DIR

# Warning With dockerfile you cannot copy file from parent directory
# This path is relative where you do the docker build command.
ARG HOST_SRC_DIR=./

WORKDIR $HOME

RUN apt update && \
	apt install -yqq --no-install-recommends \
    # General tools
    vim nano \
    # Compil tools + netcd
    libnetcdf-dev libnetcdff-dev ncview netcdf-bin nco \
    gfortran gcc cpp make \
    # GUI stuff
    libxxf86vm1 xauth && \
    # Removing the contents of the /var/lib/apt/lists directory usually saves around 30 megabytes from the final image
    rm -rf /var/lib/apt/lists/*

# Copy all files and folder inside the $SRC_DIR of the container directory
COPY $HOST_SRC_DIR $PYTHON_TOOLS_DIR

# RUN ls -lah $SRC_DIR/

# Install conda env
# use mamba for performance
RUN mamba env create --file $SRC_DIR/$CONDA_ENV_FILE && conda clean -afy

RUN mkdir $HOME/data

# Compile fortran with f2py.
RUN cd $SRC_DIR && mamba run -p /opt/conda/envs/envcroco make

# Auto activate conda env when login.
RUN conda init bash && echo "conda activate envcroco" >> $HOME/.bashrc

USER root
WORKDIR $HOME
CMD ["bash"]