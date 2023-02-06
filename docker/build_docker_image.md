# Create a docker file with python tools

You only do this if you want to create a new docker image.
It's append for exemple when you want to share new developpement of python tools to the community.

## Install new dependencie

If you install new dependencie(s), you need to update environment.yml and conda-lock.yml file.

Requirements :

- conda
- conda-lock

If you want to export all your package installed by history use :

```bash
conda env export --from-history
```

This export an environment file with only the package specifications you have explicitly requested on the command line.

It's realy better than use traditional

```bash
conda env export
```

Because last command export dependencies with hash that could not be suitable for all environnement & architecture
(linux, macos ...).

After we want to generate lock file that ensure us to install the exact same dependencies version and speed up the process.

First install conda-lock

```bash
conda install -c conda-forge conda-lock
```

Use to generate lock file for dependencies for all your env.

```bash
# To export only for linux-64 linux-arm
conda-lock -p linux-64 -p linux-aarch64 -f environment.yml
```

## Build Docker image

You can build it with buildx or normal build docker command.

### Build command

```bash
docker build -f docker/Dockerfile --compress --platform linux/amd64 --load -t croco_python_tools --progress=plain .

# Save the image to tar file
docker save croco_python_tools:latest | gzip > croco_python_tools.tar.gz
```

### Buildx command

Build the image with docker buildx and BuildKit engine enable

```bash
DOCKER_BUILDKIT=1 docker buildx build -f docker/Dockerfile -t croco_python_tools --progress=plain .
docker buildx ls
docker buildx create --use --platform linux/amd64
# For linux-64 and linux-arm
# DOCKER_BUILDKIT=1 docker buildx build -f docker/Dockerfile --compress --platform linux/arm64/v8,linux/amd64 --output type=tar,dest=croco_python_tools.tar --progress=plain .

# Save the image to tar file
docker save croco_python_tools:latest | gzip > croco_python_tools.tar.gz
```

## PS

I do not manage to compile with conda dependencies

```yaml
# environment.yml
dependencies: ...
  - gcc_linux-64
  - gxx_linux-64
  - gfortran_linux-64
  - make
  ...
```

f2py try to compile with gcc instead of gfortran.
When you use this be sure to put good path to compilator into
your make file exemple x86_64-conda_cos6-linux-gnu-gcc.
[see more info here](https://docs.conda.io/projects/conda-build/en/latest/resources/compiler-tools.html)
