# Croco python tools docker image

## Requirement

For all OS you need docker ;)

### Mac

You need [XQuartz](https://www.xquartz.org/).
Run the XQuartz application, then in it's application menu you'll find a Preferences choice.
In the Security tab check "authenticate connection" and "allow connection from network clients".
Restart XQuartz

### Windows

Never tested.
You need an X11 display server for Windows.

## Load the image

If you got the image from a tar file, first you need to load it

```bash
# Import the image from .tar to docker. This could take a moment.
docker load -i PATH_TO_THE_IMAGE
# Check if import is good. You should see the new image with the name croco_python_tools
docker image ls
```

## Start the docker

For running a docker image you have two ways. One _easier_ with a config file (see docker-compose section).
Other one with big commands (see docker run section).

NB: With docker-compose or docker run we are using an hack in order to see/edit files that are inside the docker
(exemple make_grid.py) from your OS. A better approch is going inside the container and directly modify files with "vi".
Another approch it to use codium or vscode (codium is better because is vscode compiled with privacy)
and install Docker plugin (from ms-azuretools).

### Linux

#### Docker run (use this section or Docker compose section)

```bash
# Create new directories
mkdir -p python_tools/Python_tools
cd python_tools
mkdir data
# Check if your are in directory with 2 folders (Python_tools and data) and 1 file docker-compose.linux.yml
ls -lah
# Create the volume named to expose Python_files from the docker container to your OS.
docker volume create --driver local --opt type=none --opt o=bind --opt device=$(pwd)/Python_tools python_tools_volume
# We need to export HOSTNAME for the GUI.
export HOSTNAME=$(hostname)
# Check if $HOSTNAME is not empty
echo $HOSTNAME
# Run the docker
docker run -it --name croco_python_tools -v /tmp/.X11-unix:/tmp/.X11-unix:ro -v $HOME/.Xauthority:/root/.Xauthority -v python_tools_volume:/root/Python_tools -v $(pwd)/data:/root/data -h $HOSTNAME -e DISPLAY=$DISPLAY croco_python_tools bash
```

### Docker compose (use this section or Docker run)

```bash
# Create new directories
mkdir -p python_tools/Python_tools
cd python_tools
mkdir data
# Check if your are in directory with 2 folders (Python_tools and data) and 1 file docker-compose.linux.yml
ls -lah
# We need to export HOSTNAME for the GUI.
export HOSTNAME=$(hostname)
# Check if $HOSTNAME is not empty
echo $HOSTNAME
# Start the docker
docker-compose --file docker-compose.linux.yml up -d --force-recreate
# Now you can enter in the terminal of the docker
docker-compose --file docker-compose.linux.yml exec croco_python_tools bash
```

If you do not want to write `--file docker-compose.linux.yml` for each `docker-compose` just rename
`docker-compose.linux.yml` to `docker-compose.yml`.

### Troubleshooting

If when you try to run make_grid.py. You have error like **Authorization required, but no authorization protocol specified. Unable to access the X Display, is $DISPLAY set properly?**

It's could be $DISPLAY is not correct. it should be something like : `:0`
It's could be $HOSTNAME is not correct.

You could use at your **own risk**. It's better to try to fix $DISPLAY, $HOSTNAME and Xauthority. The command `xhost +` to turns off access control to your X11 server, so that anyone who can connect
to the TCP any network port or socket that X11 listens on has full control of your display, and can record every keyboard press, or record your display, etc.

**After you use the tools do not forget to run `xhost -`**

#### Mac os

### Docker run (use this section or Docker-compose)

```bash
# Check if your are in directory with 1 folders (Python_tools) and 1 file docker-compose.mac.yml
ls -lah
# We need to find your IP on local network. Maybe you need to change en0 with your interface
export IP=$(ifconfig en0 | grep inet | grep -Eow '([0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3})')
# Check if $IP is good
echo $IP

# Allow your ip to connect to xhost
xhost +$IP
# Create the volume named to expose Python_files from the docker container to your OS.
docker volume create --driver local --opt type=none --opt o=bind --opt device=$(pwd)/Python_tools python_tools_volume
# Run the docker
docker run -it --name croco_python_tools -v /tmp/.X11-unix:/tmp/.X11-unix:ro -v $HOME/.Xauthority:/root/.Xauthority -v python_tools_volume:/root/Python_tools -v $(pwd)/data:/root/data -e DISPLAY=$IP:0 croco_python_tools bash
# When you finished remove the ip to xhost
xhost -$IP
```

#### Docker-compose (use this section or Docker run section)

```bash
# Create new directories
mkdir -p python_tools/Python_tools
cd python_tools
mkdir data
# Check if your are in directory with 2 folders (Python_tools and data) and 1 file docker-compose.linux.yml
ls -lah
# We need to find your IP on local network. Maybe you need to change en0 with your interface
export IP=$(ifconfig en0 | grep inet | grep -Eow '([0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3})')
# Check if $IP is good
echo IP=$IP
# Allow your ip to connect to xhost
xhost +$IP
# Start the docker
docker-compose --file docker-compose.mac.yml up -d --force-recreate
# Now you can enter in the terminal of the docker
docker-compose --file docker-compose.mac.yml exec croco_python_tools bash
# When you finished remove the ip to xhost
xhost -$IP
```

### Troubleshooting

If you have **"command not found: xhost"** normaly xhost should be on `opt/X11/bin/` directory.

If when you try to run make_grid.py. You have error like **Authorization required, but no authorization protocol specified. Unable to access the X Display, is $DISPLAY set properly?**

It's could be $IP is not correct. Try this command to get it. It should be something like 192.168.1.100

```bash
IP=$(ifconfig en0 | grep inet | awk '$1=="inet" {print $2}')
```

Check if XQuartz run on `:0` display. If not change the docker run command

```
display_number=`ps -ef | grep "Xquartz :\d" | grep -v xinit | awk '{ print $9; }'`
```

Maybe you have to add the `--privileged` flag to the docker run command.

Maybe you have to add the `--net host` flag to the docker run command.

If you get an error: XDG_RUNTIME_DIR not set in the environment add this to the Docker command line: -e XDG_RUNTIME_DIR=/tmp/xdgr

You can try

```
#  adds localhost to the list of authorized clients, allowing only local processes to connect
xhost +localhost
docker run -it --name croco_python_tools -v python_tools_volume:/root/Python_tools -v $(pwd)/data:/root/data -e DISPLAY=host.docker.internal:0 croco_python_tools bash
```

The xhost command configures XQuartz to allow connections from the given IP address.
You could use at your **own risk**. It's better to try to fix $DISPLAY, $HOSTNAME and Xauthority. The command `xhost +` to turns off access control to your X11 server, so that anyone who can connect
to the TCP any network port or socket that X11 listens on has full control of your display, and can record every keyboard press, or record your display, etc.

**After you use the tools do not forget to run `xhost -`**
