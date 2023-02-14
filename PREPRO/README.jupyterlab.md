
# With docker

```bash
# Import the image from .tar to docker. This could take a moment.
docker load -i croco_python_tools.tar.gz
# Check if import is good. You should see the new image with the name croco_python_tools
docker image ls
# Create the volume named to expose Python_files from the docker container to your OS.
mkdir Python_tools
docker volume create --driver local --opt type=none --opt o=bind --opt device=$(pwd)/Python_tools python_tools_volume
# Run the docker
docker run -it -p 8888:8888 --name croco_python_tools -v python_tools_volume:/root/Python_tools -v $(pwd)/data:/root/data croco_python_tools bash
# Start Jupyter lab
cd Python_tools
jupyter lab --ip=0.0.0.0 --no-browser --allow-root -e --ServerApp.token="YOUR_PASSWORD"

# Now open a browser and go to 
http://127.0.0.1:8888
```

If you are working on server you need to create a ssh tunnel.

```bash
# 8888 the first port is the port on the server opened by jupyter lab and 8080 is the port on your computer
ssh -R 8888:localhost:8080 public.example.com

# Now open a browser and go to 
http://127.0.0.1:8080
```

## Without docker

Do README.INSTALL

```bash
# Start jupyterlab
jupyter lab
```

# How to change password

```bash
jupyter lab password
```
