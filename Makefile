

all: 
	cd Modules; $(MAKE);

clean:
	cd Modules; rm -f *.f *.o *.mod *.so
	cd Modules/tools_fort_routines; rm -f *.f *.o *.mod *.so

docker_image:
	docker build -f docker/Dockerfile --compress --platform linux/amd64 --load -t croco_python_tools --progress=plain .