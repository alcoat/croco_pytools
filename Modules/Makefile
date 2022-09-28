# Makedefs used
include ../Makedefs

# compile flags
CPPFLAGS = -I$(NETCDFINC)
FCFLAGS = -c -I$(NETCDFINC) -I$(PWD)/tools_fort_routines/  --f90flags="-g -check all -CA -CB -CS -fopenmp" -lgomp

# lib flags
LDFC = $(NETCDFLIB)

# program name
PROGRAM = toolsf
MOD= $(PWD)/tools_fort_routines/*.o
# source files and objects

   SRCS = tools_fort.F

all:$(RCS)
	$(CXX) $(SRCS) tools_fort_.f
	cd tools_fort_routines; $(CC) $(CPPFLAGS) -DFLIP -c read_gshhs.c -fPIC
	cd tools_fort_routines; $(F90) -c module_file.F -fPIC;
	cd tools_fort_routines; $(F90) -c spln2d.F -fPIC;
	$(FC) $(FCFLAGS) $(MOD) -m $(PROGRAM) tools_fort_.f $(LDFC)


clean:
	cd tools_fort_routines; rm -f *.o *.mod
	rm -f *.f *.so  *.o *.mod

