#
# default makefile for ifort compiler with more or less
# appropriate options for debugging and high performance
#

%.o: %.mod

# application name
APP = m3dfft2

# list of source files
SRC = model.f90 nrtype.f90 nrutil.f90 recipes.f90 atompot.f90 mkl_dfti.f90 gfx.f90 ft_filters_3d.f90 m3dfft.f90

# list of object files
OBJ = model.o nrtype.o nrutil.o recipes.o atompot.o mkl_dfti.o gfx.o ft_filters_3d.o m3dfft.o

# list of include directory
DIR1 = /export/home/group/mkl/10.0.3.020/include
DIR2 = /export/home/group/mkl/10.0.3.020/examples/vslf/source

# define libraries needed by the linker
LIBS = -L/export/home/group/mkl/10.0.3.020/lib/em64t /export/home/group/mkl/10.0.3.020/lib/em64t/libmkl_intel_lp64.a -Wl,--start-group /export/home/group/mkl/10.0.3.020/lib/em64t/libmkl_intel_thread.a /export/home/group/mkl/10.0.3.020/lib/em64t/libmkl_core.a -Wl,--end-group /export/home/group/mkl/10.0.3.020/lib/em64t/libguide.a -lpthread

# compiler options for debugging
FC_DEBUG = ifort -g -debug all -check all -implicitnone -warn all

# compiler options for optmized running
FC_OPT = ifort -O3 -xO -ipo -no-prec-div -static


# build rules

.SUFFIXES: .f90 .o .mod
.f90.o:
	${FC_DEBUG} -c $<

debug: ${OBJ} ${SRC}
	${FC_DEBUG} -o ${APP} ${OBJ} ${LIBS}

opt: ${SRC}
	${FC_OPT} -o ${APP} ${SRC} ${LIBS}

clean:
	rm -f *.mod *.o ${APP}

