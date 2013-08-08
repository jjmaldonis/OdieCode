#
# default makefile for ifort compiler with more or less
# appropriate options for debugging and high performance
#

# application name
APP = rmc

# list of source files
SRC = ran2.f90 globals.f90 read_inputs_v2.f90 model_v2.f90 rmc_functions.f90 scattering_factors.f90 fft2dc.c slicelib.c islice.c fem1_v2.f90 gr.f90 eam.f90 rmc_v_1_0.f90 

# list of object files
OBJ = ran2.o   globals.o   read_inputs_v2.o   model_v2.o   rmc_functions.o   scattering_factors.o   fft2dc.o slicelib.o islice.o fem1_v2.o   gr.o   eam.o   rmc_v_1_0.o 

# define libraries needed by the linker
#LIBS = -L/export/apps/openmpi_intel_20130618/lib -lmpi
LIBS=
#LIBS = -L/export/apps/openmpi_intel_20130618/lib
# -rpath /state/partition1/apps/intel_20130618/composer_xe_2013.3.163/compiler/lib/

# compiler options for debugging
#FC_DEBUG = ifort -g -debug all -check all -implicitnone -warn all
#FC_DEBUG = mpif90 -g -debug all -check all -implicitnone -warn all
FC_DEBUG = mpif90 -g -debug all -check all -implicitnone -warn all -openmp -fpp
#FC_DEBUG = mpif90 -g -implicitnone -openmp -fpp

CC_DEBUG = icc -g  -Wall

# compiler options for optmized running
#FC_OPT = ifort -O3 -xO -ipo -no-prec-div -static
#FC_OPT = mpif90 -O3 -ipo -static
FC_OPT = mpif90 -O3

#CFLAGS=-I/export/apps/openmpi_intel_20130618/include/

# build rules

.SUFFIXES: .f90 .c .o
.f90.o:
	${FC_DEBUG} -c $<

.c.o:
	${CC_DEBUG} -c $<

debug: ${OBJ}
	${FC_DEBUG} -o ${APP} ${OBJ} ${LIBS}

opt: ${SRC}
	${FC_OPT} -o ${APP} ${SRC} ${LIBS}

clean:
	rm -f *.mod *.o ${APP}
