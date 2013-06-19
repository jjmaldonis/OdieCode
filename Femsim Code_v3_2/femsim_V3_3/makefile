

.SUFFIXES: .o .c
.c.o:
	$(CC) -c $(CFLAGS) $< 


# Macros

#CC = condor_compile gcc
CC = gcc
#CFLAGS = -Wall -O3
CFLAGS = -Wall -O2
LIBS = -lm
SRC=FemSim.c MathFuncs.c StartStop.c Sample.c IntensitySum.c IntensityG2.c Variance.c Pixels.c Memory.c ReadInputParameters.c Rotations.c
OBJ=FemSim.o MathFuncs.o StartStop.o Sample.o IntensitySum.o IntensityG2.o Variance.o Pixels.o Memory.o ReadInputParameters.o Rotations.o


# Explicit rule

app: $(OBJ)
	$(CC) $(CFLAGS) -o femsim $(OBJ) $(LIBS)

clean:
	rm -r *.o femsim
