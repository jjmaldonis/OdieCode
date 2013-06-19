// function headers for FemMulti

// Function headers for femsim
//
// v3.0 08-24-05 pmv

#ifndef FemSim
#define FemSim
#include <stdio.h>

// StartStop.c: initialization and shutdown functions
void initialize(char *sample, char *basename);  // reads the sample files, 
                                                // writes intial data to the log file
void getscatfunc(void);  // calculates the actual atomic scattering factors 
                         // from their parametrizations in the ScatFact array.
void output(char *basename); // write results of the calculation to various files.
void cleanup();   // free allocated memory, close files, etc.
void proglog(char *s); // write to stdout and program log file

// ReadInputParameters.c: Parse the input parameters file
int ReadParameters(FILE *f);
int GetLine(FILE *f, char *s);


// Memory.c: dynamic memory allocation functions
// allocates memory for intensity / variance calculations depending on algorithm
void allocate(void);
void alloc_sum(void);
void alloc_g2(void);
float *float1D(int n1, const char *message);
float **float2D(int n1, int n2, const char *message);
float ***float3D(int n1, int n2, int n3, const char *message);
float ****float4D(int n1, int n2, int n3, int n4, const char *message);
int *int1D( int n, const char *message );


// Sample.c: functions for dealing with sample
void ReadSampleAtoms(char sampl_path[]);
void CenterZero();
void rotate_matrix(float phi,float theta,float psi,float *a);
void rotation(float r[3][3]);
void grow(void);
void cut(void);


// Variance.c: sample rotation & variance bookkeeping
void appn_I_I2(void);
void rota_aver_I(void);
void variance(void);


// IntensitySum.c: calculating intensity by direct summation
void atomcyli(void);
void intensitysum(void);


//IntensityG2.c: calculate intensity by Keblinski g2(r, r')
void cleanG2(void);
void atomcorr(void);
void intensityG2(void);


// Pixels.c: calculate the pixel center positions by various methods
void init_pixels();
void PixelsSquareArray();
void PixelsHexagonalArray(); // not implemented
void PixelsReadFile();


// Rotations.c: calculate the model rotation angles by various methods
void init_rotations();
void RotationsUniform();


// MathFunc.s: special functions
float bessj1(float x);
float a_q(float s);
float bessj0(float x);

#endif
