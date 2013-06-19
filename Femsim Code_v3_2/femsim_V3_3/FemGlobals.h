// global variables and struct definitions for FemSim.c
// used to be called "FemData.h"
//
// changed fixed-size arrays to pointers for dynamic allocation
// removed "static" declarations on arrays 8-19-05 pmv
// adding two variables Start_Step and Stop_Step for only calculating a certain range of rotation
// by Feng Yi on 04/04/2008

#include <time.h>
#include <stdio.h>
#include <stdlib.h>

char version[64];

// simulation parameters
int Algorithm; // 1 = direct summation; 2 = g2(r, r'); 3 = multislice (not implemented)
float Q_aperture;
float K_mini, K_scal; 
int K_step;
float Delta_r;
int R_step;

// sample rotations parameters
int Do_rotations;
float Thet_ofst;
int Angl_step, phi_step, theta_step, psi_step, R_step;
int Start_Step, Stop_Step; //Begin rotation position and stop rotation position 04/04/2008
float *phi_range, *theta_range, *psi_range; // 1-D arrays to hold the range of angle values for rotations
float *phi, *theta, *psi; // 1-D arrays to hold the actual values for every orientation (NRot long)
int NRot;  // total number of rotations to be made
int Rotations;      //to count how many rotations already made

// image pixel parameters
int PixelMode; // 1 for square mesh, 2 for read a file
int Imag_mesh_x, Imag_mesh_y;
char PixelFile[512];

//define a structure for the sample atoms
struct SampAtom
{
  int AtomCode;
  int TypeCode;
  float Cord_X;
  float Cord_Y;
  float Cord_Z;
};

float BoundX, BoundY, BoundZ;
int NumElem;        // how many types of atoms in the sample 
                    // (eg. NumType=2 for SiO_2)
long int NumAtom;   // how many atoms are there in the sample
long int RNumAtom;  // how many atoms are there in the sample after a 
                    // rotation and a cut

time_t ltime1,ltime2,ltime3; //variables to ask the system time



// g2(r, r') G2Multi[MAX_NUM_ELEM][MAX_NUM_ELEM][R_STEP][IMAG_MESH_X*IMAG_MESH_Y]
float ****G2Multi;

struct SampAtom *Atom;         // for read in the original sample
struct SampAtom *GrowAtom;     // store the sample in a 3*3*3 spacial extention 
                               // with periodic boundary conditions
struct SampAtom *RotaGrowAtom; // the 3*3*3 big sample rotated
struct SampAtom *RotaAtom;     // the rotated big sample is cut into a 
                               // cube of original size (but in new position)

int AtomType[MAX_NUM_ELEM];
int AtomStat[MAX_NUM_ELEM];

//intensity[|k|][r(image)][rotations]
float ***Inte;

//(intensity[|k|][r(image)][rotations])^2
float ***Inte2;

//intensity[k][r][rotaions] averaged over image for one particular rotation
float **Inte_Imag_Aver;

//(intensity[k][r][rotaions])^2 averaged over image for one particular rotation
// Inte2_Imag_Aver[K_STEP][ANGL_STEP*ANGL_STEP];
float **Inte2_Imag_Aver;

//Vari_Imag[K_STEP][ANGL_STEP*ANGL_STEP];
float **Vari_Imag;

// all K_STEP size, 1D arrays:
float *Inte_Imag_Aver2_Aver;
float *Inte_Aver;//intensity[k][r][rotaions] averaged over image and rotation
float *Inte2_Aver;//(intensity[k][r][rotaions])^2 averaged over image and rotation
float *Vari_Aver;
float *Vari2_Aver;
float *Vari;//variation
float *Vari2;
float *VariVari;
float *Sd;

//scattering factors of 104 elements, each have 12 factors
float ScatFact[104][12];

// ScatFunc[104][K_STEP], scattering function value table ScatFact[element][|k|]
float **ScatFunc;


// define a structure for the atoms in a scattering-relevant 
// cylinder of a image point
struct CyliAtom
{
	int AtomCode;
	float d_x;//distance in x axis between sample atom and image point (after applying periodic boundary condition)
	float d_y;
	float cord_x;//coordinate x axis for the atom after applying periodic boundary condition 
	float cord_y;
	float s_l;//distance in xy plane between sample atom and image point (after applying periodic boundary condition)
};

//store atoms in the cylinder in each image point (pixel)
struct CyliAtom **AtomIn;

//count how many atoms are there in the cylinder of each image point
int *AtomCoun;

//file pointer for a log of the program run.
FILE *FpLog;
char tolog[1024];

//globals describing the pixel array
float *pix_x;   // pixel center x positions
float *pix_y;   // pixel center y positions
int npix;   // number of pixels per orientation of the model
