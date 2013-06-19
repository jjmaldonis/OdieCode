// FemSim, a program to simulation Fluctuation Electron Microscopy
//
// Computes I(k), I^2(k), and V(k) for an arbitrary atomic structure
//
// Intensities are computed  using either the direct quadruple sum 
// algorithm (Treacy and Gibson, Acta Cryst A 1996) or the g2(r, r') 
// algorithm (Dash et al J. Phys. Cond. Mat. 2002).
//
// The model can if desired be rotated into different orientations
// for better statistical sampling or to simulate an isotropic sample.
// 
// Version history:
// v1, FemMulti.c by Sanjay Khare
// v2, FemMulti.c & femsim.sh, pmv
// v3.0, FemSim.c, 8-24-05 pmv
// v3.1  2-26-06 pmv
//       change rotations from two-axis to three-axis in Sample.c
//       change orientational averaging to add 3rd axis and remove weighting factor
//       change variance calculation method
//       add error checking for negative intensities in g2
//       change input parameter from # of points in r' integral to step size
// v3.1.1 5-31-06 pmv
//       various bug fixes leading up to v3.2
// v3.2 6-4-06 pmv
//       change to rotations algorithm: moved to precalculated list & change total angle ranges
//3.3 04/04/2008 fy
// rotate to a particular rotation range

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "FemDefinitions.h"
#include "FemGlobals.h"
#include "FemSim.h"

static const char help[] = "\
Incorrect number of input parameters.  Proper usage is:\n\
\n\
     femsim <simulation parameters file> <structure .xyz> <output basename>\n\
\n\
Further help is available in the documentation or from\n\
Paul Voyles, voyles@engr.wisc.edu\n\n";

int main(int argc, char *argv[])
{
  FILE *param;
  char outname[512];

  // Set version string
  sprintf(version, "Femsim v3.2 6-4-06 pmv");
  printf("\n%s\n\n", version);

  // Check for the right number of arguments
  if(argc != 4) {
    printf("%s", help);
    exit(0);
  }

  // Open the log file.
  sprintf(outname, "%s.log", argv[3]);
  if(!(FpLog = fopen(outname, "w"))) {
    printf("Cannot open log file.  Exiting.\n");
    exit(0);
  }
  fprintf(FpLog,"Log of %s\n", version);
  time(&ltime1);
  fprintf(FpLog,"Begin running at time: %s\n\n",ctime(&ltime1));
  fflush(FpLog);

  // Open and parse the input parameters file
  if(!(param = fopen(argv[1], "r"))) {
    sprintf(tolog, "Cannot open parameters file %s Exiting.\n", argv[1]);
    proglog(tolog);
    exit(0);
  }  
  sprintf(tolog, "1. Reading simulation parameters from file %s.\n", argv[1]);
  proglog(tolog);
  ReadParameters(param);
  fclose(param);

  // read atom cordinates from input files, log simulation parameters, 
  initialize(argv[2], argv[3]);

  // calculate the scattering functions for each element 
  // (as a 2D table of atom number and k)
  getscatfunc();

  // allocate memory for intensity computation
  allocate();

  // compute the intensity at many positions on the model in many orientations.  This is
  // where most of the computational work is done.
  rota_aver_I();
  
  // Compute the variance V(k, Q) from the intensity samples computed in rota_aver_I
  variance();

  // write output files
  output(argv[3]);

  cleanup();

  fclose(FpLog);

  return 1;
}
