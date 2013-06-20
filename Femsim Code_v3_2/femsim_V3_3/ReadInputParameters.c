
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "FemDefinitions.h"
#include "FemGlobals.h"
#include "FemSim.h"
//Edited by Feng Yi on 04/04/2008
//read values for Start_step and Stop_Step


int ReadParameters(FILE *f) {

  char s[1024];

  // Gets and stores starting k, k step, and number of ks
  // from input .fem file
  if(!GetLine(f, s)) {
    proglog("Parameters file terminated prematurely.  Exiting.\n");
    exit(0);
  }
  sscanf(s, "%f %f %d", &K_mini, &K_scal, &K_step);

  // Gets and stores objective aperature, 0.61/Q, from input .fem file
  if(!GetLine(f, s)) {
    proglog("Parameters file terminated prematurely.  Exiting.\n");
    exit(0);
  }
  sscanf(s, "%f", &Q_aperture);

  // Gets and stores number of rotations, N^3/2 ??? from .fem
  if(!GetLine(f, s)) {
    proglog("Parameters file terminated prematurely.  Exiting.\n");
    exit(0);
  }
  //sscanf(s, "%d %f", &Angl_step, &Thet_ofst);
  sscanf(s, "%d %d %d", &Angl_step, &Start_Step, &Stop_Step); //Edited by Feng Yi on 04/04/2008
  Thet_ofst=0;
  if(Angl_step) {
    Do_rotations = 1;
  }
  else {
    Do_rotations = 0;
    phi_step = 1;
    theta_step = 1;
    psi_step = 1;
  }

  // Gets and stores number of x and y pixels from .fem
  // If there is an input file, read in that filename instead
  // for later calculations.
  // Both cannot be true.
  if(!GetLine(f, s)) {
    proglog("Parameters file terminated prematurely.  Exiting.\n");
    exit(0);
  }
  Imag_mesh_x = atoi(s); // check if the line starts with a number
  if(Imag_mesh_x) {
    sscanf(s, "%d %d", &Imag_mesh_x, &Imag_mesh_y);
    PixelMode = 1; // square mesh
  } else {
    sscanf(s, "%s ", PixelFile);
    PixelMode = 2; // read file
  }

  // Gets and stores which intensity algorithm to use
  // This must be LAST due to the read line in case 2:
  if(!GetLine(f, s)) {
    proglog("Parameters file terminated prematurely.  Exiting.\n");
    exit(0);
  }
  sscanf(s, "%d", &Algorithm);
  switch(Algorithm) {
  case 1: // sum
    break;
  case 2: // g2
    // Get Delta_r from the last line of .fem file
    if(!GetLine(f, s)) {
      proglog("Parameters file terminated prematurely.  Exiting.\n");
      exit(0);
    }
    sscanf(s, "%f", &Delta_r);
    R_step = (int)ceil(0.61/(Q_aperture*Delta_r));
    break;
  case 3: // multslice
    proglog("This algorithm is not yet implemented.  Exiting.\n");
    exit(0);
    break;
  default:
    sprintf(tolog, "Algorithm %d is not a choice.  Exiting.\n", Algorithm);
    proglog(tolog);
    exit(0);
  }
  
  return 1;
}

int GetLine(FILE *f, char *s) {
// GetLine skips blank lines or commented lines due to if(strlen(s)) break;
    int i=0;

    while(1) {
        // fgets reads in <arg2> characters or until EOF or EOL is reached
        fgets(s, 1024, f);
        if(!s) return 0;

        while(1) {
            if(s[i] == '#') {
                s[i] = '\0';
                break;
            }
            if(s[i] == '\0')
                break;
            i++;
        }
        if(strlen(s)) break;
    }

    return 1;
}

