// functions for calculating rotation angles for a given sample
// current does uniform rotations only, but could be expanded
// a la Pixels.c to offer other methods in the future.
//
// pmv 6/3/06
//Revised by Feng Yi on 04/04/2008

#include <stdlib.h>
#include <math.h>
#include "FemDefinitions.h"
#include "FemGlobals.h"
#include "FemSim.h"


void init_rotations() {

     RotationsUniform();

}


void RotationsUniform() {
  
  int i, j, k, count;
  float *phi_temp, *theta_temp, *psi_temp; //three temporary float arrays, 04/04/2008 Feng Yi
  float s;

  phi_step = Angl_step; // runs from 0 to 2 PI
  theta_step = (int)(ceil((float)Angl_step / 2.0)); // runs for 0 to PI, so half as many steps
  psi_step = Angl_step; // runs from 0 to 2 PI

  
  phi_range = float1D(phi_step, "Unable to allocate phi angle range array.");
  theta_range = float1D(theta_step, "Unable to allocate theta angle range array.");
  psi_range = float1D(psi_step, "Unable to allocate psi angle range array.");
  

  // first generate uniform rotation angles in phi, theta, and psi
  s = 2*PI / phi_step;
  for(i=0; i<phi_step; i++) {
    phi_range[i] = s*i;
  }

  // need to avoid theta = 0 in the Goldstein x convention.  Otherwise, the phi and psi
  // rotations are about the same axis.
  s = PI / (theta_step + 1.0);
  for(i=0; i<theta_step; i++) {
    theta_range[i] = s*(i+1);
  }

  s = 2*PI / psi_step;
  for(i=0; i<psi_step; i++) {
    psi_range[i] = s*i;
  }

  // now create and fill in rot_ang
  /*
  NRot = phi_step*theta_step*psi_step + phi_step;
  phi = float1D(NRot, "phi angle array");
  theta = float1D(NRot, "theta angle array");
  psi = float1D(NRot, "psi angle array");
*/
//Original code above

  NRot = phi_step*theta_step*psi_step + phi_step;
  phi_temp = float1D(NRot, "phi angle array");
  theta_temp = float1D(NRot, "theta angle array");
  psi_temp = float1D(NRot, "psi angle array");

 //Code written by Feng Yi
  
  // add rotations about phi with theta and psi = 0.  Also adds (0,0,0) as the first rotation
  // rotations about theta and phi and psi 0 are already captured in the list, and rotations
  // about psi with phi and theta 0 are equivalent to rotations about phi.
  count = 0;
  for(i=0; i<phi_step; i++) {
  /*
    phi[count] = phi_range[i];
    theta[count] = 0.0;
    psi[count] = 0.0;
    count++;
  */

  //Code written by Feng Yi on 04/04/2008
  phi_temp[count] = phi_range[i];
  theta_temp[count] = 0.0;
  psi_temp[count] = 0.0;
  count++;


  }

  // add rotations about all axes simultaneously
  for(i=0; i<phi_step; i++) {
    for(j=0; j<theta_step; j++) {
      for(k=0; k<psi_step; k++) {
  /*
	phi[count] = phi_range[i];
	theta[count] = theta_range[j];
	psi[count] = psi_range[k];
	count++;
  */

//code written by Feng Yi on 04/04/2008
      phi_temp[count] = phi_range[i];
	theta_temp[count] = theta_range[j];
	psi_temp[count] = psi_range[k];
	count++;

      }
    }
  }

//Code written by Feng Yi on 04/04/2008
NRot=Stop_Step-Start_Step;
count=Start_Step+phi_step;

//Assign space for three angle arrays
phi = float1D(NRot, "phi angle array");
theta = float1D(NRot, "theta angle array");
psi = float1D(NRot, "psi angle array");

for(i=0;i<NRot;i++)
{
 phi[i] = phi_temp[count] ;
 theta[i] = theta_temp[count] ;
 psi[i] = psi_temp[count] ;
 count++;

}

//Free temporary array
free(phi_temp);
free(theta_temp);
free(psi_temp);

}
