// functions to read atomic coordinates from a sample file and manipulate them
// afterwards.
//
// added ReadSampleAtoms 1-29-05 pmv
// changed rotate_matrix to Goldstein "x-convention" with Euler angles 2-24-06 pmv
// fixed error message for exceeding bounding box z to read "z" instead of "x" 3-5-08 pmv

#include <math.h>
#include <stdlib.h>
#include "FemDefinitions.h"
#include "FemGlobals.h"
#include "FemSim.h"

// reads sample in the Kirkland XYZ format
void ReadSampleAtoms(char sample_path[]) {

  FILE *sf;
  int i, z_i, n_z;
  char line[512];
  float t;

  if(!(sf = fopen(sample_path, "r"))) {
    sprintf(tolog, "Cannot open sample file %s.  Exiting.\n", sample_path);
    proglog(tolog);
    exit(0);
  }


  // read bounding box and comment line
  //count the number of atoms and the number of unique elements
  fgets(line, 512, sf);
  proglog(line);
  proglog("\n");

  if((fscanf(sf, "%g %g %g", &BoundX, &BoundY, &BoundZ)!=3)) {
    printf("Cannot read bounding box from sample.\n");
    exit(0);
  }
  // file format calls for full box size, but code wants 1/2 box size.
  BoundX /= 2.0;
  BoundY /= 2.0;
  BoundZ /= 2.0;

  NumElem=0;
  NumAtom=0;
  while(1) {
    if(!fscanf(sf, "%d %f %f %f %f %f", &z_i, &t, &t, &t, &t, &t)) {
      break;
    }
    if(z_i == -1)
      break;

    if(NumAtom==0) {
      AtomType[0] = z_i;
      NumElem=1;
    }

    n_z = 0;
    do {
      if(z_i == AtomType[n_z])
	break;
      n_z++;      
    } while(n_z < NumElem);
    if(n_z==NumElem){
      AtomType[NumElem] = z_i;
      AtomStat[z_i] = 0; // initialize AtomStat array
      NumElem++;
    }
      NumAtom++;

  }

  // Allocate memory for sample (Atom), 3x3x3 periodic extension (GrowAtom), rotated
  // 3x3x3 extension (RotaGrowAtom), and cut back down to size (RotaAtom).
  if((Atom=(struct SampAtom *)malloc(NumAtom*sizeof(struct SampAtom)))==NULL) {
    proglog("sample malloc error");
    exit(0);
  }
  if((GrowAtom=(struct SampAtom *)malloc(27*NumAtom*sizeof(struct SampAtom)))==NULL) {
    proglog("sample malloc error");
    exit(0);
  }
  if((RotaGrowAtom=(struct SampAtom *)malloc(27*NumAtom*sizeof(struct SampAtom)))==NULL) {
    proglog("sample malloc error");
    exit(0);
  }
  if((RotaAtom=(struct SampAtom *)malloc(2*NumAtom*sizeof(struct SampAtom)))==NULL) {
    proglog("sample malloc error");
    exit(0);
  }
  

  rewind(sf);
  fgets(line, 512, sf);
  fgets(line, 512, sf);
  for(i=0; i<NumAtom; i++) {
    fscanf(sf, "%d %f %f %f %f %f", &Atom[i].AtomCode, &Atom[i].Cord_X, &Atom[i].Cord_Y, &Atom[i].Cord_Z, &t, &t);
    n_z=0;
    do {
      if(Atom[i].AtomCode == AtomType[n_z]) {
	Atom[i].TypeCode = n_z;
	AtomStat[n_z]++;
	break;
      }
      n_z++;
    } while(n_z<NumElem);
  }

  fclose(sf);

  CenterZero();

}


// generate the members of a 3x3 rotation matrix in *a
// Use the Goldstein "x-convention" and Euler angles phi
// theta, psi.
void rotate_matrix(float phi,float theta,float psi,float *a) {

  *(a)   = cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi);
  *(a+1) = cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi);
  *(a+2) = sin(psi)*sin(theta);
  *(a+3) = -1.0*sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi);
  *(a+4) = -1.0*sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi);
  *(a+5) = cos(psi)*sin(theta);
  *(a+6) = sin(theta)*sin(phi);
  *(a+7) = -1.0*sin(theta)*cos(phi);
  *(a+8) = cos(theta);

}


//do the rotations of the grown sample.
void rotation(float r[3][3])
{
  int i;
  float x,y,z;
  for(i=0;i<(27*NumAtom);i++)
    {
      x=GrowAtom[i].Cord_X*r[0][0]+GrowAtom[i].Cord_Y*r[0][1]+GrowAtom[i].Cord_Z*r[0][2];
      y=GrowAtom[i].Cord_X*r[1][0]+GrowAtom[i].Cord_Y*r[1][1]+GrowAtom[i].Cord_Z*r[1][2];
      z=GrowAtom[i].Cord_X*r[2][0]+GrowAtom[i].Cord_Y*r[2][1]+GrowAtom[i].Cord_Z*r[2][2];
      RotaGrowAtom[i].AtomCode=GrowAtom[i].AtomCode;
      RotaGrowAtom[i].TypeCode=GrowAtom[i].TypeCode;
      RotaGrowAtom[i].Cord_X=x;
      RotaGrowAtom[i].Cord_Y=y;
      RotaGrowAtom[i].Cord_Z=z;
    }
}

// grow the sample to 3*3*3 times of original size using periodic boundary 
// conditions (for rotation and cut)
void grow(void)
{
  int i,j,k,a,b;
  b=0;
  for(i=-1;i<1.1;i++)
    {
      for(j=-1;j<1.1;j++)
	{
	  for(k=-1;k<1.1;k++)
	    {
	      for(a=0;a<NumAtom;a++)
		{
		  GrowAtom[b].AtomCode=Atom[a].AtomCode;
		  GrowAtom[b].TypeCode=Atom[a].TypeCode;
		  GrowAtom[b].Cord_X=Atom[a].Cord_X+i*2*BoundX;
		  GrowAtom[b].Cord_Y=Atom[a].Cord_Y+j*2*BoundY;
		  GrowAtom[b].Cord_Z=Atom[a].Cord_Z+k*2*BoundZ;
		  b++;
		}
	    }
	}
    }
}

//cut the rotated grown sample to original size in the rotated position
void cut(void)
{
  int i,j;
  RNumAtom=0;
  j=0;
  for(i=0;i<(27*NumAtom);i++)
    {
      if((fabs(RotaGrowAtom[i].Cord_X)<BoundX)&&(fabs(RotaGrowAtom[i].Cord_Y)<BoundY)&&(fabs(RotaGrowAtom[i].Cord_Z)<BoundZ))
	{
	  RotaAtom[j].AtomCode=RotaGrowAtom[i].AtomCode;
	  RotaAtom[j].TypeCode=RotaGrowAtom[i].TypeCode;
	  RotaAtom[j].Cord_X=RotaGrowAtom[i].Cord_X;
	  RotaAtom[j].Cord_Y=RotaGrowAtom[i].Cord_Y;
	  RotaAtom[j].Cord_Z=RotaGrowAtom[i].Cord_Z;
	  j++;
	}
    }
  RNumAtom=j;
  if(abs( (RNumAtom - NumAtom)/NumAtom ) > 0.01)
    proglog("WARNING: There has been a greater than 2% change in the number of atoms on rotation.\n");

}


// Center zero:
// Ensure the that model is extends from e.g. -a/2 to a/2 rather than
// 0 to a.  As a by-product checks to make sure there are no atoms outside
// the supercell box.
void CenterZero() {

  float min_x, max_x, min_y, max_y, min_z, max_z;
  float center_x, center_y, center_z;
  int i;

  min_x = Atom[0].Cord_X;
  max_x = Atom[0].Cord_X;
  min_y = Atom[0].Cord_Y;
  max_y = Atom[0].Cord_Y;
  min_z = Atom[0].Cord_Z;
  max_z = Atom[0].Cord_Z;
  
  for(i=0; i<NumAtom; i++) {
    if(Atom[i].Cord_X < min_x) min_x = Atom[i].Cord_X;
    if(Atom[i].Cord_Y < min_y) min_y = Atom[i].Cord_Y;
    if(Atom[i].Cord_Z < min_z) min_z = Atom[i].Cord_Z;
    if(Atom[i].Cord_X > max_x) max_x = Atom[i].Cord_X;
    if(Atom[i].Cord_Y > max_y) max_y = Atom[i].Cord_Y;
    if(Atom[i].Cord_Z > max_z) max_z = Atom[i].Cord_Z;
  }    

  if(max_x - min_x > 2*BoundX) {
    sprintf(tolog, "Maximum x atom separation of %f is outside the supercell size of %f.  Exiting.\n",
	   max_x-min_x, 2*BoundX);
    proglog(tolog);
    exit(0);
  }
  if(max_y - min_y > 2*BoundY) {
    sprintf(tolog, "Maximum y atom separation of %f is outside the supercell size of %f.  Exiting.\n",
	   max_y-min_y, 2*BoundY);
    proglog(tolog);
    exit(0);
  }
  if(max_z - min_z > 2*BoundZ) {
    sprintf(tolog, "Maximum z atom separation of %f is outside the supercell size of %f.  Exiting.\n",
	   max_z-min_z, 2*BoundZ);
    proglog(tolog);
    exit(0);
  }

  center_x = (min_x + max_x)/2.0; center_y = (min_y+max_y)/2.0; center_z = (min_z+max_z)/2.0;
 
  for(i=0; i<NumAtom; i++) {
    Atom[i].Cord_X -= center_x;
    Atom[i].Cord_Y -= center_y;
    Atom[i].Cord_Z -= center_z;
  }

  sprintf(tolog, "Shifting center of the model from (%g, %g, %g) to (0, 0, 0) for rotations.\n", 
	 center_x, center_y, center_z);
  proglog(tolog);
  
}
