// functions for calculating intensity from large models using the g2(r,r')
// method from Pawel Keblinski.
//
// changed atomcorr and intensity2 to use pre-calculated pixel values pmv 1-28-05
// added intensity<0 warning to intensityG2 pmv 02-25-06
// added intensity is NaN warning to intensityG2 pmb 05-31-06

#include <math.h>
#include "FemDefinitions.h"
#include "FemGlobals.h"
#include "FemSim.h"


// reset G2Multi array to all zero
void cleanG2(void) {

  int i,j,m,n;

  for(i=0;i<NumElem;i++) {
    for(j=0;j<NumElem;j++) {
      for(m=0;m<R_step;m++) {
	for(n=0;n<npix;n++) {
	    G2Multi[i][j][m][n]=0;
	}
      }
    }
  }
}


// calculate G2Multi for the sample currently in RotaAtom.
// this currently does an unrestricted double sum over all the
// atoms in the sample, which is hugely wasteful.  Only atoms within
// an Airy function width of  the center of each pixel need to be considered.
// can fix this with a cell assignment algorithm.
void atomcorr(void) {

  int m,l,j,discrete_r,count;
  float r_x,r_y,d_x1,d_y1,cord_x1,cord_y1,s_l1,d_x2,d_y2,cord_x2,cord_y2,s_l2,s_l,R_c,airy1;
  R_c=0.61/Q_aperture;
  count=0;

  //loop over image pixels
  for(m=0;m<npix;m++) {

    r_x=pix_x[m];
    r_y=pix_y[m];

    //loop over the atoms left after rotation of the sample
    for(l=0;l<RNumAtom;l++) {

      //compute displacement of atom l from the center of the pixel (r_x, r_y)
      d_x1 = RotaAtom[l].Cord_X - r_x;
      cord_x1 = RotaAtom[l].Cord_X;
      cord_y1 = RotaAtom[l].Cord_Y;
      if(d_x1 > BoundX) {
	d_x1=d_x1-2*BoundX;
	cord_x1=cord_x1-2*BoundX;
      }
      else if(d_x1 < -BoundX) {
	d_x1 = d_x1+2*BoundX;
	cord_x1 = cord_x1+2*BoundX;
      }
      
      d_y1 = RotaAtom[l].Cord_Y - r_y;
	if(d_y1 > BoundY) {
	  d_y1 = d_y1-2*BoundY;
	  cord_y1 = cord_y1-2*BoundX;
	}
	else if(d_y1 < -BoundY) {
	  d_y1 = d_y1+2*BoundY;
	  cord_y1 = cord_y1+2*BoundX;
	}
	
	s_l1=sqrt(d_x1*d_x1+d_y1*d_y1);
	
	if(s_l1<R_c) {
	  
	  airy1=a_q(s_l1);
	  
	  // loop over unique pairs of atoms - l is a counter running from 0 to RNumAtoms
	  for(j=0;j<l;j++) {
	    
	    // compute the displacement of atom j from the center of the pixel
	    d_x2=RotaAtom[j].Cord_X-r_x;cord_x2=RotaAtom[j].Cord_X;cord_y2=RotaAtom[j].Cord_Y;
	    if(d_x2>BoundX){d_x2=d_x2-2*BoundX;cord_x2=cord_x2-2*BoundX;}
	    else if(d_x2<-BoundX){d_x2=d_x2+2*BoundX;cord_x2=cord_x2+2*BoundX;}
	    d_y2=RotaAtom[j].Cord_Y-r_y;
	    if(d_y2>BoundY){d_y2=d_y2-2*BoundY;cord_y2=cord_y2-2*BoundX;}
	    else if(d_y2<-BoundY){d_y2=d_y2+2*BoundY;cord_y2=cord_y2+2*BoundX;}
	    s_l2=sqrt(d_x2*d_x2+d_y2*d_y2);
	    
	    // compute the distance between s & j
	    s_l=sqrt((cord_x1-cord_x2)*(cord_x1-cord_x2)+(cord_y1-cord_y2)*(cord_y1-cord_y2));
	    
	    // if both atoms are in the cylinder, add their contribution to G2Multi
	    if((s_l1<R_c)&&(s_l2<R_c)) {
	      if(s_l>2*R_c) {printf("sth wrong\n");continue;}
	      discrete_r=(int)(s_l*R_step/(2*R_c));
	      G2Multi[RotaAtom[l].TypeCode][RotaAtom[j].TypeCode][discrete_r][count]+= 
		2*a_q(s_l1)*a_q(s_l2);
	    }
	  }
	  // add the contributions from the self-pairs j-j
	  G2Multi[RotaAtom[l].TypeCode][RotaAtom[l].TypeCode][0][count]+=airy1*airy1;
	  
	} // end if atom in cylinder
	
    } // end loop over atoms
    count++;
  } //end loop over pixels
} //end atomcorr


// calculate intensity from G2Multi
void intensityG2(void) {

  int i,j,m,l,h;
  float k_mod,r_prime,R_c;
  R_c=0.61/Q_aperture;

  // loop over k values
  for(i=0;i<K_step;i++) {
    k_mod=i*K_scal+K_mini;

    // loop over image pixels x & y
    for(m=0;m<npix;m++) {

      // loop over r' argument in g2(r, r')
      for(j=0;j<R_step;j++) {
	r_prime=(float)(j*2.0*R_c/R_step);

	// loop over different elements, twice.
	for(l=0;l<NumElem;l++) {
	  for(h=0;h<NumElem;h++) {
	    if(G2Multi[l][h][j][m]!=0) { 
	      Inte[i][m][Rotations] += ScatFunc[AtomType[l]][i]*ScatFunc[AtomType[h]][i]*G2Multi[l][h][j][m]*bessj0(2*PI*k_mod*r_prime);
	    }
	  }
	} // end elements loop
      } // end r' loop
      if(Inte[i][m][Rotations] < 0.0) {
	sprintf(tolog, "WARNING: Intensity at k = %g, pixel %d is %g, less than zero!\n", 
		k_mod, m+1, Inte[i][m][Rotations]);
	sprintf(tolog, "%sYou may need to decrease the r' step size.\n", tolog);
	proglog(tolog);
      }
      if(isnan(Inte[i][m][Rotations])) {
	sprintf(tolog, "WARNING: Intensity at k = %g, pixel %d is NaN!\n", k_mod, m+1);
	sprintf(tolog, "%sSomething is wrong.\n", tolog);
	proglog(tolog);
      }
    } // end image pixels loop
  }
} 
