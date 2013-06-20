// functions for direct summation algorithm to calculate the intensity of a 
// small number of atoms.  Currently used for models of less than 60 atoms.
// modified atomcyli and intensity1 to use precalculated pixel positions pmv 1-28-05
// added check for negative and NaN intensities pmv 05-31-06

#include <math.h>
#include "FemDefinitions.h"
#include "FemGlobals.h"
#include "FemSim.h"

// Defines:
// void atomcyli(void)
// void intensitysum(void)

//locate the sample atoms in the scattering relevant cylinder of one image point
void atomcyli(void) {
    int m,l,j;
    float r_x,r_y,d_x,d_y,s_l,R_c,cord_x,cord_y;
    R_c=0.61/Q_aperture;

    for(m=0;m<npix;m++) {
        j=0;
        r_x=pix_x[m];
        r_y=pix_y[m];
        for(l=0;l<RNumAtom;l++) {
            if(j>=CYLI_MAX) {
                printf("Too many atoms in one pixel.  Switch to the g2(r, r') algorithm.\n");
                fprintf(FpLog, "Too many atoms in one pixel.  Switch to the g2(r, r') algorithm.\n");
                exit(0);
            }
            d_x=RotaAtom[l].Cord_X-r_x;
            if(d_x>BoundX){d_x=d_x-2*BoundX;cord_x=RotaAtom[l].Cord_X-2*BoundX;}
            else if(d_x<-BoundX){d_x=d_x+2*BoundX;cord_x=RotaAtom[l].Cord_X+2*BoundX;}
            else cord_x=RotaAtom[l].Cord_X;
            d_y=RotaAtom[l].Cord_Y-r_y;
            if(d_y>BoundY){d_y=d_y-2*BoundY;cord_y=RotaAtom[l].Cord_Y-2*BoundY;}
            else if(d_y<-BoundY){d_y=d_y+2*BoundY;cord_y=RotaAtom[l].Cord_Y+2*BoundY;}
            else cord_y=RotaAtom[l].Cord_Y;
            s_l=sqrt(d_x*d_x+d_y*d_y);
            if(s_l<R_c)
            {
                AtomIn[m][j].AtomCode=RotaAtom[l].AtomCode;
                AtomIn[m][j].d_x=d_x;
                AtomIn[m][j].d_y=d_y;
                AtomIn[m][j].cord_x=cord_x;
                AtomIn[m][j].cord_y=cord_y;
                AtomIn[m][j].s_l=s_l;
                j++;
            }
        }
        AtomCoun[m]=j;
    }
}

// calculate the intensity at one particular rotation (only need to use 
// sample atoms in the cylinder got in atomcyli())
void intensitysum(void) {

    int i,m,l,s;
    float k_mod,r_x,r_y,sigma;

    for(i=0;i<K_step;i++) {
        printf("%d\t",i);
        k_mod=i*K_scal+K_mini;

        for(m=0;m<npix;m++) {
            r_x=pix_x[m];
            r_y=pix_y[m];
            for(l=0;l<AtomCoun[m];l++) {
                for(s=0;s<AtomCoun[m];s++) {
                    sigma=sqrt((AtomIn[m][l].cord_x - AtomIn[m][s].cord_x) * (AtomIn[m][l].cord_x - AtomIn[m][s].cord_x) 
                            +(AtomIn[m][l].cord_y - AtomIn[m][s].cord_y) * (AtomIn[m][l].cord_y - AtomIn[m][s].cord_y));
                    Inte[i][m][Rotations] += ScatFunc[AtomIn[m][l].AtomCode][i]*ScatFunc[AtomIn[m][s].AtomCode][i] 
                        *a_q(AtomIn[m][l].s_l)*a_q(AtomIn[m][s].s_l) * bessj0(2*PI*k_mod*sigma);
                } // end AtomCoun 2 loop
            } // end AtomCount 1 loop

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

        } // end pixels loop
    } // end k loop
}

