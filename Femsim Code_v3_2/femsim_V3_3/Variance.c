// functions for calculating the intensity at many positions and
// hence the variance.

// 2-24-06 pmv: removed sin(theta) angular normalization from appn_I_I2

//Revised by Feng Yi on 04/04/2008

#include <math.h>
#include "FemDefinitions.h"
#include "FemGlobals.h"
#include "FemSim.h"



void appn_I_I2()
{
    int i,m;

    for(i=0;i<K_step;i++)
    {
        for(m=0;m<npix;m++)
        {
            Inte2[i][m][Rotations]=Inte[i][m][Rotations]*Inte[i][m][Rotations];
            Inte_Imag_Aver[i][Rotations]+=Inte[i][m][Rotations];
            Inte2_Imag_Aver[i][Rotations]+=Inte2[i][m][Rotations];
        }
        Inte_Imag_Aver[i][Rotations]=Inte_Imag_Aver[i][Rotations]/(npix);
        Inte2_Imag_Aver[i][Rotations]=Inte2_Imag_Aver[i][Rotations]/(npix);
        Vari_Imag[i][Rotations]=Inte2_Imag_Aver[i][Rotations]/(Inte_Imag_Aver[i][Rotations]*Inte_Imag_Aver[i][Rotations])-1;
        Inte_Aver[i]+=Inte_Imag_Aver[i][Rotations];
        Inte2_Aver[i]+=Inte2_Imag_Aver[i][Rotations];
        Inte_Imag_Aver2_Aver[i]+=Inte_Imag_Aver[i][Rotations]*Inte_Imag_Aver[i][Rotations];
        Vari_Aver[i]+=Vari_Imag[i][Rotations];
        Vari2_Aver[i]+=Vari_Imag[i][Rotations]*Vari_Imag[i][Rotations];
    }
}

void rota_aver_I(void)
{
    int j;
    float a_rotate[3][3];
    float *a_pointer;
    a_pointer=&(a_rotate[0][0]);
    grow();
    Rotations=0;

    if(Do_rotations) {
        proglog("\n6. Calculating the intensity and average over rotations of sample.\n");

        switch(Algorithm) {
            case 1: // direct summation
                proglog("\tUsing the direct summation algorithm.\n");
                break;
            case 2: // g2
                proglog("\tUsing the g2(r, r') algorithm.\n");
                sprintf(tolog, "\tDelta r' is %g Angstroms, so\n", Delta_r);
                proglog(tolog);
                sprintf(tolog, "\tr' is divided into %d steps in the cylinder.\n", R_step);
                proglog(tolog);
                break;
        }
        fflush(FpLog);

        //for(j=0;j<NRot; j++)      
        //for(j=Start_Step;j<Stop_Step; j++)      
        for(j=0;j<NRot; j++) {   //Revised here by Feng Yi on 04/04/2008
            sprintf(tolog, "\tRotation: phi=%f \t theta=%f \t psi=%f\n", phi[j], theta[j], psi[j]);
            proglog(tolog);
            rotate_matrix(phi[j],theta[j],psi[j],a_pointer);
            rotation(a_rotate);
            cut();
            switch(Algorithm) {
                case 1: // direct summation
                    atomcyli();
                    intensitysum();
                    break;
                case 2: // g2(r, r')
                    cleanG2();
                    atomcorr();
                    intensityG2();
                    break;
            }
            appn_I_I2();
            Rotations++;
            fflush(FpLog);
        }

    } // end if(Do_rotations)
    else {
        proglog("\n6. Calculating the intensities on the sample in the orientation:\n");
        sprintf(tolog, "\tphi=0.0\ttheta=0.0\tpsi=0.0\n");
        proglog(tolog);
        rotate_matrix(0.0,0.0,0.0,a_pointer);
        rotation(a_rotate);
        cut();
        switch(Algorithm) {
            case 1: // direct summation
                atomcyli();
                intensitysum();
                break;
            case 2: // g2(r, r')
                cleanG2();
                atomcorr();
                intensityG2();
                break;
        }
        appn_I_I2();//theta is not used in this case
        Rotations++;
    } // end no rotations
}


void variance(void)
{
    int i,j,m;
    double one_int, one_int2;
    proglog("\n7. Calculating the variance.\n");
    for(i=0;i<K_step;i++)
    {
        // thought at one point that scaling by a magnification factor might 
        // be necessary, but it's not.
        //Inte_Aver[i]=Inte_Aver[i]/(MAGN_FACT*MAGN_FACT);
        //Inte2_Aver[i]=Inte2_Aver[i]/(MAGN_FACT*MAGN_FACT*MAGN_FACT*MAGN_FACT);
        //Inte_Imag_Aver2_Aver[i]=Inte_Imag_Aver2_Aver[i]/(MAGN_FACT*MAGN_FACT*MAGN_FACT*MAGN_FACT);

        one_int = 0.0;
        one_int2 = 0.0;

        for(j=0;j<NRot;j++) {
            for(m=0;m<npix;m++) {
                one_int += Inte[i][m][j];
                one_int2 += Inte[i][m][j]*Inte[i][m][j];
            }
        }

        one_int /= (npix*NRot);
        one_int2 /= (npix*NRot);
        Inte_Aver[i] = one_int;
        Inte2_Aver[i] = one_int2;
        Vari[i] = one_int2 / (one_int*one_int) - 1.0;

        //      Vari[i]=Inte2_Aver[i]/(Inte_Aver[i]*Inte_Aver[i])-1.0;
        Vari2[i]=Inte2_Aver[i]/Inte_Imag_Aver2_Aver[i]-1.0;
        VariVari[i]=Vari2_Aver[i]/(Vari_Aver[i]*Vari_Aver[i])-1.0;
        Sd[i]=Inte2_Aver[i]-Inte_Aver[i]*Inte_Aver[i];
    }

}

