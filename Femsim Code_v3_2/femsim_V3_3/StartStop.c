// 

// modified initialize to use ReadSampleAtoms from Sample.c 1-29-05 pmv
//
// added dynamic memory allocation of just about everything 
// to initialize and cleanup, 8/19/05 pmv
//
// fixed bug in all intensities output when rotations are enabled 5-11-06 pmv

//Revised output considering rotation range by Feng Yi on 04/04/2008

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "FemDefinitions.h"
#include "FemGlobals.h"
#include "FemSim.h"
#include "ScatFacs.h"

void initialize(char *sample, char *basename)
{
    int i,g2length;

    // read in the sample atoms
    proglog("\n2. Reading sample atoms with description:\n");
    ReadSampleAtoms(sample);
    g2length=(int)(exp(NumElem*log(2.0)))*R_step*npix;
    proglog("Sample read finished.\n");

    // initialize the pixel array: compute (or read in) the position at which
    // to calculate the intensity
    proglog("\n3. Initializing intensity pixel positions.\n");
    init_pixels();

    // initialize form factors
    Ffactors();

    // log sample and simulation parameters.
    printf("\n4. Logging sample and simulation parameters to file.\n");
    fprintf(FpLog,"\n4. Simulation setup and parameters:\n");
    fprintf(FpLog, "Sample data:\n");
    fprintf(FpLog,"\tThe sample size is:\n\tx(%f, %f), y(%f, %f), z(%f, %f)\n", -BoundX,BoundX,-BoundY,BoundY,-BoundZ,BoundZ);
    fprintf(FpLog, "\tThe sample contains %d atoms of %d different elements.\n", (int)NumAtom, NumElem);
    for(i=0; i<NumElem; i++) {
        fprintf(FpLog, "\tThere are %d atoms with Z = %d.\n", AtomStat[i], AtomType[i]);
    }

    fprintf(FpLog, "\nImage data:\n");
    fprintf(FpLog,"\tQ_aperture=%f\n",Q_aperture);
    if(Do_rotations) {
        init_rotations();
        fprintf(FpLog,"\tRotations: %d in phi, %d in theta, and %d in psi.  Total rotations: %d\n", 
                phi_step, theta_step, psi_step, NRot);
        if(Thet_ofst) {
            fprintf(FpLog,"\t\t%g radian offset in theta, which is NOT CORRECTLY IMPLEMENTED.\n", Thet_ofst);
        }
    }
    else
        fprintf(FpLog,"\tNo rotations.\n");
    fprintf(FpLog,"\t|k|=%f->%f,has %d steps\n",K_mini,K_mini+K_step*K_scal,K_step);
    fprintf(FpLog,"\tSize of image is:\n\tx(%f, %f), y(%f, %f)\n",-BoundX,BoundX,-BoundY,BoundY);
    fprintf(FpLog,"\tNumber of pixels per image: %d\n", npix);
    fprintf(FpLog,"\n\nTotal number of intensity samples: %d\n", npix*NRot);
    fflush(FpLog);

} // end initialize


//calculate a list of scattering functions
void getscatfunc(void)
{
    int i,j;
    float Lore,Gaus;
    float k,k_squa;

    // why 104? ???
    ScatFunc = float2D(104, K_step, "ScatFunc");

    for(i=0;i<NumElem;i++) {
        for(j=0;j<K_step;j++) {
            k=j*K_scal+K_mini;
            k_squa=k*k;
            Lore = ScatFact[AtomType[i]][0]/(k_squa+ScatFact[AtomType[i]][1]) + 
                ScatFact[AtomType[i]][2]/(k_squa+ScatFact[AtomType[i]][3]) + 
                ScatFact[AtomType[i]][4]/(k_squa+ScatFact[AtomType[i]][5]);
            Gaus = ScatFact[AtomType[i]][6]*exp(-ScatFact[AtomType[i]][7]*k_squa) + 
                ScatFact[AtomType[i]][8]*exp(-ScatFact[AtomType[i]][9]*k_squa) + 
                ScatFact[AtomType[i]][10]*exp(-ScatFact[AtomType[i]][11]*k_squa);
            ScatFunc[AtomType[i]][j]=Lore+Gaus;
        }
    }
    //fprintf(FpLog,"\tScattering functions calculated.\n");
}


// write various output files.
void output(char *basename) {

    FILE *fpo1,*fpo2,*fpo3,*fpo4,*fpo5,*fpo6;
    int i,j,m,l,n;
    int count1;
    float k_mod;
    char outname[1024];

    proglog("\n8. Writing output files.\n");

    if(Do_rotations) {
        sprintf(outname, "%s_rota_stat.txt", basename);
        if((fpo1=fopen(outname,"w"))==NULL){printf("open error");goto end;}
    }
    else fpo1=NULL;
    sprintf(outname, "%s.out", basename);
    if((fpo2=fopen(outname,"w"))==NULL){printf("open error");goto end;}
    sprintf(outname, "%s_scatfunc.txt", basename);
    if((fpo3=fopen(outname,"w"))==NULL){printf("open error");goto end;}
    sprintf(outname, "%s_inte03.txt", basename);
    fpo4=fopen(outname,"w");
    sprintf(outname, "%s_inte06.txt", basename);
    fpo5=fopen(outname,"w");
    sprintf(outname, "%s_int.txt", basename);
    fpo6=fopen(outname,"w");

    count1=0;

    // write main_out.txt.  Contains in columns k, I(k), I^2(k), V(k), and std dev of I
    fprintf(fpo2, "k \t I(k) \t I^2(k) \t V(k) \t sigma(k)\n");
    for(i=0;i<K_step;i++) {
        fprintf(fpo2, "%f\t%f\t%f\t%f\t%f\n", i*K_scal+K_mini, Inte_Aver[i], Inte2_Aver[i],
                Vari[i], Sd[i]);
    }
    proglog("\tWrote *.out\n");

    // write rota_stat.txt.  Contains two blocks.  In the first block in columns, are the full range of
    // the rotation  angles phi, theta, and psi.  In the second block are the phi, theta, and psi used 
    // in order. 
    if(Do_rotations) {

        fprintf(fpo1, "Rotation angles used in this computation:\n");
        fprintf(fpo1, "phi\ttheta\tpsi\n");
        if(phi_step >= theta_step) count1 = phi_step;
        else count1 = theta_step;
        if(count1 >= psi_step) ;
        else count1 = psi_step;
        for(j=0; j<count1; j++) {
            if(j<phi_step)
                fprintf(fpo1, "%g\t", phi_range[j]);
            else
                fprintf(fpo1, "\t");
            if(j<theta_step)
                fprintf(fpo1, "%g\t", theta_range[j]);
            else
                fprintf(fpo1, "\t");
            if(j<psi_step)
                fprintf(fpo1, "%g\n", psi_range[j]);
            else
                fprintf(fpo1, "\n");
        }


        fprintf(fpo1, "\n\n\nSample orientations used in this computation:\n");
        fprintf(fpo1, "N \t phi \t theta \t psi\n");
        for(j=0; j<NRot; j++) {
            fprintf(fpo1, "%d\t%f\t%f\t%f\n", j, phi[j], theta[j], psi[j]);
        }
    }

    proglog("\tWrote *rota_stat.txt.\n");

    // write scatfunc.txt.  Contains in columns the atom type, k, and F(k) for that atom
    fprintf(fpo3, "atom type \t k \t F(k) \n");
    for(i=0;i<NumElem;i++)
    {
        m=AtomType[i];
        for(j=0;j<K_step;j++)
        {
            k_mod=j*K_scal+K_mini;
            fprintf(fpo3,"%d\t%f\t%f\n",m,k_mod,ScatFunc[m][j]);
        }
    }
    proglog("\tWrote *scatfunc.txt.\n");

    // write inte03.txt and inte06.txt.  I don't know what these are.
    for(i=0;i<npix;i++) {
        fprintf(fpo4,"%f\n",Inte[5][i][0]);
        fprintf(fpo5,"%f\n",Inte[20][i][0]);
    }
    proglog("\tWrote *inte03.txt and *inte06.txt.\n");

    // write the raw intensities for every pixel, k, and rotation to a text file
    fprintf(fpo6, "Intensity of every pixel at every k for every rotation.\n");
    if(Do_rotations) {
        fprintf(fpo6, " \t k (1/A)\n");
        fprintf(fpo6, "rotation.pixel # \t");
        for(m=0; m<K_step; m++) {
            k_mod=m*K_scal+K_mini;
            fprintf(fpo6, "%g \t", k_mod);
        }
        fprintf(fpo6, "\n");

        count1 = 0;
        for(i=0; i<NRot; i++) {
            for(l=0; l<npix; l++) {
                //fprintf(fpo6, "%d-%d \t", i, l);
                fprintf(fpo6, "%d-%d \t", i+Start_Step, l);//Added by Feng Yi on 04/04/2008
                //for actual rotation position
                for(n=0; n<K_step; n++) {
                    fprintf(fpo6, "%g \t", Inte[n][l][i]);
                }
                fprintf(fpo6, "\n");
                fflush(fpo6);
            } // end npix
        } // end NRot loop
    } // end if (Do_rotations)
    else {
        fprintf(fpo6, " \t k (1/A)\n");
        fprintf(fpo6, "pixel # \t");
        for(m=0; m<K_step; m++) {
            k_mod=m*K_scal+K_mini;
            fprintf(fpo6, "%g \t", k_mod);
        }
        fprintf(fpo6, "\n");
        count1=0;
        for(i=0; i<npix; i++) {
            fprintf(fpo6, "%d \t", i);
            for(m=0; m<K_step; m++) {
                fprintf(fpo6, "%g \t", Inte[m][i][0]);
            }
            fprintf(fpo6, "\n");
        }
    } // end write intensities
    proglog("\tWrote *int.txt.\n");


    // write some ending information to the log file
    time(&ltime2);
    sprintf(tolog,"\n\nEnd running at time: %s\n",ctime(&ltime2));
    proglog(tolog);
    ltime3=ltime2-ltime1;
    sprintf(tolog,"Total time spent: %d seconds.\n\n",(int)ltime3);
    proglog(tolog);

    // close all the files.
    if(Do_rotations)
        fclose(fpo1);
    fclose(fpo2);
    fclose(fpo3);
    fclose(fpo4);
    fclose(fpo5);
    fclose(fpo6);

    // file open error handler.
    while(0)
    {
end: printf("Error open file");
    }

}


// free allocated memory, close files, etc.
void cleanup() {

    free(pix_x);
    free(pix_y);
    free(Atom);
    free(GrowAtom);
    free(RotaGrowAtom);
    free(RotaAtom);

    free(G2Multi);
    free(Inte);
    free(Inte2);
    free(Inte_Imag_Aver);
    free(Inte2_Imag_Aver);
    free(Vari_Imag);
    free(Inte_Imag_Aver2_Aver);
    free(Inte_Aver);
    free(Inte2_Aver);
    free(Vari_Aver);
    free(Vari2_Aver);
    free(Vari);
    free(Vari2);
    free(VariVari);
    free(Sd);
    free(ScatFunc);

}


void proglog(char *s) {

    printf("%s", s);
    fprintf(FpLog, "%s", s);
    fflush(FpLog);

}
