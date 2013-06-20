// Functions for allocating and freeing memory for femsim.
// Stolen from / based on Memory.c by Earl Kirkland
// begun 02/08/05 pmv.

#include <stdlib.h>
#include <stdio.h>
#include "FemDefinitions.h"
#include "FemGlobals.h"
#include "FemSim.h"

// Allocates memory for intensity calculations
void allocate() {

    proglog("\n5. Allocating memory for computations.\n");

    //Allocate arrays common to all algorithms
    Inte = float3D(K_step, npix, NRot, "Inte");
    Inte2 = float3D(K_step, npix, NRot, "Inte2");
    Inte_Imag_Aver = float2D(K_step, NRot, "Inte_Imag_Aver");
    Inte2_Imag_Aver = float2D(K_step, NRot, "Inte2_Imag_Aver");
    Vari_Imag = float2D(K_step, NRot, "Vari_Imag");
    Inte_Imag_Aver2_Aver = float1D(K_step, "Inte_Imag_Aver2_Aver");
    Inte_Aver = float1D(K_step, "Inte_Aver");
    Inte2_Aver = float1D(K_step, "Inte2_Aver");
    Vari_Aver = float1D(K_step, "Vari_Aver");
    Vari2_Aver = float1D(K_step, "Vari2_Aver");
    Vari = float1D(K_step, "Vari");
    Vari2 = float1D(K_step, "Vari2");
    VariVari = float1D(K_step, "VariVari");
    Sd = float1D(K_step, "Sd");

    switch(Algorithm) {
        case 1: // direct summation
            alloc_sum();
            break;
        case 2: // g2(r, r')
            alloc_g2();
            break;
        case 3: // multislice
            break;
        default:
            proglog("Unknown algorithm for intensities requested.  Exiting.\n");
    }

    proglog("\tMemory allocated.\n");
}

void alloc_sum(void) {
    int i;

    AtomCoun = int1D(npix, "AtomCoun");

    //struct CyliAtom AtomIn[npix][CYLI_MAX];
    AtomIn = (struct CyliAtom**) malloc( npix*sizeof(struct CyliAtom*));
    if(AtomIn == NULL) {
        printf("Cannot allocate memory for AtomIn array of %d pixels.\n", npix);
        exit(0);
    }

    for(i=0; i<npix; i++) {
        AtomIn[i] = (struct CyliAtom*) malloc(CYLI_MAX*sizeof(struct CyliAtom));
        if(AtomIn[i] == NULL) {
            printf("Cannot finish allocating cylinder atoms array.\n");
            exit(0);
        }
    }
}

void alloc_g2(void) {

    G2Multi = float4D(NumElem, NumElem, R_step, npix, "G2Multi");


}


// Allocate a 1D float array f[n1]
float *float1D( int n1, const char *message )
{
    float *m;

    m = (float*) malloc( n1 * sizeof( float) );
    if( m == NULL ) {
        printf("float1D() cannot allocate memory size=%d: %s\n",
                n1, message);
        exit( 0 );
    }
    return( m );

}

// Allocate a 2D float array f[n1][n2]
float **float2D( int n1, int n2, const char *message ) {
    float **m;
    int i;

    m = (float**) malloc( n1 * sizeof( float* ) ); 
    if( m == NULL ) {
        printf("float2D cannot allocate pointers, size=%d: %s\n",
                n1, message );
        exit(0);
    }

    for (i=0; i<n1; i++){
        m[i] = (float *) malloc( n2 * sizeof( float ) );
        if( m[i] == NULL ){
            printf("float2D cannot allocate arrays, size=%d: %s\n",
                    n2, message );
            exit(0);
        }
    }

    return m;

}

// Allocate a 3D float array f[n1][n2][n3]
float ***float3D(int n1, int n2, int n3, const char *message) {
    float ***m;
    int i, j;

    m = (float***) malloc( n1 * sizeof( float** ));
    if( m == NULL) {
        printf("float3D cannot allocate pointers, size=%d: %s\n", n1, message);
        exit(0);
    }

    for(i=0; i<n1; i++) {
        m[i] = (float **) malloc( n2 * sizeof( float* ) );
        if( m[i] == NULL ) {
            printf("float3d cannot allocate points, size=%d: %s\n", n2, message);
            exit(0);
        }

        for(j=0; j<n2; j++) {
            m[i][j] = (float *) malloc( n3 * sizeof( float) );
            if(m[i][j] == NULL) {
                printf("float3D cannot allocate arrays, size=%d: %s\n", n3, message);
                exit(0);
            }
        }
    }

    return m;
}

// Allocate a 4D float array f[n1][n2][n3][n4]
float ****float4D(int n1, int n2, int n3, int n4, const char *message) {
    int i,j,k;
    float ****m;

    m = (float****) malloc( n1*sizeof(float***) );
    if( m == NULL) {
        printf("float4d cannot allocate needed memory.\n%s\n", message);
        exit(0);
    }

    for(i=0; i<n1; i++) {
        m[i] = (float***) malloc( n2*sizeof(float**) );
        if(m[i] == NULL) {
            printf("float4d cannot allocate needed memory.\n%s\n", message);
            exit(0);
        }

        for(j=0; j<n2; j++) {
            m[i][j] = (float**) malloc(n3*sizeof(float*) );
            if(m[i][j] == NULL) {
                printf("float4d cannot allocate needed memory.\n%s\n", message);
                exit(0);
            }

            for(k=0; k<n3; k++) {
                m[i][j][k] = (float*) malloc(n4*sizeof(float) );
                if(m[i][j][k] == NULL) {
                    printf("float4d cannot allocate needed memory.\n%s\n", message);
                    exit(0);
                }

            }
        }
    }

    return m;

}

// allocate a 1D integer array
int* int1D( int n, const char *message )
{
    int *m;

    m = (int*) malloc( n * sizeof( int ) );
    if( m == NULL ) {
        printf("int1D() cannot allocate memory size=%d: %s\n",
                n, message);
        exit( 0 );
    }
    return( m );

}  /* end int1D() */
