#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "slicelib.h"
#include <sys/time.h>

int main() {
  
  int i, j, k;
  long ltime;
  unsigned long iseed;
  struct timeval tv;
  int nseeds, nsteps, nwait;

  nseeds = 10000;
  nsteps = 10;
  nwait = 1000000;

  /*
  printf("Size of long is %d\n", sizeof(long));
  printf("2 to the 63 is %lu\n", (long)pow(2, 63));
  */
  
  for(i=0; i<nseeds; i++) {
    gettimeofday(&tv, NULL);
    ltime = (long) tv.tv_usec;
    iseed = (unsigned) ltime;
    for(j=0; j<nsteps; j++) {
      printf("%g\n", rangauss(&iseed));
    }
    for(k=0; k<nwait; k++) {}
  }

}
