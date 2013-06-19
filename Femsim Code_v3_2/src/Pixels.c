// functions for calculating pixel positions for a given sample size.
// 
// pmv, 1/28/05
// added read positions from a file option 8-24-05 pmv

#include <stdlib.h>
#include "FemDefinitions.h"
#include "FemGlobals.h"
#include "FemSim.h"

void init_pixels() {

  switch(PixelMode) {
  case 1:
    PixelsSquareArray();
    break;
  case 2:
    PixelsReadFile();
    break;
  default:
    sprintf(tolog, "Incorrect pixel mode %d in init_pixels.  Exiting.\n", PixelMode);
    proglog(tolog);
    exit(0);
  }

}

void PixelsSquareArray() {

  int i, j, n=0;
  float x, y, sx, dx, sy, dy, res;

  npix = Imag_mesh_x * Imag_mesh_y;
  pix_x = float1D(npix, "Cannot allocate memory of pix_x.");
  pix_y = float1D(npix, "Cannot allocate memory for pix_y.");

  res = 0.61/Q_aperture;
  sx = res/2.0 - BoundX;
  sy = res/2.0 - BoundY;

  dx = (2.0*BoundX - res) / (Imag_mesh_x - 1);
  dy = (2.0*BoundY - res) / (Imag_mesh_y - 1);
  
  if(dx < res/3.0) {
    sprintf(tolog, "\n\tPixel spacing of %g A in x is closer than the recommended\n", dx);
    sprintf(tolog, "%s\tminimum of the one-third the resolution, %g A.\n", tolog, res/3.0);
    proglog(tolog);
  }
  if(dy < res/3.0) {
    sprintf(tolog, "\n\tPixel spacing of %g A in y is closer than the recommended\n", dy);
    sprintf(tolog, "%s\tminimum of the one-third the resolution, %g A.\n\n", tolog, res/3.0);
    proglog(tolog);
  }

  fprintf(FpLog, "\tPixel positions on the model are:\n");

  for(i=0; i<Imag_mesh_x; i++) {
    x = sx + i*dx;
    for(j=0; j<Imag_mesh_y; j++) {
      y = sy + j*dy;
      pix_x[n] = x;
      pix_y[n] = y;
      n++;
      fprintf(FpLog, "\t%g, %g\n", x, y);
    }
  }
  fprintf(FpLog, "\n");

  proglog("\tPixel positions calculated on a square array.\n");
}

void PixelsReadFile() {

  FILE *f;
  char line[512];
  int i;

  if(!(f = fopen(PixelFile, "r"))) {
    sprintf(tolog, "Cannot open pixel positions file %s.  Exiting.\n", PixelFile);
    proglog(tolog);
    exit(0);
  }

  // Count pixels in the file
  npix = 0;
  while(fgets(line, 512, f)) npix++;
  rewind(f);

  pix_x = float1D(npix, "Cannot allocate memory of pix_x.");
  pix_y = float1D(npix, "Cannot allocate memory for pix_y.");

  fprintf(FpLog, "\tPixel positions on the model are:\n");
  i=0;
  while(fgets(line, 512, f)) {
    sscanf(line, "%f %f", &pix_x[i], &pix_y[i]);
    fprintf(FpLog, "\t%g, %g\n", pix_x[i], pix_y[i]);
    i++;
  }
  fprintf(FpLog, "\n");

  fclose(f);

  sprintf(tolog, "\tPixel positions read from file %s.\n", PixelFile);
  proglog(tolog);
}
  
