#ifndef VISUAL_H
#define VISUAL_H

#include "uvp.h"

void OUTPUTVEC(REAL **U, REAL **V, REAL **P, int imax, int jmax, REAL delx, REAL dely, int n);
void PLOT(REAL **U,int imax, int row, char*fileName, REAL Re, char*problem);

#endif // VISUAL_H
