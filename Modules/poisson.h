#ifndef POISSON_H
#define POISSON_H

#include <stdlib.h>
#include <stdio.h>
#include "real.h"
#include "fields.h"
#include "blas.h"

REAL** create2DpoissonMatrix(REAL ilength, REAL jlength,int imax, int jmax);
REAL** sampleFDgridOnCellCorners(REAL(*func)(REAL,REAL),REAL ilength,REAL jlength, int imax, int jmax);
void solveSOR(REAL** A, REAL* x, REAL* b, int rows, int cols,REAL omega, REAL epsilon, int itermax);
void applyHomogeneousNeumannBC(REAL** p, int imax, int jmax);
void applyHomogeneousDirichletBC(REAL** p, int imax, int jmax);
REAL** sampleFDgridOnCellCenters(REAL(*func)(REAL,REAL),REAL ilength, REAL jlength, int imax, int jmax);
void solveSORforPoisson(REAL** p, REAL** rhs, REAL omega, REAL epsilon, int itermax, int useNeumannBC,
                        REAL ilength, REAL jlength, int imax, int jmax);

REAL testfunc(REAL x, REAL y);
REAL function(REAL x, REAL y);
int min(int a, int b);

#endif // POISSON_H
