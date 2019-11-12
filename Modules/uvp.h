#ifndef UVP_H
#define UVP_H
#include<math.h>
#include "boundary.h"
#include "poisson.h"


void COMP_delt(REAL* delt, int imax, int jmax, REAL delx, REAL dely, REAL** U, REAL** V, REAL Re, REAL tau);

void COMP_delt_NEW(REAL* delt,int imax, int jmax, REAL delx, REAL dely, REAL** U, REAL** V, REAL Re, REAL tau);
void COMP_RHS(REAL**F, REAL**G, REAL**RHS, int imax, int jmax, REAL delt, REAL delx, REAL dely);

void COMP_FG(REAL**U, REAL**V, REAL**F, REAL**G, int imax, int jmax, REAL delt, REAL delx, REAL dely, REAL GX, REAL GY, REAL alpha, REAL Re);
void COMP_FG_NEW(REAL**U, REAL**V, REAL**F, REAL**G, int** FLAG, int imax, int jmax, REAL delt, REAL delx, REAL dely, REAL GX, REAL GY, REAL alpha, REAL Re);
void COMP_FG_NEW1(REAL**U, REAL**V, REAL**F, REAL**G, int** FLAG, int imax, int jmax, REAL delt, REAL delx, REAL dely, REAL GX, REAL GY, REAL alpha, REAL Re);

void ADAP_UV(REAL**U, REAL**V, REAL**F, REAL**G, REAL**P, int imax, int jmax, REAL delt, REAL delx, REAL dely);
void ADAP_UV_NEW(REAL**U, REAL**V, REAL**F, REAL**G, REAL**P, int**FLAG, int imax, int jmax, REAL delt, REAL delx, REAL dely);

int POISSON(REAL**P, REAL**RHS, REAL**U, REAL**V, int imax, int jmax, REAL dely, REAL delx, REAL eps, int itermax, REAL omg, REAL* res);
int POISSON_NEW(REAL**P, REAL**RHS, REAL**U, REAL**V, int** FLAG, int imax, int jmax, int flowNr, REAL delx, REAL dely, REAL eps, int itermax, REAL omg, REAL *res);

void COMP_PRESS(REAL** P, int** FLAG, int imax, int jmax, REAL delx, REAL dely);

#endif // UVP_H
