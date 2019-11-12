#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "init.h"

void SETBCOND(REAL** U, REAL** V, int imax, int jmax);
void SETBCOND_NEW(REAL**U, REAL**V, int**FLAG, int imax, int jmax, int wl, int wr, int wt, int wb);
void SETBCOND_NEW1(REAL**U, REAL**V, int** FLAG,int imax, int jmax, int wl, int wr, int wt, int wb,double t);
void SETSPECBOND(REAL** U,  REAL**V, int imax, int jmax, char* problem);
void SETBCOND_NEW_M(REAL**U, REAL**V, int**FLAG, int imax, int jmax, int wl, int wr, int wt, int wb);

#endif // BOUNDARY_H
