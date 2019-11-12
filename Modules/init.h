#ifndef INIT_H
#define INIT_H

#include "fields.h"
#include "blas.h"
#include <string.h>

void READ_PARAMETER(char* Inputfile, REAL* xlength, REAL* ylength, int* imax, int* jmax, REAL* delx,
                   REAL* dely, REAL* delt,REAL* t_end, REAL* tau, int* itermax, REAL* eps, REAL* omg,
                   REAL* alpha, REAL* Re, REAL* GX, REAL* GY, REAL* UI, REAL* VI, REAL* PI);
                   
void READ_PARAMETER_NEW(char* Inputfile, REAL* xlength, REAL* ylength, int* imax, int* jmax, 
                    REAL* delx, REAL* dely, REAL* delt, REAL* t_end, REAL* tau, int* itermax,
                    REAL* eps, REAL* omg, REAL* alpha, REAL* Re, REAL* GX, REAL* GY, REAL* UI, 
                    REAL* VI, REAL* PI, int* wl, int*wr, int* wt, int*wb);

void INIT_UVP(REAL** U, REAL** V, REAL** P, int imax,int jmax, REAL UI, REAL VI, REAL PI,char*problem);
void INIT_FLAG(char*problem, int**FLAG, int imax, int jmax);
void INIT_FLAG_NEW(char*problem, int**FLAG, int imax, int jmax);



#endif // INIT_H
