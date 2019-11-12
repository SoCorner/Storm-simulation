#ifndef BLAS_H_
#define BLAS_H_

#include <math.h>
#include "real.h"

void	axpy(REAL alpha, REAL* x, REAL* y, int len); //Computes alpha*x + y and saves result in y
REAL	dot(REAL* x, REAL* y, int len); //Computes inner produnct of x and y
REAL	nrm2(REAL* x, int len); //Computes euclidean norm of x
void	copy(REAL* x, REAL* y, int len); //Copys x onto y
void	scal(REAL alpha, REAL* x, int len); //Computes alpha*x and saves result in x

void	gemv(REAL alpha, REAL** A, REAL* x, REAL beta, REAL* y, int rows, int cols);
//Computes alpha*A*x + beta*y and saves result in y
void	scal2Dfield(REAL alpha, REAL** X, int sizeX, int sizeY); //Scales X by alpha
void	axpy2Dfield(REAL alpha, REAL** X, REAL** Y, int sizeX, int sizeY); //Computes alpha*X + Y and saves result in Y
#endif /* BLAS_H_ */
