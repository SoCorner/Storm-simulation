#include "blas.h"
void axpy(REAL alpha, REAL* x, REAL* y, int len)
{
    for(int i=0; i<len; ++i){
        y[i] = alpha*x[i] + y[i];
    }
}

REAL dot(REAL* x, REAL* y, int len)
{
    REAL res=0.0;
    for(int i=0; i<len; ++i){
        res = res + x[i]*y[i];
    }
    return res;
}

REAL nrm2(REAL* x, int len)
{
    return sqrt(dot(x,x,len));
}

void copy(REAL* x, REAL* y, int len)
{
    for(int i=0; i<len; ++i){
        y[i] = x[i];
    }
}

void scal(REAL alpha, REAL* x, int len)
{
    for(int i=0; i<len; ++i){
        x[i] = alpha*x[i];
    }
}

void gemv(REAL alpha, REAL** A, REAL* x, REAL beta, REAL* y, int rows, int cols)
{
    for(int i=0; i<rows; ++i){
        REAL sum=0.0;
        for(int j=0; j<cols; ++j){
            sum = sum + alpha*A[i][j]*x[j];
        }
        sum = sum + beta*y[i];
        y[i]=sum;
    }
}

void scal2Dfield(REAL alpha, REAL** X, int sizeX, int sizeY)
{
    for(int i=0; i<sizeX; ++i){
        for(int j=0; j<sizeY; ++j){
            X[i][j]= alpha*X[i][j];
        }
    }
}

void axpy2Dfield(REAL alpha, REAL** X, REAL** Y, int sizeX, int sizeY)
{
    for(int i=0; i<sizeX; ++i){
        axpy(alpha,X[i],Y[i],sizeY);
    }
}

