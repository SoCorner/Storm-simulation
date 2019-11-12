#include "poisson.h"
#include <math.h>
REAL** create2DpoissonMatrix(REAL ilength, REAL jlength, int imax, int jmax){
    REAL** A;
    A=createMatrix(imax*jmax,imax*jmax);
    fill2Dfield(0.0,A,imax*jmax,imax*jmax);
    REAL dx= ilength/(imax+1);
    REAL dy= jlength/(jmax+1);
    double dx2=dx*dx;
    double dy2=dy*dy;
    for(int i=0;i<imax*jmax;++i){
        for(int j=0;j<jmax*imax;++j){
            if(i==j){
                A[i][j]=2*(1/dx2+1/dy2);
            }
            if(j==jmax+i||j==i-jmax){
                    A[i][j]=-1/(dx*dx);
                }
            if((i+1)%jmax==1){
                if(j==i+1){
                    A[i][j]=-1/(dy*dy);
                }
            }
            if((i+1)%jmax==0){
                if(j==i-1){
                    A[i][j]=-1/(dy*dy);
                }
            }
            if((i+1)%jmax!=0 && (i+1)%jmax!=1){
                if(j==i+1||j==i-1){
                    A[i][j]=-1/(dx*dx);
                }
            }
        }
    }
    return A;
}
REAL** sampleFDgridOnCellCorners(REAL(*func)(REAL,REAL), REAL ilength, REAL jlength, int imax, int jmax){
    REAL** A;
    A=createMatrix(imax,jmax);
    for(int i=0;i<imax;++i){
        for(int j=0;j<jmax;++j){
            A[i][j]=func((i+1)*(ilength/(imax+1)),(j+1)*(jlength/(jmax+1)));
        }
    }
    return A;
}
int min(int a, int b){
    if(a<b){
        return a;
    }
    else return b;
}
void solveSOR(REAL** A, REAL* x, REAL* b, int rows, int cols, REAL omega, REAL epsilon, int itermax){
    REAL* y;
    y=createVector(cols);
    REAL* res;
    res=createVector(cols);
    REAL sum=0.0;
    for(int m=1;m<itermax;++m){
        sum=0.0;
        for(int k=0;k<cols;++k){
            if(k==0){
                for(int i=0;i<min(rows,cols);++i){
                y[i]=x[i];
                x[i]= (1/A[i][i])*b[i];
                }
            }
            else{
                for(int i=0;i<cols;++i){
                    sum=0.0;
                    for(int j=0;j<k-1;++j){
                    sum = sum + A[i][j]*x[j];
                    }
                    for(int l=i+1;l<cols;++l){
                    sum = sum + A[i][l]*y[l];
                    }
                }
                for(int i=0;i<min(cols,rows);++i){
                    y[i]=x[i];
                    x[i]= (1/A[i][i])*(b[i]-sum);
                    x[i]= y[i] + omega*(x[i]-y[i]);
                    res[i]=x[i]-y[i];
                }
        }

        if(nrm2(res,cols)<epsilon) {
            break;
        }
    }
    }
}
void applyHomogeneousNeumannBC(REAL** p, int imax, int jmax){
    for(int i=1;i<imax+1;++i){
        p[i][0]=p[i][1];
        p[i][jmax+1]=p[i][jmax];
    }
    for(int j=1;j<jmax+1;++j){
        p[0][j]=p[1][j];
        p[imax+1][j]=p[imax][j];
    }
}
void applyHomogeneousDirichletBC(REAL** p, int imax, int jmax){
    for(int i=1;i<imax+1;++i){
        p[i][0]=-p[i][1];
        p[i][jmax+1]=-p[i][jmax];
    }
    for(int j=1;j<jmax+1;++j){
        p[0][j]=-p[1][j];
        p[imax+1][j]=-p[imax][j];
    }
}
REAL** sampleFDgridOnCellCenters(REAL(*func)(REAL,REAL),REAL ilength, REAL jlength, int imax, int jmax){
    REAL**A;
    A=createMatrix(imax+2,jmax+2);
    for(int i=0;i<imax+2;++i){
        for(int j=0;j<jmax+2;++j){
            A[i][j]=func((i-0.5)*(ilength/(imax+2)),(j-0.5)*(jlength/(jmax+2)));
        }
    }
    return A;
}
void solveSORforPoisson(REAL** p, REAL** rhs, REAL omega, REAL epsilon, int itermax,
                        int useNeumannBC, REAL ilength, REAL jlength, int imax, int jmax){
    REAL res=0.0;
    REAL dx= ilength/(imax+2);
    REAL dy= jlength/(jmax+2);
    for(int m=1;m<itermax;++m){
        for(int i=1;i<imax+1;++i){
            if(useNeumannBC==1){
                applyHomogeneousNeumannBC(p,imax,jmax);
            }
            else{
                if(useNeumannBC==0){
                    applyHomogeneousDirichletBC(p,imax,jmax);
                }
                for(int j=1;j<jmax+1;++j){
                    p[i][j]=(1-omega)*p[i][j] + (omega/((1/(dx*dx))+(1/(dy*dy))))*((p[i+1][j]+p[i][j-1])/(dx*dx) + (p[i][j+1]+p[i][j-1])/(dy*dy) - rhs[i][j]);
                }
            }
        }
        res=0.0;
        for(int k=1;k<imax+1;++k){
            for(int l=1;l<jmax+1;++l){
                res=res + pow((p[k+1][l]-2*p[k][l]+p[k-1][l])/(dx*dx) + (p[k][l+1]-2*p[k][l]+p[k][l-1])/(dy*dy),2)/(imax*jmax);
            }
        }
        if(sqrt(res)<epsilon){
            break;
        }
    }
}



