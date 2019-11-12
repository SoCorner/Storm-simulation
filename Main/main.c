#include "visual.h"
#include "init.h"
#include "partikel.h"

double findeMax(double** mart, int imax, int jmax);
int main(int argc, char* argv[]){
    
    //==========================================| define new variables |==========================================
   
    char* path = "test.txt";
    char* filename="plot.txt";
    char* problem = "stairs";
    int n =0, flowNr=0;
    int imax, jmax, itermax, wl, wr, wt, wb;
    
    REAL t = 0.0;
    REAL s = 0.0;
    REAL eps, alpha, omg, Re, GX, GY, UI, VI, PI, delt, t_end, tau,delx, dely;
    REAL res = 0, del_vec;
    REAL xlength = 0, ylength=0;

    READ_PARAMETER(path,&xlength,&ylength,&imax,&jmax,&delx,&dely,&delt,&t_end,&tau,&itermax,&eps,&omg,&alpha,&Re,&GX,&GY,&UI,&VI,&PI);
    
    del_vec = 0.5;
    res = 1 + eps;
    
  //==========================================| define new matrix and initialize them |===========================
    REAL** U = createMatrix(imax+2,jmax+2);
    REAL** V = createMatrix(imax+2,jmax+2);
    REAL** P = createMatrix(imax+2,jmax+2);
    REAL** F = createMatrix(imax+2,jmax+2);
    REAL** G = createMatrix(imax+2,jmax+2);
    REAL** RHS = createMatrix(imax+2,jmax+2); 
    REAL** u = createMatrix(imax+2,jmax+2);
    
    int** FLAG = Create2DIntegerField(imax+2,jmax+2);
    
    fill2Dfield(0.0,F,imax+2,jmax+2);
    fill2Dfield(0.0,G,imax+2,jmax+2);
    fill2Dfield(0.0,RHS,imax+2,jmax+2);

    INIT_UVP(U,V,P,imax,jmax,UI,VI,PI,problem);
    INIT_FLAG_NEW(problem, FLAG,imax,jmax);

    int index=0;
    double test;
    
    while(t < t_end){
        printf("t is : %lf and  itermax is: %d. \n",t,index);
        COMP_delt_NEW(&delt,imax,jmax,delx,dely,U,V,Re,tau);
        SETBCOND(U, V, imax,  jmax);
        for(int i=1;i<=imax;i++){
            U[i][jmax+1]=2-U[i][jmax];
        }
        COMP_FG(U,V,F,G,imax,jmax,delt,delx,dely,GX,GY,alpha,Re);
        COMP_RHS(F,G,RHS,imax,jmax,delt,delx,dely);
        index = POISSON_NEW(P,RHS,U,V,FLAG,imax,jmax,flowNr, delx,dely,eps,itermax,omg,&res);
        ADAP_UV_NEW(U,V,F,G,P,FLAG,imax,jmax,delt,delx,dely);
        if(t - s > del_vec){
            OUTPUTVEC(U,V,P,imax,jmax,delx,dely,n);
            s = t;
            n=n+1;
        }
        t = t + delt;
        n = n + 1;
    }
    n=n+1;
    OUTPUTVEC(U,V,P,imax,jmax,delx,dely,n);
    for(int i=0;i<imax+2;++i){
        u[i][0] = U[(imax+1)/2][i];
    }
    
    writePlotfile(U, jmax+2,filename);
    char*filename1="for_plot";
    PLOT(U,imax+2, jmax+2, filename1, Re,problem);
    destroyMatrix(u,imax+2);
    destroyMatrix(U,imax+2);
    destroyMatrix(V,imax+2);
    destroyMatrix(P,imax+2);
    destroyMatrix(F,imax+2);
    destroyMatrix(G,imax+2);
    destroyMatrix(RHS,imax+2);

    printf("programm is ok!\n");
    return 0;
}

double findeMax(double** mart, int imax, int jmax){
    double umax =0;
    for(int i=0;i<imax+2;i++){
        for(int j=0;j<jmax+2;j++){
            if(fabs(mart[i][j+1]-mart[i][j])>umax){
                umax=fabs(mart[i][j]);
            }
        }
    }
    return umax;
}
