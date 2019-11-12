#include "visual.h"
#include "init.h"
#include "partikel.h"

double findeMax(double** mart, int imax, int jmax);
double findeU(struct Partikel* particles, int partcount);
double findeV(struct Partikel* particles, int partcount);
int main(int argc, char* argv[]){
    
    //==========================================| define new variables |==========================================
   
    char* path = "test_new.txt";
    char* filename="plot.txt";
    char* problem = "karman";
    int n =0, flowNr=0;
    int imax, jmax, itermax, wl, wr, wt, wb;
    
    REAL t = 0.0;
    REAL s = 0.0;
    REAL eps, alpha, omg, Re, GX, GY, UI, VI, PI, delt, t_end, tau,delx, dely;
    REAL res = 0, del_vec;
    REAL xlength = 0, ylength=0;

    
    
    READ_PARAMETER_NEW(path,&xlength,&ylength,&imax,&jmax,&delx,&dely,&delt,&t_end,
                        &tau,&itermax,&eps,&omg,&alpha,&Re,&GX,&GY,&UI,&VI,&PI,&wl,
                        &wr,&wt,&wb);
    printf("wl: %d, wr:%d, wt:%d, wb:%d\n",wl,wr,wt,wb);
    REAL** U = createMatrix(imax+2,jmax+2);
    REAL** V = createMatrix(imax+2,jmax+2);
    REAL** P = createMatrix(imax+2,jmax+2);
    REAL** F = createMatrix(imax+2,jmax+2);
    REAL** G = createMatrix(imax+2,jmax+2);
    int** FLAG = Create2DIntegerField(imax+2,jmax+2);
    REAL** RHS = createMatrix(imax+2,jmax+2); 
   
    
    del_vec = 2;
    res = 1 + eps;
    
    //================| initialize the matrix |====================
      
    INIT_UVP(U,V,P,imax,jmax,UI,VI,PI,problem);
    if(strcmp(problem,"stairs")==0){
		for(int i=0;i<imax+2;i++){
           for(int j=0;j<=jmax/2+1;j++){
               U[i][j]=0;
           }
       }
    }
    INIT_FLAG(problem, FLAG,imax,jmax);
   // print2Dfield_int(FLAG,imax+2,jmax+2);
	
    int partcount=10000, anzahl = 10000;
    struct Partikel* particles=new_particle(partcount);
    

    double posx1=0, posy1=1;
    double posx2=0, posy2=1.7;

    ParticleSeed(particles,posx1,posy1,posx2,posy2,partcount,anzahl);
    char* parfile="data//particle_0.vtk";
    writeVTKfileForParticles(parfile,partcount,particles);



    int count =0;
   


   
    //==========================================| main programm without FLAG |==========================================       
    int index;
    int nr=0;
    
    while(t < t_end){
        
        
        COMP_delt_NEW(&delt,imax,jmax,delx,dely,U,V,Re,tau);        
        SETBCOND_NEW(U, V, FLAG, imax, jmax,wl,wr,wt,wb);
        SETSPECBOND(U,V,imax,jmax,problem);
        COMP_FG_NEW(U, V, F, G,  FLAG, imax, 
                     jmax, delt, delx, dely, GX, GY, 
                    alpha, Re);
        
        COMP_RHS(F,G,RHS,imax,jmax,delt,delx,dely);
        
        index=POISSON_NEW(P,RHS,U,V,FLAG,imax,jmax,flowNr, delx,dely,eps,itermax,omg,&res);
        
        ADAP_UV_NEW(U,V,F,G,P,FLAG,imax,jmax,delt,delx,dely);

        ParticleSeed(particles,posx1,posy1,posx2,posy2,partcount,anzahl);
       printf("t is : %lf and  count  is: %d, and max distance is:%lf \n",t,count,findeU(particles,partcount));

        ParticleVelocity_NEW(U,V,particles,delx,dely,imax,jmax,partcount);
        if(count==11||count==21){
            printf("finde u!\n");
           return 1;
        }
        ParticleTransport(particles,delt,partcount,9.5,ylength,&nr); 
        count = countPar(particles, xlength, ylength, 0.3, 0.3, partcount, anzahl);

       
        if(t - s > del_vec){
            OUTPUTVEC(U,V,P,imax,jmax,delx,dely,n);
            parfile=(char*)malloc(255*sizeof(char*));
           sprintf(parfile,"data//particle_%d.vtk",n);
           writeVTKfileForParticles(parfile,partcount,particles);
            s = t;
            n=n+1;
        }
        t = t + delt;
        n = n + 1;
    }
    n=n+1;
    OUTPUTVEC(U,V,P,imax,jmax,delx,dely,n);
    
    char*filename1="for_plot";
    PLOT(U,imax+2, jmax+2, filename1, Re,problem);

    destroyMatrix(U,imax+2);
    destroyMatrix(V,imax+2);
    destroyMatrix(P,imax+2);
    destroyMatrix(F,imax+2);
    destroyMatrix(G,imax+2);
    destroyMatrix(RHS,imax+2);
    
    
    printf("programm is ok!\n");
    return 0;
}

double findeU(struct Partikel* particles, int partcount){
    double max=0;
    for(int i=0;i<partcount;i++){
        if(max<particles[i].x){
            max=particles[i].x;
        }
    }
    return max;
}

double findeV(struct Partikel* particles, int partcount){
    double max=0;
    for(int i=0;i<partcount;i++){
        //if(max>particles[i].x){
            max+=particles[i].y;
        //}
    }
    return max/partcount;
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
