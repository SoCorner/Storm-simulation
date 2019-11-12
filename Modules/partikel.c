#include "partikel.h"

struct Partikel* new_particle(int partcount){
    struct Partikel *particles=malloc(partcount*sizeof(struct Partikel));

    for(int i=0;i<partcount;i++){
        particles[i].x=-10;
        particles[i].y=-10;
        particles[i].u=0;
        particles[i].v=0;
    }

    return particles;
}


void printParticles(struct Partikel *particles,int partcount)
{
    for(int i=0;i<partcount;i++){
        printf("status of particles[%d]: x is %lf, y is %lf, u is %lf, v is %lf\n",i,particles[i].x,particles[i].y,particles[i].u,particles[i].v);
    }
}

void ParticleSeed(struct Partikel *particles, REAL posx1, REAL posy1, REAL posx2, REAL posy2,int partcount,int anzahl)
{   
    anzahl=0;
    for(int i=0;i<partcount;i++){
        if(particles[i].x<0||particles[i].y<0){
            anzahl++;
        }
    }
    if(anzahl==0){
        return;
    }else{
        REAL dx=(posx2-posx1)/anzahl, dy=(posy2-posy1)/anzahl;
        //int number=0;
        for(int i=0;i<partcount;i++){
            if(particles[i].x<0||particles[i].y<0){
                particles[i].x=posx1+dx*i;
                particles[i].y=posy1+dy*i;
            }
    }
    
    }
    
}


void ParticleTransport(struct Partikel *particles, REAL dt, int partcount, REAL xlength, REAL ylength,int* nr )
{   
   // printf("ok1\n");
   for (int k=0; k<partcount;k++){
        if(particles[k].x<(double)0|| particles[k].y<(double)0 || particles[k].x>xlength || particles[k].y>ylength){
            //printf("ok_%d\n",k);
            particles[k].x=-10;
            particles[k].y=-10;
            particles[k].u=0;
            particles[k].v=0;
            nr++;

        }
    }
    for (int k=0;k<partcount;k++){
        if(particles[k].x>=(double)0 && particles[k].y>(double)0){ 
            particles[k].x=particles[k].x+ dt*particles[k].u;
            particles[k].y=particles[k].y+ dt*particles[k].v;
        }
    }
   
}



int ParticleVelocity(REAL**U, REAL**V, struct Partikel* particles, REAL delx, REAL dely, int imax, int jmax, int partcount)
{
    for (int k=0;k<partcount;k++){
        if(particles[k].x>=-(double)1 && particles[k].y>=(double)0){
            int i=((int)(particles[k].x/ delx))+1;
            int j=(int)((particles[k].y+( dely/(double)2))/ dely)+1;
            if(i<=0||j<=0){
                return 11;
            }
            REAL x1=(i-1)* delx;
            REAL x2=i* delx;
            REAL y1=((j-1)-0.5)* dely;
            REAL y2=(j-0.5)* dely;
            REAL u1= U[i-1][j-1];
            REAL u2= U[i][j-1];
            REAL u3= U[i-1][j];
            REAL u4= U[i][j];
            particles[k].u=((x2-particles[k].x)*(y2-particles[k].y)*u1
                            +(particles[k].x-x1)*(y2-particles[k].y)*u2
                            +(x2-particles[k].x)*(particles[k].y-y1)*u3
                            +(particles[k].x-x1)*(particles[k].y-y1)*u4)/( delx* dely);
        }
    }
    for (int k=0;k<partcount;k++){
        if(particles[k].x>=-(double)1 && particles[k].y>=(double)0){
            int i=(int)((particles[k].x+( delx/(double)2))/ delx)+1;
            int j=(int)(particles[k].y/ dely)+1;
             if(i<=0||j<=0){
                return 21;
            }
            REAL x1=((i-1)-0.5)* delx;
            REAL x2=(i-0.5)* delx;
            REAL y1=(j-1)* dely;
            REAL y2=j* dely;
            REAL v1= V[i-1][j-1];
            REAL v2= V[i][j-1];
            REAL v3= V[i-1][j];
            REAL v4= V[i][j];
            particles[k].v=((x2-particles[k].x)*(y2-particles[k].y)*v1
                            +(particles[k].x-x1)*(y2-particles[k].y)*v2
                            +(x2-particles[k].x)*(particles[k].y-y1)*v3
                            +(particles[k].x-x1)*(particles[k].y-y1)*v4)/( delx* dely);
        }
    }
}
void writeVTKfileForParticles(const char *fileName, int partcount, struct Partikel *particles)
{
    FILE* vtkFile = fopen(fileName, "w");
    if (vtkFile == NULL)
    {
        printf("Fehler beim Oeffnen der Datei.");
    }
    fprintf(vtkFile, "# vtk DataFile Version 3.0\n");
    fprintf(vtkFile,"Partikel\n");
    fprintf(vtkFile,"ASCII\n");
    fprintf(vtkFile,"DATASET POLYDATA\n");
    fprintf(vtkFile,"POINTS  %d  double\n",partcount);
    for(int k=0;k<partcount;k++){
        fprintf(vtkFile,"%e %e 0.0\n",particles[k].x,particles[k].y);
    }
    fclose(vtkFile);

}
int countPar(struct Partikel *particles, double xlength, double ylength,double delx, double dely, int partcount, int anzahl)
{   
    anzahl =0;
    for(int i=0;i<partcount;i++){
        if(particles[i].x>=xlength-delx||particles[i].y>=2-dely){
            anzahl++;
        }
    }
    return anzahl;
}