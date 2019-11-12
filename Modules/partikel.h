#ifndef PARTIKEL_H_
#define PARTIKEL_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "real.h"

struct Partikel{
    REAL x,y;
    REAL u,v;
};

struct Partikel* new_particle(int partcount);
void printParticles(struct Partikel *particles,int partcount);
void ParticleSeed(struct Partikel *particles, REAL posx1, REAL posy1, REAL posx2, REAL posy2, int partcount, int anzahl);
int countPar(struct Partikel *particles, double xlength, double ylength, double delx, double dely,int partcount, int anzahl);

int ParticleVelocity(REAL**U, REAL**V, struct Partikel* particles, REAL delx, REAL dely, int imax, int jmax, int partcount);
void ParticleTransport(struct Partikel *particles, REAL dt, int partcount, REAL xlength, REAL ylength,int* nr );
void writeVTKfileForParticles(const char* fileName, int partcount, struct Partikel *particles);
#endif
