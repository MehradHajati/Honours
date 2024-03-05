#ifndef AMOEBA_H
#define AMOEBA_H


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define DIMENSIONS 2 // number of parameters
#define NMAX 500
#define ALPHA 1.0
#define BETA 0.5
#define GAMMA 2.0

#define TINY 1.0e-10
#define LAMBDA 1.2
#define FTOL 1e-4

// Particle structure
typedef struct {
    double position[DIMENSIONS];
    double velocity[DIMENSIONS];
    double best_position[DIMENSIONS];
    double best_value;
} Particle;

#define GET_PSUM for (j=1;j<=ndim;j++) { for (i=1,sum=0.0;i<=mpts;i++) sum += p[i][j]; psum[j]=sum;}

void main();

// void nrerror(char *error_text);

// double *vector(int nl, int nh);

// void free_vector(double *v, int nl, int nh);

// double **matrix(int nrl, int nrh, int ncl, int nch);

// void free_matrix(double **m, int nrl, int nrh, int ncl, int nch);

void getNeighbor(double *current_solution, double *new_solution, double bounds[][2]);

double simulatedAnnealing(double cooling_rate, double bounds[][2]);

double objectiveFunction(double *current_solution);

void checkBounds(double *new_soliution, double bounds[][2]);

double randBounds(double bounds[2]);

double particleSwarm(double bounds[][2]);

void moveParticle(Particle *particle, double *global_best_position, double bounds[][2]);

void createParticle(Particle *particle, double bounds[][2]);

#endif // AMOEBA_H