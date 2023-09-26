#ifndef AMOEBA_H
#define AMOEBA_H

#include "BandContrastAFMMapper.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define NMAX 1000
#define ALPHA 1.0
#define BETA 0.5
#define GAMMA 2.0

#define TINY 1.0e-10
#define FTOL 1e-4

#define STRETCH 1.2
#define PUSH 0.1
#define TWIST 0.0001

#define GET_PSUM for (j=1;j<=ndim;j++) { for (i=1,sum=0.0;i<=mpts;i++) sum += p[i][j]; psum[j]=sum;}

void runAmoeba(BandContrast *bcMeasured, AFMData afm, BandContrast *bcTilted, BandContrastAFMMapper *bcAFMmOut, double mStdDev, double simStdDev, double *asbs[], int fitLevel);

double amoeba_chisq(BandContrast *bcMeasured, AFMData afm, BandContrast *bcTilted, BandContrastAFMMapper *bcAFMmOut, double mStdDev, double simStdDev, double *params, double *asbs[], int ndim);

int amoeba(double **p, double *y, int ndim, double ftol,double (*funk)(), int *nfunk, BandContrast *bcMeasured, AFMData afm, BandContrast *bcTilted, BandContrastAFMMapper *bcAFMmOut, double mStdDev, double simStdDev, double *asbs[]);

double amotry(double **p, double *y, double *psum,int ndim, double (*funk)(), int ihi, int *nfunk,double fac, BandContrast *bcMeasured, AFMData afm, BandContrast *bcTilted, BandContrastAFMMapper *bcAFMmOut, double mStdDev, double simStdDev, double *asbs[]);

void nrerror(char *error_text);

double *vector(int nl, int nh);

void free_vector(double *v, int nl, int nh);

double **matrix(int nrl, int nrh, int ncl, int nch);

void printSimplex(double **p, double *y, int ndim);

double **constructP(double *asbs[], int fitLevel);

#endif // AMOEBA_H
