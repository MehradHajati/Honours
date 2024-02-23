#include "Amoeba.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define DIMENSIONS 12 // number of parameters
#define UpperBound 1000
#define LowerBound -1000

void runAmoeba(BandContrast *bcMeasured, AFMData afm, BandContrast *bcTilted, BandContrastAFMMapper *bcAFMmOut, double mStdDev, double simStdDev, double *asbs[], int fitLevel) {
    double **p, *y;  // initial gueess as simplex corners
    int nfunk, row, col, ndim;
    double best[DIMENSIONS];  

    // making sure that the input of user is within reason
    if(fitLevel < 1) fitLevel = 1;
    else if(fitLevel > 6) fitLevel = 6;
    
    // creating the array which will hold the values, i think
    p = constructP(asbs, fitLevel);
    ndim = fitLevel * 2;
    y = (double*)malloc(sizeof(double) * (ndim));

    printf("calling annealing\n");

    simulated_annealing(bcMeasured, afm, bcTilted, bcAFMmOut, mStdDev, simStdDev, p[0], asbs, ndim, 0.01, 0.01, 100, best);

    

    y[0] = amoeba_chisq(bcMeasured, afm, bcTilted, bcAFMmOut, mStdDev, simStdDev, p[0], asbs, ndim);

    for(row = 1; row < ndim + 1; row++) {
        for(col = 0; col < ndim; col++) {
            p[row][col] = p[0][col];
        }
        p[row][row-1] *= LAMBDA;
        y[row] = amoeba_chisq(bcMeasured, afm, bcTilted, bcAFMmOut, mStdDev, simStdDev, p[row], asbs, ndim);
    }

    printSimplex(p, y, ndim);

    int ilo = amoeba(p,y,ndim,FTOL,amoeba_chisq,&nfunk, bcMeasured, afm, bcTilted, bcAFMmOut, mStdDev, simStdDev, asbs);

    // update asbs
    for(row = 0; row < fitLevel; row++){
        *(asbs[row]) = p[ilo][2 * row];

        if(row == 0 || row == 5){
            *(asbs[row + 6]) = p[ilo][2 * row + 1];
        }
        else{
            *(asbs[row + 6 + ((row & 1) == 1 ? 1 : -1)]) = p[ilo][2 * row + 1];
        }
    }

    for(row = 0; row < ndim + 1; row++){
        free(p[row]);
    }
    free(p);
    free(y);
}

double amoeba_chisq(BandContrast *bcMeasured, AFMData afm, BandContrast *bcTilted, BandContrastAFMMapper *bcAFMmOut, double mStdDev, double simStdDev, double *simplexCorner, double *asbs[], int ndim){
    double chiSquared = 0.0, a0, a1, a2, a3, a4, a5, b0, b1, b2, b3, b4, b5;
    int row, col, overlappingPoints = 0, fitLevel;
    fitLevel = ndim / 2;

    // Should default as and bs be current as and bs even if not fitting those parameters. YES!
    a0 = simplexCorner[0];
    a1 = (fitLevel >= 2 ? simplexCorner[2] : *(asbs[1]));
    a2 = (fitLevel >= 3 ? simplexCorner[4] : *(asbs[2]));
    a3 = (fitLevel >= 4 ? simplexCorner[6] : *(asbs[3]));
    a4 = (fitLevel >= 5 ? simplexCorner[8] : *(asbs[4]));
    a5 = (fitLevel == 6 ? simplexCorner[10] : *(asbs[5]));
    
    b0 = simplexCorner[1];
    b1 = (fitLevel >= 3 ? simplexCorner[5] : *(asbs[7]));
    b2 = (fitLevel >= 2 ? simplexCorner[3] : *(asbs[8]));
    b3 = (fitLevel >= 5 ? simplexCorner[9] : *(asbs[9]));
    b4 = (fitLevel >= 4 ? simplexCorner[7] : *(asbs[10]));
    b5 = (fitLevel == 6 ? simplexCorner[11] : *(asbs[11]));

    // here it is creating the mapping to be compared and find the difference
    *bcAFMmOut = bandContrastAFMMapper_map(bcMeasured, afm, a0, a1, a2, a3, a4, a5, b0, b1, b2, b3, b4, b5); // Scales set to 1
    // what to about the scaling here?
    bandContrast_scaleTo255(&bcAFMmOut->map[GREYSCALE_LAYER], bcAFMmOut->nrow, bcAFMmOut->ncol);

    // objective function is here
    for(row = 0; row < bcAFMmOut->nrow; row++){
        for(col = 0; col < bcAFMmOut->ncol; col++){
            // what to do about the transparency here
            if(bcAFMmOut->map[GREYSCALE_LAYER][row][col] < GREYSCALE_DEFAULT * 255.0){ // Transparency
                // Chi Squared and difference here?
                chiSquared += (bcAFMmOut->map[GREYSCALE_LAYER][row][col] - bcTilted->greyScale[row][col]) * (bcAFMmOut->map[GREYSCALE_LAYER][row][col] - bcTilted->greyScale[row][col]) / sqrt((mStdDev*mStdDev)*(simStdDev*simStdDev));
                overlappingPoints++;
            }
        }
    }

    if(overlappingPoints != 0) chiSquared /= (double)overlappingPoints;
    else chiSquared = 99999999;
    bandContrastAFMMapper_free(bcAFMmOut);
    return chiSquared;

}

int amoeba(double **p, double *y, int ndim, double ftol,double (*funk)(), int *nfunk, BandContrast *bcMeasured, AFMData afm, BandContrast *bcTilted, BandContrastAFMMapper *bcAFMmOut, double mStdDev, double simStdDev, double *asbs[]){
	int i,j,ilo,ihi,inhi,mpts=ndim+1;
	double ytry,ysave,sum,rtol,*psum;

    //printf("ftol = %f\n",ftol);
    //printf("ndim = %d\n",ndim);
    //printf("mpts = %d \n", mpts);
    fflush(stdout);
	psum = (double*)malloc(sizeof(double) * ndim);
	*nfunk=0;

    // GET_PSUM
    for (j = 0; j < ndim; j++) { 
        sum = 0.0;
        for (i = 0; i < mpts; i++) {
            sum += p[i][j]; 
            psum[j] = sum;
        }
    }
	for (;;) {
		ilo=0;
		ihi = y[0] > y[1] ? (inhi = 1,0) : (inhi = 0,1);
		for (i = 0; i < mpts; i++) {
			if (y[i] < y[ilo]) ilo=i;
			if (y[i] > y[ihi]) {
				inhi=ihi;
				ihi=i;
			} else if (y[i] > y[inhi] && i != ihi){
                inhi = i;
            } 
		}
        printSimplex(p, y, ndim);
		rtol = 2.0 * fabs(y[ihi] - y[ilo]) / (fabs(y[ihi]) + fabs(y[ilo]) + TINY);
        printf("rtol: %f\n", rtol);
		if (rtol < ftol) return ilo;
		if (*nfunk >= NMAX) {
			printf("nfunk = %d\n",*nfunk);
            return ilo;
		}
        *nfunk += 2;
		ytry = amotry(p, y, psum, ndim, funk, ihi, nfunk, -ALPHA, bcMeasured, afm, bcTilted, bcAFMmOut, mStdDev, simStdDev, asbs);
		if (ytry <= y[ilo]){
            ytry = amotry(p, y, psum, ndim, funk, ihi, nfunk, GAMMA, bcMeasured, afm, bcTilted, bcAFMmOut, mStdDev, simStdDev, asbs);
        }
		else if (ytry >= y[inhi]) {
			ysave = y[ihi];
			ytry = amotry(p, y, psum, ndim, funk, ihi, nfunk, BETA, bcMeasured, afm, bcTilted, bcAFMmOut, mStdDev, simStdDev, asbs);
			if (ytry >= ysave) {
				for (i = 0; i < mpts; i++) {
					if (i != ilo) {
						for (j = 0; j < ndim; j++) {
							psum[j] = 0.5 * (p[i][j] + p[ilo][j]);
							p[i][j] = psum[j];
						}
						y[i]=(*funk)(bcMeasured, afm, bcTilted, bcAFMmOut, mStdDev, simStdDev, psum, asbs, ndim);
					}
				}
				*nfunk += ndim;
                // GET_PSUM
				for (j = 0; j < ndim; j++) { 
                    sum = 0.0;
                    for (i = 0; i < mpts; i++) {
                        sum += p[i][j]; 
                        psum[j] = sum;
                    }
                }
			}
		}
        else --(*nfunk);
	}
    free(psum);
}

double amotry(double **p, double *y, double *psum,int ndim, double (*funk)(), int ihi, int *nfunk,double fac, BandContrast *bcMeasured, AFMData afm, BandContrast *bcTilted, BandContrastAFMMapper *bcAFMmOut, double mStdDev, double simStdDev, double *asbs[]){
	int col;
	double fac1,fac2,ytry,*ptry;

	ptry = (double*)malloc(sizeof(double) * ndim);
	fac1 = (1.0 - fac) / ndim;
	fac2 = fac1 - fac;
	for (col = 0; col < ndim; col++) {
        ptry[col] = psum[col] * fac1 - p[ihi][col] * fac2;
    }
	ytry = (*funk)(bcMeasured, afm, bcTilted, bcAFMmOut, mStdDev, simStdDev, ptry, asbs, ndim);
	if (ytry < y[ihi]) {
		y[ihi]=ytry;
		for (col = 0; col < ndim; col++) {
			psum[col] += ptry[col] - p[ihi][col];
			p[ihi][col] = ptry[col];
		}
	}
	free(ptry);
	return ytry;

}

/**
 * @brief 
 * 
 * @param p 
 * @param y 
 * @param ndim 
 */
void printSimplex(double **p, double *y, int ndim){    
    int row, col;

    // show the optimized simplex
    for(row = 0; row < ndim + 1; row++) {
        printf("row=%d: ", row );
        for(col = 0; col < ndim; col++) {
            printf("%f ", p[row][col]);
        }
        printf("value = %f\n", y[row]);
    }
}

/**
 * @brief not sure what it is doing here 
 * 
 * @param asbs 
 * @param fitLevel 
 * @return double** 
 */
double **constructP(double *asbs[], int fitLevel){
    int row, col;
    int ndim = fitLevel * 2;

    // creating a 2d double array
    double **p = (double**)malloc(sizeof(double) * (ndim + 1));

    // creating 1d arrays in the indices of the 2d array created above
    for(row = 0; row < ndim + 1; row++){
        p[row] = (double*)malloc(sizeof(double) * ndim);
    }

    for(col = 0; col < fitLevel; col++){
        // set a_i
        p[0][2 * col] = *(asbs[col]);
        // set b_i
        if(col == 0 || col == 5){
            p[0][2 * col + 1] = *(asbs[col + 6]);
        }
        else{
            p[0][2 * col + 1] = *(asbs[col + 6 + ((col & 1) == 1 ? 1 : -1)]);
        }
    }

    return p;
}

// Generate a neighbor solution
void get_neighbor(double *current_solution[], double *new_solution[], double temp) {
    for (int i = 0; i < DIMENSIONS; i++) {
        double adjustment = (rand() / (double)RAND_MAX) - 0.5;
        double fresh = *(current_solution[i]) + adjustment * temp;
        new_solution[i] = &fresh;
        double lowerBound = -1000;
        double upperBound = 1000;
        if (*(new_solution[i]) < LowerBound) new_solution[i] = &lowerBound;
        if (*(new_solution[i]) > UpperBound) new_solution[i] = &upperBound;
    }
}

// Simulated annealing algorithm
void simulated_annealing(BandContrast *bcMeasured, AFMData afm, BandContrast *bcTilted, BandContrastAFMMapper *bcAFMmOut, double mStdDev, double simStdDev, double *simplexCorner, double *asbs[], int ndim, double start_temp, double end_temp, double cooling_rate, double *best_solution) {
    srand(time(NULL)); // Seed the random number generator

    double temp = start_temp;
    double s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12;
    double n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12;
    double lowerBound = -1000;
    double upperBound = 1000;
    double *current_solution[DIMENSIONS] = {&s1, &s2, &s3, &s4, &s5, &s6, &s7, &s8, &s9, &s10, &s11, &s12};
    double *new_solution[DIMENSIONS] = {&n1, &n2, &n3, &n4, &n5, &n6, &n7, &n8, &n9, &n10, &n11, &n12};

    // Initialize current solution with random values between -10 and 10
    for (int i = 0; i < DIMENSIONS; i++) {
        double adj = (rand() / (double)RAND_MAX) * 20 - 10;
        current_solution[i] = &adj;
    }

    printf("made it to before the call");

    double current_energy = amoeba_chisq(bcMeasured, afm, bcTilted, bcAFMmOut, mStdDev, simStdDev, simplexCorner, current_solution, ndim);
    printf("made it to after the call");
    while (temp > end_temp) {
        get_neighbor(current_solution, new_solution, temp);
        double new_energy = amoeba_chisq(bcMeasured, afm, bcTilted, bcAFMmOut, mStdDev, simStdDev, simplexCorner, current_solution, ndim);

        // Calculate change in energy
        double energy_change = new_energy - current_energy;

        // Decide if we should accept the new solution
        if (energy_change < 0) {
            // New solution is better, accept it
            for (int i = 0; i < DIMENSIONS; i++) {
                current_solution[i] = new_solution[i];
            }
            current_energy = new_energy;
        } else if (exp(-energy_change / temp) > (rand() / (double)RAND_MAX)) {
            // New solution is worse, accept it with a probability decreasing with temperature
            for (int i = 0; i < DIMENSIONS; i++) {
                current_solution[i] = new_solution[i];
            }
            current_energy = new_energy;
        }

        // Cool down
        temp *= 1 - cooling_rate;
    }

    // Copy the best solution
    for (int i = 0; i < DIMENSIONS; i++) {
        best_solution[i] = *(current_solution[i]);
    }
    printf("Best solution:\n");
    for (int i = 0; i < DIMENSIONS; i++) {
        printf("x[%d] = %f\n", i, best_solution[i]);
    }
}

/*int main() {
    double start_temp = 100.0;
    double end_temp = 0.01;
    double cooling_rate = 0.01;
    double best_solution[DIMENSIONS];

    simulated_annealing(start_temp, end_temp, cooling_rate, best_solution);

    printf("Best solution:\n");
    for (int i = 0; i < DIMENSIONS; i++) {
        printf("x[%d] = %f\n", i, best_solution[i]);
    }
    printf("Objective function value: %f\n", objective_function(best_solution));

    return 0;
}*/