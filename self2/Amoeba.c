#include "Amoeba.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define DIMENSIONS 12 // number of parameters
#define SWARM_SIZE 20 // number of particles in the swarm

typedef struct {
    double position[DIMENSIONS];
    double velocity[DIMENSIONS];
    double best_position[DIMENSIONS];
    double best_value;
} Particle;

Particle swarm[SWARM_SIZE];
double global_best_position[DIMENSIONS];
double global_best_value;

void runAmoeba(BandContrast *bcMeasured, AFMData afm, BandContrast *bcTilted, BandContrastAFMMapper *bcAFMmOut, double mStdDev, double simStdDev, double *asbs[], int fitLevel) {
    double **p, *y;  // initial gueess as simplex corners
    int nfunk, row, col, ndim;
    double best[DIMENSIONS];  

    // creating the upper and lower bounds
    double bounds[DIMENSIONS][2] = {{600, 640}, {3, 3.5}, {0, 0.003}, {0, 5e-6}, {0, 5e-6}, {0, 5e-6}, {630, 650}, {0, 0.003}, {2, 2.8}, {0, 5e-6}, {0, 5e-6}, {0, 5e-6}};

    // making sure that the input of user is within reason
    if(fitLevel < 1) fitLevel = 1;
    else if(fitLevel > 6) fitLevel = 6;
    
    // creating the array which will hold the values, i think
    p = constructP(asbs, fitLevel);
    ndim = fitLevel * 2;
    y = (double*)malloc(sizeof(double) * (ndim));

    printf("calling annealing\n");

    simulated_annealing(bcMeasured, afm, bcTilted, bcAFMmOut, mStdDev, simStdDev, 100, 1, 0.01, bounds);

    
    /*
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
    free(y);*/
}

double amoeba_chisq(BandContrast *bcMeasured, AFMData afm, BandContrast *bcTilted, BandContrastAFMMapper *bcAFMmOut, double mStdDev, double simStdDev, double *simplexCorner, double asbs[], int ndim){
    double chiSquared = 0.0, a0, a1, a2, a3, a4, a5, b0, b1, b2, b3, b4, b5;
    int row, col, overlappingPoints = 0, fitLevel;
    fitLevel = ndim / 2;

    // Should default as and bs be current as and bs even if not fitting those parameters. YES!
    a0 = simplexCorner[0];
    a1 = (fitLevel >= 2 ? simplexCorner[2] : (asbs[1]));
    a2 = (fitLevel >= 3 ? simplexCorner[4] : (asbs[2]));
    a3 = (fitLevel >= 4 ? simplexCorner[6] : (asbs[3]));
    a4 = (fitLevel >= 5 ? simplexCorner[8] : (asbs[4]));
    a5 = (fitLevel == 6 ? simplexCorner[10] : (asbs[5]));
    
    b0 = simplexCorner[1];
    b1 = (fitLevel >= 3 ? simplexCorner[5] : (asbs[7]));
    b2 = (fitLevel >= 2 ? simplexCorner[3] : (asbs[8]));
    b3 = (fitLevel >= 5 ? simplexCorner[9] : (asbs[9]));
    b4 = (fitLevel >= 4 ? simplexCorner[7] : (asbs[10]));
    b5 = (fitLevel == 6 ? simplexCorner[11] : (asbs[11]));

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

/*
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

*/
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


/**
 * @brief This function is used to randomly generate a neighbor of the current solution
 * 
 * @param current_solution The current solution
 * @param new_solution The newly randomly generated solution which will be close to the current solution
 * @param bounds The bounds of the of the search space
 */
void get_neighbor(double *current_solution, double *new_solution, double bounds[][]) {
    for (int i = 0; i < DIMENSIONS; i++) {
        double r = (double)rand() / (double)RAND_MAX;
        double adjustment;
        if(i == 0 || i == 6){
            adjustment = 10;
        }
        else if(i == 1 || i == 8){
            adjustment = 0.1;
        }
        else if(i == 2 || i == 7){
            adjustment = 0.0001;
        }
        else{
            adjustment = 0.0000001;
        }
        new_solution[i] = current_solution[i] + adjustment * (r-0.5);
    }
    checkBounds(new_solution, bounds);
}

/**
 * @brief Method to check if the newly/randomly generated solution is within the bounds
 * 
 * @param new_solution The solution given by the get_neighbor method wchich is partially random and depends on the last solution
 * @param bounds A 2d array of bounds, 12 by 2, the first column defines the lower bound and the second column defines the upper bound
 */
void checkBounds(double* new_solution, double *bounds[]){
    for(int i = 0; i < DIMENSIONS; i++){
        if(new_solution[i] < bounds[i][0]){
            printf("This index %d has hit the lower bound\n",i);
            new_solution[i] = bounds[i][0];
        }
        else if(new_solution[i] > bounds[i][1]){
            printf("This index %d has hit the upper bound\n",i);
            new_solution[i] = bounds[i][1];
        }
    }
}

// Simulated annealing algorithm
void simulated_annealing(BandContrast *bcMeasured, AFMData afm, BandContrast *bcTilted, BandContrastAFMMapper *bcAFMmOut, double mStdDev, double simStdDev, double start_temp, double end_temp, double cooling_rate, double *bounds[]){
    srand(time(NULL)); // Seed the random number generator

    double temp = start_temp;
    double current_solution[DIMENSIONS] = {640.0000, 3.2201, 0.0016, 2.44141e-06, 2.44141e-06, 2.44141e-06, 640.0000, 0.0016, 2.3988, 2.44141e-06, 2.44141e-06, 2.44141e-06};
    double new_solution[DIMENSIONS];
    int iter = 1;

    double current_energy = objective_function(bcMeasured, afm, bcTilted, bcAFMmOut, mStdDev, simStdDev, current_solution);

    while (temp > end_temp) {
        printf("Current Temp is: %f\n and number of iteration is: %d", temp, iter);
        get_neighbor(current_solution, new_solution, bounds);
        double new_energy = objective_function(bcMeasured, afm, bcTilted, bcAFMmOut, mStdDev, simStdDev, new_solution);

        // Calculate change in energy
        double energy_change = new_energy - current_energy;

        // Decide if we should accept the new solution
        if (energy_change < 0) {
            // New solution is better, accept it
            for (int i = 0; i < DIMENSIONS; i++) {
                current_solution[i] = new_solution[i];
            }
            current_energy = new_energy;
        } else{
            double r = (double)rand() / (double)RAND_MAX;
            if (exp(-energy_change / temp) > r) {
                // New solution is worse, accept it with a probability decreasing with temperature
                for (int i = 0; i < DIMENSIONS; i++) {
                    current_solution[i] = new_solution[i];
                }
                current_energy = new_energy;
            } 
        }

        // Cool down
        temp = cooling_rate / log(iter + 2);
        iter++;
    }

    printf("Best solution:\n");
    for (int i = 0; i < DIMENSIONS; i++){
        printf("x[%d] = %f\n", i, current_solution[i]);
    }
    printf("with chi squared = %f\n", current_energy);
}


/**
 * @brief Function that mimicks ameoba_chisq function but is more streamlined for our purposes
 * 
 * @param bcMeasured Band contrast data
 * @param afm AFM data
 * @param bcTilted Simulated Band Contrast data
 * @param bcAFMmOut original band contrast mapped
 * @param mStdDev Standard deviation of the measured data
 * @param simStdDev standard deviation of the simulated data
 * @param solution variables used for the mapping
 * @return double chi squared value with the current solution vector
 */
double objective_function(BandContrast *bcMeasured, AFMData afm, BandContrast *bcTilted, BandContrastAFMMapper *bcAFMmOut, double mStdDev, double simStdDev, double *solution){

    double chiSquared = 0.0;
    int row, col, overlappingPoints = 0;

    // here it is creating the mapping to be compared and find the difference
    *bcAFMmOut = bandContrastAFMMapper_map(bcMeasured, afm, solution[0], solution[1], solution[2], solution[3], solution[4], solution[5], solution[6], solution[7], solution[8], solution[9], solution[10], solution[11]); // Scales set to 1
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

    if(overlappingPoints != 0) {
        chiSquared /= (double)overlappingPoints;
    }
    else {
        chiSquared = 99999999;
    }
    bandContrastAFMMapper_free(bcAFMmOut);
    printf("chi-squared is: %f\n", chiSquared);

    return chiSquared;
}



/* Initialize particles in the swarm
void initialize_swarm() {
    for (int i = 0; i < SWARM_SIZE; i++) {
        for (int j = 0; j < DIMENSIONS; j++) {
            swarm[i].position[j] = ((double)rand() / RAND_MAX) * 20.0 - 10.0; // Random initial position between -10 and 10
            swarm[i].velocity[j] = ((double)rand() / RAND_MAX) - 0.5; // Random initial velocity
            swarm[i].best_position[j] = swarm[i].position[j]; // Initial best position is the starting position
        }
        swarm[i].best_value = objective_function(swarm[i].position);
        if (i == 0 || swarm[i].best_value < global_best_value) {
            global_best_value = swarm[i].best_value;
            for (int j = 0; j < DIMENSIONS; j++) {
                global_best_position[j] = swarm[i].position[j];
            }
        }
    }
}

// Update velocity and position of each particle
void update_particles() {
    double w = INERTIA_WEIGHT;
    for (int i = 0; i < SWARM_SIZE; i++) {
        for (int j = 0; j < DIMENSIONS; j++) {
            double r_p = (double)rand() / RAND_MAX;
            double r_g = (double)rand() / RAND_MAX;
            swarm[i].velocity[j] = w * swarm[i].velocity[j] + 
                                   PHI_P * r_p * (swarm[i].best_position[j] - swarm[i].position[j]) + 
                                   PHI_G * r_g * (global_best_position[j] - swarm[i].position[j]);
            swarm[i].position[j] += swarm[i].velocity[j];
        }
        double value = objective_function(swarm[i].position);
        if (value < swarm[i].best_value) {
            swarm[i].best_value = value;
            for (int j = 0; j < DIMENSIONS; j++) {
                swarm[i].best_position[j] = swarm[i].position[j];
            }
            if (value < global_best_value) {
                global_best_value = value;
                for (int j = 0; j < DIMENSIONS; j++) {
                    global_best_position[j] = swarm[i].position[j];
                }
            }
        }
    }
}

/*int main() {
    srand(time(NULL));
    initialize_swarm();

    for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
        update_particles();
        // Optionally, decay the inertia weight to balance exploration and exploitation
        INERTIA_WEIGHT *= INERTIA_WEIGHT_DECAY;
    }

    printf("Global best position:\n");
    for (int i = 0; i < DIMENSIONS; i++) {
        printf("x[%d] = %f\n", i, global_best_position[i]);
    }
    printf("Objective function value: %f\n", global_best_value);

    return 0;
}*/

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