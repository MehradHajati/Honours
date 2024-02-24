#include "Amoeba.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define DIMENSIONS 12 // number of parameters
#define SWARM_SIZE 30 // Number of particles in the swarm
#define MAX_ITERATIONS 200 // Maximum number of iterations
#define PHI_P 2.05 // ratio for personal particle component
#define PHI_S 2.05 // ratio for Social swarm component

// Particle structure
typedef struct {
    double position[DIMENSIONS];
    double velocity[DIMENSIONS];
    double best_position[DIMENSIONS];
    double best_value;
} Particle;


void runAmoeba(BandContrast *bcMeasured, AFMData afm, BandContrast *bcTilted, BandContrastAFMMapper *bcAFMmOut, double mStdDev, double simStdDev, double *asbs[], int fitLevel) {
    double **p, *y;  // initial gueess as simplex corners
    int nfunk, row, col, ndim;
    double best[DIMENSIONS];  

    // creating the upper and lower bounds
    double bounds[DIMENSIONS][2] = {{600, 700}, {3.2, 3.3}, {-0.2, 0.2}, {-1e-4, 1e-4}, {-1e-4, 1e-4}, {-1e-4, 1e-4}, {610, 680}, {-0.2, 0.2}, {1, 1.3}, {-1e-4, 1e-4}, {-1e-4, 1e-4}, {-1e-4, 1e-4}};
    //double upperBounds[DIMENSIONS] = {640, 3.5, 0.003, 5e-6, 5e-6, 5e-6, 650, 0.003, 2.8, 5e-6, 5e-6, 5e-6};
    //double lowerBounds[DIMENSIONS] = {600, 3, 0, 0, 0, 0, 630, 0, 2, 0, 0, 0};

    // making sure that the input of user is within reason
    if(fitLevel < 1) fitLevel = 1;
    else if(fitLevel > 6) fitLevel = 6;
    
    // creating the array which will hold the values, i think
    p = constructP(asbs, fitLevel);
    ndim = fitLevel * 2;
    y = (double*)malloc(sizeof(double) * (ndim));

    printf("calling annealing\n");

    simulatedAnnealing(bcMeasured, afm, bcTilted, bcAFMmOut, mStdDev, simStdDev, 20, bounds);

    
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
void getNeighbor(double *current_solution, double *new_solution) {
    for (int i = 0; i < DIMENSIONS; i++) {
        double r = (double)rand() / (double)RAND_MAX;
        double adjustment;
        if(i == 0){
            adjustment = 10;
        }
        else if(i == 6){
            adjustment = 0.7;
        }
        else if(i == 1){
            adjustment = 0.01;
        }
        else if(i == 8){
            adjustment = 0.03;
        }
        else if(i == 2 || i == 7){
            adjustment = 0.04;
        }
        else{
            adjustment = 0.00002;
        }
        new_solution[i] = current_solution[i] + adjustment * (r-0.5);
    }
}

/**
 * @brief Method to check if the newly/randomly generated solution is within the bounds
 * 
 * @param new_solution The solution given by the get_neighbor method wchich is partially random and depends on the last solution
 * @param bounds A 2d array of bounds, 12 by 2, the first column defines the lower bound and the second column defines the upper bound
 */
void checkBounds(double *new_solution, double bounds[][2]) {
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
void simulatedAnnealing(BandContrast *bcMeasured, AFMData afm, BandContrast *bcTilted, BandContrastAFMMapper *bcAFMmOut, double mStdDev, double simStdDev, double cooling_rate, double bounds[][2]){
    
    // Seed the random number generator
    srand(time(NULL)); 

    double temp = 1;
    double current_solution[DIMENSIONS] = {640.0000, 3.2201, 0.0016, 2.44141e-06, 2.44141e-06, 2.44141e-06, 640.0000, 0.0016, 2.3988, 2.44141e-06, 2.44141e-06, 2.44141e-06};
    double new_solution[DIMENSIONS];
    int iter;

    double current_energy = objectiveFunction(bcMeasured, afm, bcTilted, bcAFMmOut, mStdDev, simStdDev, current_solution);

    for(iter = 1; iter < MAX_ITERATIONS; iter++) {
        printf("Current Temp is: %f and number of iteration is: %d\n", temp, iter);
        getNeighbor(current_solution, new_solution);
        checkBounds(new_solution, bounds);
        double new_energy = objectiveFunction(bcMeasured, afm, bcTilted, bcAFMmOut, mStdDev, simStdDev, new_solution);

        // Decide if we should accept the new solution
        if ( new_energy < current_energy) {
            // New solution is better, accept it
            for (int i = 0; i < DIMENSIONS; i++) {
                current_solution[i] = new_solution[i];
            }
            current_energy = new_energy;
        } else{ // New solution is worse, accept it with a probability decreasing with temp
            double r = (double)rand() / (double)RAND_MAX;
            double e = exp( (current_energy-new_energy) / temp );
            if(e > r) {
                printf("Accepting the worst solution\n");
                for (int i = 0; i < DIMENSIONS; i++) {
                    current_solution[i] = new_solution[i];
                }
                current_energy = new_energy;
            } 
        }

        // Cool down
        temp = cooling_rate / log(iter + 2);
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
double objectiveFunction(BandContrast *bcMeasured, AFMData afm, BandContrast *bcTilted, BandContrastAFMMapper *bcAFMmOut, double mStdDev, double simStdDev, double *solution){

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

    if(overlappingPoints != 0) chiSquared /= (double)overlappingPoints;
    else chiSquared = 99999999;
    bandContrastAFMMapper_free(bcAFMmOut);
    printf("chi-squared is: %f\n", chiSquared);
    return chiSquared;
}


/**
 * @brief Function to return random number between the upper and lower bounds
 * 
 * @param bounds array of size the first one is the lower bound and the second one is the upper bound
 * @return double random value between the upper and lower bounds
 */
double randBounds(double bounds[2]){
    double range = bounds[0] - bounds[1];
    double div = RAND_MAX / range;
    return bounds[1] + ( (double)rand() / div);
}


/**
 * @brief Create a Particle object
 * 
 * @param particle Memory locatio for the particle
 * @param bounds upper and lower bounds of the search space
 */
void createParticle(Particle *particle, double bounds[][2]){
    for (int i = 0; i < DIMENSIONS; i++) {
        // randomly assign a value for the position of particle between the lower and upper bounds
        particle->position[i] = randBounds(bounds[i]);
        // randomly assign a value between -0.5 and 0.5 for the velocity of the particle
        particle->velocity[i] = ( (double)rand() / (double)RAND_MAX) - 0.5;
        particle->best_position[i] = particle->position[i];
    }
}

// Update velocity and position of the particle
void moveParticle(Particle *particle, double *global_best_position) {
    for (int i = 0; i < DIMENSIONS; i++) {
        double r = (double)rand() / (double)RAND_MAX;
        double adjustment;
        if(i == 0){
            adjustment = 10;
        }
        else if(i == 6){
            adjustment = 0.7;
        }
        else if(i == 1){
            adjustment = 0.01;
        }
        else if(i == 8){
            adjustment = 0.03;
        }
        else if(i == 2 || i == 7){
            adjustment = 0.04;
        }
        else{
            adjustment = 0.00002;
        }

        double personal_velocity = PHI_P * adjustment * r * (particle->best_position[i] - particle->position[i]);
        double social_velocity = PHI_S * adjustment * r * (global_best_position[i] - particle->position[i]);
        particle->velocity[i] = personal_velocity + social_velocity;
        particle->position[i] += particle->velocity[i];
    }
}

// Particle Swarm Optimization algorithm
void particleSwarm(BandContrast *bcMeasured, AFMData afm, BandContrast *bcTilted, BandContrastAFMMapper *bcAFMmOut, double mStdDev, double simStdDev, double bounds[][2]) {
    
    Particle swarm[SWARM_SIZE];
    double global_best_value = INFINITY;
    double global_best_position[DIMENSIONS];

    // Create the particles
    for (int i = 0; i < SWARM_SIZE; i++) {
        createParticle(&swarm[i], bounds);
        swarm[i].best_value = objectiveFunction(bcMeasured, afm, bcTilted, bcAFMmOut, mStdDev, simStdDev, swarm[i].position);
        if (swarm[i].best_value < global_best_value) {
            global_best_value = swarm[i].best_value;
            for (int j = 0; j < DIMENSIONS; j++) {
                global_best_position[j] = swarm[i].position[j];
            }
        }
    }

    // Looping until MAX Iterations is reached
    for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
        for (int i = 0; i < SWARM_SIZE; i++) {

            // Moving the current particle
            moveParticle(&swarm[i], global_best_position);

            // checking if the particle is still within the search space
            checkBounds(swarm[i].position, bounds);

            // calculating the new objective function value for the particle after moving
            double current_value = objectiveFunction(bcMeasured, afm, bcTilted, bcAFMmOut, mStdDev, simStdDev, swarm[i].position);

            // Checking if the current objective function value is better than particle's best, if yes update best value and best position
            if (current_value < swarm[i].best_value) {
                swarm[i].best_value = current_value;
                for (int j = 0; j < DIMENSIONS; j++) {
                    swarm[i].best_position[j] = swarm[i].position[j];
                }
            }

            // Checking if the current objective function value is better than globla best, if yes update best value and best position
            if (current_value < global_best_value) {
                global_best_value = current_value;
                for (int j = 0; j < DIMENSIONS; j++) {
                    global_best_position[j] = swarm[i].position[j];
                }
            }
        }
    }

    printf("Best solution:\n");
    for (int i = 0; i < DIMENSIONS; i++){
        printf("x[%d] = %f\n", i, global_best_position[i]);
    }
    printf("with chi squared = %f\n", global_best_value);

}