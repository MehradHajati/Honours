#include "Amoeba.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define DIMENSIONS 12 // number of parameters
#define SWARM_SIZE 10 // Number of particles in the swarm
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

    // creating the upper and lower bounds
    double bounds[DIMENSIONS][2] = {{600, 700}, {3.2, 3.3}, {-0.2, 0.2}, {-1e-4, 1e-4}, {-1e-4, 1e-4}, {-1e-4, 1e-4}, {610, 680}, {-0.2, 0.2}, {1, 1.3}, {-1e-4, 1e-4}, {-1e-4, 1e-4}, {-1e-4, 1e-4}};
    //double upperBounds[DIMENSIONS] = {640, 3.5, 0.003, 5e-6, 5e-6, 5e-6, 650, 0.003, 2.8, 5e-6, 5e-6, 5e-6};
    //double lowerBounds[DIMENSIONS] = {600, 3, 0, 0, 0, 0, 630, 0, 2, 0, 0, 0};


    //simulatedAnnealing(bcMeasured, afm, bcTilted, bcAFMmOut, mStdDev, simStdDev, 20, bounds);

    particleSwarm(bcMeasured, afm, bcTilted, bcAFMmOut, mStdDev, simStdDev, bounds);
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
        // particle is not moving to velocity is 0
        particle->velocity[i] = 0;
        particle->best_position[i] = particle->position[i];
    }
}

/**
 * @brief Moves the particle to a new position
 * 
 * @param particle Particle to be moved
 * @param global_best_position Best position found so far
 */
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
        // Calculating the personal velocity based on the best position for the particle
        double personal_velocity = PHI_P * adjustment * r * (particle->best_position[i] - particle->position[i]);
        // Calculating the social velocity based on the global best position
        double social_velocity = PHI_S * adjustment * r * (global_best_position[i] - particle->position[i]);
        // changing the veloctity of the particle
        particle->velocity[i] = personal_velocity + social_velocity;
        // moving the particle
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
        printf("Iteration %d with chi squared = %f\n", iter, global_best_value);
    }

    printf("Best solution:\n");
    for (int i = 0; i < DIMENSIONS; i++){
        printf("x[%d] = %f\n", i, global_best_position[i]);
    }
    printf("with chi squared = %f\n", global_best_value);

}