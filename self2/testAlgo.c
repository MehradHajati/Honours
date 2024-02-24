#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "testAlgo.h"

#define DIMENSIONS 2 // number of parameters
#define SWARM_SIZE 10 // Number of particles in the swarm
#define MAX_ITERATIONS 10000 // Maximum number of iterations
#define PHI_P 0.5 // ratio for personal particle component
#define PHI_S 0.5 // ratio for Social swarm component

// Particle structure
typedef struct {
    double position[DIMENSIONS];
    double velocity[DIMENSIONS];
    double best_position[DIMENSIONS];
    double best_value;
} Particle;


void main(int argc, char *argv[]) { 

    // creating the upper and lower bounds
    double bounds[DIMENSIONS][2] = {{-1000, 1000}, {-1000, 1000}};
    //double upperBounds[DIMENSIONS] = {640, 3.5, 0.003, 5e-6, 5e-6, 5e-6, 650, 0.003, 2.8, 5e-6, 5e-6, 5e-6};
    //double lowerBounds[DIMENSIONS] = {600, 3, 0, 0, 0, 0, 630, 0, 2, 0, 0, 0};


    //simulatedAnnealing(5, bounds);

    particleSwarm(bounds);
}

/**
 * @brief This function is used to randomly generate a neighbor of the current solution
 * 
 * @param current_solution The current solution
 * @param new_solution The newly randomly generated solution which will be close to the current solution
 * @param bounds The bounds of the of the search space
 */
void getNeighbor(double *current_solution, double *new_solution) {
    // Seed the random number generator
    srand(time(NULL)); 
    for (int i = 0; i < DIMENSIONS; i++) {
        double r = (double)rand() / (double)RAND_MAX;
        double adjustment = 10;
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
            printf("This index %d has hit the lower bound with value %f\n",i, new_solution[i]);
            new_solution[i] = bounds[i][0];
        }
        else if(new_solution[i] > bounds[i][1]){
            printf("This index %d has hit the upper bound with value %f\n", i, new_solution[i]);
            new_solution[i] = bounds[i][1];
        }
    }
}

// Simulated annealing algorithm
void simulatedAnnealing(double cooling_rate, double bounds[][2]){

    double temp = 1;
    double current_solution[DIMENSIONS] = {0, 0};
    double new_solution[DIMENSIONS];
    int iter;

    double current_energy = objectiveFunction(current_solution);

    for(iter = 1; iter < MAX_ITERATIONS; iter++) {
        printf("Current Temp is: %f and number of iteration is: %d with chi-squared: %f\n", temp, iter, current_energy);
        getNeighbor(current_solution, new_solution);
        checkBounds(new_solution, bounds);
        double new_energy = objectiveFunction(new_solution);

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

double objectiveFunction(double *solution){

    return 0.5 * (pow(solution[0], 4) - 16*pow(solution[0], 2) + 5*solution[0] + pow(solution[1], 4) - 16*pow(solution[1], 2) + 5*solution[1]);
}


/**
 * @brief Function to return random number between the upper and lower bounds
 * 
 * @param bounds array of size the first one is the lower bound and the second one is the upper bound
 * @return double random value between the upper and lower bounds
 */
double randBounds(double bounds[2]){
    // Seed the random number generator
    srand(time(NULL)); 
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
    // Seed the random number generator
    srand(time(NULL)); 
    for (int i = 0; i < DIMENSIONS; i++) {
        // randomly assign a value for the position of particle between the lower and upper bounds
        particle->position[i] = randBounds(bounds[i]);
        // particle is not moving to velocity is 0
        particle->velocity[i] = ( (double)rand() / (double)RAND_MAX) - 0.5;
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
    // Seed the random number generator
    srand(time(NULL)); 
    for (int i = 0; i < DIMENSIONS; i++) {
        double r = (double)rand() / (double)RAND_MAX;
        double adjustment = 1;
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
void particleSwarm(double bounds[][2]) {
    
    Particle swarm[SWARM_SIZE];
    double global_best_value = INFINITY;
    double global_best_position[DIMENSIONS];

    // Create the particles
    for (int i = 0; i < SWARM_SIZE; i++) {
        createParticle(&swarm[i], bounds);
        swarm[i].best_value = objectiveFunction(swarm[i].position);
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
            double current_value = objectiveFunction(swarm[i].position);

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