#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "testAlgo.h"

#define SWARM_SIZE 100 // Number of particles in the swarm
#define MAX_ITERATIONS_SA 1000 // Maximum number of iterations for simulated annealing
#define MAX_ITERATIONS_PS 100 // Maximum number of iterations for particle swarm
#define PHI_P 0.5 // ratio for personal particle component
#define PHI_S 0.5 // ratio for Social swarm component
#define PHI_C 0.5 // ratio for current velocity
#define stepCoeff 10 // Coefficient for the step size


void main(int argc, char *argv[]) { 

    // Seed the random number generator
    srand(time(NULL));
    
    // creating the upper and lower bounds
    double bounds[DIMENSIONS][2] = {{-200, 200}, {-0.25, 0.25}, {-0.3, 0.3}, {-2e-4, 2e-4}, {-2e-4, 2e-4}, {-2e-4, 2e-4}, {-40, 40}, {-0.3, 0.3}, {-0.35, 0.35}, {-2e-4, 2e-4}, {-2e-4, 2e-4}, {-2e-4, 2e-4}};
    //double upperBounds[DIMENSIONS] = {640, 3.5, 0.003, 5e-6, 5e-6, 5e-6, 650, 0.003, 2.8, 5e-6, 5e-6, 5e-6};
    //double lowerBounds[DIMENSIONS] = {600, 3, 0, 0, 0, 0, 630, 0, 2, 0, 0, 0};

    // variable to keep track of the time
    clock_t start, end;
    double cpu_time_used;

    start = clock();
    simulatedAnnealing(0.1, bounds);
    end = clock();

    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("SA CPU time: %f\n", cpu_time_used);


    /*start = clock();
    particleSwarm(bounds);
    end = clock();

    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("PSO CPU time: %f\n", cpu_time_used);*/
}

/**
 * @brief This function is used to randomly generate a neighbor of the current solution
 * 
 * @param current_solution The current solution
 * @param new_solution The newly randomly generated solution which will be close to the current solution
 * @param bounds The bounds of the of the search space
 */
void getNeighbor(double *current_solution, double *new_solution, double bounds[][2]) {
    for (int i = 0; i < DIMENSIONS; i++) {
        double r = (double)rand() / (double)RAND_MAX;
        double stepSize = (bounds[i][1] - bounds[i][0]) / stepCoeff;
        new_solution[i] = current_solution[i] + stepSize * (r-0.5);
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
            // bring it back within bounds and take one step towards the upper bound
            new_solution[i] = bounds[i][0] + ((bounds[i][1] - bounds[i][0]) / stepCoeff);
        }
        else if(new_solution[i] > bounds[i][1]){
            printf("This index %d has hit the upper bound with value %f\n", i, new_solution[i]);
            // bring it back within bounds and take one step towards the lower bound
            new_solution[i] = bounds[i][1] - ((bounds[i][1] - bounds[i][0]) / stepCoeff);
        }
    }
}

// Simulated annealing algorithm
void simulatedAnnealing(double cooling_rate, double bounds[][2]){

    double temp = 1;
    // creating a random solution
    double current_solution[DIMENSIONS] = {randBounds(bounds[0]), randBounds(bounds[1]), randBounds(bounds[2]), randBounds(bounds[3]), randBounds(bounds[4]), randBounds(bounds[5]), randBounds(bounds[6]), randBounds(bounds[7]), randBounds(bounds[8]), randBounds(bounds[9]), randBounds(bounds[10]), randBounds(bounds[11])};
    double new_solution[DIMENSIONS];
    int iter;

    double current_energy = objectiveFunction(current_solution);

    double best = current_energy;
    double* best_solution = current_solution;

    for(iter = 0; iter < MAX_ITERATIONS_SA; iter++) {

        // Cool down
        temp = cooling_rate / log(iter + 2);

        if(iter % 100 == 0){
            printf("Current Temp is: %f and number of iteration is: %d with current chi-sqr %f and best chi-sqr: %f\n", temp, iter, current_energy, best);
        }
        
        getNeighbor(current_solution, new_solution, bounds);
        checkBounds(new_solution, bounds);
        double new_energy = objectiveFunction(new_solution);

        // rplacing the best energy 
        if(best > new_energy) {
            best = new_energy;
            best_solution = new_solution;
        }

        // Decide if we should accept the new solution
        if ( new_energy < current_energy) {
            // New solution is better, accept it
            for (int i = 0; i < DIMENSIONS; i++) {
                current_solution[i] = new_solution[i];
            }
            current_energy = new_energy;
        } else{ // New solution is worse, accept it with a probability decreasing with temp
            double r = ((double)rand() / (double)RAND_MAX);
            double e = exp( (current_energy-new_energy) / temp );
            if(e > r) {
                for (int i = 0; i < DIMENSIONS; i++) {
                    current_solution[i] = new_solution[i];
                }
                current_energy = new_energy;
            } 
        }
    }

    printf("Best solution:\n");
    for (int i = 0; i < DIMENSIONS; i++){
        printf("x[%d] = %f\n", i, best_solution[i]);
    }
    printf("with chi squared = %f\n", best);
}

double objectiveFunction(double *solution){
    //return 100 * pow( (solution[1] - pow(solution[0], 2) ), 2 ) + pow( (1-solution[0]), 2 );
    //return 0.5 * (pow(solution[0], 4) - 16*pow(solution[0], 2) + 5*solution[0] + pow(solution[1], 4) - 16*pow(solution[1], 2) + 5*solution[1]);
    return pow(solution[0], 2) + pow(solution[1], 2) + pow(solution[2], 2) + pow(solution[3], 2) + pow(solution[4], 2) + pow(solution[5], 2) + pow(solution[6], 2) + pow(solution[7], 2) + pow(solution[8], 2) + pow(solution[9], 2) + pow(solution[10], 2) + pow(solution[11], 2);

}


/**
 * @brief Function to return random number between the upper and lower bounds
 * 
 * @param bounds array of size the first one is the lower bound and the second one is the upper bound
 * @return double random value between the upper and lower bounds
 */
double randBounds(double bounds[2]){
    double range = bounds[1] - bounds[0];
    double div = RAND_MAX / range;
    return bounds[0] + ( (double)rand() / div);
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
        // randomly aassign velocity for the particle
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
void moveParticle(Particle *particle, double *global_best_position, double bounds[][2]) {
    for (int i = 0; i < DIMENSIONS; i++) {
        double r = (double)rand() / (double)RAND_MAX;
        double stepSize = (bounds[i][1] - bounds[i][0]) / stepCoeff;
        // Calculating the personal velocity based on the best position for the particle
        double personal_velocity = PHI_P * stepSize * r * (particle->best_position[i] - particle->position[i]);
        // Calculating the social velocity based on the global best position
        double social_velocity = PHI_S * stepSize * r * (global_best_position[i] - particle->position[i]);
        // changing the veloctity of the particle
        particle->velocity[i] = PHI_C * particle->velocity[i] + personal_velocity + social_velocity;
        // moving the particle
        particle->position[i] = particle->position[i] + particle->velocity[i];
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
    for (int iter = 0; iter < MAX_ITERATIONS_PS; iter++) {
        for (int i = 0; i < SWARM_SIZE; i++) {

            // Moving the current particle
            moveParticle(&swarm[i], global_best_position, bounds);

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
        printf("Iteration %d with current chi squared = %f\n", iter, global_best_value);
    }

    printf("Best solution:\n");
    for (int i = 0; i < DIMENSIONS; i++){
        printf("x[%d] = %f\n", i, global_best_position[i]);
    }
    printf("with chi squared = %f\n", global_best_value);

}