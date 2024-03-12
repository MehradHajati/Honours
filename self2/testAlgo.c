#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include "testAlgo.h"

#define SWARM_SIZE 10 // Number of particles in the swarm  which is ten times the number of dimensions
#define MAX_ITERATIONS_SA 1000 // Maximum number of iterations for simulated annealing
#define MAX_ITERATIONS_PS 100 // Maximum number of iterations for particle swarm
#define PHI_P 0.9 // ratio for personal particle component
#define PHI_S 0.9 // ratio for Social swarm component
#define PHI_C 0.9 // ratio for current velocity
#define stepCoeff 3 // Coefficient for the step size


void main(int argc, char *argv[]) { 

    // Seed the random number generator
    srand(time(NULL));
    
    // creating the upper and lower bounds
    //double bounds[DIMENSIONS][2] = {{-200, 200}, {-0.25, 0.25}, {-0.3, 0.3}, {-2e-4, 2e-4}, {-2e-4, 2e-4}, {-2e-4, 2e-4}, {-40, 40}, {-0.3, 0.3}, {-0.35, 0.35}, {-2e-4, 2e-4}, {-2e-4, 2e-4}, {-2e-4, 2e-4}};
    //double bounds[DIMENSIONS][2] = {{-5, 5}, {-5, 5}, {-5, 5}, {-5, 5}, {-5, 5}, {-5, 5}, {-5, 5}, {-5, 5}, {-5, 5}, {-5, 5}, {-5, 5}, {-5, 5}};
    double bounds[DIMENSIONS][2] = {{-5, 5}, {-5, 5}};
    //double upperBounds[DIMENSIONS] = {640, 3.5, 0.003, 5e-6, 5e-6, 5e-6, 650, 0.003, 2.8, 5e-6, 5e-6, 5e-6};
    //double lowerBounds[DIMENSIONS] = {600, 3, 0, 0, 0, 0, 630, 0, 2, 0, 0, 0};

    // variable to keep track of the time
    clock_t start, end;
    double sum_time_SA = 0;
    double sum_time_PS = 0;
    double sum_SA = 0;
    double sum_PS = 0;


    /*for(int i = 0; i < 1; i++){

    
        start = clock();
        sum_SA += simulatedAnnealing(0.1, bounds);
        end = clock();

        sum_time_SA += ((double) (end - start)) / CLOCKS_PER_SEC;
    }

    printf("For SA the total chi is: %f and time used is: %f\n", sum_SA/1, sum_time_SA/1);*/

    for(int i = 0; i < 20; i++){
        start = clock();
        sum_PS += particleSwarm(bounds);
        end = clock();

        sum_time_PS += ((double) (end - start)) / CLOCKS_PER_SEC;
    }

    
    printf("For PS the total chi is: %f and time used is: %f\n", sum_PS/1, sum_time_PS/1);

    //double eg[DIMENSIONS] = {0, 0};
    //double eg[DIMENSIONS] = {-2.903, -2.903, -2.903, -2.903, -2.903, -2.903, -2.903, -2.903, -2.903, -2.903, -2.903, -2.903};
    /*double eg[DIMENSIONS] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double meanSums[DIMENSIONS];
    double StdSums[DIMENSIONS];
    double output[DIMENSIONS*2];

    for(int i = 0; i < 10000; i++){
        MetropolisHasting(bounds, eg, output);
        for(int j = 0; j < DIMENSIONS; j++){
            meanSums[j] += output[j];
            StdSums[j] += output[j+DIMENSIONS];
        }

    }
    
    for(int i = 0; i < DIMENSIONS; i++){
        printf("Index %d with Mean %f and std is: %f\n", i, meanSums[i]/10000, StdSums[i]/10000);
    }*/
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
        double r = ((double)rand() / (double)RAND_MAX) - 0.5;
        double stepSize = (bounds[i][1] - bounds[i][0]) / stepCoeff;
        new_solution[i] = current_solution[i] + stepSize * r;
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
            //printf("This index %d has hit the lower bound with value %f\n",i, new_solution[i]);
            // bring it back within bounds
            new_solution[i] = bounds[i][0];
        }
        else if(new_solution[i] > bounds[i][1]){
            //printf("This index %d has hit the upper bound with value %f\n", i, new_solution[i]);
            // bring it back within bounds
            new_solution[i] = bounds[i][1];
        }
    }
}

// Simulated annealing algorithm
double simulatedAnnealing(double cooling_rate, double bounds[][2]){

    double temp = 1;
    // creating a random solution
    double current_solution[DIMENSIONS] = {randBounds(bounds[0]), randBounds(bounds[1])};//, randBounds(bounds[2]), randBounds(bounds[3]), randBounds(bounds[4]), randBounds(bounds[5]), randBounds(bounds[6]), randBounds(bounds[7]), randBounds(bounds[8]), randBounds(bounds[9]), randBounds(bounds[10]), randBounds(bounds[11])};
    double new_solution[DIMENSIONS];
    int iter;

    double current_energy = objectiveFunction(current_solution);

    double best = current_energy;
    double best_solution[DIMENSIONS];

    for(int i = 0; i < DIMENSIONS; i++){
        best_solution[i] = current_solution[i];
    }

    //FILE* file = fopen("SA_one_best.txt", "w");

    for(iter = 0; iter < MAX_ITERATIONS_SA; iter++) {

        // Cool down
        temp = cooling_rate / log(iter + 2);

        /*if(iter % 100 == 0){
            printf("Current Temp is: %f and number of iteration is: %d with current chi-sqr %f and best chi-sqr: %f\n", temp, iter, current_energy, best);
        }*/
        
        getNeighbor(current_solution, new_solution, bounds);
        checkBounds(new_solution, bounds);
        double new_energy = objectiveFunction(new_solution);

        // rplacing the best energy 
        if(best > new_energy) {
            best = new_energy;
            for (int i = 0; i < DIMENSIONS; i++) {
                best_solution[i] = new_solution[i];
                //fprintf(file, "%f ", best_solution[i]);
            }
            //fprintf(file, "\n");
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

    /*FILE* file = fopen("SA_output_best.txt", "a");

    for (int j = 0; j < DIMENSIONS; j++) {
        fprintf(file, "%lf ", best_solution[j]);
    }
    fprintf(file, "\n");

    fclose(file);*/

    /*printf("Best solution:\n");
    for (int i = 0; i < DIMENSIONS; i++){
        printf("x[%d] = %f\n", i, best_solution[i]);
    }
    printf("with chi squared = %f\n", best);*/
    return best;
}

double objectiveFunction(double *solution){
    // Squaring Function
    //return (pow(solution[0], 2)) + pow(solution[1], 2) + pow(solution[2], 2) + pow(solution[3], 2) + pow(solution[4], 2) + pow(solution[5], 2) + pow(solution[6], 2) + pow(solution[7], 2) + pow(solution[8], 2) + pow(solution[9], 2) + pow(solution[10], 2) + pow(solution[11], 2);
    
    /*Styblinski Tang Function
    double sum = 0.0;
    for (int i = 0; i < DIMENSIONS; i++) {
        sum += solution[i] * solution[i] * solution[i] * solution[i] - 16 * solution[i] * solution[i] + 5 * solution[i];
    }
    return sum / 2.0;*/

    //Three-Hump Camel Function
    double x = solution[0];
    double y = solution[1];
    return 2*x*x - 1.05*x*x*x*x + (x*x*x*x*x*x)/6 + x*y + y*y;
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
        // the particle is not moving at the begginning so velocity is 0
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
void moveParticle(Particle *particle, double *global_best_position, double bounds[][2]) {
    for (int i = 0; i < DIMENSIONS; i++) {
        double r = (double)rand() / (double)RAND_MAX;
        // Calculating the personal velocity based on the best position for the particle
        double personal_velocity = PHI_P * r * (particle->best_position[i] - particle->position[i]);
        // Calculating the social velocity based on the global best position
        double social_velocity = PHI_S * r * (global_best_position[i] - particle->position[i]);
        // changing the veloctity of the particle
        particle->velocity[i] = (PHI_C * particle->velocity[i]) + personal_velocity + social_velocity;
        // moving the particle
        particle->position[i] = particle->position[i] + particle->velocity[i];
    }
}

// Particle Swarm Optimization algorithm
double particleSwarm(double bounds[][2]) {
    
    Particle swarm[SWARM_SIZE];
    double global_best_value = INFINITY;
    double global_best_position[DIMENSIONS];
    bool flag = true;

    /*FILE* file1 = fopen("PS_one_particle1.txt", "w");
    FILE* file2 = fopen("PS_one_particle2.txt", "w");
    FILE* file3 = fopen("PS_one_particle3.txt", "w");
    FILE* file4 = fopen("PS_one_particle4.txt", "w");
    FILE* file5 = fopen("PS_one_particle5.txt", "w");

    FILE* file6 = fopen("PS_one_best_particle1.txt", "w");
    FILE* file7 = fopen("PS_one_best_particle2.txt", "w");
    FILE* file8 = fopen("PS_one_best_particle3.txt", "w");
    FILE* file9 = fopen("PS_one_best_particle4.txt", "w");
    FILE* file10 = fopen("PS_one_best_particle5.txt", "w");

    FILE* file11 = fopen("PS_group_best.txt", "w");

    FILE* file12 = fopen("PS_output.txt", "a");*/

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
                    //fprintf(file11, "%lf ", global_best_position[j]);
                }
                //fprintf(file11, "\n");
            }
            
            // checking when the global best position is within the 0.25 square
            if(abs(global_best_position[0]) < 0.25 && abs(global_best_position[1]) < 0.25 && flag) {
                printf("iteration number is: %d\n", iter);
                flag = false;
            }


            /*if(i == 1){
                for (int j = 0; j < DIMENSIONS; j++) {
                    fprintf(file1, "%lf ", swarm[i].position[j]);
                    fprintf(file6, "%lf ", swarm[i].best_position[j]);
                }
                fprintf(file1, "\n");
                fprintf(file6, "\n");
            }
            if(i == 2){
                for (int j = 0; j < DIMENSIONS; j++) {
                    fprintf(file2, "%lf ", swarm[i].position[j]);
                    fprintf(file7, "%lf ", swarm[i].best_position[j]);
                }
                fprintf(file2, "\n");
                fprintf(file7, "\n");
            }
            if(i == 3){
                for (int j = 0; j < DIMENSIONS; j++) {
                    fprintf(file3, "%lf ", swarm[i].position[j]);
                    fprintf(file8, "%lf ", swarm[i].best_position[j]);
                }
                fprintf(file3, "\n");
                fprintf(file8, "\n");
            }
            if(i == 4){
                for (int j = 0; j < DIMENSIONS; j++) {
                    fprintf(file4, "%lf ", swarm[i].position[j]);
                    fprintf(file9, "%lf ", swarm[i].best_position[j]);
                }
                fprintf(file4, "\n");
                fprintf(file9, "\n");
            }
            if(i == 5){
                for (int j = 0; j < DIMENSIONS; j++) {
                    fprintf(file5, "%lf ", swarm[i].position[j]);
                    fprintf(file10, "%lf ", swarm[i].best_position[j]);
                }
                fprintf(file5, "\n");
                fprintf(file10, "\n");
            }*/
            
        }
        //printf("Iteration %d with current chi squared = %f\n", iter, global_best_value);
    }

    /*for(int i = 0; i < DIMENSIONS; i++){
        fprintf(file12, "%lf ", global_best_position[i]);
    }
    fprintf(file12, "\n");*/

    /*printf("Best solution:\n");
    for (int i = 0; i < DIMENSIONS; i++){
        printf("x[%d] = %f\n", i, global_best_position[i]);
    }
    printf("with chi squared = %f\n", global_best_value);*/
    /*fclose(file1);
    fclose(file2);
    fclose(file3);
    fclose(file4);
    fclose(file5);
    fclose(file6);
    fclose(file7);
    fclose(file8);
    fclose(file9);
    fclose(file10);
    fclose(file11);
    fclose(file12);*/
    
    return global_best_value;
}


// Metropolis-Hasting algorithm
void MetropolisHasting(double bounds[][2], double* solution, double* output){

    double temp = 1;
    double acceptCounter = 0;

    // setting the given solution as the current solution
    double* current_solution = solution;
    double new_solution[DIMENSIONS];
    int iter;

    // getting the current objective function value
    double current_energy = objectiveFunction(current_solution);

    // creating the arrays that will hold all the solutions value parameters and setting its first entry
    double all_solutions[MAX_ITERATIONS_MTH+1][DIMENSIONS];
    for(int i = 0; i < DIMENSIONS; i++){
        all_solutions[0][i] = current_solution[i];
    }

    for(iter = 0; iter < MAX_ITERATIONS_MTH; iter++) {

        /*if(iter % 100 == 0){
            printf("number of iteration is: %d with current chi-sqr %f\n", iter, current_energy);
        }*/
        
        getNeighbor(current_solution, new_solution, bounds);
        checkBounds(new_solution, bounds);
        double new_energy = objectiveFunction(new_solution);

        // Decide if we should accept the new solution
        if ( new_energy < current_energy) {
            // New solution is better, accept it
            // incrementing counter
            acceptCounter++;
            for (int i = 0; i < DIMENSIONS; i++) {
                current_solution[i] = new_solution[i];
            }
            current_energy = new_energy;
            

        } 
        else{ // New solution is worse, accept it with a probability decreasing with temp
            double r = ((double)rand() / (double)RAND_MAX);
            double e = exp( (current_energy-new_energy) / temp );
            if(e > r) {
                // incrementing counter
                acceptCounter++;
                for (int i = 0; i < DIMENSIONS; i++) {
                    current_solution[i] = new_solution[i];
                }
                current_energy = new_energy;
            } 
        }

        // placing the current solution in the all_solutions array
        for(int i = 0; i < DIMENSIONS; i++) {
            all_solutions[iter+1][i] = current_solution[i];
        }
    }    

    //printf("the acceptance rate is: %f\n", 100 * acceptCounter / MAX_ITERATIONS_MTH);

    double xbar[DIMENSIONS] = {0, 0};//, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double xsig[DIMENSIONS] = {0, 0};//, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    // calculating the mean and stds here
    for(int i = 0; i < (MAX_ITERATIONS_MTH+1); i++){
        double d1[DIMENSIONS];
        double d2[DIMENSIONS];
        double n2 = i + 2;

        // d1 = x - xbar
        for(int j = 0; j < DIMENSIONS; j++){
            d1[j] = all_solutions[i][j] - xbar[j];
        }

        // xbar = xbar + (d1/n2)
        for(int j = 0; j < DIMENSIONS; j++){
            xbar[j] = xbar[j] + (d1[j] / n2);
        }

        // d2 = x - xbar
        for(int j = 0; j < DIMENSIONS; j++){
            d2[j] = all_solutions[i][j] - xbar[j];
        }

        // xsig = xsig + (d1*d2)
        for(int j = 0; j < DIMENSIONS; j++){
            xsig[j] = xsig[j] + (d1[j] * d2[j]);
        }

    }

    // adding the mean and stds array to the output array
    for (int i = 0; i < DIMENSIONS; i++){
        output[i] = xbar[i];
        output[i+DIMENSIONS] = sqrt(xsig[i]/MAX_ITERATIONS_MTH);
    }
}