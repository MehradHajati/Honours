#include "FacetsReader.h"

/**
 * @brief Creates a facets struct form the given theta and phi map files.
 * 
 * @param thetaFn The path to the theta map. 
 * @param phiFn The path to the phi map.
 * @return A facets struct. 
 */
Facets facets_readFromFiles(char *thetaFn, char *phiFn){
    int row, col, nrows, ncols;
    char *line, *tok;
    Facets facets;
    FILE *fTheta = fopen(thetaFn, "r");
    FILE *fPhi   = fopen(phiFn, "r");

    printf("Determining Resolution...\n");
    fflush(stdout);

    nrows = 0;
    ncols = 0;
    // while the length of the read in line is greater than 0, the loop will continue
    while(strlen(line = readLine(&fTheta)) > 0){
        if(nrows == 0){
            // counting the number of columns and rows in the files
            tok = strtok(line, " ");
            while(tok != NULL){
                ncols++;
                tok = strtok(NULL, " ");
            }
        }
        free(line);
        nrows++;
    }

    fclose(fTheta);

    printf("Number of rows in %s: %d\n", thetaFn, nrows);
    printf("Number of cols in %s: %d\n", thetaFn, ncols);
    fflush(stdout);
    facets = facets_new(nrows, ncols);

    printf("Reading Facets data...\n");
    fflush(stdout);

    fTheta = fopen(thetaFn, "r");
    // reading the data here
    for(col = 0; col < ncols; col++){
        line = readLine(&fTheta);
        // using strtok and space to get each of the individual numbers
        tok = strtok(line, " ");
        for(row = 0; row < nrows; row++){
            // using atof to convert the string (tok) into an integer
            facets.thetaMap[row][col] = atof(tok);
            tok = strtok(NULL, " ");
        }
        free(line);

        line = readLine(&fPhi);
        // using strtok and space to get each of the individual numbers
        tok = strtok(line, " ");
        for(row = 0; row < nrows; row++){
            // using atof to convert the string (tok) into an integer
            facets.phiMap[row][col] = atof(tok);
            tok = strtok(NULL, " ");
        }
        free(line);
    }

    fclose(fTheta);
    fclose(fPhi);
    return facets;
}
