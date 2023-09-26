#include "FacetsReader.h"

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

    while(strlen(line = readLine(&fTheta)) > 0){
        if(nrows == 0){
            // count cols
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

    for(col = 0; col < ncols; col++){
        line = readLine(&fTheta);
        tok = strtok(line, " ");
        for(row = 0; row < nrows; row++){
            facets.thetaMap[row][col] = atof(tok);
            tok = strtok(NULL, " ");
        }
        free(line);

        line = readLine(&fPhi);
        tok = strtok(line, " ");
        for(row = 0; row < nrows; row++){
            facets.phiMap[row][col] = atof(tok);
            tok = strtok(NULL, " ");
        }
        free(line);
    }

    fclose(fTheta);
    fclose(fPhi);
    return facets;
}
