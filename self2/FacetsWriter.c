#include "FacetsWriter.h"

/**
 * @brief Writes the theta and phi maps from the given facets struct.
 * 
 * @param facets The facets from which to write theta and phi maps.
 * @param thetaPath The path to the theta map file.
 * @param phiPath The path to the phi map file.
 */
void facets_write(Facets facets, char *thetaPath, char *phiPath){
    int row, col;
    FILE *thetaFile, *phiFile;
    // opening the files to write to them
    thetaFile = fopen(thetaPath, "w");
    phiFile = fopen(phiPath, "w");

    printf("Writing Facets...\n");
    fflush(stdout);
    // for loop to write the data of the matrices to files
    for(col = 0; col < facets.ncol; col++){
        for(row = 0; row < facets.nrow; row++){
            fprintf(thetaFile, "%.5lf ", facets.thetaMap[row][col]);
            fprintf(phiFile, "%.5lf ", facets.phiMap[row][col]);
        }
        fprintf(thetaFile, "\n");
        fprintf(phiFile, "\n");
    }

    fclose(thetaFile);
    fclose(phiFile);
}
