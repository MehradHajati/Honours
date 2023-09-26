#include "FacetsWriter.h"

void facets_write(Facets facets, char *thetaPath, char *phiPath){
    int row, col;
    FILE *thetaFile, *phiFile;
    thetaFile = fopen(thetaPath, "w");
    phiFile = fopen(phiPath, "w");

    printf("Writing Facets...\n");
    fflush(stdout);
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
