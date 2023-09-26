#include "AFMDataReader.h"
#define MAXPATH 260

AFMData afmData_readFromFile(char *fileName){
    int row, col;
    char *line, *tok;
    char *binnedName;

    printf("File Name: %s\n", fileName);
    fflush(stdout);

    int xCapacity = 128, yCapacity = 128;
    double **tmpp, *tmp;
    AFMData afmData;
    afmData.xResolution = 0;
    afmData.yResolution = 0;

    afmData.zValues = (double**)malloc(sizeof(double*) * xCapacity);

    // Read file again for values
    FILE *file = fopen(fileName, "r");

    printf("Reading AFM data...\n");
    fflush(stdout);

    while(strlen(line = readLine(&file)) > 0){
        if(afmData.xResolution == xCapacity){
            xCapacity *= 2;
            tmpp = (double **)realloc(afmData.zValues, sizeof(double*) * xCapacity);
            if(tmpp == NULL){
                printf("Failed to reallocate space for line.");
                exit(1);
            }
            afmData.zValues = tmpp;
        }
        afmData.zValues[afmData.xResolution] = (double*)calloc(yCapacity, sizeof(double));

        if(afmData.xResolution == 0){
            tok = strtok(line, "\t");
            while(tok != NULL){
                if(afmData.yResolution == yCapacity){
                    yCapacity *= 2;
                    tmp = (double*)realloc(afmData.zValues[afmData.xResolution], sizeof(double) * yCapacity);
                    if(tmp == NULL){
                        printf("Failed to reallocate space for line.");
                        exit(1);
                    }
                    afmData.zValues[afmData.xResolution] = tmp;
                }
                afmData.zValues[afmData.xResolution][afmData.yResolution] = atof(tok);
                tok = strtok(NULL, "\t");
                afmData.yResolution++;
            }
        }
        else{
            col = 0;
            tok = strtok(line, "\t");
            while(tok != NULL){
                afmData.zValues[afmData.xResolution][col] = atof(tok);
                tok = strtok(NULL, "\t");
                col++;
            }
        }
        free(line);
        afmData.xResolution++;
    }

    fclose(file);
    return afmData;
}

AFMData afmData_readBinnedFromFile(char *fileName){
    int row, col, xCapacity, yCapacity;
    char line[100];
    AFMData afmData;

    printf("File Name: %s\n", fileName);
    fflush(stdout);

    FILE *file;
    
    file = fopen(fileName, "r");
    fgets(line, 100, file);
    xCapacity = atoi(line);
    printf("xCapacity: %d\n", xCapacity);
    fgets(line, 100, file);
    yCapacity = atoi(line);
    printf("yCapacity: %d\n", yCapacity);

    afmData = afmData_new(xCapacity, yCapacity);

    for(row = 0; row < xCapacity; row++){
        for(col = 0; col < yCapacity; col++){
            fgets(line, 100, file);
            afmData.zValues[row][col] = atof(line);
        }
    }
    fclose(file);
    return afmData;
}
