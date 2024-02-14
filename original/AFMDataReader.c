#include "AFMDataReader.h"
#define MAXPATH 260


/**
 * @brief Creates an AFMData struct from the data in the given file that has not been binned. Meaning from the Atotech file
 * 
 * @param fileName The name of the file to read.
 * @return AFMData An AFMData struct.
 */
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

    // allocating the memory for the outer matrix of the 2d matrix
    afmData.zValues = (double**)malloc(sizeof(double*) * xCapacity);

    // Read file again for values
    FILE *file = fopen(fileName, "r");

    printf("Reading AFM data...\n");
    fflush(stdout);

    // while the length of the next line is not zero continue reading
    while(strlen(line = readLine(&file)) > 0){
        // working with the lines and allocating space for them here
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
            // seperating the line by tab and taking each of the numbers one by one
            tok = strtok(line, "\t");
            // while the tok is not null continue reading and allocating space for the numbers, if tok is that means we at the end of the line
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

/**
 * @brief Creates an AFMData struct from the data in the given file that has been binned.
 * 
 * @param fileName The name of the file to read.
 * @return AFMData An AFMData struct.
 */
AFMData afmData_readBinnedFromFile(char *fileName){
    int row, col, xCapacity, yCapacity;
    char line[100];
    AFMData afmData;

    printf("File Name: %s\n", fileName);
    fflush(stdout);

    FILE *file;
    // opening the file to read it
    file = fopen(fileName, "r");
    fgets(line, 100, file);
    // getting the x and y capacity
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
