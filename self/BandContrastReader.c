#include "BandContrastReader.h"
#define MAXPATH 260

/**
 * @brief Creates a band contrast struct from the data in the given file.
 * 
 * @param path The name of the file to read.
 * @return A band contrast struct. 
 */
BandContrast bandContrast_readFromFile(char *path){
    int row, col, nrows, ncols;
    char *line, *tok;
    FILE *file;
    BandContrast bc;
    char *fn;


    nrows = 0;
    ncols = 0;

    fn = (char *)malloc(sizeof(char) * MAXPATH);
    strcpy(fn, path);
    strcat(fn, "BandContrast.txt");
    file = fopen(fn, "r");
    printf("Determining Resolution for file %s\n", fn);
    fflush(stdout);
    //counting the number of columns and rows in this while loop
    while(strlen(line = readLine(&file)) > 0){
        if(nrows == 0){
            // count cols
            tok = strtok(line, " ");
            while(tok != NULL){
                ncols++;
                tok = strtok(NULL, " ");
            }
        }
        nrows++;
        free(line);
    }
    fclose(file);
    printf("Number of rows in %s: %d\n", fn, nrows);
    printf("Number of cols in %s: %d\n", fn, ncols);
    fflush(stdout);

    bc = bandContrast_new(nrows, ncols);

    printf("Reading Band Contrast data...\n");
    file = fopen(fn, "r");
    free(fn);
    // for loop to read each line of the file for the band contrast data
    for(row = 0; row < nrows; row++){
        line = readLine(&file);
        tok = strtok(line, " ");
        for(col = 0; col < ncols; col++){
            bc.greyScale[row][col] = atof(tok) / 255.0;
            tok = strtok(NULL, " ");
        }
        free(line);
    }
    fclose(file);

    printf("Reading Euler data...\n");

    // reading the red part of the euler data with the for loop to go over every line in the file and the strtok will seperate them based on space
    fn = (char *)malloc(sizeof(char) * MAXPATH);
    strcpy(fn, path);
    strcat(fn, "EBSDred.txt");
    printf("File is %s\n",fn);
    file = fopen(fn, "r");
    free(fn);
    for(row = 0; row < nrows; row++){
        line = readLine(&file);
        tok = strtok(line, " ");
        for(col = 0; col < ncols; col++){
            bc.EBSDred[row][col] = atoi(tok);
            tok = strtok(NULL, " ");
        }
        free(line);
    }
    fclose(file);

    // reading the green part of the euler data with the for loop to go over every line in the file and the strtok will seperate them based on space
    fn = (char *)malloc(sizeof(char) * MAXPATH);
    strcpy(fn, path);
    strcat(fn, "EBSDgreen.txt");
    printf("File is %s\n",fn);
    file = fopen(fn, "r");
    free(fn);
    for(row = 0; row < nrows; row++){
        line = readLine(&file);
        tok = strtok(line, " ");
        for(col = 0; col < ncols; col++){
            bc.EBSDgreen[row][col] = atoi(tok);
            tok = strtok(NULL, " ");
        }
        free(line);
    }
    fclose(file);

    // reading the blue part of the euler data with the for loop to go over every line in the file and the strtok will seperate them based on space
    fn = (char *)malloc(sizeof(char) * MAXPATH);
    strcpy(fn, path);
    strcat(fn, "EBSDblue.txt");
    printf("File is %s\n",fn);
    file = fopen(fn, "r");
    free(fn);
    for(row = 0; row < nrows; row++){
        line = readLine(&file);
        tok = strtok(line, " ");
        for(col = 0; col < ncols; col++){
            bc.EBSDblue[row][col] = atoi(tok);
            tok = strtok(NULL, " ");
        }
        free(line);
    }
    fclose(file);

    return(bc);
}
