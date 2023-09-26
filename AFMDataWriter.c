#include "AFMDataWriter.h"

void afmData_write(AFMData *afmData, char *path, int rowStart, int rowEnd){
    printf("AFM Path: %s\n", path);
    int row, col;
    FILE *out;
    out = fopen(path, "w");

    for(row = rowStart; row < rowEnd; row++){
        for(col = 0; col < afmData->yResolution; col++){
            fprintf(out, "%.5lf\t", afmData->zValues[row][col]);
        }
        fprintf(out, "\n");
    }

    fclose(out);
}