#include "AFM.h"


/**
 * @brief Allocates memory for AFMData with the given resolution.
 * 
 * @param xResolution The number of data points in the x-direction the probe will store.
 * @param yResolution The number of data points in the y-direction the probe will store.
 * @return AFMData with the given resolution.
 */
AFMData afmData_new(int xResolution, int yResolution){
    int i;
    AFMData afmData;
    afmData.xResolution = xResolution;
    afmData.yResolution = yResolution;
    // allocating the needed space for the outer matrix
    afmData.zValues = (double**)malloc(sizeof(double*) * xResolution);
    // for loop to allocated the needed space for the inner matrices at each index, creating a 2d matrix
    for(i = 0; i < xResolution; i++){
        afmData.zValues[i] = (double*)calloc(yResolution, sizeof(double));
    }

    return afmData;
}


/**
 * @brief Frees memory allocated within the AFMData.
 * 
 * @param afmData The AFMData to free the memory of.
 */
void afmData_free(AFMData* afmData){
    int i;
    for(i = 0; i < afmData->xResolution; i++){
        free(afmData->zValues[i]);
    }
    free(afmData->zValues);
}


/**
 * @brief Bins the given AFM data by a factor of 2.
 * 
 * @param afm The AFMData to bin.
 */
void afmData_binBy2(AFMData *afm){
    double **zVals;
    int newXRes, newYRes, oldRow, oldCol, newRow, newCol, dir, numDirs = 4;
    int dirs[4][2] = {
        {0,0},
        {1,0},
        {0,1},
        {1,1}
    };
    // setting the new resolutions
    newXRes = afm->xResolution / 2;
    newYRes = afm->yResolution / 2;
    // allocating new memory for the data
    zVals = (double **)malloc(sizeof(double*) * newXRes);
    // filling in indices of zVals with matrices to make it a 2d matrix
    for(oldRow = 0; oldRow < newXRes; oldRow++){
        zVals[oldRow] = (double *)calloc(newYRes, sizeof(double));
    }
    // binning the data by a factor of 2
    for(newRow = 0; newRow < newXRes; newRow++){
        for(newCol = 0; newCol < newYRes; newCol++){
            zVals[newRow][newCol] = 0;
            oldRow = newRow * 2;
            oldCol = newCol * 2;
            for(dir = 0; dir < numDirs; dir++){
                zVals[newRow][newCol] += afm->zValues[oldRow + dirs[dir][0]][oldCol + dirs[dir][1]];
            }
            zVals[newRow][newCol] /= 8.0; // 4x2 for averaging and new pixel size.
        }
    }
    // Freeing the memory after we are done with it
    for(oldRow = 0; oldRow < afm->xResolution; oldRow++){
        free(afm->zValues[oldRow]);
    }
    free(afm->zValues);
    afm->zValues = NULL;
    afm->zValues = zVals;
    afm->xResolution = newXRes;
    afm->yResolution = newYRes;
}


/**
 * @brief Flips the given AFM data across the middle row.
 * 
 * @param afm The AFM data to flip.
 */
void afmData_flipY(AFMData *afm){
    double *tmpRow;
    int row;
    // for loop to flip the data
    for(row = 0; row < afm->xResolution / 2; row++){
        tmpRow = afm->zValues[row];
        afm->zValues[row] = afm->zValues[afm->xResolution - 1 - row];
        afm->zValues[afm->xResolution - 1 - row] = tmpRow;
    }
}


/**
 * @brief Copies AFM data.
 * 
 * @param afm The AFM data to copy.
 * @param afmCopy The copied version
 */
void afmData_copy(AFMData *afm, AFMData *afmCopy){
    int row, col;
    // using the new afm data method to create the copied version
    *afmCopy = afmData_new(afm->xResolution, afm->yResolution);

    // Copying the resolutions
    afmCopy->xResolution = afm->xResolution;
    afmCopy->yResolution = afm->yResolution;
    // for loops to copy the data over
    for(row = 0; row < afm->xResolution; row++){
        for(col = 0; col < afm->yResolution; col++){
            afmCopy->zValues[row][col] = afm->zValues[row][col];
        }
    }
}
