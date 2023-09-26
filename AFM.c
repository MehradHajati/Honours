#include "AFM.h"

AFMData afmData_new(int xResolution, int yResolution){
    int i;
    AFMData afmData;
    afmData.xResolution = xResolution;
    afmData.yResolution = yResolution;

    afmData.zValues = (double**)malloc(sizeof(double*) * xResolution);
    for(i = 0; i < xResolution; i++){
        afmData.zValues[i] = (double*)calloc(yResolution, sizeof(double));
    }

    return afmData;
}

void afmData_free(AFMData* afmData){
    int i;
    for(i = 0; i < afmData->xResolution; i++){
        free(afmData->zValues[i]);
    }
    free(afmData->zValues);
}

void afmData_binBy2(AFMData *afm){
    double **zVals;
    int newXRes, newYRes, oldRow, oldCol, newRow, newCol, dir, numDirs = 4;
    int dirs[4][2] = {
        {0,0},
        {1,0},
        {0,1},
        {1,1}
    };
    newXRes = afm->xResolution / 2;
    newYRes = afm->yResolution / 2;
    zVals = (double **)malloc(sizeof(double*) * newXRes);
    for(oldRow = 0; oldRow < newXRes; oldRow++){
        zVals[oldRow] = (double *)calloc(newYRes, sizeof(double));
    }

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

    for(oldRow = 0; oldRow < afm->xResolution; oldRow++){
        free(afm->zValues[oldRow]);
    }
    free(afm->zValues);
    afm->zValues = NULL;
    afm->zValues = zVals;
    afm->xResolution = newXRes;
    afm->yResolution = newYRes;
}

void afmData_flipY(AFMData *afm){
    double *tmpRow;
    int row;
    for(row = 0; row < afm->xResolution / 2; row++){
        tmpRow = afm->zValues[row];
        afm->zValues[row] = afm->zValues[afm->xResolution - 1 - row];
        afm->zValues[afm->xResolution - 1 - row] = tmpRow;
    }
}

void afmData_copy(AFMData *afm, AFMData *afmCopy){
    int row, col;

    *afmCopy = afmData_new(afm->xResolution, afm->yResolution);

    afmCopy->xResolution = afm->xResolution;
    afmCopy->yResolution = afm->yResolution;
    for(row = 0; row < afm->xResolution; row++){
        for(col = 0; col < afm->yResolution; col++){
            afmCopy->zValues[row][col] = afm->zValues[row][col];
        }
    }
}
