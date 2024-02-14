#include "Tilt.h"

Sample sample_new(AFMData *afm, BandContrast *bc){
    int row, col;
    Sample smpl;
    smpl.numRows = afm->xResolution;
    smpl.numCols = afm->yResolution;

    smpl.data = (double ***)malloc(smpl.numRows * sizeof(double));
    smpl.tilt = (double ***)malloc(smpl.numRows * sizeof(double));
    smpl.assigned = (int **)malloc(smpl.numRows * sizeof(int *));
    smpl.origCol = (int **)malloc(smpl.numRows * sizeof(int *));
    smpl.origRow = (int **)malloc(smpl.numRows * sizeof(int *));
    smpl.bestZ = (double **)malloc(smpl.numRows * sizeof(double *));
    // allocate origng stuff 

    for(row = 0; row < smpl.numRows; row++){
        smpl.data[row] = (double **)malloc(smpl.numCols * sizeof(double));
        smpl.tilt[row] = (double **)malloc(smpl.numCols * sizeof(double));
        smpl.assigned[row] = (int *)calloc(smpl.numCols, sizeof(int));
        smpl.origCol[row] = (int *)calloc(smpl.numCols, sizeof(int));
        smpl.origRow[row] = (int *)calloc(smpl.numCols, sizeof(int));
        smpl.bestZ[row] = (double *)calloc(smpl.numCols, sizeof(double));

        for(col = 0; col < smpl.numCols; col++){
            smpl.data[row][col] = (double *)malloc(SAMPLE_DEPTH * sizeof(double));
            smpl.tilt[row][col] = (double *)malloc(SAMPLE_DEPTH * sizeof(double));
            
            smpl.data[row][col][Z_VALUES] = afm->zValues[row][col];
            // smpl.data[row][col][EULER] = (double)rgb->values[row][col];
            smpl.data[row][col][BAND_CONTRAST] = (double)bc->greyScale[row][col];

            smpl.tilt[row][col][Z_VALUES] = Z_DEFAULT;
            // smpl.tilt[row][col][EULER] = 0;
            smpl.tilt[row][col][BAND_CONTRAST] = 0;
        }
    }
    return smpl;
}

void sample_free(Sample *smpl){
    int row, col;
    for(row = 0; row < smpl->numRows; row++){
        for(col = 0; col < smpl->numCols; col++){
            free(smpl->data[row][col]);
            free(smpl->tilt[row][col]);
        }
        free(smpl->data[row]);
        free(smpl->tilt[row]);
        free(smpl->assigned[row]);
        free(smpl->origCol[row]);
        free(smpl->origRow[row]);
        free(smpl->bestZ[row]);
        // free origin stuff
    }
    free(smpl->data);
    free(smpl->tilt);
    free(smpl->assigned);
    free(smpl->origRow);
    free(smpl->origCol);
    free(smpl->bestZ);
}

void sample_tilt(Sample *smpl, double phi){
// last argument is always 1.0
    int row, col, isMappable, tmpRow, zCount = 0, lastMappedRow, newCol, newRow;
    double minZ = 1000000.0, avgZ = 0.0;
    Vector3 pos, newPos;

    // Find minZ and avgZ
    for(row = 0; row < smpl->numRows; row++){
        for(col = 0; col < smpl->numCols; col++){
            zCount++;
            avgZ += smpl->data[row][col][Z_VALUES];
            if(minZ > smpl->data[row][col][Z_VALUES]){
                minZ = smpl->data[row][col][Z_VALUES];
            }
        }
    }
    avgZ /= (double)zCount;
    //minZ -= avgZ; // from Sawyer, seems to be not needed 

    // Tilt everything
    for(col = 0; col < smpl->numCols; col++){
        isMappable = 0;
        lastMappedRow = -1;
        for(row = 0; row < smpl->numRows; row++){
            pos = vector3_new(col, row, (smpl->data[row][col][Z_VALUES] - minZ - avgZ));
            newPos = pos;
            newPos.x -= (double)smpl->numCols * 0.5;  // center everything for the tilt
            newPos.y -= (double)smpl->numRows * 0.5;
            quaternion_rotateVector3AxisAngle(&newPos, -1.0, 0, 0, phi);            
            newPos.x += (double)smpl->numCols * 0.5;  // push everything back
            newPos.y += (double)smpl->numRows * 0.5;

            newCol = floor(newPos.x + 0.5);
            newRow = floor(newPos.y + 0.5);
            isMappable = (newRow > lastMappedRow);  // this sets the isMappable to 1 or 0

            if( isMappable &&
                sample_mapOldPointToTilt( smpl, (int)pos.y, (int)pos.x, newRow, newCol,
                smpl->data[(int)pos.y][(int)pos.x][Z_VALUES] - minZ, smpl->data[row][col][Z_VALUES] - minZ - avgZ))
            {
                lastMappedRow = newRow; 
            }   
        }
    }
}

int sample_mapOldPointToTilt(Sample *smpl, int row, int col, int newRow, int newCol, double newZ, double oldZ){
    if(0 <= newRow && newRow < smpl->numRows && 0 <= newCol && newCol < smpl->numCols){
        smpl->tilt[newRow][newCol][Z_VALUES] = oldZ;
        smpl->tilt[newRow][newCol][BAND_CONTRAST] = smpl->data[row][col][BAND_CONTRAST];
        smpl->bestZ[newRow][newCol] = newZ;
        smpl->assigned[newRow][newCol] = 1;
        smpl->origRow[newRow][newCol] = row;
        smpl->origCol[newRow][newCol] = col;
        return 1;
    }
    return 0;
}

AFMData sample_getAFM(int numRows, int numCols, double ****data){
    int row, col;
    AFMData afmData = afmData_new(numRows, numCols);

    for(row = 0; row < numRows; row++){
        for(col = 0; col < numCols; col++){
            afmData.zValues[row][col] = (*data)[row][col][Z_VALUES];
        }
    }

    return afmData;
}

BandContrast sample_getBandContrast(int numRows, int numCols, double ****data){
    int row, col;
    BandContrast bc = bandContrast_new(numRows, numCols);

    for(row = 0; row < numRows; row++){
        for(col = 0; col < numCols; col++){
            bc.greyScale[row][col] = (*data)[row][col][BAND_CONTRAST];
        }
    }

    return bc;
}

AFMData sample_getTiltedAFMFromSample(Sample *smpl){
    return sample_getAFM(smpl->numRows, smpl->numCols, &smpl->tilt);
}

AFMData sample_getAFMFromSample(Sample *smpl){
    return sample_getAFM(smpl->numRows, smpl->numCols, &smpl->data);
}

BandContrast sample_getTiltedBandContrastFromSample(Sample *smpl){
    return sample_getBandContrast(smpl->numRows, smpl->numCols, &smpl->tilt);
}

BandContrast sample_getBandContrastFromSample(Sample *smpl){
    return sample_getBandContrast(smpl->numRows, smpl->numCols, &smpl->data);
}
