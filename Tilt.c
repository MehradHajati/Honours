#include "Tilt.h"
 
/**
 * @brief Creates a new sample struct that can be tilted.
 * 
 * @param afm The afm data to tilt with the sample.
 * @param bc The band contrast data to tilt with the sample.
 * @return The new sample (untilted). 
 */
Sample sample_new(AFMData *afm, BandContrast *bc){
    int row, col;
    Sample smpl;
    smpl.numRows = afm->xResolution;
    smpl.numCols = afm->yResolution;

    // allocating all the first layer of memory required
    smpl.data = (double ***)malloc(smpl.numRows * sizeof(double));
    smpl.tilt = (double ***)malloc(smpl.numRows * sizeof(double));
    smpl.assigned = (int **)malloc(smpl.numRows * sizeof(int *));
    smpl.origCol = (int **)malloc(smpl.numRows * sizeof(int *));
    smpl.origRow = (int **)malloc(smpl.numRows * sizeof(int *));
    smpl.bestZ = (double **)malloc(smpl.numRows * sizeof(double *));
    // allocate origng stuff 

    for(row = 0; row < smpl.numRows; row++){
        // allocating the memory for the second layer of memory to create the matrices
        smpl.data[row] = (double **)malloc(smpl.numCols * sizeof(double));
        smpl.tilt[row] = (double **)malloc(smpl.numCols * sizeof(double));
        smpl.assigned[row] = (int *)calloc(smpl.numCols, sizeof(int));
        smpl.origCol[row] = (int *)calloc(smpl.numCols, sizeof(int));
        smpl.origRow[row] = (int *)calloc(smpl.numCols, sizeof(int));
        smpl.bestZ[row] = (double *)calloc(smpl.numCols, sizeof(double));

        for(col = 0; col < smpl.numCols; col++){
            // allocating the third layer for 3d matrices
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

/**
 * @brief Frees memory allocated by the given sample.
 * 
 * @param smpl The sample to free.
 */
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

/**
 * @brief Tilts the given sample.
 * 
 * @param smpl The sample to tilt.
 * @param phi The angle by which to tilt the sample.
 */
void sample_tilt(Sample *smpl, double phi){
// last argument is always 1.0
    int row, col, tmpRow, zCount = 0, lastMappedRow, newCol, newRow;
    double minZ = 1000000.0, avgZ = 0.0;
    double newZ, oldZ;
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
        lastMappedRow = -1;
        for(row = 0; row < smpl->numRows; row++){
            pos = vector3_new(col, row, (smpl->data[row][col][Z_VALUES] - minZ - avgZ));
            newPos = pos;
            // center everything for the tilt
            newPos.x -= (double)smpl->numCols * 0.5;  
            newPos.y -= (double)smpl->numRows * 0.5;
            quaternion_rotateVector3AxisAngle(&newPos, -1.0, 0, 0, phi);  
            // push everything back          
            newPos.x += (double)smpl->numCols * 0.5;  
            newPos.y += (double)smpl->numRows * 0.5;
            newCol = floor(newPos.x + 0.5);
            newRow = floor(newPos.y + 0.5);
            newZ = smpl->data[row][col][Z_VALUES] - minZ;
            oldZ = newZ -avgZ;

            //mapOldPointToTilt
            if(newRow>lastMappedRow && newRow>=0 && newRow<smpl->numRows && newCol>=0 && newCol<smpl->numCols) { 
                lastMappedRow = newRow; 
                smpl->tilt[newRow][newCol][Z_VALUES] = oldZ;
                smpl->tilt[newRow][newCol][BAND_CONTRAST] = smpl->data[row][col][BAND_CONTRAST];
                smpl->bestZ[newRow][newCol] = newZ;
                smpl->assigned[newRow][newCol] = 1;
                smpl->origRow[newRow][newCol] = row;
                smpl->origCol[newRow][newCol] = col;
            }
        }
    }
}
 
/**
 * @brief Creates an AFMData struct from the AFM data stored in data.
 * 
 * @param numRows The number of rows in the data. 
 * @param numCols The number of columnc in the data.
 * @param data The data containing the afm data (from a sample struct).
 * @return An AFMData struct. 
 */
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

/**
 * @brief Creates a BandContrast struct from the band contrast data stored in data.
 * 
 * @param numRows The number of rows in the data.
 * @param numCols The number of columnc in the data.
 * @param data The data containing the band contrast data (from a sample struct).
 * @return A BandContrast struct. 
 */
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

 
/**
 * @brief Creates an AFMData struct from the tilted AFM data in the given sample.
 * 
 * @param smpl The sample from which to get the AFM data.
 * @return An AFMData struct.
 */
AFMData sample_getTiltedAFMFromSample(Sample *smpl){
    return sample_getAFM(smpl->numRows, smpl->numCols, &smpl->tilt);
}
 
/**
 * @brief Creates a BandContrast struct from the tilted band contrast data in the given sample.
 * 
 * @param smpl The sample from which to get the band contrast data.
 * @return A BandContrast struct. 
 */
AFMData sample_getAFMFromSample(Sample *smpl){
    return sample_getAFM(smpl->numRows, smpl->numCols, &smpl->data);
}

/**
 * @brief Creates an AFMData struct from the untilted AFM data in the given sample.
 * 
 * @param smpl The sample from which to get the AFM data.
 * @return An AFMData struct.
 */
BandContrast sample_getTiltedBandContrastFromSample(Sample *smpl){
    return sample_getBandContrast(smpl->numRows, smpl->numCols, &smpl->tilt);
}

/**
 * @brief Creates a BandContrast struct from the untilted band contrast data in the given sample.
 * 
 * @param smpl The sample from which to get the band contrast data.
 * @return A BandContrast struct. 
 */
BandContrast sample_getBandContrastFromSample(Sample *smpl){
    return sample_getBandContrast(smpl->numRows, smpl->numCols, &smpl->data);
}
