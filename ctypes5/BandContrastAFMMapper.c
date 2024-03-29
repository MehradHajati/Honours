#include "BandContrastAFMMapper.h"

 
/**
 * @brief Creates a new band contrast AFM mapper.
 * 
 * @param nrow The number of rows.
 * @param ncol The number of columns.
 * @return A band contrast AFM mapper with the given dimensions.
 */
BandContrastAFMMapper bandContrastAFMMapper_new(int nrow, int ncol){
    int row, col, layer;
    BandContrastAFMMapper bcAFMm;
    bcAFMm.nrow = nrow;
    bcAFMm.ncol = ncol;

    //allocating the space for the 3d matrix/map within the band contrast mapper 
    bcAFMm.map = (double ***)malloc(sizeof(double **) * NUMBER_OF_LAYERS_IN_BCAFMM);

    for(layer = 0; layer < NUMBER_OF_LAYERS_IN_BCAFMM; layer++){
        bcAFMm.map[layer] = (double **)malloc(sizeof(double *) * bcAFMm.nrow);
        for(row = 0; row < bcAFMm.nrow; row++){
            bcAFMm.map[layer][row] = (double *)malloc(sizeof(double) * bcAFMm.ncol);
        }
    }

    // setting all the deafult values for the 3d matrix/map
    for(row = 0; row < bcAFMm.nrow; row++){
        for(col = 0; col < bcAFMm.ncol; col++){
            bcAFMm.map[GREYSCALE_LAYER][row][col] = GREYSCALE_DEFAULT;
            bcAFMm.map[OLD_ROW_LAYER][row][col] = OLD_POSITION_DEFAULT;
            bcAFMm.map[OLD_COL_LAYER][row][col] = OLD_POSITION_DEFAULT;
        }
    }

    return bcAFMm;
}

/**
 * @brief Frees memory allocated by the given band contrast AFM mapper.
 * 
 * @param bcAFMm The band contrast AFM mapper to free.
 */
void bandContrastAFMMapper_free(BandContrastAFMMapper *bcAFMm){
    int row, col, layer;
    // freeing the memory for every 3rd level entry in the 3d matrix
    for(layer = 0; layer < NUMBER_OF_LAYERS_IN_BCAFMM; layer++){
        for(row = 0; row < bcAFMm->nrow; row++){
            free(bcAFMm->map[layer][row]);
        }
        // freeing the 2nd level entries in the 3d matrix
        free(bcAFMm->map[layer]);
    }
    // freeing the memory for the map
    free(bcAFMm->map);
}


/**
 * @brief Maps the given band contrast to the given AFM according to the given parameters.
 * 
 * @param bcMeasured The band contrast to map.
 * @param afmTilted The AFM upon which to map the band contrast.
 * @param a0 x pos
 * @param a1 x scale
 * @param a2 x skew
 * @param a3 x stretch
 * @param a4 x curve
 * @param a5 x star?
 * @param b0 y pos.
 * @param b1 y skew.
 * @param b2 y scale.
 * @param b3 y curve.
 * @param b4 y stretch.
 * @param b5 y star?
 * @return The measured band contrast mapped onto the given AFM. 
 */
BandContrastAFMMapper bandContrastAFMMapper_map(BandContrast bcMeasured, AFMData afmTilted, double a0, double a1, double a2, double a3, double a4, double a5, double b0, double b1, double b2, double b3, double b4, double b5){
    //printf("Mapping band contrast onto AFM...\n");
    fflush(stdout);
    int row, col, i, j, newRow, newCol, midRow, midCol, fracIters = 11;
    double fraction = 0.1, colDiff, rowDiff;

    midRow = bcMeasured.nrow / 2; // Y
    midCol = bcMeasured.ncol / 2; // X
    BandContrastAFMMapper bcAFMm = bandContrastAFMMapper_new(afmTilted.xResolution, afmTilted.yResolution);

// NOTE: if value is GREYSCALE_DEFAULT then that pixel has not been mapped
    for(row = 0; row < bcMeasured.nrow; row++){
        rowDiff = row - midRow; // y - Y
        for(col = 0; col < bcMeasured.ncol; col++){
            colDiff = -(col - midCol); // x - X
            
            newCol = floor(a0 + a1*rowDiff + a2*colDiff + a3*rowDiff*rowDiff + a4*colDiff*colDiff + a5*rowDiff*colDiff + 0.5); //x from X, Y
            newRow = floor(b0 + b1*rowDiff + b2*colDiff + b3*rowDiff*rowDiff + b4*colDiff*colDiff + b5*rowDiff*colDiff + 0.5); //y from X, Y

            // If within bounds and still default value:
            if(newRow >= 0 && newRow < bcAFMm.nrow && newCol >= 0 && newCol < bcAFMm.ncol && bcAFMm.map[GREYSCALE_LAYER][newRow][newCol] == GREYSCALE_DEFAULT){
                bcAFMm.map[GREYSCALE_LAYER][newRow][newCol] = bcMeasured.greyScale[row][col];
                bcAFMm.map[OLD_ROW_LAYER][newRow][newCol] = row;
                bcAFMm.map[OLD_COL_LAYER][newRow][newCol] = col;
            }
        }
    }

    // Fill fractional positions:
    for(row = 0; row < bcMeasured.nrow; row++){
        for(i = 0; i < fracIters; i++){
            rowDiff = (row + fraction * i - midRow); // y - Y
            for(col = 0; col < bcMeasured.ncol; col++){
                for(j = 0; j < fracIters; j++){
                    colDiff = -(col + fraction * j - midCol); // x - X

                    newCol = floor(a0 + a1*rowDiff + a2*colDiff + a3*rowDiff*rowDiff + a4*colDiff*colDiff + a5*rowDiff*colDiff + 0.5); //x from X, Y
                    newRow = floor(b0 + b1*rowDiff + b2*colDiff + b3*rowDiff*rowDiff + b4*colDiff*colDiff + b5*rowDiff*colDiff + 0.5); //y from X, Y

                    // If within bounds and still default value:
                    if(newRow >= 0 && newRow < bcAFMm.nrow && newCol >= 0 && newCol < bcAFMm.ncol && bcAFMm.map[GREYSCALE_LAYER][newRow][newCol] == GREYSCALE_DEFAULT){
                        bcAFMm.map[GREYSCALE_LAYER][newRow][newCol] = bcMeasured.greyScale[row][col];
                        bcAFMm.map[OLD_ROW_LAYER][newRow][newCol] = row;
                        bcAFMm.map[OLD_COL_LAYER][newRow][newCol] = col;
                    }
                }                
            }
        }
    }

    return bcAFMm;
}

/**
 * @brief Computes chi squared.
 * 
 * @param bcAFMm The band contrast AFM Mapper of which to get chi squared.
 * @param bcTilted The band contrast to which to compare the mapping.
 * @param mStdDev The standard deviation in the mapped area.
 * @param simStdDev The standard deviation in the band contrast.
 * @return Chi squared of the mapping.
 */
double bandContrastAFMMapper_chiSquared(BandContrastAFMMapper *bcAFMm, BandContrast *bcTilted, double mStdDev, double simStdDev){
    double chiSquared = 0.0;
    int row, col, overlappingPoints = 0;
    double diff;

    for(row = 0; row < bcAFMm->nrow; row++){
        for(col = 0; col < bcAFMm->ncol; col++){
            if(bcAFMm->map[GREYSCALE_LAYER][row][col] < GREYSCALE_DEFAULT * 255.0){ // Transparency
                // Chi Squared and difference here?
                diff = bcAFMm->map[GREYSCALE_LAYER][row][col] - bcTilted->greyScale[row][col];
                chiSquared += diff*diff;
                overlappingPoints++;
            }
        }
    }
    if(overlappingPoints != 0) chiSquared /= ((double)overlappingPoints * sqrt((mStdDev*mStdDev)*(simStdDev*simStdDev)));
    else chiSquared = -99999999;
    printf("Chi Squared: %f\n", chiSquared);
    return chiSquared;
}
