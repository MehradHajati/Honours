#include "BandContrast.h"

 
/**
 * @brief Creates a new band contrast with the given dimensions and allocated the memory space for each component.
 * 
 * @param nrow The number of rows. 
 * @param ncol The number of columns.
 * @return BandContrast A band contrast struct.
 */
BandContrast bandContrast_new(int nrow, int ncol){
    int row;
    BandContrast bc;

    bc.nrow = nrow;
    bc.ncol = ncol;
    // allocating the required space for the pictures of EBDS, one for grey, one for red, blue and green
    bc.greyScale = (double**)malloc(sizeof(double*) * nrow);
    bc.EBSDred = (int**)malloc(sizeof(int*) * nrow);
    bc.EBSDgreen = (int**)malloc(sizeof(int*) * nrow);
    bc.EBSDblue = (int**)malloc(sizeof(int*) * nrow);
    for(row = 0; row < bc.nrow; row++){
        // bc.greyScale[row] = (double*)malloc(sizeof(double) * ncol);
        bc.greyScale[row] = (double *)calloc(ncol, sizeof(double));
        bc.EBSDred[row] = (int *)calloc(ncol, sizeof(int));
        bc.EBSDgreen[row] = (int *)calloc(ncol, sizeof(int));
        bc.EBSDblue[row] = (int *)calloc(ncol, sizeof(int));
    }

    return bc;
}
 
/**
 * @brief Frees memory allocate by the given band contrast struct.
 * 
 * @param bc The band contrast to free. 
 */
void bandContrast_free(BandContrast *bc){
    int row;
    for(row = 0; row < bc->nrow; row++){
        free(bc->greyScale[row]);
        free(bc->EBSDred[row]);
        free(bc->EBSDgreen[row]);
        free(bc->EBSDblue[row]);
    }
    free(bc->greyScale);
    free(bc->EBSDred);
    free(bc->EBSDgreen);
    free(bc->EBSDblue);
}


/**
 * @brief Creates a band contrast based on the given light source.
 * 
 * @param facets The facets of the sample.
 * @param light Vector in the direction of the illumination.
 * @return BandContrast A light-based band contrast.
 */
BandContrast bandContrast_light(Facets *facets, Vector3 light){
    return bandContrast_fromPhiThetaMaps(facets->phiMap, facets->thetaMap, facets->nrow, facets->ncol, light);
}

/**
 * @brief A utility used by bandContrast_light.
 * 
 * @param phimap The Phimap data.
 * @param thetamap The Thetamap data.
 * @param nrow Number of rows.
 * @param ncol Number of cols.
 * @param light Vector in the direction of the illumination.
 * @return BandContrast A light-based band contrast.
 */
BandContrast bandContrast_fromPhiThetaMaps(double **phimap, double **thetamap, int nrow, int ncol, Vector3 light){
    int row, col;
    double cos2, angle;
    Vector3 n, negLight;
    // allocating the memory for the new bandcontrast object
    BandContrast bc = bandContrast_new(nrow, ncol);

    // converting all the componenet in the vector to negative version of themselves
    negLight = vector3_scale(light, -1.0);
    
    for(row = 0; row < nrow; row++){
        for(col = 0; col < ncol; col++){
            // Don't have info so setting the greyscale equal to zero
            if(thetamap[row][col] == THETA_DEFAULT){
                bc.greyScale[row][col] = 0;
                continue;
            }

            n = vector3_fromThetaPhi(thetamap[row][col], phimap[row][col]);
            angle = vector3_angleBetween(n, negLight);
            // printf("Angle between @ [%d,%d]: %f\n", row, col, angle * 180.0 / M_PI);
            // If angle is sufficiently small
            if(angle < M_PI_2){
                cos2 = cos(angle);
                //cos2 *= cos2;
                bc.greyScale[row][col] = cos2;
            }
            // Tilted away from light
            else bc.greyScale[row][col] = 0;
        }
    }

    return bc;
}

/**
 * @brief Create a band contrast based on a detector view. Monte Carlo approach.
 * 
 * @param facets The facets of the sample.
 * @param afm The height information of the sample.
 * @return BandContrast A detector-based band contrast.
 */
BandContrast bandContrast_detector(Facets *facets, AFMData *afm){
    int i, j, thetaIters, phiIters, randomIters;
    double theta, phi, deltaTheta, deltaPhi, avgZ, thetaMax, phiMax;
    Vector3 afmCOM, randomVec;
    Quaternion qTilt;
    BandContrast bc = bandContrast_new(afm->xResolution, afm->yResolution);

    thetaIters = 10;
    phiIters = 5;
    randomIters = 100;
    thetaMax = 40.0 * M_PI / 180.0;
    phiMax = 2.0 * M_PI;

    deltaTheta = thetaMax / (double)thetaIters;
    deltaPhi = phiMax / (double)phiIters;

    bandContrast_prepAFMForDetector(afm, &avgZ, &afmCOM);

    for(i = 1; i <= randomIters; i++){
        if(i % 10 == 0 || i == 1){
            printf("Tilt %d/%d\n", i, randomIters);
            fflush(stdout);
        }
        theta = randomPercentage() * thetaMax;
        phi = randomPercentage() * phiMax;
        qTilt = bandContrast_tiltFromThetaPhi(theta, phi);
        bandContrast_tiltForDetector(&bc, afm, qTilt, avgZ);
    }

    // Normalize greyscale
    for(i = 0; i < bc.nrow; i++){
        for(j = 0; j < bc.ncol; j++){
            // bc.greyScale[i][j] /= (double)(phiIters * thetaIters - phiIters + 1);
            bc.greyScale[i][j] /= (double)randomIters;
        }
    }

    return bc;
}

/**
 * @brief Creates a list of band contrast structs (light, detector, and a linear combination of the two) scaled to match the average greyscale of the measured band contrast.
 * 
 * @param facets The facets of the smaple.
 * @param afm The height information of the sample.
 * @param light Vector in the direction of the illumination.
 * @param alpha The fraction of the linear combination given to the light-based band contrast, rest is detector-based
 * @param measured The measured band contrast.
 * @return A list of band contrasts.
 */
BandContrast *simulateBandContrast(Facets *facets, AFMData *afm, Vector3 light, double alpha, BandContrast *measured){
    int i, j, row, col, randomIters, numBCs, iLight = 0, iVis = 1, iLinComboA = 2;
    double theta, phi, avgZ, avgVis, avgLight, avgLinCombo, maxVis, maxResult, thetaMax, phiMax, iVisScaleFactor, iLightScaleFactor, iLinComboScaleFactor, avgMeasured;
    Vector3 afmCOM, randomVec;
    Quaternion qTilt;
    numBCs = 3;
    BandContrast *bcs;

    // the random iters are for the monte carlo method and are set automatically at 100, for more iterations this number would have to be increased
    randomIters = 100;
    // the theta and phi max are set as 70 degrees and 2 degrees, these are preset values and are converted to radians. 
    thetaMax = 70.0 * M_PI / 180.0;
    phiMax = 2.0 * M_PI;
    
    bcs = (BandContrast *)calloc(numBCs, sizeof(BandContrast));
    // calling the bandContrast_light method which will create the light view image which tell us what the light can see based on the given light angle
    bcs[iLight] = bandContrast_light(facets, light);
    // this for loop depends on the preset value of iVis which is set to 1 so it only creates one detector view picture
    // in order for more detector pictures to be created the preset iVis value needs to be increased. 
    for(i = iVis; i < numBCs; i++){
        bcs[i] = bandContrast_new(afm->xResolution, afm->yResolution);
    }

    // find the minimum and the average and then reduces all z data by the minimum and the average aswell
    // This way the smallest z value is at 0
    bandContrast_prepAFMForDetector(afm, &avgZ, &afmCOM);

    // Create iVis
    // this for loop depends on the randomIters value which is preset to 100 so a hundred iterations
    for(i = 1; i <= randomIters; i++){
        if(i % 10 == 0){
            printf("iVis tilt %d/%d...\n", i, randomIters);
            fflush(stdout);
        }
        // doing the monte carlo method here by creating random theta and phi values
        theta = randomPercentage() * thetaMax;
        phi = randomPercentage() * phiMax;
        qTilt = bandContrast_tiltFromThetaPhi(theta, phi);
        bandContrast_tiltForDetector(&bcs[iVis], afm, qTilt, avgZ);
    }

    // Normalize iVis greyscale
    printf("normalizing iVis...\n");
    fflush(stdout);
    for(i = 0; i < bcs[iVis].nrow; i++){
        for(j = 0; j < bcs[iVis].ncol; j++){
            bcs[iVis].greyScale[i][j] = (bcs[iVis].greyScale[i][j]*bcs[iVis].greyScale[i][j]) / (double)(randomIters * randomIters);
        }
    }

    // Scale iVis
    printf("Scaling iVis based on measured...\n");
    // mAvg / iVis avg
    // double iVisMax = -1;
    double iVisAvg = bandContrast_averageGrey(&bcs[iVis], 300, bcs[iVis].nrow - 300, 150, bcs[iVis].ncol - 150);
    printf("iVisAvg Before: %f\n", iVisAvg);
    avgLight = bandContrast_averageGrey(&bcs[iLight], 300, bcs[iLight].nrow - 300, 150, bcs[iLight].ncol - 150);
    printf("iLightAvg Before: %f\n", avgLight);
    avgMeasured = bandContrast_averageGrey(measured, 0, measured->nrow, 0, measured->ncol);
    iVisScaleFactor = avgMeasured / iVisAvg;
    iLightScaleFactor = avgMeasured / avgLight;

    printf("iVisScaleFactor: %f\n", iVisScaleFactor);
    printf("iLightScaleFactor: %f\n", iLightScaleFactor);
    for(row = 0; row < bcs[iVis].nrow; row++){
        for(col = 0; col < bcs[iVis].ncol; col++){
            bcs[iVis].greyScale[row][col] *= iVisScaleFactor;
            bcs[iLight].greyScale[row][col] *= iLightScaleFactor;
        }
    }
    iVisAvg = bandContrast_averageGrey(&bcs[iVis], 300, bcs[iVis].nrow - 300, 150, bcs[iVis].ncol - 150);
    avgLight = bandContrast_averageGrey(&bcs[iLight], 300, bcs[iLight].nrow - 300, 150, bcs[iLight].ncol - 150);
    printf("iVisAvg After: %f\n", iVisAvg);
    printf("iLightAvg After: %f\n", avgLight);

    printf("Creating linCombo...\n\n");
    fflush(stdout);

    for(i = 0; i < bcs[iVis].nrow; i++){
        for(j = 0; j < bcs[iVis].ncol; j++){
            // iLinComboA
            bcs[iLinComboA].greyScale[i][j] = alpha * bcs[iLight].greyScale[i][j] + (1-alpha) * bcs[iVis].greyScale[i][j];
        }
    }

    avgLinCombo = bandContrast_averageGrey(&bcs[iLinComboA], 300, bcs[iLinComboA].nrow - 300, 150, bcs[iLinComboA].ncol - 150);
    iLinComboScaleFactor = avgMeasured / avgLinCombo;
    for(row = 0; row < bcs[iLinComboA].nrow; row++){
        for(col = 0; col < bcs[iLinComboA].ncol; col++){
            bcs[iLinComboA].greyScale[row][col] *= iLinComboScaleFactor;
        }
    }

    return bcs;
}


/**
 * @brief Does tilting for the detector-based band contrast.
 * 
 * @param bc The detector-based band contrast.
 * @param afm The height information of the sample.
 * @param qTilt The quaternion by which to tilt the sample.
 * @param avgZ The average height of the sample (used for tilting about COM).
 */
void bandContrast_tiltForDetector(BandContrast *bc, AFMData *afm, Quaternion qTilt, double avgZ){
    int row, col, newRow, newCol, isVisible, *lastVisibleRows, visRowWidth;
    Vector3 pos, newPos;

    visRowWidth = (int)(sqrt((double)(afm->yResolution * afm->yResolution + afm->xResolution * afm->xResolution) * 3));

    lastVisibleRows = (int *)malloc(sizeof(int) * visRowWidth);
    for(col = 0; col < visRowWidth; col++){
        lastVisibleRows[col] = -9999999;
    }

    for(row = 0; row < afm->xResolution; row++){
        for(col = 0; col < afm->yResolution; col++){
            pos = vector3_new(col, row, (afm->zValues[row][col] - avgZ));
            newPos = pos;
            newPos.x -= (double)afm->yResolution * 0.5;
            newPos.y -= (double)afm->xResolution * 0.5;
            quaternion_rotateVector3ByQuat(&newPos, qTilt);
            newPos.x += (double)afm->yResolution * 0.5;
            newPos.y += (double)afm->xResolution * 0.5;

            newRow = floor(newPos.y + 0.5);
            newCol = floor(newPos.x + 0.5) + (visRowWidth / 3.0);

            // // Check out of bounds
            if(newCol < 0 || newCol >= visRowWidth) {
                continue;
            }

            // Check if visible
            isVisible = newRow > lastVisibleRows[newCol];

            if(isVisible){
                bc->greyScale[row][col] += 1.0;
                lastVisibleRows[newCol] = newRow;
            }
        }
    }

    free(lastVisibleRows);
}


/**
 * @brief Ensures the minimum of the height values rests at 0.
 * 
 * @param afm The AFM data to prep.
 * @param avgZ Stores the average height in this variable.
 * @param afmCOM Stores the COM of the sample in this variable.
 */
void bandContrast_prepAFMForDetector(AFMData *afm, double *avgZ, Vector3 *afmCOM){
    int row, col;
    // preset value to find the smallest z value in the data, so that atleast one z value will be smaller than this value
    double minZ = 1000000;
    // setting this value to zero since it will be used for the total of the z values first and then to hold the average of the z values
    *avgZ = 0;
    // Prep afm
    // going through the data and finding the min z value and keeping track of the total of the z values in the data
    for(row = 0; row < afm->xResolution; row++){
        for(col = 0; col < afm->yResolution; col++){
            // adding for the calculating of the average later
            *avgZ += afm->zValues[row][col];
            // comparing the current value to the minZ and replacing it if it smaller
            minZ = minZ < afm->zValues[row][col] ? minZ : afm->zValues[row][col];
        }
    }
    // calculating the average of the z values by diving the total number of pixels by the total
    *avgZ /= (double)(afm->xResolution * afm->yResolution);
    // when we decrease everything by minz, we also need to adjust the average since it will also decrease by minZ
    *avgZ -= minZ;
    // reducing all the values of the data by the minimum that they have, so that the lowest value in the data is zero
    for(row = 0; row < afm->xResolution; row++){
        for(col = 0; col < afm->yResolution; col++){
            afm->zValues[row][col] -= minZ;
        }
    }
    // creating a new vector with the x being half of x-resolution and the y being half of the y-resolution and the z component being the average of the z values
    *afmCOM = vector3_new(afm->xResolution / 2.0, afm->yResolution / 2.0, *avgZ);
}

 
/**
 * @brief Creates a rotation quaternion from the given theta and phi angles.
 * 
 * @param theta Radians down from detector normal.
 * @param phi Radians around detector normal.
 * @return A rotation quaternion.
 */
Quaternion bandContrast_tiltFromThetaPhi(double theta, double phi){
    double st, yTheta, xTheta;
    Vector3 negXAxis;
    Quaternion q70, qY, qX, qTilt;

    // If the phi angle is more than or equal to pi, the theta value is negated
    xTheta = (phi >= M_PI) ? -theta : theta;

    negXAxis = vector3_new(-1.0,0.0,0.0);

    // Rotates by 70 degrees about the x-axis, aligning sample with center of detector.
    q70 = quaternion_new(negXAxis.x * (35.0 * M_PI / 180.0), 0, 0, cos(35.0 * M_PI / 180.0));

    // Rotation about Y axis
    while(phi >= M_PI) phi -= M_PI;
    yTheta = M_PI_2 - phi;
    st = sin(yTheta / 2.0);
    qY = quaternion_new(0, st, 0, cos(yTheta / 2.0));

    // Rotation about X axis
    st = sin(xTheta / 2.0);
    qX = quaternion_new(st, 0, 0, cos(xTheta / 2.0));

    // Combine rotations into one
    qTilt = quaternion_multiply(qX, qY);
    quaternion_normalize(&qTilt);
    qTilt = quaternion_multiply(qTilt, q70);
    quaternion_normalize(&qTilt);

    return qTilt;
}


/**
 * @brief Computes the average grey value in the given band contrast within the given bounds.
 * 
 * @param bc The band contrast from which to get the average grey.
 * @param startRow The starting row.
 * @param endRow The ending row.
 * @param startCol The starting column.
 * @param endCol The ending column.
 * @return The average grey.
 */
double bandContrast_averageGrey(BandContrast *bc, int startRow, int endRow, int startCol, int endCol){
    int row, col, nrow, ncol;
    double avg = 0.0;

    nrow = endRow - startRow;
    ncol = endCol - startCol;

    for(row = startRow; row < endRow; row++){
        for(col = startCol; col < endCol; col++){
            avg += bc->greyScale[row][col];
        }
    }

    avg /= (double)(nrow * ncol);
    return avg;
}


/**
 * @brief Computes the standard deviation in the given band contrast within the given bounds.
 * 
 * @param bc The band contrast from which to get the standard deviation.
 * @param startRow The starting row.
 * @param endRow The ending row.
 * @param startCol The starting column.
 * @param endCol The ending column.
 * @return The standard deviation. 
 */
double bandContrast_stdDev(BandContrast *bc, int startRow, int endRow, int startCol, int endCol){
    int row, col, nrow, ncol;
    double xi2 = 0, xi = 0;

    nrow = endRow - startRow;
    ncol = endCol - startCol;

    for(row = startRow; row < endRow; row++){
        for(col = startCol; col < endCol; col++){
            xi2 += bc->greyScale[row][col] * bc->greyScale[row][col];
            xi += bc->greyScale[row][col];
        }
    }

    double deviation = sqrt((xi2 - xi*xi/(double)(nrow * ncol))/(double)(nrow * ncol));
    return deviation;
}


/**
 * @brief Multiplies all values by 255.
 * 
 * @param vals The values to scale.
 * @param nrow The number of rows to scale.
 * @param ncol The number of columns to scale.
 */
void bandContrast_scaleTo255(double ***vals, int nrow, int ncol){
    int row, col;

    for(row = 0; row < nrow; row++){
        for(col = 0; col < ncol; col++){
            (*vals)[row][col] = floor((*vals)[row][col] * 255 + 0.5);
        }
    }
}


/**
 * @brief Flood-fills 0 values in the given band contrast.
 * 
 * @param bc The band contrast to fill.
 */
void bandContrast_fillGaps(BandContrast *bc){
    int dirs[4][2] = {
        { 1, 0 },
        { 0, 1 },
        { -1, 0 },
        { 0, -1 }
    };
    int min = 0, max = 1, found;
    int row, col, dir, newRow, newCol, numDirs = 4, avgCount, **minMax;
    double avg;

    minMax = (int **)malloc(bc->ncol * sizeof(int *));
    for(col = 0; col < bc->ncol; col++){
        minMax[col] = (int *)malloc(sizeof(int) * 2);
        minMax[col][min] = 0;
        minMax[col][max] = bc->ncol;

        found = 0;
        row = 0;
        while(row < bc->nrow && !found){
            if((int)bc->greyScale[row][col] != 0){
                minMax[col][min] = row;
                found = 1;
            }
            row++;
        }

        found = 0;
        row = bc->nrow-1;
        while(row >= 0 && !found){
            if((int)bc->greyScale[row][col] != 0){
                minMax[col][max] = row;
                found = 1;
            }
            row--;
        }
    }



    for(row = 0; row < bc->nrow; row++){
        for(col = 0; col < bc->ncol; col++){
            if(row >= minMax[col][min] && row <= minMax[col][max] && (int)bc->greyScale[row][col] == 0){
                avg = 0;
                avgCount = 0;
                // fill hole
                for(dir = 0; dir < numDirs; dir++){
                    newRow = row + dirs[dir][0];
                    newCol = col + dirs[dir][1];
                    if(newRow >= 0 && newRow < bc->nrow && newCol >= 0 && newCol < bc->ncol && (int)bc->greyScale[newRow][newCol] != 0){
                        avg += bc->greyScale[newRow][newCol];
                        avgCount++;
                    }
                }
                if(avgCount != 0){
                    bc->greyScale[row][col] = avg / (double)avgCount;
                }
            }
        }
    }

    for(col = 0; col < bc->ncol; col++){
        free(minMax[col]);
    }
    free(minMax);
}
