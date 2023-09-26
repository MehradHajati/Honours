#include "Facets.h"

Facets facets_new(int nrow, int ncol){
    int row, col;
    Facets facets;
    facets.nrow = nrow;
    facets.ncol = ncol;
    facets.thetaMap = (double **)malloc(sizeof(double*) * facets.nrow);
    facets.phiMap = (double **)malloc(sizeof(double*) * facets.nrow);
    for(row = 0; row < facets.nrow; row++){
        facets.thetaMap[row] = (double*)malloc(sizeof(double) * facets.ncol);
        facets.phiMap[row] = (double*)malloc(sizeof(double) * facets.ncol);
        for(col = 0; col < facets.ncol; col++){
            facets.thetaMap[row][col] = THETA_DEFAULT;
            facets.phiMap[row][col] = PHI_DEFAULT;
        }
    }
    return facets;
}

Face face_new(int row, int col, int size, double rssAvg, double theta, double phi){
    Face face;
    face.row = row;
    face.col = col;
    face.size = size;
    face.rssAvg = rssAvg;
    face.theta = theta;
    face.phi = phi;
    return face;
}

Facets facets_compute(AFMData *afmData, int minFacetSize, int maxFacetSize, double maxRSSAvg, int binIters){
    AFMData afmBinned = *afmData;
    //afmData_flipY(&afmBinned);
    int binIter;
    for(binIter = 0; binIter < binIters; binIter++){
        printf("Binning...\n");
        fflush(stdout);
        afmData_binBy2(&afmBinned);
    }
    Facets facets = facets_new(afmBinned.xResolution, afmBinned.yResolution);
    int tileSize = maxFacetSize;

    while(tileSize >= minFacetSize){
        facets_tileArea(&afmBinned, &facets, tileSize, maxRSSAvg);
        tileSize -= 2;
    }

    afmData_free(&afmBinned);
    return facets;
}

void facets_tileArea(AFMData *afmData, Facets *facets, int tileSize, double maxRSSAvg){
    time_t startTime = time(0);
    int validFacesSize, i, row, col;
    Face *validFaces = facets_findValidFaces(facets, afmData, tileSize, &validFacesSize, maxRSSAvg);
    Face cur;
    qsort(validFaces, validFacesSize, sizeof(Face), facets_cmpFacesByRSS);

    for(i = 0; i < validFacesSize; i++){
        cur = validFaces[i];
        // printf("RSS_avg: %f\n", cur.rssAvg);
        for(row = cur.row; row < cur.row + cur.size; row++){
            for(col = cur.col; col < cur.col + cur.size; col++){
                if(facets->thetaMap[row][col] == THETA_DEFAULT){
                    facets->thetaMap[row][col] = cur.theta;
                    facets->phiMap[row][col] = cur.phi;
                }
            }
        }
    }

    int area = tileSize * tileSize * validFacesSize;
    printf("patch: %d valid: %d area: %d\n", tileSize, validFacesSize, area);
    free(validFaces);
    time_t endTime = time(0);
    printf("Time to tile size %d: %ld seconds\n", tileSize, endTime - startTime);
    fflush(stdout);
}

Face *facets_findValidFaces(Facets *facets, AFMData *afmData, int tileSize, int *outArraySize, double maxRSSAvg){
    int row, col, tileRow, tileCol, count = 0, capacity = 16, all = 0;
    int tileMoveDist = tileSize / 2;
    Face tmpFace;
    Face *faces = NULL;

    // Prep beta Matrix
    struct Matrix *X = matrix_new(tileSize * tileSize, 3);

    for(row = 0; row < X->nrow; row++){
        X->vals[row][0] = 1;
        X->vals[row][1] = (row / 3) - 1;
        X->vals[row][2] = (row % 3) - 1;
    }

    struct Matrix *T = matrix_copy(X);
    matrix_transpose(&T);
    struct Matrix *I = matrix_multiply(T, X);
    matrix_inverse(&I);
    struct Matrix *betaMatrix = matrix_multiply(I, T);

    for(tileRow = 0; tileRow + tileSize <= afmData->xResolution; tileRow += tileMoveDist){
        for(tileCol = 0; tileCol + tileSize <= afmData->yResolution; tileCol += tileMoveDist){
            if(facets->thetaMap[tileRow][tileCol] != THETA_DEFAULT) continue;
            all++;
            if(facets_findFace(afmData, &tmpFace, betaMatrix, X, tileRow, tileCol, tileSize, maxRSSAvg)){
                // try to add tmpFace to faces
                facets_addToFaces(&faces, tmpFace, &count, &capacity);
            }            
        }
    }
    printf("\nall: %d\n", all);
    *outArraySize = count;
    matrix_free(betaMatrix);
    matrix_free(I);
    matrix_free(T);
    matrix_free(X);
    return faces;
}

int facets_findFace(AFMData *afmData, Face *tmpFace, struct Matrix *betaMatrix, struct Matrix *X, int tileRow, int tileCol, int tileSize, double maxRSSAvg)
{
    int row, col;
    double theta, phi, aa, bb, cc, coeffMagnitude, nMagnitude;
    struct Matrix *zVals, *betaHat, *yHat, *e, *et, *RSS;
    Vector3 a, b, n;

    zVals = matrix_new(tileSize*tileSize, 1);
    for(row = tileRow; row < tileRow + tileSize; row++){
        for(col = tileCol; col < tileCol + tileSize; col++){
            zVals->vals[(row - tileRow) * tileSize + (col - tileCol)][0] = afmData->zValues[row][col];
        }
    }
    betaHat = matrix_multiply(betaMatrix, zVals);
    yHat = matrix_multiply(X, betaHat);
    matrix_scalarMultiply(yHat, -1.0);
    e = matrix_add(zVals, yHat);
    matrix_scalarMultiply(yHat, -1.0);
    et = matrix_copy(e);
    matrix_transpose(&et);
    RSS = matrix_multiply(e, et);

    a = vector3_new(1, 0, betaHat->vals[1][0]);
    b = vector3_new(0, 1, betaHat->vals[2][0]);
    n = vector3_new(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y-a.y*b.x);

    nMagnitude = sqrt(n.x*n.x + n.y*n.y + n.z*n.z);
    n = vector3_scale(n, 1.0/nMagnitude);

    // printf("Normal at {%d,%d}: <%f,%f,%f>\n", tileRow, tileCol, n.x, n.y, n.z);

    theta = acos(n.z) * 180.0 / M_PI;
    phi = atan2(n.x,n.y) * 180.0 / M_PI;

    // create face
    //printf("Calculated RSS: %f\n", RSS->vals[0][0]);
    *tmpFace = face_new(tileRow, tileCol, tileSize, RSS->vals[0][0] / (((double)(tileSize * tileSize)) - 3), theta, phi);

    // if(tileRow == 0 && tileCol == 0 && tileSize == 5){
    //    printf("betaMatrix:\n%s\n", matrix_toString(betaMatrix));
    //     printf("Y:\n%s\n", matrix_toString(zVals));
    //     printf("betaHat:\n%s\n", matrix_toString(betaHat));
    //     printf("yHat:\n%s\n", matrix_toString(yHat));
    //     printf("e:\n%s\n", matrix_toString(e));
    //     printf("et:\n%s\n", matrix_toString(et));
    //     printf("RSS:\n%f\n", RSS->vals[0][0]);
    //     printf("RSS_avg:\n%f\n", tmpFace->rssAvg);
    //     printf("theta:\n%f\n", theta);
    //     printf("phi:\n%f\n", phi); 
    //     fflush(stdin);
    // }

    // computeRSS(tmpFace, afmData, coeff);

    // clean up
    matrix_free(zVals);
    matrix_free(betaHat);
    matrix_free(yHat);
    matrix_free(e);
    matrix_free(et);
    matrix_free(RSS);
    return tmpFace->rssAvg < maxRSSAvg;
}

void facets_addToFaces(Face **faces, Face toAdd, int *count, int *capacity){
    if(*faces == NULL) *faces = (Face *)calloc(sizeof(Face), *capacity);
    
    (*faces)[(*count)++] = toAdd;
    if(*count == *capacity){
        *capacity *= 2;
        Face *tmp = realloc(*faces, sizeof(Face) * (*capacity));
        if(!tmp){
            printf("Failed to reallocate memory for Faces.");
            exit(-1);
        }
        *faces = tmp;
    }
}

int facets_cmpFacesByRSS(const void *f1, const void *f2){
    double diff = ((Face*)f1)->rssAvg - ((Face*)f2)->rssAvg;

    if(diff < 0) return -1;
    else if(diff == 0) return 0;
    else return 1;
}
