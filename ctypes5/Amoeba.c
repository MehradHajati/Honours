#include "Amoeba.h"

BandContrast bcMeasuredA;
AFMData afmA;
BandContrast bcTiltedA;
BandContrastAFMMapper bcAFMmOutA;

BandContrast readBandContrast(FILE* file);
BandContrastAFMMapper readBandContrastAFMMapper(FILE* file);
AFMData readAFMData(FILE* file);

void runAmoeba(double mStdDev, double simStdDev, double *asbs[], int fitLevel) {

    FILE *bcmFileA = fopen("bcMeasured.txt", "r");
    FILE *afmFileA = fopen("afm.txt", "r");
    FILE *bcTiltedFileA = fopen("bcTilted.txt", "r");
    FILE *bcAFMmFileA = fopen("bcAFMm.txt", "r");

    bcMeasuredA = readBandContrast(bcmFileA);
    afmA = readAFMData(afmFileA);
    bcTiltedA = readBandContrast(bcTiltedFileA);
    bcAFMmOutA =  readBandContrastAFMMapper(bcAFMmFileA);

    double **p, *y;  // initial gueess as simplex corners
    int nfunk, row, col, ndim;

    // making sure that the input of user is within reason
    if(fitLevel < 1) fitLevel = 1;
    else if(fitLevel > 6) fitLevel = 6;

    // creating the array which will hold the values, i think
    p = constructP(asbs, fitLevel);
    ndim = fitLevel * 2;
    y = (double*)malloc(sizeof(double) * (ndim));
    
    y[0] = amoeba_chisq(mStdDev, simStdDev, p[0], ndim);

    for(row = 1; row < ndim + 1; row++) {
        for(col = 0; col < ndim; col++) {
            p[row][col] = p[0][col];
        }
        p[row][row-1] *= LAMBDA;
        y[row] = amoeba_chisq(mStdDev, simStdDev, p[row], ndim);
    }

    printSimplex(p, y, ndim);

    int ilo = amoeba(p,y,ndim,FTOL,amoeba_chisq,&nfunk, mStdDev, simStdDev, asbs);

    // update asbs
    for(row = 0; row < fitLevel; row++){
        *(asbs[row]) = p[ilo][2 * row];

        if(row == 0 || row == 5){
            *(asbs[row + 6]) = p[ilo][2 * row + 1];
        }
        else{
            *(asbs[row + 6 + ((row & 1) == 1 ? 1 : -1)]) = p[ilo][2 * row + 1];
        }
    }

    for(row = 0; row < ndim + 1; row++){
        free(p[row]);
    }
    free(p);
    free(y);
}

double amoeba_chisq(double mStdDev, int simStdDev, double asbs[], int ndim){
    printf("Printing the good stuff here\n");
    printf("%d\n", mStdDev);
    printf("%d\n", simStdDev);
    printf("%d\n", ndim);
    printf("the first entry is: %d\n", asbs[0]);
    int i = 5;
    printf("i is set to: %d\n", i);
    FILE *bcmFileA = fopen("bcMeasured.txt", "r");
    printf("hello\n");
    if(bcmFileA == NULL){
        printf("file is null\n");
    }
    else{
        printf("file is not null\n");
    }
    fscanf(bcmFileA, "%d", i);
    printf("%d\n", i);
    FILE *afmFileA = fopen("afm.txt", "r");
    FILE *bcTiltedFileA = fopen("bcTilted.txt", "r");
    FILE *bcAFMmFileA = fopen("bcAFMm.txt", "r");

    printf("here\n");
    bcMeasuredA = readBandContrast(bcmFileA);
    printf("here 2\n");
    afmA = readAFMData(afmFileA);
    bcTiltedA = readBandContrast(bcTiltedFileA);
    bcAFMmOutA =  readBandContrastAFMMapper(bcAFMmFileA);

    
    double chiSquared = 0.0, a0, a1, a2, a3, a4, a5, b0, b1, b2, b3, b4, b5;
    int row, col, overlappingPoints = 0, fitLevel;
    fitLevel = ndim / 2;

    // Should default as and bs be current as and bs even if not fitting those parameters. YES!
    a0 = asbs[0];
    a1 = asbs[1];
    a2 = asbs[2];
    a3 = asbs[3];
    a4 = asbs[4];
    a5 = asbs[5];
    
    b0 = asbs[1];
    b1 = asbs[7];
    b2 = asbs[8];
    b3 = asbs[9];
    b4 = asbs[10];
    b5 = asbs[11];

    // here it is creating the mapping to be compared and find the difference
    bcAFMmOutA = bandContrastAFMMapper_map(bcMeasuredA, afmA, a0, a1, a2, a3, a4, a5, b0, b1, b2, b3, b4, b5); // Scales set to 1
    // what to about the scaling here?
    bandContrast_scaleTo255(bcAFMmOutA.map[GREYSCALE_LAYER], bcAFMmOutA.nrow, bcAFMmOutA.ncol);

    // objective function is here
    for(row = 0; row < bcAFMmOutA.nrow; row++){
        for(col = 0; col < bcAFMmOutA.ncol; col++){
            // what to do about the transparency here
            if(bcAFMmOutA.map[GREYSCALE_LAYER][row][col] < GREYSCALE_DEFAULT * 255.0){ // Transparency
                // Chi Squared and difference here?
                chiSquared += (bcAFMmOutA.map[GREYSCALE_LAYER][row][col] - bcTiltedA.greyScale[row][col]) * (bcAFMmOutA.map[GREYSCALE_LAYER][row][col] - bcTiltedA.greyScale[row][col]) / sqrt((mStdDev*mStdDev)*(simStdDev*simStdDev));
                overlappingPoints++;
            }
        }
    }

    if(overlappingPoints != 0) chiSquared /= (double)overlappingPoints;
    else chiSquared = 99999999;
    bandContrastAFMMapper_free(&bcAFMmOutA);
    return chiSquared;

}

int amoeba(double **p, double *y, int ndim, double ftol,double (*funk)(), int *nfunk, double mStdDev, double simStdDev, double *asbs[]){
	int i,j,ilo,ihi,inhi,mpts=ndim+1;
	double ytry,ysave,sum,rtol,*psum;

    //printf("ftol = %f\n",ftol);
    //printf("ndim = %d\n",ndim);
    //printf("mpts = %d \n", mpts);
    fflush(stdout);
	psum = (double*)malloc(sizeof(double) * ndim);
	*nfunk=0;

    // GET_PSUM
    for (j = 0; j < ndim; j++) { 
        sum = 0.0;
        for (i = 0; i < mpts; i++) {
            sum += p[i][j]; 
            psum[j] = sum;
        }
    }
	for (;;) {
		ilo=0;
		ihi = y[0] > y[1] ? (inhi = 1,0) : (inhi = 0,1);
		for (i = 0; i < mpts; i++) {
			if (y[i] < y[ilo]) ilo=i;
			if (y[i] > y[ihi]) {
				inhi=ihi;
				ihi=i;
			} else if (y[i] > y[inhi] && i != ihi){
                inhi = i;
            } 
		}
        printSimplex(p, y, ndim);
		rtol = 2.0 * fabs(y[ihi] - y[ilo]) / (fabs(y[ihi]) + fabs(y[ilo]) + TINY);
        printf("rtol: %f\n", rtol);
		if (rtol < ftol) return ilo;
		if (*nfunk >= NMAX) {
			printf("nfunk = %d\n",*nfunk);
            return ilo;
		}
        *nfunk += 2;
		ytry = amotry(p, y, psum, ndim, funk, ihi, nfunk, -ALPHA, mStdDev, simStdDev, asbs);
		if (ytry <= y[ilo]){
            ytry = amotry(p, y, psum, ndim, funk, ihi, nfunk, GAMMA, mStdDev, simStdDev, asbs);
        }
		else if (ytry >= y[inhi]) {
			ysave = y[ihi];
			ytry = amotry(p, y, psum, ndim, funk, ihi, nfunk, BETA, mStdDev, simStdDev, asbs);
			if (ytry >= ysave) {
				for (i = 0; i < mpts; i++) {
					if (i != ilo) {
						for (j = 0; j < ndim; j++) {
							psum[j] = 0.5 * (p[i][j] + p[ilo][j]);
							p[i][j] = psum[j];
						}
						y[i]=(*funk)(bcMeasuredA, afmA, bcTiltedA, bcAFMmOutA, mStdDev, simStdDev, psum, asbs, ndim);
					}
				}
				*nfunk += ndim;
                // GET_PSUM
				for (j = 0; j < ndim; j++) { 
                    sum = 0.0;
                    for (i = 0; i < mpts; i++) {
                        sum += p[i][j]; 
                        psum[j] = sum;
                    }
                }
			}
		}
        else --(*nfunk);
	}
    free(psum);
}

double amotry(double **p, double *y, double *psum,int ndim, double (*funk)(), int ihi, int *nfunk,double fac, double mStdDev, double simStdDev, double *asbs[]){
	int col;
	double fac1,fac2,ytry,*ptry;

	ptry = (double*)malloc(sizeof(double) * ndim);
	fac1 = (1.0 - fac) / ndim;
	fac2 = fac1 - fac;
	for (col = 0; col < ndim; col++) {
        ptry[col] = psum[col] * fac1 - p[ihi][col] * fac2;
    }
	ytry = (*funk)(bcMeasuredA, afmA, bcTiltedA, bcAFMmOutA, mStdDev, simStdDev, ptry, asbs, ndim);
	if (ytry < y[ihi]) {
		y[ihi]=ytry;
		for (col = 0; col < ndim; col++) {
			psum[col] += ptry[col] - p[ihi][col];
			p[ihi][col] = ptry[col];
		}
	}
	free(ptry);
	return ytry;

}

/**
 * @brief 
 * 
 * @param p 
 * @param y 
 * @param ndim 
 */
void printSimplex(double **p, double *y, int ndim){    
    int row, col;

    // show the optimized simplex
    for(row = 0; row < ndim + 1; row++) {
        printf("row=%d: ", row );
        for(col = 0; col < ndim; col++) {
            printf("%f ", p[row][col]);
        }
        printf("value = %f\n", y[row]);
    }
}

/**
 * @brief not sure what it is doing here 
 * 
 * @param asbs 
 * @param fitLevel 
 * @return double** 
 */
double **constructP(double *asbs[], int fitLevel){
    int row, col;
    int ndim = fitLevel * 2;

    // creating a 2d double array
    double **p = (double**)malloc(sizeof(double) * (ndim + 1));

    // creating 1d arrays in the indices of the 2d array created above
    for(row = 0; row < ndim + 1; row++){
        p[row] = (double*)malloc(sizeof(double) * ndim);
    }

    for(col = 0; col < fitLevel; col++){
        // set a_i
        p[0][2 * col] = *(asbs[col]);
        // set b_i
        if(col == 0 || col == 5){
            p[0][2 * col + 1] = *(asbs[col + 6]);
        }
        else{
            p[0][2 * col + 1] = *(asbs[col + 6 + ((col & 1) == 1 ? 1 : -1)]);
        }
    }

    return p;
}

BandContrast readBandContrast(FILE* file) {

    if(file == NULL){
        printf("file is null\n");
    }
    printf("start\n");
    int nrow, ncol;
    
    fscanf(file, "%d", nrow);
    fscanf(file, "%d", ncol);
    printf("opened files\n");
    BandContrast newBC = bandContrast_new(nrow, ncol);
    for(int i = 0; i < newBC.nrow; i++) {
        for(int j = 0; j < newBC.ncol; j++) {
            fscanf(file, "%lf", &newBC.greyScale[i][j]);
        }
    }

    return newBC;
}

BandContrastAFMMapper readBandContrastAFMMapper(FILE* file) {
    int nrow, ncol;
    
    fscanf(file, "%d", nrow);
    fscanf(file, "%d", ncol);

    BandContrastAFMMapper newMap = bandContrastAFMMapper_new(nrow, ncol);
    for(int i = 0; i < newMap.nrow; i++) {
        for(int j = 0; j < newMap.ncol; j++) {
            fscanf(file, "%lf", &newMap.map[GREYSCALE_LAYER][i][j]);
        }
    }

    return newMap;
}

AFMData readAFMData(FILE* file) {
    int xres, yres;
    
    fscanf(file, "%d", xres);
    fscanf(file, "%d", yres);

    AFMData newafm = afmData_new(xres, yres);
    for(int i = 0; i < newafm.xResolution; i++) {
        for(int j = 0; j < newafm.yResolution; j++) {
            fscanf(file, "%lf", &newafm.zValues[i][j]);
        }
    }

    return newafm;
}

