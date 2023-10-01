// with functions from nrutil.c
#include "Amoeba.h"

void runAmoeba(BandContrast *bcMeasured, AFMData afm, BandContrast *bcTilted, BandContrastAFMMapper *bcAFMmOut, double mStdDev, double simStdDev, double *asbs[], int fitLevel) {
    double **p, *y;  // initial gueess as simplex corners
    int nfunk, i, j, ndim;
    int index;

    if(fitLevel < 1) fitLevel = 1;
    if(fitLevel > 6) fitLevel = 6;

    ndim = fitLevel * 2;

    p = constructP(asbs, fitLevel);
    y = (double*)malloc(sizeof(double) * (ndim+1));
    
    printf("fit level: %d\n", fitLevel);
    y[1] = amoeba_chisq(bcMeasured, afm, bcTilted, bcAFMmOut, mStdDev, simStdDev, p[1], asbs, ndim);

    for(i=2; i<= ndim+1; i++) {
        for(j=1; j<= ndim; j++) {
            p[i][j] = p[1][j];
        }
    }

    p[2][1] *= STRETCH;       // always fit a0 and b0
    p[3][2] *= STRETCH;
    
    if (fitLevel > 1) {
        p[4][3] *= STRETCH;   // doing a1 and b2
        p[5][4] *= STRETCH;
    }

    if (fitLevel > 2) {
        p[6][5] += PUSH;   // move a2 and b1 away from zero
        p[7][6] += PUSH;
    }

    if (fitLevel > 4 && fitLevel <6 ) {  // bad news
        printf("invalid level ...\n");
    }


    if (fitLevel == 6) {
        p[8][7] += TWIST;   // move quadratic coefficients away from zero
        p[9][8] += TWIST;
        p[10][9] += TWIST;
        p[11][10] += TWIST;
        p[12][11] += TWIST;
        p[13][12] += TWIST;
    }

    for(i=2; i<= ndim+1; i++) {
        y[i] = amoeba_chisq(bcMeasured, afm, bcTilted, bcAFMmOut, mStdDev, simStdDev, p[i], asbs, ndim);
    }

    int ilo = amoeba(p,y,ndim,FTOL,amoeba_chisq,&nfunk,bcMeasured,afm,bcTilted,bcAFMmOut,mStdDev,simStdDev,asbs);

    // update asbs
    for(i = 1; i <= fitLevel; i++){
        *(asbs[i-1]) = p[ilo][2 * (i - 1) + 1];

        if(i == 1 || i == 6){
            *(asbs[i - 1 + 6]) = p[ilo][2 * i];
        }
        else{
            *(asbs[i - 1 + 6 + ((i & 1) == 1 ? -1 : 1)]) = p[ilo][2 * i];
        }
    }

    // free p
    free(y);
}

double **constructP(double *asbs[], int fitLevel){
    int i;
    int ndim = fitLevel * 2;

    double **p = matrix(1,ndim+1,1,ndim);

    //   <a0 a1 a2 a3 a4 a5 b0 b1 b2 b3 b4 b5>
    //   level 1 fit a0, b0
    //   level 2 fit a0, b0, a1, b2
    //   level 3 fit a0, b0, a1, b2, a2, b1
    //   level 4 fit a0, b0, a1, b2, a2, b1, a3, ...

    for(i = 1; i <= fitLevel; i++){
        // set a_i
        p[1][2*i-1] = *(asbs[i - 1]);
        // set b_i
        if(i == 1 || i == 6){
            p[1][2*i] = *(asbs[6+i-1]);
        }
        else{
            p[1][2*i] = *(asbs[6+i-1 + ((i & 1) == 1 ? -1 : 1)]);
        }
    }

    return p;
}

double amoeba_chisq(BandContrast *bcMeasured, AFMData afm, BandContrast *bcTilted, BandContrastAFMMapper *bcAFMmOut, double mStdDev, double simStdDev, double *simplexCorner, double *asbs[], int ndim){
    double chiSquared = 0.0, a0, a1, a2, a3, a4, a5, b0, b1, b2, b3, b4, b5;
    int row, col, overlappingPoints = 0, fitLevel;
    fitLevel = ndim / 2;

    // Should default as and bs be current as and bs even if not fitting those parameters. YES!
    a0 =                  simplexCorner[1];
    a1 = (fitLevel >= 2 ? simplexCorner[3]  : *(asbs[1]));
    a2 = (fitLevel >= 3 ? simplexCorner[5]  : *(asbs[2]));
    a3 = (fitLevel >= 4 ? simplexCorner[7]  : *(asbs[3]));
    a4 = (fitLevel >= 5 ? simplexCorner[9]  : *(asbs[4]));
    a5 = (fitLevel == 6 ? simplexCorner[11] : *(asbs[5]));
    
    b0 =                  simplexCorner[2];
    b1 = (fitLevel >= 3 ? simplexCorner[6]  : *(asbs[7]));
    b2 = (fitLevel >= 2 ? simplexCorner[4]  : *(asbs[8]));
    b3 = (fitLevel >= 5 ? simplexCorner[10] : *(asbs[9]));
    b4 = (fitLevel >= 4 ? simplexCorner[8]  : *(asbs[10]));
    b5 = (fitLevel == 6 ? simplexCorner[12] : *(asbs[11]));

    *bcAFMmOut = bandContrastAFMMapper_map(bcMeasured, afm, a0, a1, a2, a3, a4, a5, b0, b1, b2, b3, b4, b5); // Scales set to 1
    bandContrast_scaleTo255(&bcAFMmOut->map[GREYSCALE_LAYER], bcAFMmOut->nrow, bcAFMmOut->ncol);

    for(row = 0; row < bcAFMmOut->nrow; row++){
        for(col = 0; col < bcAFMmOut->ncol; col++){
            if(bcAFMmOut->map[GREYSCALE_LAYER][row][col] < GREYSCALE_DEFAULT * 255.0){ // Transparency
                // Chi Squared and difference here?
                chiSquared += (bcAFMmOut->map[GREYSCALE_LAYER][row][col] - bcTilted->greyScale[row][col]) * (bcAFMmOut->map[GREYSCALE_LAYER][row][col] - bcTilted->greyScale[row][col]) / sqrt((mStdDev*mStdDev)*(simStdDev*simStdDev));
                overlappingPoints++;
            }
        }
    }

    if(overlappingPoints != 0) chiSquared /= (double)overlappingPoints;
    else chiSquared = 99999999;
    bandContrastAFMMapper_free(bcAFMmOut);
    return chiSquared;
}

int amoeba(double **p, double *y, int ndim, double ftol,double (*funk)(), int *nfunk, BandContrast *bcMeasured, AFMData afm, BandContrast *bcTilted, BandContrastAFMMapper *bcAFMmOut, double mStdDev, double simStdDev, double *asbs[])
{
	int i,j,ilo,ihi,inhi,mpts=ndim+1;
	double ytry,ysave,sum,rtol,*psum;

    //printf("ftol = %f\n",ftol);
    //printf("ndim = %d\n",ndim);
    //printf("mpts = %d \n", mpts);
    fflush(stdout);
	psum=vector(1,ndim);
	*nfunk=0;

    fflush(stdout);
	GET_PSUM
	for (;;) {
		ilo=1;
		ihi = y[1]>y[2] ? (inhi=2,1) : (inhi=1,2);
		for (i=1;i<=mpts;i++) {
			if (y[i] < y[ilo]) ilo=i;
			if (y[i] > y[ihi]) {
				inhi=ihi;
				ihi=i;
			} else if (y[i] > y[inhi] && i != ihi) inhi=i;
		}
        printSimplex(p, y, ndim);
		rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+TINY);
        printf("chisq = %f, rtol = %f\n", y[ilo], rtol);
        fflush(stdout);
		if (rtol < ftol) return ilo;
		if (*nfunk >= NMAX) {
			printf("nfunk = %d\n",*nfunk);
			nrerror("Too many iterations in AMOEBA\n");
		}
        *nfunk += 2;
		ytry=amotry(p,y,psum,ndim,funk,ihi,nfunk,-ALPHA, bcMeasured, afm, bcTilted, bcAFMmOut, mStdDev, simStdDev, asbs);
		if (ytry <= y[ilo])
			ytry=amotry(p,y,psum,ndim,funk,ihi,nfunk,GAMMA, bcMeasured, afm, bcTilted, bcAFMmOut, mStdDev, simStdDev, asbs);
		else if (ytry >= y[inhi]) {
			ysave=y[ihi];
			ytry=amotry(p,y,psum,ndim,funk,ihi,nfunk,BETA, bcMeasured, afm, bcTilted, bcAFMmOut, mStdDev, simStdDev, asbs);
			if (ytry >= ysave) {
				for (i=1;i<=mpts;i++) {
					if (i != ilo) {
						for (j=1;j<=ndim;j++) {
							psum[j]=0.5*(p[i][j]+p[ilo][j]);
							p[i][j]=psum[j];
						}
						y[i]=(*funk)(bcMeasured, afm, bcTilted, bcAFMmOut, mStdDev, simStdDev, psum, asbs, ndim);
					}
				}
				*nfunk += ndim;
				GET_PSUM
			}
		}
        else --(*nfunk);
	}
	free_vector(psum,1,ndim);
}

double amotry(double **p, double *y, double *psum,int ndim, double (*funk)(), int ihi, int *nfunk,double fac, BandContrast *bcMeasured, AFMData afm, BandContrast *bcTilted, BandContrastAFMMapper *bcAFMmOut, double mStdDev, double simStdDev, double *asbs[])
{
	int j;
	double fac1,fac2,ytry,*ptry,*vector();
	void nrerror(),free_vector();

	ptry=vector(1,ndim);
	fac1=(1.0-fac)/ndim;
	fac2=fac1-fac;
	for (j=1;j<=ndim;j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
	ytry=(*funk)(bcMeasured, afm, bcTilted, bcAFMmOut, mStdDev, simStdDev, ptry, asbs, ndim);
	// ++(*nfunk);
	if (ytry < y[ihi]) {
		y[ihi]=ytry;
		for (j=1;j<=ndim;j++) {
			psum[j] += ptry[j]-p[ihi][j];
			p[ihi][j]=ptry[j];
		}
	}
	free_vector(ptry,1,ndim);
	return ytry;
}

void printSimplex(double **p, double *y, int ndim){
    int i, j;

  // show the optimized simplex
  for(i=1; i<=ndim+1; i++) {
    printf("i=%d: ", i );
    for(j=1; j<= ndim; j++) {
      printf("%g ", p[i][j]);
    }
    printf("chi2 = %f\n", y[i]);
  }
}

double **matrix(int nrl, int nrh, int ncl, int nch)
{
        int i;
        double **m;

        m = (double **) malloc( (unsigned)(nrh - nrl + 1) * sizeof(double *) );
        if (!m)
                nrerror("allocation failure 1 in matrix()");
        m -= nrl;

        for (i = nrl; i <= nrh; i++)
        {
                m[i] = (double *) malloc( (unsigned)(nch - ncl + 1)*sizeof(double) );
                if (!m[i])
                        nrerror("allocation failure 2 in matrix()");
                m[i] -= ncl;
        }
        return m;
}

/* free_vector() := Frees a double vector allocated by vector() */
void free_vector(double *v, int nl, int nh)
{       
        free(v + nl);
}  

/* nrerror() := Numerical Recipes standard error handler */
void nrerror(char *error_text)
{       
        fprintf(stderr, "Numerical Recipes run-time error...\n");
        fprintf(stderr, "%s\n", error_text);
        fprintf(stderr, "...now exiting to system...\n");
        exit(1);        
}               
        
/* *vector() := Allocates a double vector with range [nl..nh] */
double *vector(int nl, int nh)
{
    //printf("Entered vector with nl, nh: %d %d\n", nl, nh);
    fflush(stdout);
    double *v;

    fflush(stdout);
    v = (double *)malloc((nh-nl+1) * sizeof(double));
    fflush(stdout);
    if (!v){
            printf("ERROR!");
            fflush(stdout);
            nrerror("\nallocation failure in vector()");
    }
    fflush(stdout);
    return v-nl;
}
