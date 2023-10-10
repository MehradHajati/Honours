#define true 1
#define false 0
#define bool int
#define MAXPATH 260
#define NUMBCS 3

#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "AFMDataReader.h"
#include "AFMDataWriter.h"
#include "Facets.h"
#include "FacetsWriter.h"
#include "FacetsReader.h"
#include "BandContrast.h"
#include "BandContrastWriter.h"
#include "BandContrastReader.h"
#include "BandContrastAFMMapper.h"
#include "Tilt.h"
#include "SampleStretcher.h"
#include "BitmapWriter.h"
#include "Amoeba.h"

#if (defined(_WIN32) || defined(__WIN32__))
#define mkdir(A, B) mkdir(A)
#endif

const char *rootDir = "/Facets/samples/";
const int requiredArgc = 7;

void runBandContrastSim(char *fileName, char *simOutputDir, char *measuredFN, char *thetaFN, char *phiFN, char *sampleDir, double phi, double lightDeg, double width, double alpha, bool readFacets);

int main(int argc, char *argv[]){

    if(argc != requiredArgc){
        printf("%d arguments are required:\nSampleName, tiltAngle, lightAngle, widthInUM, alpha (light fraction for linCombo [0-1]), readFacets (0 to calculate facets, 1 to read from files)\n", requiredArgc - 1);
        return 0;
    }

    time_t startTime = time(0);
    // Prep cmd line arguments
    bool readFacets = atoi(argv[6]);
    double alpha;
    alpha = atof(argv[5]);
    // Prep paths
    char *sampleDir, *simOutputDir, *afmFile, *measuredBCpath, *thetaFN, *phiFN;
    sampleDir = (char *)malloc(sizeof(char) * MAXPATH);
    // setting up the location of the sample directories here for easy reading and saving
    strcpy(sampleDir, "/Users/mhaja");
    strcat(sampleDir, rootDir);
    strcat(sampleDir, argv[1]);
    strcat(sampleDir, "/");
    
    // Make the simOutput directory if it does not exist.
    simOutputDir = (char *)malloc(sizeof(char) * MAXPATH);
    strcpy(simOutputDir, sampleDir);
    strcat(simOutputDir, "simOutput");
    int mkdirResult = mkdir(simOutputDir, 0755);

    // afm
    afmFile = (char *)malloc(sizeof(char) * MAXPATH);
    strcpy(afmFile, sampleDir);
    strcat(afmFile, "Input/");
    strcat(afmFile, argv[1]);
    strcat(afmFile, ".txt");
    // Measured BC
    measuredBCpath = (char *)malloc(sizeof(char) * MAXPATH);
    strcpy(measuredBCpath, sampleDir);
    strcat(measuredBCpath, "Output/");
    // Theta
    thetaFN = (char *)malloc(sizeof(char) * MAXPATH);
    strcpy(thetaFN, sampleDir);
    strcat(thetaFN, "Output/thetamap_0.001.txt");
    // Phi
    phiFN = (char *)malloc(sizeof(char) * MAXPATH);
    strcpy(phiFN, sampleDir);
    strcat(phiFN, "Output/phimap_0.001.txt");
    
    // this method is the magic method and everything happens in it
    runBandContrastSim(afmFile, simOutputDir, measuredBCpath, thetaFN, phiFN, sampleDir, atof(argv[2]), atof(argv[3]), atof(argv[4]), alpha, readFacets);

    // after the run we are freeing the used up memory from the method
    free(sampleDir);
    free(simOutputDir);
    free(afmFile);
    free(measuredBCpath);
    free(thetaFN);
    free(phiFN);

    // printing out the time that the program was in use
    time_t endTime = time(0);
    printf("Time to complete: %ld seconds\n", endTime - startTime);
}

void runBandContrastSim(char *fileName, char *simOutputDir, char *measuredFN, char *thetaFN, char *phiFN, char *sampleDir, double phi, double lightDeg, double width, double alpha, bool readFacets){

    int row, col, i, j, iLight = 0, iVis = 1, iLinComboA = 2; 
    //double minZ = 10000000, maxZ = -1000000, measuredOpacity = 0.75, a0, a1, a2, a3, a4, a5, b0, b1, b2, b3, b4, b5;
    double minZ = 10000000, maxZ = -1000000, measuredOpacity = 0.50, a0, a1, a2, a3, a4, a5, b0, b1, b2, b3, b4, b5;
    FILE *binned;
    char *binnedName;
    int binnedThere;
    AFMData afm;
    int facetTypeFlag;

    // Read or compute the binned AFM data
    // allocating the memory 
    binnedName = (char *)malloc(sizeof(char) * MAXPATH);
    strcpy(binnedName, fileName);
    strcat(binnedName, "Binned");
    printf("Looking for binned data at %s.\n", binnedName);
    if( NULL == (binned = fopen(binnedName, "r"))) {
        binnedThere = 0;
        printf("No binned data found.\n");
        afm = afmData_readFromFile(fileName);
        // flipping the afm data across the middle row
        afmData_flipY(&afm);

        // Get min and max Z
        for(row = 0; row < afm.xResolution; row++){
            for(col = 0; col < afm.yResolution; col++){
                minZ = minZ < afm.zValues[row][col] ? minZ : afm.zValues[row][col];
                maxZ = maxZ > afm.zValues[row][col] ? maxZ : afm.zValues[row][col];
            }
        }
        // Shift and scale Z values to pixels
        for(row = 0; row < afm.xResolution; row++){
            for(col = 0; col < afm.yResolution; col++){
                afm.zValues[row][col] -= minZ;
                afm.zValues[row][col] *= afm.xResolution / (double)width;
            }
        }
        printf("MaxZ: %f\n", maxZ);
        printf("MinZ: %f\n", minZ);
        printf("MaxZ - MinZ: %f\n", maxZ - minZ);
        printf("AFM data size is %d\n", afm.xResolution);
        // In this while loop we are binning the afm data by a factor of 2 until its resolution is 1280 by 1280
        while(afm.xResolution > 1280) {
            afmData_binBy2(&afm);
            printf("AFM data size is %d\n", afm.xResolution);
        }
        printf("AFM: [%d,%d]\n", afm.xResolution, afm.yResolution);
        // now write the result to use it next time
        fclose(binned);
        binned = fopen(binnedName, "w");
        // for loops to write the binned data to files
        for(row = 0; row < afm.xResolution; row++) {
            fprintf(binned, "%.5lf", afm.zValues[row][0]);
            for(col = 1; col < afm.yResolution; col++) {
                fprintf(binned, "\t%.5lf", afm.zValues[row][col]);
            }
            fprintf(binned,"\n");
        }
        fclose(binned);
    }
    else {
        binnedThere = 1;
        printf("binned is there\n");
        afm = afmData_readFromFile(binnedName);
    }
    fclose(binned);

    // Read or compute the facets
    Facets facets;
    if(readFacets) facets = facets_readFromFiles(thetaFN, phiFN);
    // the numbers in the method below are given to us by Prof Bruening, which he found after many summers of testing and playing around with them
    // do not change them
    else facets = facets_compute(&afm, 3, 21, 0.001, 0);
    printf("Facets dimensions: %dx%d\n", facets.nrow, facets.ncol);

    // Read measured band contrast
    BandContrast bcMeasured;
    bcMeasured = bandContrast_readFromFile(measuredFN);
    bandContrast_scaleTo255(&bcMeasured.greyScale, bcMeasured.nrow, bcMeasured.ncol);
    double mAvg, mStdDev, mMin = DBL_MAX, mMax = DBL_MIN;
    mAvg = bandContrast_averageGrey(&bcMeasured, 0, bcMeasured.nrow, 0, bcMeasured.ncol);
    mStdDev = bandContrast_stdDev(&bcMeasured, 0, bcMeasured.nrow, 0, bcMeasured.ncol);
    // here we are finding the lowest and highest values in the greyscale picture of the bandcontrast
    for(row = 0; row < bcMeasured.nrow; row++){
        for(col = 0; col < bcMeasured.ncol; col++){
            // if the current pixel in the greyscale is lower than the min then replace it, if not keep mMin as mMin
            mMin = (mMin < bcMeasured.greyScale[row][col]) ? mMin : bcMeasured.greyScale[row][col];
            // if the current pixel in the greyscale is higher than the max then replace it, if not keep mMax as mMax
            mMax = (mMax < bcMeasured.greyScale[row][col]) ? bcMeasured.greyScale[row][col]: mMax;
        }
    }
    printf("Measured band contrast statistics: min, mean, max, stddev: %d, %f, %d, %f\n", (int)mMin, mAvg, (int)mMax, mStdDev);
    // divding all the values of the greyscale picture in the bandcontrast objest by 255
    for(row = 0; row < bcMeasured.nrow; row++){
        for(col = 0; col < bcMeasured.ncol; col++){
            bcMeasured.greyScale[row][col] = bcMeasured.greyScale[row][col] / 255.0;
        }
    }

    // Simulate the band constrast (called "Sample")
    printf("Starting simulation\n");
    BandContrastAFMMapper bcAFMm;
    double rad = (180.0 - lightDeg) * M_PI / 180.0;
    printf("the value of M_pi is: %d\n", M_PI);
    Vector3 light; 
    light.x = 0; light.y = cos(rad); light.z = -sin(rad);
    light = vector3_scale(light, 1 / vector3_magnitude(light));
    // this is a function that does the heavy lifting
    BandContrast *bcs = simulateBandContrast(&facets, &afm, light, alpha, &bcMeasured);  
    printf("Simulated BC: %dx%d\n", bcs[iVis].ncol, bcs[iVis].nrow);
    // creating the memory and places and freeing the memory for the three pictures of the simulated band contrast
    // the first pis is the light picture which shows you what is illuminated based on the angle of the light
    // the second pictures is the detector view which used monte carlo and tells you if the detector was a globe what would it see.
    // the third picture is a linear combination of the two of these pictures
    Sample smpls[NUMBCS] = {sample_new(&afm, &bcs[iLight]), sample_new(&afm, &bcs[iVis]), sample_new(&afm, &bcs[iLinComboA])};
    for(i = 0; i < NUMBCS; i++){
        bandContrast_free(&bcs[i]);
    }
    free(bcs);
    char *bcTiltFNs[NUMBCS] = {
        "/light.txt",
        "/detector.txt",
        "/BCsimulated.txt"
    };
    char *bcTiltPGMFNs[NUMBCS] = {
        "/light.pgm",
        "/detector.pgm",
        "/BCsimulated.pgm"
    };
    int toTilt[NUMBCS] = {
        iLight,
        iVis,
        iLinComboA
    };

    // Tilt the simulated band contrast
    BandContrast bcTilted, bcDifference, bcSolidOverlap, bcTransOverlap;    
    double simAvg, simStdDev, simMin = DBL_MAX, simMax = DBL_MIN, val;
    // this for loop is used to tilt, then stretch and scale all three pictures created above
    for(i = 0; i < NUMBCS; i++){  // NUMBCS = 3 right now
        printf("Tilting sample %d/%d...\n", i+1, NUMBCS);
        // tilting each of the pictures here
        sample_tilt(&smpls[toTilt[i]], phi);
        printf("Stretching...\n");
        // stretching the simulated sample to fill the same area as the original data
        sampleStretcher_stretch(&smpls[toTilt[i]],phi);
        // creates a band contrast construct for the simulated band contrast
        bcTilted = sample_getTiltedBandContrastFromSample(&smpls[toTilt[i]]);
        // scaling the simulated data to 255
        bandContrast_scaleTo255(&bcTilted.greyScale, bcTilted.nrow, bcTilted.ncol);
        if(i != iLight) bandContrast_fillGaps(&bcTilted);
        // calculating the average grey and standrad deviation of the simulated band contrast after it got tilted
        simAvg = bandContrast_averageGrey(&bcTilted, 300, bcTilted.nrow - 300, 150, bcTilted.ncol - 150);
        simStdDev = bandContrast_stdDev(&bcTilted, 300, bcTilted.nrow - 300, 150, bcTilted.ncol - 150);
        for(row = 300; row < bcTilted.nrow - 300; row++){
            for(col = 150; col < bcTilted.ncol - 150; col++){
                // keeping track of the min and max and replacing them if a new one is found
                simMin = (simMin < bcTilted.greyScale[row][col]) ? simMin : bcTilted.greyScale[row][col];
                simMax = (simMax < bcTilted.greyScale[row][col]) ? bcTilted.greyScale[row][col] : simMax;
            }
        }
        printf("%s BEFORE: min, mean, max, stddev: %d, %f, %d, %f\n",
            bcTiltFNs[toTilt[i]], (int)simMin, simAvg, (int)simMax, simStdDev);

        char *bcFNpgm = (char *)malloc(sizeof(char) * MAXPATH);
        strcpy(bcFNpgm, simOutputDir);
        strcat(bcFNpgm, bcTiltPGMFNs[toTilt[i]]);
        // here we write     ../light.pgm, ../detector.pgm" and ../BCsimulated.pgm" just before deleting it all ...
        bitmapWriter_writePGM(bcFNpgm, bcTilted.greyScale, 0, bcTilted.nrow, 0, bcTilted.ncol, true);
        free(bcFNpgm);
        bandContrast_free(&bcTilted);
    }

    char *input = (char *)malloc(sizeof(char) * MAXPATH);
    char *tok;
    int tokCount, fitLevel = 1;
    bool validInput;
    input[0] = '\0';

    bcTilted = sample_getTiltedBandContrastFromSample(&smpls[iLinComboA]);
    
    bandContrast_scaleTo255(&bcTilted.greyScale, bcTilted.nrow, bcTilted.ncol);
    bandContrast_fillGaps(&bcTilted);
    // Simulated
    simAvg = bandContrast_averageGrey(&bcTilted, 300, bcTilted.nrow - 300, 150, bcTilted.ncol - 150);
    simStdDev = bandContrast_stdDev(&bcTilted, 300, bcTilted.nrow - 300, 150, bcTilted.ncol - 150);
    for(row = 300; row < bcTilted.nrow - 300; row++){
        for(col = 150; col < bcTilted.ncol - 150; col++){
            simMin = (simMin < bcTilted.greyScale[row][col]) ? simMin : bcTilted.greyScale[row][col];
            simMax = (simMax < bcTilted.greyScale[row][col]) ? bcTilted.greyScale[row][col] : simMax;
        }
    }

    a0 = afm.xResolution / 2;
    b0 = afm.yResolution / 2;
    a1 = 0.8*afm.xResolution/bcMeasured.nrow;
    b2 = 1.2*afm.yResolution * cos(phi * M_PI / 180.0) / bcMeasured.ncol;
    a2 = 1.0 / a0;
    b1 = 1.0 / b0;
    a3 = 1.0 / (a0*a0);
    b4 = 1.0 / (b0*b0);
    a4 = 1.0 / (a0*a0);
    b3 = 1.0 / (b0*b0);
    a5 = 1.0 / (a0*b0);
    b5 = 1.0 / (a0*b0);

    double *asbs[12] = {
        &a0, &a1, &a2, &a3, &a4, &a5, &b0, &b1, &b2, &b3, &b4, &b5
    };

    char *coeffNChiFN = (char *)malloc(sizeof(char) * MAXPATH);
    strcpy(coeffNChiFN, simOutputDir);
    strcat(coeffNChiFN, "/coeffNChi.txt");
    bool  firstInteraction = true;

    while(strcmp(input, "q") != 0){
        printf("\n------------------------------------------\n");
        // provide input options
        if(firstInteraction){
            printf("\nInput new coefficients in format: a0 a1 a2 a3 a4 a5 b0 b1 b2 b3 b4 b5\n");
            printf("Entering \'.\' for a coefficient will use previous value.\n");
            printf("Or type \'a #\', with # being fitLevel in {1, 2, 3, 6}\' to run amoeba from the current paramters.\n");
            printf("\'f\' to calculate the facets types.\n");
            printf("Or type \'q\' to quit.\n");
            firstInteraction = !firstInteraction;
        }
        // print current values
        printf("\nCurrent coefficients:\n<a0 a1 a2 a3 a4 a5 b0 b1 b2 b3 b4 b5>\n<%.4f %.4f %.4f %g %g %g %.4f %.4f %.4f %g %g %g>\n",
            a0, a1, a2, a3, a4, a5, b0, b1, b2, b3, b4, b5);
        printf("\n> ");
        // validate input
        fgets(input, MAXPATH, stdin); // includes newline
        input[strlen(input)-1] = '\0'; // trim newline
        if(strcmp(input, "q") == 0) continue;
        facetTypeFlag = 0;
        if(strcmp(input, "f") == 0) facetTypeFlag = 1;
        if(input[0] == 'a'){
            tok = strtok(input, " ");
            tok = strtok(NULL, " ");
            if(tok != NULL){
                fitLevel = atoi(tok);
            }
            runAmoeba(&bcMeasured, afm, &bcTilted, &bcAFMm, mStdDev, simStdDev, asbs, fitLevel);
        }
        else if(facetTypeFlag == 0) {
            // parse input
            tok = strtok(input, " ");
            tokCount = 0;
            while(tok != NULL){                
                if(strcmp(tok, ".") != 0){
                    *asbs[tokCount] = atof(tok);
                }
                tokCount++;
                tok = strtok(NULL, " ");
            }
        }
        // run new mapping!
        bcAFMm = bandContrastAFMMapper_map(&bcMeasured, afm, a0, a1, a2, a3, a4, a5, b0, b1, b2, b3, b4, b5);
        bandContrast_scaleTo255(&bcAFMm.map[GREYSCALE_LAYER], bcAFMm.nrow, bcAFMm.ncol);

        char * bcAFMmFN = (char *)malloc(sizeof(char) * MAXPATH);
        strcpy(bcAFMmFN, simOutputDir);
        strcat(bcAFMmFN, "/measured_mapped.pgm");
        bitmapWriter_writePGM(bcAFMmFN, bcAFMm.map[GREYSCALE_LAYER], 0, bcAFMm.nrow, 0, bcAFMm.ncol, true);
        free(bcAFMmFN);

        // option to produce a check image
        // char *bcFNpgm = (char *)malloc(sizeof(char) * MAXPATH);
        // strcpy(bcFNpgm, simOutputDir);
        // strcat(bcFNpgm, "/BCsimulated1.pgm");
        // bitmapWriter_writePGM(bcFNpgm, bcTilted.greyScale, 0, bcTilted.nrow, 0, bcTilted.ncol, true);
        // free(bcFNpgm);

        // Overlap
        bcSolidOverlap  = bandContrast_new(bcTilted.nrow, bcTilted.ncol);
        bcTransOverlap  = bandContrast_new(bcTilted.nrow, bcTilted.ncol);
        bcDifference    = bandContrast_new(bcTilted.nrow, bcTilted.ncol);

        double chiSquared = bandContrastAFMMapper_chiSquared(&bcAFMm, &bcTilted, mStdDev, simStdDev);

        for(row = 0; row < bcAFMm.nrow; row++){
            for(col = 0; col < bcAFMm.ncol; col++){
                //printf("%d ", (int)bcAFMm.map[GREYSCALE_LAYER][row][col]);
                if(bcAFMm.map[GREYSCALE_LAYER][row][col] >= GREYSCALE_DEFAULT * 255.0){  // is white, not mapped, show background

                    bcSolidOverlap.greyScale[row][col]    = bcTilted.greyScale[row][col];
                    bcTransOverlap.greyScale[row][col]    = bcTilted.greyScale[row][col];
                    bcDifference.greyScale[row][col]      = 127.5;  //255.0/2.0
                }
                else{ // mapped, show data
                    bcSolidOverlap.greyScale[row][col] = bcAFMm.map[GREYSCALE_LAYER][row][col];
                    bcTransOverlap.greyScale[row][col] = (measuredOpacity)*bcAFMm.map[GREYSCALE_LAYER][row][col] + (1.0-measuredOpacity)*bcTilted.greyScale[row][col];
                    bcDifference.greyScale[row][col]   = 0.5*((bcTilted.greyScale[row][col] - bcAFMm.map[GREYSCALE_LAYER][row][col]) + 255.0);
                }
            }
        }

        // Append coeffs and chi squared to file
        FILE *coeffNChi = fopen(coeffNChiFN, "a");
        for(i = 0; i < 12; i++){
            fprintf(coeffNChi, "%.4f ", *(asbs[i]));
        }
        fprintf(coeffNChi, "%.4f\n", chiSquared);
        fclose(coeffNChi);

        char *bcDiffFN = (char *)malloc(sizeof(char) * MAXPATH);
        strcpy(bcDiffFN, simOutputDir);
        strcat(bcDiffFN, "/diff.pgm");
        bitmapWriter_writePGM(bcDiffFN, bcDifference.greyScale, 0, bcDifference.nrow, 0, bcDifference.ncol, true);
        free(bcDiffFN);

        // measured solid on top of simulated
        char *bcSolidOverlapFN = (char *)malloc(sizeof(char) * MAXPATH);
        strcpy(bcSolidOverlapFN, simOutputDir);
        strcat(bcSolidOverlapFN, "/BCsolid.pgm");
        bitmapWriter_writePGM(bcSolidOverlapFN, bcSolidOverlap.greyScale, 0, bcSolidOverlap.nrow, 0, bcSolidOverlap.ncol, true);
        free(bcSolidOverlapFN);

        // measured with measuredOpacity on top of simulated
        char * bcAFMOverlapFN = (char *)malloc(sizeof(char) * MAXPATH);
        strcpy(bcAFMOverlapFN, simOutputDir);
        strcat(bcAFMOverlapFN, "/BCfitted.pgm");
        bitmapWriter_writePGM(bcAFMOverlapFN, bcTransOverlap.greyScale, 0, bcAFMm.nrow, 0, bcAFMm.ncol, true);
        free(bcAFMOverlapFN);

        bandContrast_free(&bcDifference);
        
        if(facetTypeFlag) {

            // calculate Euler red on top of simulated and plot it
            BandContrast bcRedOverlap;    
            bcRedOverlap = bandContrast_new(bcTilted.nrow, bcTilted.ncol);
            printf("size tilted %d %d\n", bcTilted.nrow, bcTilted.ncol);
            printf("size bcAFM %d %d\n", bcAFMm.nrow, bcAFMm.ncol);
            for(row = 0; row < bcAFMm.nrow; row++) {
                for(col = 0; col < bcAFMm.ncol; col++) {
                    if(bcAFMm.map[GREYSCALE_LAYER][row][col] >= GREYSCALE_DEFAULT * 255.0) {  //white, unassigned, show background
                        bcRedOverlap.greyScale[row][col] = bcTilted.greyScale[row][col];
                    }
                    else {  // overlay the Euler data
                        bcRedOverlap.greyScale[row][col] =
                            bcMeasured.EBSDred[ (int)bcAFMm.map[OLD_ROW_LAYER][row][col] ][ (int)bcAFMm.map[OLD_COL_LAYER][row][col] ];
                    }
                }
            }
            // plot red on top of simulated
            char *bcRedOverlapFN = (char *)malloc(sizeof(char) * MAXPATH);
            strcpy(bcRedOverlapFN, simOutputDir);
            strcat(bcRedOverlapFN, "/BCred.pgm");
            bitmapWriter_writePGM(bcRedOverlapFN, bcRedOverlap.greyScale, 0, bcRedOverlap.nrow, 0, bcRedOverlap.ncol, true);
            free(bcRedOverlapFN);
            //bandContrast_free(&bcRedOverlap);

            // plot EBSD assigned on top of original AFM grid
            AFMData afm_copy = afmData_new(afm.xResolution, afm.yResolution);
            //afmData_copy(&afm, &afm_copy);
            for(row = 0; row < afm_copy.yResolution; row++) {
                for(col = 0; col < afm_copy.xResolution; col++) {
                    afm_copy.zValues[row][col] = 0.0;
                }
            }
            int orow, ocol, row1, col1, countwhite = 0;
            for(row = 0; row < afm_copy.yResolution; row++) {
                for(col = 0; col < afm_copy.xResolution; col++) {

                    orow = smpls[iLinComboA].origRow[row][col];
                    ocol = smpls[iLinComboA].origCol[row][col];

                    row1 = (int)bcAFMm.map[OLD_ROW_LAYER][row][col];
                    col1 = (int)bcAFMm.map[OLD_COL_LAYER][row][col];
                    //printf("%d %d - %d %d\n", row, col, orow, ocol);
                    if(afm_copy.zValues[orow][ocol] != 255.0) countwhite ++;
                    afm_copy.zValues[orow][ocol] = 255.0;

                    //afm_copy.zValues[orow][ocol] = bcRedOverlap.greyScale[row][col]; //far off
                    //   afm_copy.zValues[orow][ocol] = bcRedOverlap.greyScale[row][col]; //far off
                    if(row1 > -1 && col1 > -1) {
                        afm_copy.zValues[orow][ocol] = bcMeasured.EBSDred[row1][col1];  //far off
                    //  //afm_copy.zValues[orow][ocol] = bcRedOverlap.greyScale[row1][col1]; //does not work at all
                    }

                }
            }
            printf("Countwhite is %d\n", countwhite);
            char *assignedMapFN = (char *)malloc(sizeof(char) * MAXPATH);
            strcpy(assignedMapFN, simOutputDir);
            strcat(assignedMapFN, "/mapAssigned.pgm");
            bitmapWriter_writePGM(assignedMapFN, afm_copy.zValues, 0, afm.yResolution, 0, afm.yResolution, true);
            printf("hello\n");
            free(assignedMapFN);


//if(row == col) printf("row col %d %d: %d %d\n", row, col, smpls[iLinComboA].origRow[row][col], smpls[iLinComboA].origCol[row][col]);
            // plot red on top of simulated
            BandContrast eulerOnTheta, thetaCheck;    
            eulerOnTheta = bandContrast_new(facets.nrow, facets.ncol);
            thetaCheck = bandContrast_new(facets.nrow, facets.ncol);
            for(row = 0; row < bcAFMm.nrow; row++) {
                for(col = 0; col < bcAFMm.ncol; col++) {
                    if(facets.thetaMap[row][col] > THETA_DEFAULT)
                        eulerOnTheta.greyScale[row][col] = 255*facets.thetaMap[row][col]/90.0;
                    else
                        eulerOnTheta.greyScale[row][col] = 0.0; // unassigned, invalid
                    //if(row == col) printf("%d %f %f\n",row, facets.thetaMap[row][col], eulerOnTheta.greyScale[row][col]);
                }
            }
            char *eulerOnThetaFN = (char *)malloc(sizeof(char) * MAXPATH);
            strcpy(eulerOnThetaFN, simOutputDir);
            strcat(eulerOnThetaFN, "/eulerOnTheta.pgm");
            bitmapWriter_writePGM(eulerOnThetaFN, eulerOnTheta.greyScale, 0, facets.nrow, 0, facets.ncol, true);
            free(eulerOnThetaFN);
            bandContrast_free(&eulerOnTheta);
            afmData_free(&afm_copy);
            printf("still alive\n");
        } // end if(facetTypeFlag)
    }
    bandContrast_free(&bcSolidOverlap);
    free(input);
    free(coeffNChiFN);
    bandContrastAFMMapper_free(&bcAFMm);

    if(!readFacets) facets_write(facets, thetaFN, phiFN);

    afmData_free(&afm);
    bandContrast_free(&bcMeasured);
    bandContrast_free(&bcTilted);
    for(i = 0; i < NUMBCS; i++){
        sample_free(&smpls[i]);
    }
}
