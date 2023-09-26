#include <stdio.h>
#include "SampleStretcher.h"

void sampleStretcher_stretch(Sample *smpl, double phi){
    int row, col, newRow;
    int min = 1000000;
    int max = -1;
    double **tmpZVals, **tmpBandContrastVals;
    
    tmpZVals = (double **)malloc(sizeof(double *) * smpl->numRows);
    tmpBandContrastVals = (double **)malloc(sizeof(double *) * smpl->numRows);

    // Get min and max row vals
    for(row = 0; row < smpl->numRows; row++){
        for(col = 0; col < smpl->numCols; col++){
	    if(smpl->tilt[row][col][Z_VALUES] != Z_DEFAULT){
		min = ((min < row) ? min : row);
		max = ((max > row) ? max : row);
	    }
	}
     }

     printf("stretcher %d %d %g\n", min, max, phi);
     // Build stretched array
     for(row = 0; row < smpl->numRows; row++){
         tmpZVals[row]            = (double *)calloc(smpl->numCols, sizeof(double));
         tmpBandContrastVals[row] = (double *)calloc(smpl->numCols, sizeof(double));
         // this makes a copy of the data in the tmp arrays
	 for(col = 0; col < smpl->numCols; col++){
             newRow = min + ((double)row / (double)smpl->numRows) * (max - min);
             //do nothing: looks kind of right
             //smpl->origRow[newRow][col] = row;  // looks wrong, seems to undo move
             //smpl->origRow[row][col] = newRow;  // this gives a narrow strip
	     tmpZVals[row][col]            = smpl->tilt[newRow][col][Z_VALUES];
	     tmpBandContrastVals[row][col] = smpl->tilt[newRow][col][BAND_CONTRAST];
	 }
     }

     // Here the origianl arrays get overwritten with the calculations in the tmp arrays
     for(row = 0; row < smpl->numRows; row++){
	 for(col = 0; col < smpl->numCols; col++){
             smpl->tilt[row][col][Z_VALUES]      = tmpZVals[row][col];	
             smpl->tilt[row][col][BAND_CONTRAST] = tmpBandContrastVals[row][col];		
 	 }
     }

     // Free memory
     for(row = 0; row < smpl->numRows; row++){
	 free(tmpZVals[row]);
         free(tmpBandContrastVals[row]);
     }
     free(tmpZVals);
     free(tmpBandContrastVals);
}
