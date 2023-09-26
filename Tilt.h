#ifndef TILT_H
#define TILT_H

#define SAMPLE_DEPTH 2
//#define SAMPLE_DEPTH 5
#define Z_VALUES 0
#define BAND_CONTRAST 1
#define NONVALID -99
#define Z_DEFAULT -1

#include "AFM.h"
#include "BandContrast.h"

typedef struct{

    int numRows, numCols;
    double ***data;

    int **assigned;
    double **bestZ;
    double ***tilt;
    int **origCol;  // is along x
    int **origRow;  // is along y

} Sample;

/// @brief Creates a new sample struct that can be tilted.
/// @param afm The afm data to tilt with the sample.
/// @param bc The band contrast data to tilt with the sample.
/// @return The new sample (untilted). 
Sample sample_new(AFMData *afm, BandContrast *bc);

/// @brief Frees memory allocated by the given sample.
/// @param smpl The sample to free.
void sample_free(Sample *smpl);

/// @brief Tilts the given sample.
/// @param smpl The sample to tilt.
/// @param phi The angle by which to tilt the sample.
/// @param zToPxScale Scales the z afm values to pixels during tilting.
void sample_tilt(Sample *smpl, double phi);

/// @brief Creates an AFMData struct from the AFM data stored in data.
/// @param numRows The number of rows in the data. 
/// @param numCols The number of columnc in the data.
/// @param data The data containing the afm data (from a sample struct).
/// @return An AFMData struct.
AFMData sample_getAFM(int numRows, int numCols, double ****data);

/// @brief Creates a BandContrast struct from the band contrast data stored in data.
/// @param numRows The number of rows in the data. 
/// @param numCols The number of columnc in the data.
/// @param data The data containing the band contrast data (from a sample struct).
/// @return A BandContrast struct.
BandContrast sample_getBandContrast(int numRows, int numCols, double ****data);

/// @brief Creates an AFMData struct from the tilted AFM data in the given sample.
/// @param smpl The sample from which to get the AFM data.
/// @return An AFMData struct.
AFMData sample_getTiltedAFMFromSample(Sample *smpl);

/// @brief Creates a BandContrast struct from the tilted band contrast data in the given sample.
/// @param smpl The sample from which to get the band contrast data.
/// @return A BandContrast struct.
BandContrast sample_getTiltedBandContrastFromSample(Sample *smpl);

/// @brief Creates an AFMData struct from the untilted AFM data in the given sample.
/// @param smpl The sample from which to get the AFM data.
/// @return An AFMData struct.
AFMData sample_getAFMFromSample(Sample *smpl);

/// @brief Creates a BandContrast struct from the untilted band contrast data in the given sample.
/// @param smpl The sample from which to get the band contrast data.
/// @return A BandContrast struct.
BandContrast sample_getBandContrastFromSample(Sample *smpl);

#endif // TILT_H
