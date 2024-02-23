#ifndef BANDCONTRAST_H
#define BANDCONTRAST_H

#include "Facets.h"
#include "Quaternion.h"
#include "Vector3.h"

typedef struct {

    int nrow, ncol;
    double **greyScale;
    int **EBSDred;
    int **EBSDgreen;
    int **EBSDblue;

} BandContrast;

/// @brief Creates a new band contrast with the given dimensions.
/// @param nrow The number of rows. 
/// @param ncol The number of columns.
/// @return A band contrast struct.
BandContrast bandContrast_new(int nrow, int ncol);

/// @brief Frees memory allocate by the given band contrast struct.
/// @param bc The band contrast to free. 
void bandContrast_free(BandContrast *bc);

/// @brief Creates a band contrast based on the given light source.
/// @param facets The facets of the sample.
/// @param light Vector in the direction of the illumination.
/// @return A light-based band contrast.
BandContrast bandContrast_light(Facets *facets, Vector3 light);

/// @brief A utility used by bandContrast_light.
BandContrast bandContrast_fromPhiThetaMaps(double **phimap, double **thetamap, int nrow, int ncol, Vector3 light);

/// @brief Create a band contrast based on a detector view. Monte Carlo approach.
/// @param facets The facets of the sample.
/// @param afm The height information of the sample.
/// @return A detector-based band contrast.
BandContrast bandContrast_detector(Facets *facets, AFMData *afm);

/// @brief Creates a list of band contrast structs (light, detector, and a linear combination of the two) scaled to match the average greyscale of the measured band contrast.
/// @param facets The facets of the smaple.
/// @param afm The height information of the sample.
/// @param light Vector in the direction of the illumination.
/// @param alpha The fraction of the linear combination given to the light-based band contrast, rest is detector-based
/// @param measured The measured band contrast.
/// @return A list of band contrasts.
BandContrast *simulateBandContrast(Facets *facets, AFMData *afm, Vector3 light, double alpha, BandContrast *measured);

/// @brief Does tilting for the detector-based band contrast.
/// @param bc The detector-based band contrast.
/// @param afm The height information of the sample.
/// @param qTilt The quaternion by which to tilt the sample.
/// @param avgZ The average height of the sample (used for tilting about COM).
void bandContrast_tiltForDetector(BandContrast *bc, AFMData *afm, Quaternion qTilt, double avgZ);

/// @brief Ensures the minimum of the height values rests at 0.
/// @param afm The AFM data to prep.
/// @param avgZ Stores the average height in this variable.
/// @param afmCOM Stores the COM of the sample in this variable.
void bandContrast_prepAFMForDetector(AFMData *afm, double *avgZ, Vector3 *afmCOM);

/// @brief Creates a rotation quaternion from the given theta and phi angles.
/// @param theta Radians down from detector normal. 
/// @param phi Radians around detector normal.
/// @return A rotation quaternion.
Quaternion bandContrast_tiltFromThetaPhi(double theta, double phi);

/// @brief Computes the average grey value in the given band contrast within the given bounds.
/// @param bc The band contrast from which to get the average grey.
/// @param startRow The starting row.
/// @param endRow The ending row.
/// @param startCol The starting column.
/// @param endCol The ending column.
/// @return The average grey.
double bandContrast_averageGrey(BandContrast *bc, int startRow, int endRow, int startCol, int endCol);

/// @brief Computes the standard deviation in the given band contrast within the given bounds.
/// @param bc The band contrast from which to get the standard deviation.
/// @param startRow The starting row.
/// @param endRow The ending row.
/// @param startCol The starting column.
/// @param endCol The ending column.
/// @return The standard deviation.
double bandContrast_stdDev(BandContrast *bc, int startRow, int endRow, int startCol, int endCol);

/// @brief Multiplies all values by 255.
/// @param vals The values to scale.
/// @param nrow The number of rows to scale.
/// @param ncol The number of columns to scale.
void bandContrast_scaleTo255(double **vals, int nrow, int ncol);

/// @brief Flood-fills 0 values in the given band contrast.
/// @param bc The band contraat to fill.
void bandContrast_fillGaps(BandContrast *bc);

#endif // BANDCONTRAST_H
