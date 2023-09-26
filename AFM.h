#ifndef PROBE_H
#define PROBE_H
/// The number of atoms the probe must contact to consider the z-position being on the
/// surface. Must be non-negative.
#define ATOM_COUNT 3
#define RANDOM_PROBES_PER_BIN 10
/// The precision of the probe.
#define TOLERANCE 0.000001

#include "Utility.h"

typedef struct
{
    int xResolution;
    int yResolution;
    double** zValues;
} AFMData;

/// @brief Allocates memory for AFMData with the given resolution.
/// @param xResolution The number of data points in the x-direction the probe will 
/// store.
/// @param yResolution The number of data points in the y-direction the probe will 
/// store. 
/// @return AFMData with the given resolution.
AFMData afmData_new(int xResolution, int yResolution);

/// @brief Frees memory allocated within the AFMData.
/// @param afmData The AFMData to free the memory of. 
void afmData_free(AFMData* afmData);

/// @brief Bins the given AFM data by a factor of 2.
/// @param afm The AFMData to bin.
void afmData_binBy2(AFMData *afm);

/// @brief Flips the given AFM data across the middle row.
/// @param afm The AFM data to flip.
void afmData_flipY(AFMData *afm);

/// @brief Copies AFM data.
/// @param afm The AFM data to copy.
void afmData_copy(AFMData *afm, AFMData *afmCopy);

#endif // PROBE_H
