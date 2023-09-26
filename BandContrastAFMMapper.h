// Ralf: note that the old_row and old_colums layers keept track of from where the mapped BC data come from, so the Euler colors can be accessed through them


#ifndef BANDCONTRASTAFMMAPPER_H
#define BANDCONTRASTAFMMAPPER_H

#include "BandContrast.h"
#include "AFM.h"

#define GREYSCALE_DEFAULT 1
#define OLD_POSITION_DEFAULT -1
#define GREYSCALE_LAYER 0
#define OLD_ROW_LAYER 1
#define OLD_COL_LAYER 2
#define NUMBER_OF_LAYERS_IN_BCAFMM 3

typedef struct {

    int nrow, ncol;
    double ***map;

} BandContrastAFMMapper;

/// @brief Creates a new band contrast AFM mapper.
/// @param nrow The number of rows.
/// @param ncol The number of columns.
/// @return A band contrast AFM mapper with the given dimensions.
BandContrastAFMMapper bandContrastAFMMapper_new(int nrow, int ncol);

/// @brief Frees memory allocated by the given band contrast AFM mapper.
/// @param bcAFMm The band contrast AFM mapper to free.
void bandContrastAFMMapper_free(BandContrastAFMMapper *bcAFMm);

/// @brief Maps the given band contrast to the given AFM according to the given parameters.
/// @param bcMeasured The band contrast to map.
/// @param afmTilted The AFM upon which to map the band contrast.
/// @param a0 y pos.
/// @param a1 x scale.
/// @param a2 x skew.
/// @param a3 x stretch.
/// @param a4 x curve.
/// @param a5 x star?
/// @param b0 y pos.
/// @param b1 y skew.
/// @param b2 y scale.
/// @param b3 y curve.
/// @param b4 y stretch.
/// @param b5 y star?
/// @return The measured band contrast mapped onto the given AFM.
BandContrastAFMMapper bandContrastAFMMapper_map(BandContrast *bcMeasured, AFMData afmTilted, double a0, double a1, double a2, double a3, double a4, double a5, double b0, double b1, double b2, double b3, double b4, double b5);

/// @brief Computes chi squared.
/// @param bcAFMm The band contrast AFM Mapper of which to get chi squared.
/// @param bcTilted The band contrast to which to compare the mapping. 
/// @param mStdDev The standard deviation in the mapped area.
/// @param simStdDev The standard deviation in the band contrast.
/// @return Chi squared of the mapping.
double bandContrastAFMMapper_chiSquared(BandContrastAFMMapper *bcAFMm, BandContrast *bcTilted, double mStdDev, double simStdDev);

#endif // BANDCONTRASTAFMMAPPER_H
