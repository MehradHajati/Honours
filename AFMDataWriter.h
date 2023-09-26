#ifndef AFMDATAWRITER_H
#define AFMDATAWRITER_H

#include "AFM.h"

/// @brief Writes the afmData to the given path.
/// @param afmData The data to be writen.
/// @param path The file path to write the values to.
void afmData_write(AFMData *afmData, char *path, int rowStart, int rowEnd);

#endif // AFMDATAWRITER_H