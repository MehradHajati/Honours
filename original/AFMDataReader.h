#ifndef AFMDATAREADER_H
#define AFMDATAREADER_H
#define MAX_AFM_LINE_WIDTH 100000

#include <stdio.h>
#include "AFM.h"
#include "Utility.h"

/// @brief Creates an AFMData struct from the data in the given file.
/// @param fileName The name of the file to read.
/// @return An AFMData struct.
AFMData afmData_readFromFile(char *fileName);

/// @brief Creates an AFMData struct from the data in the given file.
/// @param fileName The name of the file to read.
/// @return An AFMData struct.
AFMData afmData_readBinnedFromFile(char *fileName);

#endif // AFMDATAREADER_H