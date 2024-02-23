#ifndef BANDCONTRASTREADER_H
#define BANDCONTRASTREADER_H

#include <stdio.h>
#include "BandContrast.h"

/// @brief Creates a band contrast struct from the data in the given file.
/// @param fn The name of the file to read.
/// @return A band contrast struct.
BandContrast bandContrast_readFromFile(char *fn);

#endif // BANDCONTRASTREADER_H
