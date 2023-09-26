#ifndef BANDCONTRASTWRITER_H
#define BANDCONTRASTWRITER_H

#include <stdio.h>
#include "BandContrast.h"

/// @brief Writes the given band contrast struct to the given path.
/// @param bc The band contrast to write.
/// @param path The path to the file to which to write.
/// @param startRow The row at which to start writing.
/// @param endRow The row at which to stop writing.
/// @param startCol The column at which to start writing.
/// @param endCol The column at which to stop writing.
void bandContrast_write(BandContrast *bc, char *path, int rowStart, int rowEnd, int colStart, int colEnd);

#endif // BANDCONTRASTWRITER_H