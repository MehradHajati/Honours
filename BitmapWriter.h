#ifndef BITMAPWRITER_H
#define BITMAPWRITER_H

#include <stdio.h>

/// @brief Writes a PGM formatted file from the given values.
/// @param fn The name of the output file (extension included).
/// @param vals The values to write [0-255].
/// @param startRow The row at which to start writing.
/// @param endRow The row at which to stop writing.
/// @param startCol The column at which to start writing.
/// @param endCol The column at which to stop writing.
/// @param yFlip 1 to write from bottom. 0 to write from top.
void bitmapWriter_writePGM(char *fn, double **vals, int startRow, int endRow, int startCol, int endCol, int yFlip);


#endif // BITMAPWRITER_H
