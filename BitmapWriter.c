#include "BitmapWriter.h"


/**
 * @brief Writes a PGM formatted file from the given values.
 * 
 * @param fn The name of the output file (extension included).
 * @param vals The values to write [0-255].
 * @param startRow The row at which to start writing.
 * @param endRow The row at which to stop writing.
 * @param startCol The column at which to start writing.
 * @param endCol The column at which to stop writing.
 * @param yFlip 1 to write from bottom. 0 to write from top.
 */
void bitmapWriter_writePGM(char *fn, double **vals, int startRow, int endRow, int startCol, int endCol, int yFlip){
    int row, col, nrows, ncols;
    nrows = endRow - startRow;
    ncols = endCol - startCol;
    FILE *pgm = fopen(fn, "wb");
    printf(".pgm filename: %s\n", fn);

    // Magic Number for writing greyscale 255 in ASCII
    fprintf(pgm, "P2\n");

    // Width and height
    fprintf(pgm, "%d %d\n", ncols, nrows);

    // Max value
    fprintf(pgm, "%d\n", 255);

    // Write array to pgm.
    for(row = (yFlip ? endRow-1 : startRow); row != (yFlip ? startRow-1 : endRow); row += (yFlip ? -1 : 1)){
        for(col = startCol; col < endCol; col++){
            //if(row == col) printf("writer %d %d\n", row, (int)vals[row][col]);
            fprintf(pgm, "%d ", ((int)vals[row][col] < 255 ? (int)vals[row][col] : 255));
        }
        fprintf(pgm, "\n");
    }
    fclose(pgm);
}
