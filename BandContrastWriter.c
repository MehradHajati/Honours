#include "BandContrastWriter.h"

 
/**
 * @brief Writes the given band contrast struct to the given path.
 * 
 * @param bc The band contrast to write.
 * @param path The path to the file to which to write.
 * @param rowStart The row at which to start writing.
 * @param rowEnd The row at which to stop writing.
 * @param colStart The column at which to start writing.
 * @param colEnd The column at which to stop writing.
 */
void bandContrast_write(BandContrast *bc, char *path, int rowStart, int rowEnd, int colStart, int colEnd){
    printf("Band Contrast Path: %s\n", path);
    int row, col, i;
    FILE *out;
    // opening the file to write
    out = fopen(path, "w");
    // for loops to write each row and col
    for(row = rowStart; row < rowEnd; row++){
        for(col = colStart; col < colEnd; col++){
            fprintf(out, "%d ", ((int)bc->greyScale[row][col] > 255 ? 255 : (int)bc->greyScale[row][col]));
        }
        fprintf(out, "\n");
    }
    fclose(out);
}
