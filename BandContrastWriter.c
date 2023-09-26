#include "BandContrastWriter.h"

void bandContrast_write(BandContrast *bc, char *path, int rowStart, int rowEnd, int colStart, int colEnd){
    printf("Band Contrast Path: %s\n", path);
    int row, col, i;
    FILE *out;

    out = fopen(path, "w");
    for(row = rowStart; row < rowEnd; row++){
        for(col = colStart; col < colEnd; col++){
            fprintf(out, "%d ", ((int)bc->greyScale[row][col] > 255 ? 255 : (int)bc->greyScale[row][col]));
        }
        fprintf(out, "\n");
    }
    fclose(out);
}
