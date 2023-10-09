#include "Matrix.h"

/**
 * @brief Creates a new Matrix with the given number of rows and columns.
 * 
 * @param nrow The number of rows.
 * @param ncol The number of columns.
 * @return A new Matrix with the given number of rows and columns.
 */
struct Matrix *matrix_new(int nrow, int ncol){
    struct Matrix *m = (struct Matrix *)malloc(sizeof(struct Matrix));
    m->vals = (double **)malloc(sizeof(double *) * nrow);
    int i;
    for(i = 0; i < nrow; i++){
        m->vals[i] = (double *)calloc(ncol, sizeof(double));
    }
    m->nrow = nrow;
    m->ncol = ncol;
    return m;
}
 
/**
 * @brief Creates a new Matrix with the same values as the given Matrix.
 * 
 * @param M The matrix to copy.
 * @return A copy of the given matrix. 
 */
struct Matrix *matrix_copy(struct Matrix *M){
    struct Matrix *cpy = matrix_new(M->nrow, M->ncol);
    int row, col;
    for(row = 0; row < M->nrow; row++){
        for(col = 0; col < M->ncol; col++){
            cpy->vals[row][col] = M->vals[row][col];
        }
    }
    return cpy;
}
 
/**
 * @brief Frees memory allocated for the given Matrix.
 * 
 * @param M The Matrix to free.
 */
void matrix_free(struct Matrix *M){
    int i;
    for(i = 0; i < M->nrow; i++){
        free(M->vals[i]);
    }
    free(M->vals);
    free(M);
}

/**
 * @brief Multiplies M1 by M2 using matrix multiplication. Requires the number of columns in M1 to equal the number of rows in M2.
 * 
 * @param M1 The first Matrix.
 * @param M2 The second Matrix.
 * @return The product of M1 and M2 using matrix multiplication.
 */
struct Matrix *matrix_multiply(struct Matrix *M1, struct Matrix *M2){
    if(M1->ncol != M2->nrow){
        printf("The given matrices cannot be multiplied together.\n A:%dx%d\nB:%dx%d\n", M1->nrow, M1->ncol, M2->nrow, M2->ncol);
        exit(1);
    }
    int row, col;
    struct Matrix *product = matrix_new(M1->nrow, M2->ncol);

    for(row = 0; row < product->nrow; row++){
        for(col = 0; col < product->ncol; col++){
            product->vals[row][col] = matrix_rowColDotProduct(M1, M2, row, col, M1->ncol);
        }
    }
    return product;
}

/**
 * @brief Scales the given Matrix by the given scalar.
 * 
 * @param M The Matrix to scale.
 * @param scalar he scalar by which to multiply M.
 */
void matrix_scalarMultiply(struct Matrix *M, double scalar){
    int row, col;
    for(row = 0; row < M->nrow; row++){
        for(col = 0; col < M->ncol; col++){
            M->vals[row][col] *= scalar;
        }
    }
}

/**
 * @brief Gets the sum of the two given matrices.
 * 
 * @param M1 The first matrix.
 * @param M2 The second Matrix.
 * @return The result of the addition. 
 */
struct Matrix * matrix_add(struct Matrix *M1, struct Matrix *M2){
    struct Matrix *sum = matrix_new(M1->ncol, M1->nrow);
    int row, col;
    for(row = 0; row < M1->nrow; row++){
        for(col = 0; col < M2->ncol; col++){
            sum->vals[col][row] = M1->vals[row][col] + M2->vals[row][col];
        }
    }
    return sum;
}


/**
 * @brief Gets the transpose of the given Matrix.
 * 
 * @param M The matrix to transpose.
 */
void matrix_transpose(struct Matrix **M){
    struct Matrix *T = matrix_new((*M)->ncol, (*M)->nrow);
    int row, col;
    for(row = 0; row < (*M)->nrow; row++){
        for(col = 0; col < (*M)->ncol; col++){
            T->vals[col][row] = (*M)->vals[row][col];
        }
    }
    matrix_free(*M);
    *M = NULL;
    *M = T;
}

/**
 * @brief Gets the inverse of the given matrix if one exists.
 * 
 * @param M The matrix to invert.
 */
void matrix_inverse(struct Matrix **M){
    int row, col;
    int size = (*M)->nrow;
    struct Matrix *augmented, *inv;
    augmented = matrix_new(size, size*size);
    inv = matrix_new(size, size);
    
    for(row = 0; row < size; row++){
        for(col = 0; col < size; col++){
            augmented->vals[row][col] = (*M)->vals[row][col];
        }
        augmented->vals[row][row+size] = 1.0;
    }
    matrix_rref(augmented);
    for(row = 0; row < size; row++){
        for(col = 0; col < size; col++){
            inv->vals[row][col] = augmented->vals[row][col + size];
        }
    }
    matrix_free(augmented);
    matrix_free(*M);
    *M = NULL;
    *M = inv;
}

/**
 * @brief Gets the left inverse (Moore Penrose) of the given Matrix: (M^T * M)^(-1) * M^T
 * 
 * @param M The matrix to get the left inverse of.
 * @return The left inverse of M.
 */
struct Matrix *matrix_leftInverse(struct Matrix *M){
    struct Matrix *T = matrix_copy(M), *TM, *leftInv;

    matrix_transpose(&T);
    TM = matrix_multiply(T, M);
    matrix_inverse(&TM);
    leftInv = matrix_multiply(TM, T); // TM has been inverted

    matrix_free(T);
    matrix_free(TM);
    return leftInv;
}

/**
 * @brief Gets the reduced row echelon form of the given matrix
 * 
 * @param M The matrix to reduce.
 */
void matrix_rref(struct Matrix *M){
    int nrow, ncol, lead, row, col, i;
    double *tmpRow, lv;
    nrow = M->nrow;
    ncol = M->ncol;
    lead = 0;    

    for(row = 0; row < nrow; row++){
        if(lead >= ncol) break;
        i = row;
        while(fabs(M->vals[i][lead]) < EPS){
            if(++i == nrow){
                i = row;
                if(++lead == ncol) return;
            }
        }
        tmpRow = M->vals[row];
        M->vals[row] = M->vals[i];
        M->vals[i] = tmpRow;
        lv = M->vals[row][lead];
        for(col = 0; col < ncol; col++){
            M->vals[row][col] /= lv;
        }
        for(i = 0; i < nrow; i++){
            if(i != row){
                lv = M->vals[i][lead];
                for(col = 0; col < ncol; col++){
                    M->vals[i][col] -= lv * M->vals[row][col];
                }
            }
        }
        lead++;
    }
}


/**
 * @brief Computes the dot product of the given arrays.
 * 
 * @param M1 The first array.
 * @param M2 The second array.
 * @param row The row in question.
 * @param col The column in question.
 * @param size The size of the arrays.
 * @return The dot product of row and col.
 */
double matrix_rowColDotProduct(struct Matrix *M1, struct Matrix *M2, int row, int col, int size){
    int i;
    double total = 0.0;
    for(i = 0; i < size; i++){
        total += M1->vals[row][i] * M2->vals[i][col];
    }
    return total;
}

/**
 * @brief  Creates a string representation of the given matrix.
 * 
 * @param M The matrix to get a string representation of.
 * @return A string representing the given matrix.
 */
char *matrix_toString(struct Matrix *M){
    static int BASE_STR_LEN = 16;
    int count = 0, capacity = BASE_STR_LEN;
    char *str = (char *)malloc(sizeof(char) * capacity);
    str[0] = '\0';
    char buffer[100], *tmp;
    int row, col;
    for(row = 0; row < M->nrow; row++){
        for(col = 0; col < M->ncol; col++){
            sprintf(buffer, "%f", M->vals[row][col]);
            appendToStr(&str, buffer, &count, &capacity);
            appendToStr(&str, " ", &count, &capacity);
        }
        appendToStr(&str, "\n", &count, &capacity);
    }
    return str;
}


/**
 * @brief Creates an rotation matrix that rotates degrees degrees about the x axis.
 * 
 * @param degrees degrees to rotate by
 * @return struct Matrix 
 */
struct Matrix *matrix_xRot(double degrees){
    double rads = degrees * PI / 180.0;
    struct Matrix *rotMatrix = matrix_new(3,3);

    rotMatrix->vals[0][0] = 1; rotMatrix->vals[0][1] = 0; rotMatrix->vals[0][2] = 0;    
    rotMatrix->vals[1][0] = 0; rotMatrix->vals[1][1] = cos(rads); rotMatrix->vals[1][2] = -sin(rads);    
    rotMatrix->vals[2][0] = 0; rotMatrix->vals[2][1] = sin(rads); rotMatrix->vals[2][2] = cos(rads);
    return rotMatrix;
}

/**
 * @brief Creates an rotation matrix that rotates degrees degrees about the y axis.
 * 
 * @param degrees degrees to rotate by
 * @return struct Matrix
 */
struct Matrix *matrix_yRot(double degrees){
    double rads = degrees * PI / 180.0;
    struct Matrix *rotMatrix = matrix_new(3,3);

    rotMatrix->vals[0][0] = cos(rads); rotMatrix->vals[0][1] = 0; rotMatrix->vals[0][2] = sin(rads);    
    rotMatrix->vals[1][0] = 0; rotMatrix->vals[1][1] = 1; rotMatrix->vals[1][2] = 0;    
    rotMatrix->vals[2][0] = -sin(rads); rotMatrix->vals[2][1] = 0; rotMatrix->vals[2][2] = cos(rads);
    return rotMatrix;
}

/**
 * @brief Creates an rotation matrix that rotates degrees degrees about the z axis.
 * 
 * @param degrees degrees to rotate by
 * @return struct Matrix
 */
struct Matrix *matrix_zRot(double degrees){
    double rads = degrees * PI / 180.0;
    struct Matrix *rotMatrix = matrix_new(3,3);

    rotMatrix->vals[0][0] = cos(rads); rotMatrix->vals[0][1] = -sin(rads); rotMatrix->vals[0][2] = 0;    
    rotMatrix->vals[1][0] = sin(rads); rotMatrix->vals[1][1] = cos(rads); rotMatrix->vals[1][2] = 0;    
    rotMatrix->vals[2][0] = 0; rotMatrix->vals[2][1] = 0; rotMatrix->vals[2][2] = 1;
    return rotMatrix;
}


