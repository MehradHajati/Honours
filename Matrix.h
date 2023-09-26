#ifndef MATRIX_H
#define MATRIX_H
#define EPS 0.00000001
#define PI 3.14159265359

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "Utility.h"

struct Matrix{
    int nrow;
    int ncol;
    double **vals;
};

/// @brief Creates a new Matrix with the given number of rows and columns.
/// @param nrow The number of rows.
/// @param ncol The number of columns.
/// @return A new Matrix with the given number of rows and columns.
struct Matrix *matrix_new(int nrow, int ncol);

/// @brief Creates a new Matrix with the same values as the given Matrix.
/// @param M The matrix to copy.
/// @return A copy of the given matrix.
struct Matrix *matrix_copy(struct Matrix *M);

/// @brief Frees memory allocated for the given Matrix.
/// @param M The Matrix to free.
void matrix_free(struct Matrix *M);

/// @brief Multiplies M1 by M2 using matrix multiplication. Requires the number of columns in M1 to equal the number of rows in M2.
/// @param M1 The first Matrix.
/// @param M2 The second Matrix.
/// @return The product of M1 and M2 using matrix multiplication.
struct Matrix *matrix_multiply(struct Matrix *M1, struct Matrix *M2);

/// @brief Scales the given Matrix by the given scalar.
/// @param M The Matrix to scale.
/// @param scalar The scalar by which to multiply M.
/// @return A new Matrix equal to the M*scalar.
void matrix_scalarMultiply(struct Matrix *M, double scalar);

/// @brief Gets the sum of the two given matrices.
struct Matrix * matrix_add(struct Matrix *M1, struct Matrix *M2);

/// @brief Gets the transpose of the given Matrix.
/// @param M The matrix to transpose.
/// @return The transpose of the given matrix.
void matrix_transpose(struct Matrix **M);

/// @brief Gets the inverse of the given matrix if one exists.
/// @param M The matrix to invert.
/// @return The inverse of te given matrix.
void matrix_inverse(struct Matrix **M);

/// @brief Gets the left inverse (Moore Penrose) of the given Matrix: (M^T * M)^(-1) * M^T
/// @param M The matrix to get the left inverse of.
/// @return The left inverse of M.
struct Matrix *matrix_leftInverse(struct Matrix *M);

/// @brief Gets the reduced row echelon form of the given matrix
/// @param M The matrix to reduce.
/// @return A matrix equal to the reduced row echilon form of M.
void matrix_rref(struct Matrix *M);

/// @brief Computes the dot product of the given arrays.
/// @param row The first array.
/// @param col The second array.
/// @param size The size of the arrays.
/// @return The dot product of row and col.
double matrix_rowColDotProduct(struct Matrix *M1, struct Matrix *M2, int row, int col, int size);

/// @brief Creates a string representation of the given matrix.
/// @param M The matrix to get a string representation of.
/// @return A string representing the given matrix.
char *matrix_toString(struct Matrix *M);

/// @brief Creates an rotation matrix that rotates degrees degrees about the x axis.
struct Matrix *matrix_xRot(double degrees);

/// @brief Creates an rotation matrix that rotates degrees degrees about the y axis.
struct Matrix *matrix_yRot(double degrees);

/// @brief Creates an rotation matrix that rotates degrees degrees about the z axis.
struct Matrix *matrix_zRot(double degrees);

#endif // MATRIX_H