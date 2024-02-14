#ifndef FACETS_H
#define FACETS_H
#define _USE_MATH_DEFINES

#define MIN_FACET_SIZE 3
#define MAX_FACET_SIZE 131 // Ideal value??
#define THETA_DEFAULT -5
#define PHI_DEFAULT -200
#define MAX_RSS_AVG 0.001 // Ideal value??

#include <time.h>
#include <math.h>
#include "AFM.h"
#include "Matrix.h"
#include "Vector3.h"

typedef struct {

    int nrow;
    int ncol;
    double **thetaMap;
    double **phiMap;

} Facets;

typedef struct {

    int row;
    int col;
    int size;
    double rssAvg;
    double theta;
    double phi;

} Face;

/// @brief Creates a new Facets struct with the given number of rows and columns for storing theta and phi values.
/// @param nrow The number of rows.
/// @param ncol The number of columns.
/// @return A new Facets struct.
Facets facets_new(int nrow, int ncol);

/// @brief Creates a new face for determining facets.
/// @param row Row position of the face.
/// @param col Column position of the face.
/// @param size The size of the face.
/// @param rssAvg The reduced sum of the squares average.
/// @param theta The theta of the face normal.
/// @param phi The phi of the face normal.
/// @return A face struct.
Face face_new(int row, int col, int size, double rssAvg, double theta, double phi);

/// @brief Computes the facets from the given AFM data.
/// @param afmData The AFM data for which to get facet information.
/// @param minFacetSize The smallest facet size (ODD).
/// @param maxFacetSize The largest facets size (ODD).
/// @param maxRSSAvg The largest reduced sum of the squares allowed.
/// @param binIters The number of time to bin the given AFM data before computing the facets.  
/// @return The facets for the given AFM data.
Facets facets_compute(AFMData *afmData, int minFacetSize, int maxFacetSize, double maxRSSAvg, int binIters);

/// @brief Repeatedly places "tiles" of the given size on the afmData grid and attempts to find sections with sufficiently low RSS averages.
/// If found, theta and phi angles are computed and stored in the given Facets struct where available.
/// @param afmData The data to place tiles on.
/// @param facets The Facets struct to store theta and phi information in.
/// @param facetSize The size of each tile to place.
void facets_tileArea(AFMData *afmData, Facets *facets, int facetSize, double maxRSSAvg);

/// @brief Creates an array of faces that have an RSS average within the required limits.
/// @param afmData The data to create faces from.
/// @param tileSize The size of a face.
/// @param outArraySize The size of the Face array created.
/// @return An array of faces.
Face *facets_findValidFaces(Facets *facets, AFMData *afmData, int tileSize, int *outArraySize, double maxRSSAvg);

/// @brief Attempts to find a face at the given location in the afmData. If found, the face is stored in tmpFace and a non-zero integer is returned.
/// @param afmData The afmData to find a face in.
/// @param tmpFace The face to store the found face in.
/// @param tileRow The row of the upper left corner of the face.
/// @param tileCol The column of the upper left corner of the face.
/// @param tileSize The size of the face to make.
/// @return 0 if a face is not found. Non-zero otherwise.
int facets_findFace(AFMData *afmData, Face *tmpFace, struct Matrix *betaMatrix, struct Matrix *X, int tileRow, int tileCol, int tileSize, double maxRSSAvg);

/// @brief Adds toAdd to faces.
/// @param faces A pointer to the array of faces to add to.
/// @param toAdd The face to add.
/// @param count The number of faces currently in faces.
/// @param capacity The capacity of faces.
void facets_addToFaces(Face **faces, Face toAdd, int *count, int *capacity);

/// @brief A comparator for sorting Face structs by RSS.
/// @param f1 Face 1
/// @param f2 Face 2
/// @return The priority of f1 where smaller numbers have higher priority.
int facets_cmpFacesByRSS(const void *f1, const void *f2);

#endif // FACETS_H