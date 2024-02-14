#ifndef FACETSWRITER_H
#define FACETSWRITER_H

#include "Facets.h"

/// @brief Writes the theta and phi maps from the given facets struct.
/// @param facets The facets from which to write theta and phi maps. 
/// @param thetaPath The path to the theta map file.
/// @param phiPath The path to the phi map file.
void facets_write(Facets facets, char *thetaPath, char *phiPath);

#endif // FACETSWRITER_H