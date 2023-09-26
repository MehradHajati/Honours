#ifndef FACETSREADER_H
#define FACETSREADER_H

#include <stdio.h>
#include "Facets.h"

/// @brief Creates a facets struct form the given theta and phi map files.
/// @param thetaFn The path to the theta map. 
/// @param phiFn The path to the phi map.
/// @return A facets struct.
Facets facets_readFromFiles(char *thetaFn, char *phiFn);

#endif // FACETSREADER_H