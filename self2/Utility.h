#ifndef UTILITY_H
#define UTILITY_H

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "Vector3.h"

/// @brief Generates a random number between 0 and 1.
/// @return A random number between 0 and 1.
double randomPercentage();

/// @brief Appends toAppend to the end of str.
/// @param str The string to append to.
/// @param toAppend The string to append to the end of str.
/// @param count The number of characters currently in str.
/// @param capacity The capacity of str.
void appendToStr(char **str, char *toAppend, int *count, int *capacity);

/// @brief Reads the next line from the given file.
/// @param file The file to read from.
/// @return The line that was read from the file.
char * readLine(FILE **file);

#endif // UTILITY_H