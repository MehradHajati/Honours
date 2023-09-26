#ifndef VECTOR_H
#define VECTOR_H
#define _USE_MATH_DEFINES

#include <stdlib.h>
#include <math.h>
#include <float.h>

typedef struct
{
    double x, y, z;
} Vector3;

/// @brief Creates a new Vector3 with the given x, y, z.
/// @param x The x-component of this Vector3.
/// @param y The y-component of this Vector3.
/// @param z The z-component of this Vector3.
/// @return A new Vector3 with the given x, y, z.
Vector3 vector3_new(double x, double y, double z);

/// @brief Gets the euclidean distance between the tips of the Vector3's v1 and v2.
/// @param v1 The first vector.
/// @param v2 The second vector.
/// @return The euclidean distance between the tips of v1 and v2.
double vector3_Distance3(Vector3 v1, Vector3 v2);

/// @brief Gets the euclidean distance between the tips of the Vector3's v1 and v2 IGNORING THE Z COMPONENT.
/// @param v1 The first vector.
/// @param v2 The second vector.
/// @return The euclidean distance between the tips of v1 and v2 IGNORING THE Z COMPONENT.
double vector3_Distance2(Vector3 v1, Vector3 v2);

/// @brief Computes the shortest distance between two Vecter3's under periodic 
/// boundary conditions. Ignores the Z-Component.
/// @param v1 The first vector.
/// @param v2 The second vector.
/// @param loBounds The low bounds on the periodic boundary.
/// @param hiBounds The high bounds on the periodic boundary.
/// @return The shortest distance between two Vecter3's under periodic 
/// boundary conditions. Ignores the Z-Component.
double vector3_shortestPeriodicBoundaryDistance2D(Vector3 v1, Vector3 v2, Vector3 loBounds, Vector3 hiBounds);

/// @brief Computes the shortest distance between two Vecter3's under periodic 
/// boundary conditions.
/// @param v1 The first vector.
/// @param v2 The second vector.
/// @param loBounds The low bounds on the periodic boundary.
/// @param hiBounds The high bounds on the periodic boundary.
/// @return The shortest distance between two Vecter3's under periodic 
/// boundary conditions. Ignores the Z-Component.
double vector3_shortestPeriodicBoundaryDistance3D(Vector3 v1, Vector3 v2, Vector3 loBounds, Vector3 hiBounds);

/// @brief Scales the given vector by the given scalar.
/// @param v The vector to scale.
/// @param scalar The scalar by which to scale the vector.
/// @return The scaled vector.
Vector3 vector3_scale(Vector3 v, double scalar);

/// @brief Computes the angle between two vectors in radians.
/// @param v1 The first vector.
/// @param v2 The second vector.
/// @return The angle between the two vectors in radians.
double vector3_angleBetween(Vector3 v1, Vector3 v2);

/// @brief Computes the dot product of the two vectors.
/// @param v1 The first vector.
/// @param v2 The second vector.
/// @return The dot product of the two given vectors.
double vector3_dot(Vector3 v1, Vector3 v2);

/// @brief Computes the cross product of the two vectors.
/// @param v1 The first vector. 
/// @param v2 The second vector.
/// @return The cross product between the two given vectors.
Vector3 vector3_cross(Vector3 v1, Vector3 v2);

/// @brief Gets the magnitude/length of the given vector.
/// @param v The vector of which to get the magnitude. 
/// @return The magnitude/length of the given vector.
double vector3_magnitude(Vector3 v);

/// @brief Adds the given Vector3s together and returns the result.
/// @param v1 The first Vector3.
/// @param v2 The second Vector3.
/// @return The sum of the twon Vector3s.
Vector3 vector3_add(Vector3 v1, Vector3 v2);

/// @brief Generates a random Vector3 with a size between 0 and maxRadius (no z component).
/// @param maxRadius The largest size of the Vector allowable.
/// @return A random Vector within a circle of radius maxRadius.
Vector3 vector3_random(double maxRadius);

/// @brief Create a vector from the given theta and phi angles.
/// @param theta The angle down in degrees. 
/// @param phi The angle around in degrees.
/// @return A vector created from theta and phi.
Vector3 vector3_fromThetaPhi(double theta, double phi);

#endif // VECTOR_H