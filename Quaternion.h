#ifndef QUATERNION_H
#define QUATERNION_H
#define _USE_MATH_DEFINES

#include <math.h>
#include <stdio.h>
#include "Vector3.h"

typedef struct {

    double x;
    double y;
    double z;
    double w;

} Quaternion;

/// @brief Creates a Quaternion with the given values.
Quaternion quaternion_new(double x, double y, double z, double w);

/// @brief Creates a Quaternion where all of its values are 0. Used as NULL.
Quaternion quaternion_zero();

/// @brief Checks if the given quaternions are equivalent.
/// @param q1 The first quaternion.
/// @param q2 The second quaternion
/// @return True if they are equal, false otherwise.
int quaternion_equals(Quaternion q1, Quaternion q2);

/// @brief The Hamiltonian Product of two quaternions. If q1 and q2 are rotation quaternions (magnitudes of 1), this 
/// returns the result of rotating q2 by q1. 
/// @param q1 The first quaternion.
/// @param q2 The second quaternion.
/// @return The hamiltonian product of the two given quaternions.
Quaternion quaternion_multiply(Quaternion q1, Quaternion q2);

/// @brief Normalizes this quaternion to have a magnitude of 1, making it a rotation quaternion.
/// @param quat The quaternion to normalize.
void quaternion_normalize(Quaternion *quat);

/// @brief Converts the given euler angles to a rotation quaternion. Assumes the euler angles are intrinsic and about the Z-X-Z axes.
/// @param euler The euler angles to convert.
/// @return The rotation quaternion representation of the given euler angles.
Quaternion quaternion_fromEulerAngles(Vector3 euler);

/// @brief Converts the given quaternion to intrinsic proper euler angles in radians in the Z-X-Z 
/// sequence following the generalized converstion algorithm found here: 
/// https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9648712/
/// @param quat The quaternion to convert.
/// @return A Vector3 containing the euler angles in the Z-X-Z order.
Vector3 quaternion_toEulerAngles(Quaternion quat);

/// @brief Converts the given quaternion to intrinsic proper euler angles in degrees in the Z-X-Z 
/// sequence following the generalized converstion algorithm found here: 
/// https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9648712/
/// @param quat The quaternion to convert.
/// @return A Vector3 containing the euler angles in the Z-X-Z order.
Vector3 quaternion_toEulerAnglesDegrees(Quaternion quat);

/// @brief Rotates the given quaternion about degrees many degrees about the given axis.
/// @param quat The quaternion to rotate.
/// @param degrees The number of degrees by which to rotate the quaternion.
/// @param axis The axis to rotate the quaternion around ('x','y', or 'z').
void quaternion_rotateAboutAxisByDegrees(Quaternion *quat, double degrees, char axis);

/// @brief Rotates the given quaternion intrinsically by p1, P, and p2.
/// @param quat The quaternion to rotate
/// @param p1 The angle by which to rotate around the Z axis the first time.
/// @param P The angle by which to rotate around the X axis.
/// @param p2 The angle by which to rotate around the Z axis the second time.
void quaternion_applyIntrinsicZXZRotation(Quaternion *quat, double p1, double P, double p2);

/// @brief Rotates the given vector about an axis made from x, y, z by deg degrees.
/// @param v The vector to rotate
/// @param x The x component of the axis.
/// @param y The y component of the axis.
/// @param z The z component of the axis.
/// @param deg Degrees by which to rotate.
void quaternion_rotateVector3AxisAngle(Vector3 *v, double x, double y, double z, double deg);

/// @brief Rotates the given vector by the given quaternion.
/// @param v The vector to rotate. 
/// @param q The quaternion by which to rotate the vector.
void quaternion_rotateVector3ByQuat(Vector3 *v, Quaternion q);

#endif // QUATERNION_H