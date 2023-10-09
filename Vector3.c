#include "Vector3.h"

 
/**
 * @brief Creates a new Vector3 with the given x, y, z.
 * 
 * @param x The x-component of this Vector3.
 * @param y The y-component of this Vector3.
 * @param z The z-component of this Vector3.
 * @return A new Vector3 with the given x, y, z.
 */
Vector3 vector3_new(double x, double y, double z){
    Vector3 v;
    v.x = x;
    v.y = y;
    v.z = z;
    return v;
}


/**
 * @brief Gets the euclidean distance between the tips of the Vector3's v1 and v2.
 * 
 * @param v1 The first vector. 
 * @param v2 The second vector.
 * @return The euclidean distance between the tips of v1 and v2.
 */
double vector3_Distance3(Vector3 v1, Vector3 v2){
    double xDiff = v1.x - v2.x;
    double yDiff = v1.y - v2.y;
    double zDiff = v1.z - v2.z;
    return sqrt((xDiff*xDiff) + (yDiff*yDiff) + (zDiff*zDiff));
}


/**
 * @brief Gets the euclidean distance between the tips of the Vector3's v1 and v2 IGNORING THE Z COMPONENT.
 * 
 * @param v1 The first vector.
 * @param v2 The second vector.
 * @return The euclidean distance between the tips of v1 and v2 IGNORING THE Z COMPONENT.
 */
double vector3_Distance2(Vector3 v1, Vector3 v2){
    double xDiff = v1.x - v2.x;
    double yDiff = v1.y - v2.y;
    return sqrt((xDiff*xDiff) + (yDiff*yDiff));
}


/**
 * @brief Computes the shortest distance between two Vecter3's under periodic boundary conditions. Ignores the Z-Component.
 * 
 * @param v1 The first vector.
 * @param v2 The second vector.
 * @param loBound The low bounds on the periodic boundary.
 * @param hiBound The high bounds on the periodic boundary.
 * @return The shortest distance between two Vecter3's under periodic boundary conditions. Ignores the Z-Component.
 */
double vector3_shortestPeriodicBoundaryDistance2D(Vector3 v1, Vector3 v2, Vector3 loBound, Vector3 hiBound){
    double xDiff = fabs(v1.x - v2.x);
    double yDiff = fabs(v1.y - v2.y);

    double v1ToLoBoundX = fabs(v1.x - loBound.x); 
    double v1ToHiBoundX = fabs(v1.x - hiBound.x);
    double v1ToLoBoundY = fabs(v1.y - loBound.y);
    double v1ToHiBoundY = fabs(v1.y - hiBound.y);

    double v2ToLoBoundX = fabs(v2.x - loBound.x);
    double v2ToHiBoundX = fabs(v2.x - hiBound.x);
    double v2ToLoBoundY = fabs(v2.y - loBound.y);
    double v2ToHiBoundY = fabs(v2.y - hiBound.y);

    double v1SmallestX = v1ToLoBoundX < v1ToHiBoundX ? v1ToLoBoundX : v1ToHiBoundX;
    double v1SmallestY = v1ToLoBoundY < v1ToHiBoundY ? v1ToLoBoundY : v1ToHiBoundY;
    double v2SmallestX = v2ToLoBoundX < v2ToHiBoundX ? v2ToLoBoundX : v2ToHiBoundX;
    double v2SmallestY = v2ToLoBoundY < v2ToHiBoundY ? v2ToLoBoundY : v2ToHiBoundY;

    double boundaryX = v1SmallestX + v2SmallestX;
    double boundaryY = v1SmallestY + v2SmallestY;

    xDiff = xDiff < boundaryX ? xDiff : boundaryX;
    yDiff = yDiff < boundaryY ? yDiff : boundaryY;

    return sqrt((xDiff*xDiff) + (yDiff*yDiff));
}

 
/**
 * @brief Computes the shortest distance between two Vecter3's under periodic boundary conditions.
 * 
 * @param v1 The first vector.
 * @param v2 The second vector.
 * @param loBound The low bounds on the periodic boundary.
 * @param hiBound The high bounds on the periodic boundary.
 * @return The shortest distance between two Vecter3's under periodic boundary conditions. Ignores the Z-Component.
 */
double vector3_shortestPeriodicBoundaryDistance3D(Vector3 v1, Vector3 v2, Vector3 loBound, Vector3 hiBound){
    double xDiff = fabs(v1.x - v2.x);
    double yDiff = fabs(v1.y - v2.y);
    double zDiff = fabs(v1.z - v2.z);

    double v1ToLoBoundX = fabs(v1.x - loBound.x); 
    double v1ToHiBoundX = fabs(v1.x - hiBound.x);
    double v1ToLoBoundY = fabs(v1.y - loBound.y);
    double v1ToHiBoundY = fabs(v1.y - hiBound.y);

    double v2ToLoBoundX = fabs(v2.x - loBound.x);
    double v2ToHiBoundX = fabs(v2.x - hiBound.x);
    double v2ToLoBoundY = fabs(v2.y - loBound.y);
    double v2ToHiBoundY = fabs(v2.y - hiBound.y);

    double v1SmallestX = v1ToLoBoundX < v1ToHiBoundX ? v1ToLoBoundX : v1ToHiBoundX;
    double v1SmallestY = v1ToLoBoundY < v1ToHiBoundY ? v1ToLoBoundY : v1ToHiBoundY;
    double v2SmallestX = v2ToLoBoundX < v2ToHiBoundX ? v2ToLoBoundX : v2ToHiBoundX;
    double v2SmallestY = v2ToLoBoundY < v2ToHiBoundY ? v2ToLoBoundY : v2ToHiBoundY;

    double boundaryX = v1SmallestX + v2SmallestX;
    double boundaryY = v1SmallestY + v2SmallestY;

    xDiff = xDiff < boundaryX ? xDiff : boundaryX;
    yDiff = yDiff < boundaryY ? yDiff : boundaryY;

    return sqrt((xDiff*xDiff) + (yDiff*yDiff) + (zDiff*zDiff));
}


/**
 * @brief Adds the given Vector3s together and returns the result.
 * 
 * @param v1 The first Vector3.
 * @param v2 The second Vector3.
 * @return The sum of the twon Vector3s.
 */
Vector3 vector3_add(Vector3 v1, Vector3 v2){
    return vector3_new(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}


/**
 * @brief Scales the given vector by the given scalar.
 * 
 * @param v The vector to scale.
 * @param scalar The scalar by which to scale the vector.
 * @return The scaled vector. 
 */
Vector3 vector3_scale(Vector3 v, double scalar){
    return vector3_new(v.x * scalar, v.y * scalar, v.z * scalar);
}


/**
 * @brief Computes the angle between two vectors in radians.
 * 
 * @param v1 The first vector.
 * @param v2 The second vector.
 * @return The angle between the two vectors in radians.
 */
double vector3_angleBetween(Vector3 v1, Vector3 v2){
    return acos(vector3_dot(v1, v2) / (vector3_magnitude(v1) * vector3_magnitude(v2)));
}

 
/**
 * @brief Computes the dot product of the two vectors.
 * 
 * @param v1 The first vector.
 * @param v2 The second vector.
 * @return The dot product of the two given vectors.
 */
double vector3_dot(Vector3 v1, Vector3 v2){
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}
 
/**
 * @brief Computes the cross product of the two vectors.
 * 
 * @param v1 The first vector.
 * @param v2 The second vector.
 * @return The cross product between the two given vectors.
 */
Vector3 vector3_cross(Vector3 v1, Vector3 v2){
    return vector3_new(
        v1.y*v2.z - v1.z*v2.y, 
        v1.z*v2.x - v1.x*v2.z, 
        v1.x*v2.y - v1.y*v2.x);
}


/**
 * @brief Gets the magnitude/length of the given vector.
 * 
 * @param v The vector of which to get the magnitude.
 * @return The magnitude/length of the given vector.
 */
double vector3_magnitude(Vector3 v){
    return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

/**
 * @brief Generates a random Vector3 with a size between 0 and maxRadius (no z component).
 * 
 * @param maxRadius The largest size of the Vector allowable.
 * @return A random Vector within a circle of radius maxRadius.

 */
Vector3 vector3_random(double maxRadius){
    double theta = 360.0 * ((double)rand() / (double)RAND_MAX); // Random theta
    double r = maxRadius * ((double)rand() / (double)RAND_MAX); // Random radius between 0 and maxRadius
    Vector3 noise = vector3_new(r * cos(theta), r * sin(theta), 0);
    return noise;
}


/**
 * @brief Create a vector from the given theta and phi angles.
 * 
 * @param theta The angle down in degrees. 
 * @param phi The angle around in degrees.
 * @return A vector created from theta and phi.
 */
Vector3 vector3_fromThetaPhi(double theta, double phi){
    theta *= M_PI / 180.0;
    phi *= M_PI / 180.0;
    Vector3 n;
    double st = sin(theta);
    n.x = st * cos(phi);
    n.y = st * sin(phi);
    n.z = cos(theta);    
    return n;
}
