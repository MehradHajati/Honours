#include "Vector3.h"

Vector3 vector3_new(double x, double y, double z){
    Vector3 v;
    v.x = x;
    v.y = y;
    v.z = z;
    return v;
}

double vector3_Distance3(Vector3 v1, Vector3 v2){
    double xDiff = v1.x - v2.x;
    double yDiff = v1.y - v2.y;
    double zDiff = v1.z - v2.z;
    return sqrt((xDiff*xDiff) + (yDiff*yDiff) + (zDiff*zDiff));
}

double vector3_Distance2(Vector3 v1, Vector3 v2){
    double xDiff = v1.x - v2.x;
    double yDiff = v1.y - v2.y;
    return sqrt((xDiff*xDiff) + (yDiff*yDiff));
}

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

Vector3 vector3_add(Vector3 v1, Vector3 v2){
    return vector3_new(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

Vector3 vector3_scale(Vector3 v, double scalar){
    return vector3_new(v.x * scalar, v.y * scalar, v.z * scalar);
}

double vector3_angleBetween(Vector3 v1, Vector3 v2){
    return acos(vector3_dot(v1, v2) / (vector3_magnitude(v1) * vector3_magnitude(v2)));
}

double vector3_dot(Vector3 v1, Vector3 v2){
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

Vector3 vector3_cross(Vector3 v1, Vector3 v2){
    return vector3_new(
        v1.y*v2.z - v1.z*v2.y, 
        v1.z*v2.x - v1.x*v2.z, 
        v1.x*v2.y - v1.y*v2.x);
}

double vector3_magnitude(Vector3 v){
    return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

Vector3 vector3_random(double maxRadius){
    double theta = 360.0 * ((double)rand() / (double)RAND_MAX); // Random theta
    double r = maxRadius * ((double)rand() / (double)RAND_MAX); // Random radius between 0 and maxRadius
    Vector3 noise = vector3_new(r * cos(theta), r * sin(theta), 0);
    return noise;
}

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
