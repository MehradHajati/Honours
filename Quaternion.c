#include "Quaternion.h"

Quaternion quaternion_new(double x, double y, double z, double w){
    Quaternion quat;
    quat.x = x;
    quat.y = y;
    quat.z = z;
    quat.w = w;
    return quat;
}

Quaternion quaternion_zero(){
    return quaternion_new(0,0,0,0);
}

int quaternion_equals(Quaternion q1, Quaternion q2){
    return q1.w == q2.w && q1.x == q2.x && q1.y == q2.y && q1.z == q2.z;
}

Quaternion quaternion_multiply(Quaternion q1, Quaternion q2){
    Quaternion quat;

    quat.w = q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z;
    quat.x = q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y;
    quat.y = q1.w * q2.y - q1.x * q2.z + q1.y * q2.w + q1.z * q2.x;
    quat.z = q1.w * q2.z + q1.x * q2.y - q1.y * q2.x + q1.z * q2.w;

    return quat;
}

void quaternion_normalize(Quaternion *quat){    
    double magnitude = sqrt(quat->w*quat->w + quat->x*quat->x + quat->y*quat->y + quat->z*quat->z);
    quat->w /= magnitude;
    quat->x /= magnitude;
    quat->y /= magnitude;
    quat->z /= magnitude;
}

Quaternion quaternion_fromEulerAngles(Vector3 euler){
    euler = vector3_scale(euler, M_PI / 180.0);

    double thetaPlus, thetaMinus, cP, ttP, ttM;
    Quaternion quat = quaternion_new(0,0,0,0);
    thetaPlus = (euler.x + euler.z)/2.0;
    thetaMinus = (euler.x - euler.z)/2.0;
    cP = cos(euler.y);
    ttP = tan(thetaPlus);
    ttM = tan(thetaMinus);

    quat.w = sqrt((cP + 1.0) / (2.0 * (1.0 + ttP * ttP)));
    quat.x = sqrt((1.0 - cP) / (2.0 * (1.0 + ttM * ttM)));
    quat.y = quat.x * ttM;
    quat.z = quat.w * ttP;
    quaternion_normalize(&quat);
    return quat;
}

Vector3 quaternion_toEulerAngles(Quaternion quat){
    if(quaternion_equals(quat, quaternion_zero())) return vector3_new(0,0,0);
    double theta1, theta2, theta3, thetaPlus, thetaMinus, a, b, c, d, epsilon;
    int i, j, k;
    i = 3;
    j = 1;
    k = 2;

    epsilon = (double)((i - j) * (j - k) * (k - i))/2.0; // equivalent to Levi-Civita (-1 for counter-clockwise)
    a = quat.w;
    b = quat.z;
    c = quat.x;
    d = quat.y * epsilon;

    theta2 = acos(2.0 * ((a*a+b*b)/(a*a+b*b+c*c+d*d)) - 1.0);
    thetaPlus = atan2(b, a);
    thetaMinus = atan2(d, c);

    if(theta2 == 0){
        theta1 = 0;
        theta3 = 2 * thetaPlus;
    }
    else if(theta2 == M_PI_2){
        theta1 = 0;
        theta3 = 2 * thetaMinus;
    }
    else{
        theta1 = thetaPlus - thetaMinus;
        theta3 = thetaPlus + thetaMinus;
    }

    // Wrap:
    while(theta1 < 0) theta1 += 2 * M_PI;
    while(theta1 >= 2 * M_PI) theta1 -= 2 * M_PI;

    while(theta2 < 0) theta2 += M_PI;
    while(theta2 >= M_PI) theta2 -= M_PI;

    while(theta3 < 0) theta3 += 2 * M_PI;
    while(theta3 >= 2 * M_PI) theta3 -= 2 * M_PI;

    // These are being returned in an intrinsic (crystal) order.
    // Reverse order for extrinsic (lab).
    return vector3_new(theta3, theta2, theta1);
}

Vector3 quaternion_toEulerAnglesDegrees(Quaternion quat){
    return vector3_scale(quaternion_toEulerAngles(quat), 180.0 / M_PI);
}

void quaternion_rotateAboutAxisByDegrees(Quaternion *quat, double degrees, char axis){
    Quaternion rotateBy;
    double theta = degrees * M_PI / 180.0;
    if(axis == 'x') rotateBy = quaternion_new(sin(theta/2.0),0,0, cos(theta/2.0));
    else if(axis == 'y') rotateBy = quaternion_new(0,sin(theta/2.0),0, cos(theta/2.0));
    else if(axis == 'z') rotateBy = quaternion_new(0,0,sin(theta/2.0), cos(theta/2.0));
    else return;
    *quat = quaternion_multiply(rotateBy, *quat);
    quaternion_normalize(quat);
}

void quaternion_applyIntrinsicZXZRotation(Quaternion *quat, double p1, double P, double p2){
    quaternion_rotateAboutAxisByDegrees(quat, p2, 'z');
    quaternion_rotateAboutAxisByDegrees(quat, P, 'x');
    quaternion_rotateAboutAxisByDegrees(quat, p1, 'z');
}

void quaternion_rotateVector3AxisAngle(Vector3 *v, double x, double y, double z, double deg){
    double theta, st, axisMagnitude;
    Vector3 axis = vector3_new(x, y, z);
    Quaternion rot, rotNeg, vec;

    theta = deg * M_PI / 180.0;
    st = sin(theta / 2.0);
    axisMagnitude = sqrt(x*x + y*y + z*z);
    axis = vector3_scale(axis, 1.0/axisMagnitude);

    rot = quaternion_new(axis.x * st, axis.y * st, axis.z * st, cos(theta / 2.0));
    quaternion_normalize(&rot);
    vec = quaternion_new(v->x, v->y, v->z, 0);
    rotNeg = quaternion_new(-rot.x, -rot.y, -rot.z, rot.w);

    vec = quaternion_multiply(rot, vec);
    vec = quaternion_multiply(vec, rotNeg);

    v->x = vec.x;
    v->y = vec.y;
    v->z = vec.z;
}

void quaternion_rotateVector3ByQuat(Vector3 *v, Quaternion q){
    Quaternion conj, vecQuat;
    conj = quaternion_new(-q.x, -q.y, -q.z, q.w);
    vecQuat = quaternion_new(v->x, v->y, v->z, 0);

    vecQuat = quaternion_multiply(q, vecQuat);
    vecQuat = quaternion_multiply(vecQuat, conj);

    v->x = vecQuat.x;
    v->y = vecQuat.y;
    v->z = vecQuat.z;
}
