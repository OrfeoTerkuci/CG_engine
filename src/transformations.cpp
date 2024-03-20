#include "transformations.h"

Matrix scaleFigure(const double scale) {
    Matrix s;
    s(1, 1) = scale;
    s(2, 2) = scale;
    s(3, 3) = scale;
    s(4, 4) = 1;
    return s;
}

Matrix rotateX(const double angle) {
    Matrix mX;
    mX(2, 2) = cos(angle);
    mX(2, 3) = sin(angle);
    mX(3, 2) = -sin(angle);
    mX(3, 3) = cos(angle);
    return mX;
}

Matrix rotateY(const double angle) {
    Matrix mY;
    mY(1, 1) = cos(angle);
    mY(1, 3) = -sin(angle);
    mY(3, 1) = sin(angle);
    mY(3, 3) = cos(angle);
    return mY;
}

Matrix rotateZ(const double angle) {
    Matrix mZ;
    mZ(1, 1) = cos(angle);
    mZ(1, 2) = sin(angle);
    mZ(2, 1) = -sin(angle);
    mZ(2, 2) = cos(angle);
    return mZ;
}

Matrix translate(const Vector3D &vector) {
    Matrix t;
    t(4, 1) = vector.x;
    t(4, 2) = vector.y;
    t(4, 3) = vector.z;
    return t;
}

void toPolar(const Vector3D &point, double &theta, double &phi, double &r) {
    // Get the points
    double x = point.x;
    double y = point.y;
    double z = point.z;
    // Calculate r
    r = sqrt(x * x + y * y + z * z);
    // Calculate theta
    theta = atan2(y, x);
    // Calculate phi
    if (r != 0) {
        phi = acos(z / r);
    } else {
        phi = 0;
    }
}

Matrix eyePointTrans(const Vector3D &eyePoint, const Vector3D &viewDir, double &theta,
                     double &phi, double &r) {
    // Make vector
    double dR;
    toPolar(eyePoint, theta, phi, r);
    toPolar(-viewDir, theta, phi, dR);

    Matrix m;
    m(1, 1) = -sin(theta);
    m(1, 2) = -cos(theta) * cos(phi);
    m(1, 3) = cos(theta) * sin(phi);
    m(2, 1) = cos(theta);
    m(2, 2) = -sin(theta) * cos(phi);
    m(2, 3) = sin(theta) * sin(phi);
    m(3, 2) = sin(phi);
    m(3, 3) = cos(phi);
    m(4, 3) = -r;

    return m;
}

point2D doProjection(const Vector3D &point, const double d) {
    double x1 = (d * point.x) / -point.z;
    double y1 = (d * point.y) / -point.z;
    return {x1, y1};
}

void getLinePointIndex(face &face, figure *&f, lines2D &lines, color &lineColor,
                       const double d) {
    // Create new variables
    int bIndex, eIndex;

    for (unsigned long i = 0; i < face.pointIndexes.size(); ++i) {
        // Get begin and end index
        bIndex = face.pointIndexes.at(i);
        eIndex = face.pointIndexes.at(i == face.pointIndexes.size() - 1 ? 0 : i + 1);
        // Get points
        Vector3D beginP = f->points[bIndex];
        Vector3D endP = f->points[eIndex];
        // Convert Vector3D to point2D
        point2D newBeginP = doProjection(beginP, d);
        point2D newEndP = doProjection(endP, d);
        // Create new line
        line2D newLine(newBeginP, newEndP, lineColor, beginP.z, endP.z);
        //newLine.z1 = beginP.z;
        //newLine.z2 = endP.z;
        lines.push_back(newLine);
    }
}

lines2D doProjection(figures3D &figs, const double d) {
    lines2D lines;
    color color;
    for (figure *f: figs) {
        color = f->ambientReflection + f->diffuseReflection + f->specularReflection;
        for (face face: f->faces) {
            // Get points index - loop through
            getLinePointIndex(face, f, lines, color, d);
        }
    }
    return lines;
}
