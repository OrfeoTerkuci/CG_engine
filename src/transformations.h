#ifndef ENGINE_TRANSFORMATIONS_H
#define ENGINE_TRANSFORMATIONS_H

#include "geometry/shapes.h"
#include "vector/vector3d.h"
#include <cmath>

Matrix scaleFigure(double scale);

Matrix rotateX(double angle);

Matrix rotateY(double angle);

Matrix rotateZ(double angle);

Matrix translate(const Vector3D &vector);

void toPolar(const Vector3D &point, double &theta, double &phi, double &r);

Matrix eyePointTrans(const Vector3D &eyePoint, const Vector3D &viewDir, double &theta,
                     double &phi, double &r);

point2D doProjection(const Vector3D &point, double d);

void getLinePointIndex(face &face, figure *&f, lines2D &lines, color &lineColor,
                       double d);

lines2D doProjection(figures3D &figs, double d);

#endif // ENGINE_TRANSFORMATIONS_H
