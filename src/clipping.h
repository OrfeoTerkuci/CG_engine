#ifndef ENGINE_CLIPPING_H
#define ENGINE_CLIPPING_H

#include "geometry/shapes.h"
#include "vector/vector3d.h"
#include <vector>

void clipTriangleOneOutNear(Vector3D &a, Vector3D &b, Vector3D &c, int &indB, int &indC,
                            double &dVal, std::vector<Vector3D> &newPoints,
                            std::vector<face> &newFaces);

void clipTriangleTwoOutNear(Vector3D &a, Vector3D &b, Vector3D &c, int &indA, double &dVal,
                            std::vector<Vector3D> &newPoints, std::vector<face> &newFaces);

void clipTriangleOneOutFar(Vector3D &a, Vector3D &b, Vector3D &c, int &indB, int &indC,
                           double &dVal, std::vector<Vector3D> &newPoints,
                           std::vector<face> &newFaces);

void clipTriangleTwoOutFar(Vector3D &a, Vector3D &b, Vector3D &c, int &indA, double &dVal,
                           std::vector<Vector3D> &newPoints, std::vector<face> &newFaces);

void clipNear(figure *&originalFigure, std::vector<face> &newFaces, double &dNear);

void clipFar(figure *&originalFigure, std::vector<face> &newFaces, double &dFar);

void clipTriangleOneOutLeft(Vector3D &a, Vector3D &b, Vector3D &c, int &indB, int &indC,
                            double &dVal, double &dNear,
                            std::vector<Vector3D> &newPoints,
                            std::vector<face> &newFaces);

void clipTriangleTwoOutLeft(Vector3D &a, Vector3D &b, Vector3D &c, int &indA, double &dVal,
                            double &dNear, std::vector<Vector3D> &newPoints,
                            std::vector<face> &newFaces);

void clipLeft(figure *&originalFigure, std::vector<face> &newFaces, double &left,
              double &dNear);

void clipRight(figure *&originalFigure, std::vector<face> &newFaces, double &right,
               double &dNear);

void clipTriangleOneOutTop(Vector3D &a, Vector3D &b, Vector3D &c, int &indB, int &indC,
                           double &dVal, double &dNear,
                           std::vector<Vector3D> &newPoints,
                           std::vector<face> &newFaces);

void clipTriangleTwoOutTop(Vector3D &a, Vector3D &b, Vector3D &c, int &indA, double &dVal,
                           double &dNear, std::vector<Vector3D> &newPoints,
                           std::vector<face> &newFaces);

void clipTop(figure *&originalFigure, std::vector<face> &newFaces, double &top,
             double &dNear);

void clipBottom(figure *&originalFigure, std::vector<face> &newFaces, double &bottom,
                double &dNear);

void triangulate(face &originalFace, std::vector<face> &newFaces);

void clipView(figures3D &originalFigures, double &dNear, double &dFar, double &hFov,
              double &aspectRatio);

#endif // ENGINE_CLIPPING_H
