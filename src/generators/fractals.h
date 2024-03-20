#ifndef ENGINE_FRACTALS_H
#define ENGINE_FRACTALS_H

#include <vector>
#include <map>
#include "platonicBodies.h"

void
insertToPentagon(int indexA, int indexD, int indexI,
                 std::map<int, std::vector<int> > &pentagons);

int isElementOf(std::vector<Vector3D> &points, Vector3D &p);

void splitTriangleHexagon(face &originalTriangle, figure *&originalFigure,
                          std::vector<face> &newFaces,
                          std::map<int, std::vector<int> > &pentagons,
                          std::vector<Vector3D> &newPoints);

figure *
createBuckyBall(std::vector<double> &ambientCoefficient,
                std::vector<double> &diffuseCoefficient,
                std::vector<double> &specularCoefficient,
                double &reflectionCoefficient);

figures3D
createMengerSponge(int nrIterations, Matrix &m,
                   std::vector<double> &ambientCoefficient,
                   std::vector<double> &diffuseCoefficient,
                   std::vector<double> &specularCoefficient,
                   double &reflectionCoefficient);

void generateFractal(figures3D &fractal, int nrIterations, double scale);

#endif //ENGINE_FRACTALS_H
