#ifndef ENGINE_PLATONICBODIES_H
#define ENGINE_PLATONICBODIES_H

#include <vector>
#include "../transformations.h"


figure *createCube(std::vector<double> &ambientCoefficient,
                   std::vector<double> &diffuseCoefficient,
                   std::vector<double> &specularCoefficient,
                   double &reflectionCoefficient);

figure *
createTetrahedron(std::vector<double> &ambientCoefficient,
                  std::vector<double> &diffuseCoefficient,
                  std::vector<double> &specularCoefficient,
                  double &reflectionCoefficient);

figure *
createOctahedron(std::vector<double> &ambientCoefficient,
                 std::vector<double> &diffuseCoefficient,
                 std::vector<double> &specularCoefficient,
                 double &reflectionCoefficient);

figure *
createIcosahedron(std::vector<double> &ambientCoefficient,
                  std::vector<double> &diffuseCoefficient,
                  std::vector<double> &specularCoefficient,
                  double &reflectionCoefficient);

figure *
createDodecahedron(std::vector<double> &ambientCoefficient,
                   std::vector<double> &diffuseCoefficient,
                   std::vector<double> &specularCoefficient,
                   double &reflectionCoefficient);

void splitTriangle(face &originalTriangle, figure *&originalFigure,
                   std::vector<face> &newFaces);

figure *
createSphere(__attribute__((unused)) const double &radius, const int &n,
             std::vector<double> &ambientCoefficient,
             std::vector<double> &diffuseCoefficient,
             std::vector<double> &specularCoefficient, double &reflectionCoefficient);

figure *createCone(const int &n, const double &height,
                   std::vector<double> &ambientCoefficient,
                   std::vector<double> &diffuseCoefficient,
                   std::vector<double> &specularCoefficient,
                   double &reflectionCoefficient);

figure *createCylinder(const int &n, const double &height,
                       std::vector<double> &ambientCoefficient,
                       std::vector<double> &diffuseCoefficient,
                       std::vector<double> &specularCoefficient,
                       double &reflectionCoefficient);

figure *createTorus(const double &r, const double &rBig, const int &n, const int &m,
                    std::vector<double> &ambientCoefficient,
                    std::vector<double> &diffuseCoefficient,
                    std::vector<double> &specularCoefficient,
                    double &reflectionCoefficient);

#endif //ENGINE_PLATONICBODIES_H
