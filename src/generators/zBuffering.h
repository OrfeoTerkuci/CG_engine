#ifndef ENGINE_ZBUFFERING_H
#define ENGINE_ZBUFFERING_H

#include "../clipping.h"
#include "../easy_image.h"
#include "../lighting/light.h"
#include "../transformations.h"
#include <algorithm>
#include <cassert>
#include <limits>
#include <vector>

const double posInf = std::numeric_limits<double>::infinity();
const double negInf = -std::numeric_limits<double>::infinity();

class zBuffer : public std::vector<std::vector<double>>
{
public:
    zBuffer(int width, int height);
};

void getProjectedPoints(figures3D &figs, double &d, double &xRange, double &yRange,
                        double &xMin, double &yMin, double &xMax, double &yMax);

double
calculateInvZ(unsigned int i, unsigned int iMin, unsigned int iMax, const double &z0,
              const double &z1);

void drawZBufferLine(zBuffer &zBuffer, img::EasyImage &image,
                     unsigned int x0, unsigned int y0, double &z0,
                     unsigned int x1, unsigned int y1, double &z1,
                     const color &lineColor);

img::EasyImage draw2DzBufferLines(const lines2D &lines, int size,
                                  std::vector<double> &backgroundColor);

double calculateIntersection(const int &yI, point2D const &p, point2D const &q);

double getAngle(const Vector3D &a, const Vector3D &b, const Vector3D &c,
                const Vector3D &direction);

void drawZBufferTriangle(zBuffer &buffer, img::EasyImage &image,
                         Vector3D const &a, Vector3D const &b, Vector3D const &c,
                         double &d, double &dx, double &dy,
                         color &ambientReflection, color &diffuseReflection,
                         color &specularReflection,
                         double &reflectionCoefficient,
                         lights3D &lights);

img::EasyImage
draw2DzBufferTriangle(const int &size, std::vector<double> &backgroundColor,
                      figures3D &figures,
                      lights3D &lights);

#endif // ENGINE_ZBUFFERING_H
