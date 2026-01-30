#ifndef ENGINE_WIREFRAMES_H
#define ENGINE_WIREFRAMES_H

#include "../clipping.h"
#include "../geometry/shapes.h"
#include "../ini_configuration.h"
#include "../l_parser/l_parser.h"
#include "../lighting/light.h"
#include "../transformations.h"
#include "draw3d.h"
#include "fractals.h"
#include "lSystems3D.h"
#include "platonicBodies.h"
#include "wireframes.h"
#include <iostream>
#include <vector>

#include <fstream>

figures3D drawWireframe(std::vector<double> &eye, int &nrFigures,
                        const ini::Configuration &configuration, lights3D &lights);

void getViewDir(bool viewCone, const Vector3D &eyePoint,
                const ini::Configuration &configuration, Vector3D &viewDir,
                double &dNear, double &dFar, double &hFov,
                double &aspectRatio);

#endif // ENGINE_WIREFRAMES_H
