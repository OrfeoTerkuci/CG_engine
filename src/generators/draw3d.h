#ifndef ENGINE_DRAW3D_H
#define ENGINE_DRAW3D_H

#include "../geometry/shapes.h"
#include "../ini_configuration.h"
#include "../transformations.h"
#include "../vector/vector3d.h"
#include <vector>

figure *
drawLineDrawing(double &scale, double &rotX, double &rotY, double &rotZ, int &nrPoints,
                int &nrLines,
                const ini::Configuration &configuration, std::vector<double> &lineColor,
                std::vector<double> &center,
                Matrix &mEye, int &i);

#endif // ENGINE_DRAW3D_H
