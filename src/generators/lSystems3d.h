#ifndef ENGINE_LSYSTEMS3D_H
#define ENGINE_LSYSTEMS3D_H

#include "../geometry/shapes.h"
#include "../l_parser/l_parser.h"
#include "../vector/vector3d.h"
#include <cmath>
#include <vector>

std::string
getEndString3D(const LParser::LSystem3D &lSystem, std::string &startingString,
               std::string &endingString);

figure *
createSystemLines3D(const LParser::LSystem3D &lSystem, std::string &startingString,
                    double &angle, std::vector<double> &lineColor,
                    Vector3D &currentPoint, int currentC);

figure *
drawSystem3D(const LParser::LSystem3D &lSystem, std::vector<double> &lineColor);

#endif // ENGINE_LSYSTEMS3D_H
