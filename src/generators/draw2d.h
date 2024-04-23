#ifndef ENGINE_DRAW2D_H
#define ENGINE_DRAW2D_H

#include "../color/color.h"
#include "../easy_image.h"
#include "../geometry/shapes.h"
#include "../l_parser/l_parser.h"
#include <cmath>
#include <vector>

img::EasyImage draw2DLines(const lines2D &lines, const int size, std::vector<double> &backgroundColor);

std::string getEndString(const LParser::LSystem2D &l_system, std::string &startingString, std::string &endingString);

lines2D
createSystemLines(const LParser::LSystem2D &l_system, lines2D &lines, std::string &startingString, std::string &endingString,
                  double &startingAngle, double &angle, std::vector<double> &lineColor, point2D &currentPoint,
                  int current_c);

lines2D drawSystem2D(const LParser::LSystem2D &l_system, const int &size, std::vector<double> &backgroundColor,
                     std::vector<double> &lineColor);

#endif // ENGINE_DRAW2D_H
