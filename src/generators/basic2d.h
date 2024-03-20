#ifndef ENGINE_BASIC2D_H
#define ENGINE_BASIC2D_H

#include <vector>
#include <cmath>
#include "../easy_image.h"
#include "../color/color.h"

img::EasyImage createColorRectangle(int &width, int &height);

img::EasyImage
createBlocks(int &imageWidth, int &imageHeight, int &blocksInX, int &blocksInY,
             std::vector<double> &colorWhite,
             std::vector<double> &colorBlack, bool &invertColors);

img::EasyImage createQuarterCircle(int &imageWidth, int &imageHeight, int &linesNumber,
                                   std::vector<double> &backgroundColor,
                                   std::vector<double> &lineColor);

img::EasyImage createEye(int &imageWidth, int &imageHeight, int &linesNumber,
                         std::vector<double> &backgroundColor,
                         std::vector<double> &lineColor);

img::EasyImage createDiamond(int &imageWidth, int &imageHeight, int &linesNumber,
                             std::vector<double> &backgroundColor,
                             std::vector<double> &lineColor);

#endif //ENGINE_BASIC2D_H
