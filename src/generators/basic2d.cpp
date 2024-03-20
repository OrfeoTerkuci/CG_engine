#include "basic2d.h"


img::EasyImage createColorRectangle(int &width, int &height) {
    img::EasyImage image(width, height);
    for (unsigned int i = 0; i < width; i++) {
        for (unsigned int j = 0; j < height; j++) {
            image(i, j).red = i;
            image(i, j).green = j;
            image(i, j).blue = (i + j) % width;
        }
    }
    return image;
}

img::EasyImage
createBlocks(int &imageWidth, int &imageHeight, int &blocksInX, int &blocksInY,
             std::vector<double> &colorWhite,
             std::vector<double> &colorBlack, bool &invertColors) {
    // Get block width
    int blockWidth = imageWidth / blocksInX;
    // Get block height
    int blockHeight = imageHeight / blocksInY;

    // Scale colors
    std::vector<unsigned int> newColorWhite = scaleColors(colorWhite);
    std::vector<unsigned int> newColorBlack = scaleColors(colorBlack);
    // Create image
    img::EasyImage image(imageWidth, imageHeight);
    // Loop through each pixel
    for (unsigned int i = 0; i < imageWidth; i++) {
        for (unsigned int j = 0; j < imageHeight; j++) {
            // Calculate block in which the pixel sits
            int coordX = int(i / blockWidth);
            int coordY = int(j / blockHeight);

            // Fill block for white
            if ((coordX + coordY) % 2 == invertColors) {
                image(i, j).red = newColorBlack.at(0);
                image(i, j).green = newColorBlack.at(1);
                image(i, j).blue = newColorBlack.at(2);
            }
                // Fill block for black
            else {
                image(i, j).red = newColorWhite.at(0);
                image(i, j).green = newColorWhite.at(1);
                image(i, j).blue = newColorWhite.at(2);
            }

        }
    }
    return image;
}

img::EasyImage createQuarterCircle(int &imageWidth, int &imageHeight, int &linesNumber,
                                   std::vector<double> &backgroundColor,
                                   std::vector<double> &lineColor) {
    // Scale height
    int heightScale = imageHeight / (linesNumber - 1);
    // Create image
    img::EasyImage image(imageWidth, imageHeight);
    // Scale colors
    std::vector<unsigned int> newBgColor = scaleColors(backgroundColor);
    std::vector<unsigned int> newLineColor = scaleColors(lineColor);
    // Set background color
    image.clear(img::Color(newBgColor.at(0), newBgColor.at(1), newBgColor.at(2)));
    // Loop through image height
    for (unsigned int i = 0; i <= imageHeight; i += heightScale) {
        if (i == imageHeight) {
            image.draw_line(0, imageHeight - 1, imageWidth - 1, imageHeight - 1,
                            img::Color(newLineColor.at(0), newLineColor.at(1),
                                       newLineColor.at(2)));
        } else {
            image.draw_line(0, i, i, imageHeight - 1,
                            img::Color(newLineColor.at(0), newLineColor.at(1),
                                       newLineColor.at(2)));
        }
    }
    return image;
}

img::EasyImage createEye(int &imageWidth, int &imageHeight, int &linesNumber,
                         std::vector<double> &backgroundColor,
                         std::vector<double> &lineColor) {
    // Scale height
    int heightScale = imageHeight / (linesNumber - 1);
    // Create image
    img::EasyImage image(imageWidth, imageHeight);
    // Scale colors
    std::vector<unsigned int> newBgColor = scaleColors(backgroundColor);
    std::vector<unsigned int> newLineColor = scaleColors(lineColor);
    // Set background color
    image.clear(img::Color(newBgColor.at(0), newBgColor.at(1), newBgColor.at(2)));

    for (unsigned int i = 0; i <= imageHeight; i += heightScale) {
        if (i == imageHeight) {
            image.draw_line(0, i - 1, i - 1, i - 1,
                            img::Color(newLineColor.at(0), newLineColor.at(1),
                                       newLineColor.at(2)));
        }
        if (i == imageWidth) {
            image.draw_line(i - 1, 0, i - 1, imageHeight - 1,
                            img::Color(newLineColor.at(0), newLineColor.at(1),
                                       newLineColor.at(2)));
        } else {
            image.draw_line(0, i, i, imageWidth - 1,
                            img::Color(newLineColor.at(0), newLineColor.at(1),
                                       newLineColor.at(2)));
            image.draw_line(i, 0, imageWidth - 1, i,
                            img::Color(newLineColor.at(0), newLineColor.at(1),
                                       newLineColor.at(2)));
        }
    }

    return image;
}

img::EasyImage createDiamond(int &imageWidth, int &imageHeight, int &linesNumber,
                             std::vector<double> &backgroundColor,
                             std::vector<double> &lineColor) {
    // Scale width
    int widthScale = imageWidth / (linesNumber - 1);
    // Scale height
    int heightScale = imageHeight / (linesNumber - 1);
    // Find middle point
    int midX = imageWidth / 2;
    int midY = imageHeight / 2;
    // Create image
    img::EasyImage image(imageWidth, imageHeight);
    // Scale colors
    std::vector<unsigned int> newBgColor = scaleColors(backgroundColor);
    std::vector<unsigned int> newLineColor = scaleColors(lineColor);
    // Set background color
    image.clear(img::Color(newBgColor.at(0), newBgColor.at(1), newBgColor.at(2)));
    // Loop through image height
    for (unsigned int i = 0; i <= midY; i += heightScale / 2) {
        image.draw_line(midX, i, midX - 1 + i, midX - 1,
                        img::Color(newLineColor.at(0), newLineColor.at(1),
                                   newLineColor.at(2)));
        image.draw_line(midX, i, midX - i, midX - 1,
                        img::Color(newLineColor.at(0), newLineColor.at(1),
                                   newLineColor.at(2)));
    }
    // Loop through image width
    for (unsigned int i = 0; i <= midX; i += widthScale / 2) {
        image.draw_line(i, midY, midX - 1, midY - 1 + i,
                        img::Color(newLineColor.at(0), newLineColor.at(1),
                                   newLineColor.at(2)));
        image.draw_line(imageWidth - 1 - i, midY, midX - 1, midY - 1 + i,
                        img::Color(newLineColor.at(0), newLineColor.at(1),
                                   newLineColor.at(2)));
    }

    return image;
}