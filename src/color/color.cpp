#include "color.h"

color color::operator+(const color &ref) const {
    color col;
    // Sum up the Color components
    col.red = this->red + ref.red;
    col.green = this->green + ref.green;
    col.blue = this->blue + ref.blue;
    return col;
}

color color::operator*(const color &ref) const {
    color col;
    col.red = this->red * ref.red;
    col.green = this->green * ref.green;
    col.blue = this->blue * ref.blue;
    return col;
}

color color::operator*(const double d) const {
    color col;
    col.red = this->red * d;
    col.green = this->green * d;
    col.blue = this->blue * d;
    return col;
}

std::vector<unsigned int> scaleColors(std::vector<double> &originalColor) {
    std::vector<unsigned int> newColors;
    newColors.reserve(originalColor.size());
    for (double c: originalColor) {
        unsigned int newColor = lround(c * 255);
        newColors.push_back(newColor);
    }
    return newColors;
}