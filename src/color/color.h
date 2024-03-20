#ifndef ENGINE_COLOR_H
#define ENGINE_COLOR_H

#include <vector>
#include <cmath>

class color {
public:
    double red;
    double green;
    double blue;

    /**
     * @brief Constructs a color type object
     * @param red Value of red component
     * @param green Value of green component
     * @param blue Value of blue component
     */
    color(double red, double green, double blue) : red(red), green(green),
                                                   blue(blue) {};

    explicit color(std::vector<double> newColor) : red(newColor.at(0)), green(newColor.at(1)), blue(newColor.at(2)) {};

    color() : red(0), green(0), blue(0) {}

    color operator+(const color &ref) const;

    color operator*(const color &ref) const;

    color operator*(double d) const;

    virtual ~color() = default;
};

std::vector<unsigned int> scaleColors(std::vector<double> &originalColor);

#endif //ENGINE_COLOR_H
