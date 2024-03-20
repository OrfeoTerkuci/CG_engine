#include "draw2d.h"

img::EasyImage draw2DLines(const lines2D &lines, const int size, std::vector<double> &backgroundColor) {
    // Declare colors vector
    std::vector<double> originalColors;
    std::vector<unsigned int> newColors;
    // Scale colors
    std::vector<unsigned int> bgColor = scaleColors(backgroundColor);
    // Determine xMin , yMin , xMax , yMax
    double xMin = 0, yMin = 0, xMax = 0, yMax = 0;
    // Loop through all the lines
    for (const line2D& l: lines) {
        // Assign the points to temp variables
        double p1_x = l.p1.x;
        double p1_y = l.p1.y;
        double p2_x = l.p2.x;
        double p2_y = l.p2.y;
        // Check coordinates of the first point
        if (xMax <= p1_x) {
            xMax = p1_x;
        }
        if (xMin >= p1_x) {
            xMin = p1_x;
        }
        if (yMax <= p1_y) {
            yMax = p1_y;
        }
        if (yMin >= p1_y) {
            yMin = p1_y;
        }
        // Check coordinates of the second point
        if (xMax <= p2_x) {
            xMax = p2_x;
        }
        if (xMin >= p2_x) {
            xMin = p2_x;
        }
        if (yMax <= p2_y) {
            yMax = p2_y;
        }
        if (yMin >= p2_y) {
            yMin = p2_y;
        }
    }

    // Calculate range along the x-axis and y-axis
    double x_range = xMax - xMin;
    double y_range = yMax - yMin;
    // Calculate image dimensions
    double imageWidth = size * (x_range / std::max(x_range, y_range));
    double imageHeight = size * (y_range / std::max(x_range, y_range));
    // Create the image file
    img::EasyImage image(lround(imageWidth), lround(imageHeight));
    image.clear(img::Color(bgColor.at(0), bgColor.at(1), bgColor.at(2)));
    // Determine the scaling factor
    double d = 0.95 * (imageWidth / x_range);
    // Calculate for x and y
    double DC_x = d * ((xMin + xMax) / 2);
    double DC_y = d * ((yMin + yMax) / 2);
    double dx = imageWidth / 2 - DC_x;
    double dy = imageHeight / 2 - DC_y;
    // Declare temporary variables
    //int x1 , x2 , y1 , y2;
    // Loop through the lines again
    for (line2D l: lines) {
        // Multiply all the points by the scaling factor
        l.p1.x *= d;
        l.p1.y *= d;
        l.p2.x *= d;
        l.p2.y *= d;
        // Add dx and dy to each point's coordinate
        l.p1.x += dx;
        l.p1.y += dy;
        l.p2.x += dx;
        l.p2.y += dy;
        // Round the coordinates
        int x1 = lround(l.p1.x);
        int x2 = lround(l.p2.x);
        int y1 = lround(l.p1.y);
        int y2 = lround(l.p2.y);
        // Fetch and rescale the colors
        originalColors = {l.lineColor.red, l.lineColor.green, l.lineColor.blue};
        newColors = scaleColors(originalColors);
        // Draw the lines
        image.draw_line(x1, y1, x2, y2, img::Color(newColors.at(0), newColors.at(1), newColors.at(2)));
    }
    return image;
}

std::string getEndString(const LParser::LSystem2D &l_system, std::string &startingString, std::string &endingString) {
    // Replace symbols
    for (char c: startingString) {
        // Add the operators
        if (c == '-' || c == '+' || c == '(' || c == ')') {
            endingString += c;
        }
            // Replace the string
        else {
            endingString += l_system.get_replacement(c);
        }
    }
    // Update startingString
    startingString = endingString;
    endingString = "";
    return startingString;
}

lines2D
createSystemLines(const LParser::LSystem2D &l_system, lines2D &lines, std::string &startingString, std::string &endingString,
                  double &startingAngle, double &angle, std::vector<double> &lineColor, point2D &currentPoint,
                  int current_c) {
    // Create variables to store the current x and y coordinate
    std::vector<double> current_x;
    std::vector<double> current_y;
    std::vector<double> current_angle;
    // Loop through characters in initiating string
    for (int i = current_c; i < startingString.length(); i++) {
        // If angle must increase
        if (startingString[i] == '+') {
            startingAngle += angle;
        }
            // If angle must decrease
        else if (startingString[i] == '-') {
            startingAngle -= angle;
        } else if (startingString[i] == '(') {
            // Save coordinates, draw everything within the bracket, return to saved coordinates
            current_x.push_back(currentPoint.x);
            current_y.push_back(currentPoint.y);
            current_angle.push_back(startingAngle);
        } else if (startingString[i] == ')') {
            currentPoint.x = current_x[current_x.size() - 1];
            current_x.pop_back();
            currentPoint.y = current_y[current_y.size() - 1];
            current_y.pop_back();
            startingAngle = current_angle[current_angle.size() - 1];
            current_angle.pop_back();
        }
            // If we must draw
        else if (l_system.draw(startingString[i])) {
            line2D line(point2D(currentPoint.x, currentPoint.y),
                        point2D(currentPoint.x + cos(startingAngle), currentPoint.y + sin(startingAngle)),
                        color(lineColor[0], lineColor[1], lineColor[2])
            );
            lines.push_back(line);
            currentPoint.x += cos(startingAngle);
            currentPoint.y += sin(startingAngle);
        }
            // If we mustn't draw
        else {
            currentPoint.x += cos(startingAngle);
            currentPoint.y += sin(startingAngle);
        }
    }
    return lines;
}

lines2D drawSystem2D(const LParser::LSystem2D &l_system, const int &size, std::vector<double> &backgroundColor,
                     std::vector<double> &lineColor) {
    // Create the list of lines
    lines2D lines;

    // Get all the components of the LSystem
    double angle = l_system.get_angle();
    double startingAngle = l_system.get_starting_angle();

    // Convert angles to radians
    angle = (angle * M_PI) / 180;
    startingAngle = (startingAngle * M_PI) / 180;

    // Get iterations
    unsigned int nr_iterations = l_system.get_nr_iterations();

    // Get initiating string
    std::string startingString = l_system.get_initiator();
    std::string endingString;
    // Replace symbols
    for (int i = 0; i < nr_iterations; i++) {
        startingString = getEndString(l_system, startingString, endingString);
    }
    point2D currentPoint(0, 0);
    return createSystemLines(l_system, lines, startingString, endingString, startingAngle, angle, lineColor,
                             currentPoint, 0);
}