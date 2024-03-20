#include "lSystems3d.h"

std::string
getEndString3D(const LParser::LSystem3D &lSystem, std::string &startingString,
               std::string &endingString) {
    // Replace symbols
    for (char c: startingString) {
        // Add the operators
        if (c == '-' || c == '+' || c == '(' || c == ')' || c == '^' || c == '&' ||
            c == '/' || c == '\\' || c == '|') {
            endingString += c;
        }
            // Replace the string
        else {
            endingString += lSystem.get_replacement(c);
        }
    }
    // Update startingString
    startingString = endingString;
    endingString = "";
    return startingString;
}

figure *
createSystemLines3D(const LParser::LSystem3D &lSystem, std::string &startingString,
                    double &angle, std::vector<double> &lineColor,
                    Vector3D &currentPoint, int currentC) {
    // Create variables to store the current x and y coordinate
    std::vector<double> currentX;
    std::vector<double> currentY;
    std::vector<double> currentZ;
    // Create new figure
    figure *newFigure;
    // Create points and faces vectors
    std::vector<Vector3D> points;
    std::vector<face> faces;
    // Add the origin to the points vector
    points.push_back(Vector3D::point(currentPoint));
    // Create orientation vectors
    Vector3D h = Vector3D::point(1, 0, 0);
    Vector3D l = Vector3D::point(0, 1, 0);
    Vector3D u = Vector3D::point(0, 0, 1);
    // Create vectors to update the orientations
    Vector3D hNew;
    Vector3D lNew;
    Vector3D uNew;
    // Create containers for orientation vectors
    std::vector<Vector3D> currentH;
    std::vector<Vector3D> currentL;
    std::vector<Vector3D> currentU;
    std::vector<int> currentIndex;
    // Loop through characters in initiating string
    int count = 0;
    for (int i = currentC; i < startingString.length(); i++) {
        // Rudder left
        if (startingString[i] == '+') {
            // Calculate the new orientations
            hNew = (h * cos(angle) + l * sin(angle));
            lNew = (l * cos(angle) - h * sin(angle));
            // Update the orientations
            h = hNew;
            l = lNew;
        }
            // Rudder right
        else if (startingString[i] == '-') {
            // Calculate the new orientations
            hNew = (h * cos(-angle) + l * sin(-angle));
            lNew = (l * cos(-angle) - h * sin(-angle));
            // Update the orientations
            h = hNew;
            l = lNew;
        }
            // Pitch up
        else if (startingString[i] == '^') {
            // Calculate the new orientations
            hNew = (h * cos(angle) + u * sin(angle));
            uNew = (u * cos(angle) - h * sin(angle));
            // Update the orientations
            h = hNew;
            u = uNew;
        }
            // Pitch down
        else if (startingString[i] == '&') {
            // Calculate the new orientations
            hNew = (h * cos(-angle) + u * sin(-angle));
            uNew = (u * cos(-angle) - h * sin(-angle));
            // Update the orientations
            h = hNew;
            u = uNew;
        }
            // Roll left
        else if (startingString[i] == '\\') {
            // Calculate the new orientations
            lNew = (l * cos(angle) - u * sin(angle));
            uNew = (l * sin(angle) + u * cos(angle));
            // Update the orientations
            l = lNew;
            u = uNew;
        }
            // Roll right
        else if (startingString[i] == '/') {
            // Calculate the new orientations
            lNew = (l * cos(-angle) - u * sin(-angle));
            uNew = (l * sin(-angle) + u * cos(-angle));
            // Update the orientations
            l = lNew;
            u = uNew;
        }
            // Turn around
        else if (startingString[i] == '|') {
            h = -h;
            l = -l;
        } else if (startingString[i] == '(') {
            // Save coordinates, draw everything within the bracket, return to saved coordinates
            currentX.push_back(currentPoint.x);
            currentY.push_back(currentPoint.y);
            currentZ.push_back(currentPoint.z);
            currentH.push_back(h);
            currentL.push_back(l);
            currentU.push_back(u);
            currentIndex.push_back(count);
        } else if (startingString[i] == ')') {

            currentPoint.x = currentX.back();
            currentX.pop_back();

            currentPoint.y = currentY.back();
            currentY.pop_back();

            currentPoint.z = currentZ.back();
            currentZ.pop_back();

            h = currentH.back();
            currentH.pop_back();

            l = currentL.back();
            currentL.pop_back();

            u = currentU.back();
            currentU.pop_back();

            count = currentIndex.back();
            currentIndex.pop_back();
        }
            // If we must draw
        else if (lSystem.draw(startingString[i])) {
            // Create new points
            currentPoint += h;
            Vector3D newPoint = Vector3D::point(currentPoint);
            points.push_back(newPoint);
            // Create new faces (the lines)
            faces.push_back(face({count, (int) points.size() - 1}));
            count = (int) points.size() - 1;
        }
            // If we mustn't draw
        else {
            currentPoint += h;
            Vector3D newPoint = Vector3D::point(currentPoint);
            points.push_back(currentPoint);
        }

    }
    // Create the new figure
    newFigure = new figure(points, faces,
                           color(lineColor.at(0), lineColor.at(1), lineColor.at(2)));
    return newFigure;
}

figure *
drawSystem3D(const LParser::LSystem3D &lSystem, std::vector<double> &lineColor) {
    // Get all the components of the LSystem3D
    double angle = lSystem.get_angle();

    // Convert angles to radians
    angle = (angle * M_PI) / 180;

    // Get iterations
    unsigned int nrIterations = lSystem.get_nr_iterations();

    // Get initiating string
    std::string startingString = lSystem.get_initiator();
    std::string endingString;
    // Replace symbols
    for (int i = 0; i < nrIterations; i++) {
        startingString = getEndString3D(lSystem, startingString, endingString);
    }
    Vector3D currentPoint = Vector3D::point(0, 0, 0);
    return createSystemLines3D(lSystem, startingString, angle, lineColor, currentPoint,
                               0);
}
