#pragma clang diagnostic push
#pragma ide diagnostic ignored "bugprone-narrowing-conversions"
#include "easy_image.h"
#include "ini_configuration.h"
#include "l_parser/l_parser.h"
#include "vector/vector3d.h"
#include <fstream>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <list>
#include <algorithm>
// For VS Code
/*
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
*/
using namespace std;

// Declaring data structures
class Color {
public:
    double red;
    double green;
    double blue;

    Color(double red, double green, double blue) : red(red), green(green), blue(blue) {}

    virtual ~Color() {

    }
};

class Point2D {
public:
    double x;
    double y;

    Point2D(double x, double y) : x(x), y(y) {}

    virtual ~Point2D() {

    }
};

class Point3D {

public:
    double x;
    double y;
    double z;

    Point3D(double x, double y, double z) : x(x), y(y), z(z) {}

    virtual ~Point3D() {

    }
};

class Line2D {
public:
    Point2D p1;
    Point2D p2;
    Color color;

    Line2D(const Point2D &p1, const Point2D &p2, const Color &color) : p1(p1), p2(p2), color(color) {}

    virtual ~Line2D() {

    }
};

class Face {

public:
    // These indexes refer to points in the 'points' vector of the Figure-class
    vector<int> point_indexes;

    Face(const vector<int> &point_indexes) : point_indexes(point_indexes) {}

    virtual ~Face() {

    }
};

class Figure {

public:
    vector<Vector3D> points;
    vector<Face> faces;
    Color color;

    Figure(const vector<Vector3D> &points, const vector<Face> &faces, const Color &color) : points(points),
                                                                                            faces(faces),
                                                                                            color(color) {}
    void applyTransformation(const Matrix &m){
        // Multiply each vector with the matrix
        for(Vector3D &v : points){
            v *= m;
        }
    }

    virtual ~Figure() {

    }
};

// Declare new types

typedef vector<Line2D> Lines2D;

typedef vector<Figure*> Figures3D;

// Main functionality

// Session 0: Optional

vector<unsigned int> scaleColors(vector<double> &originalColor){
    vector<unsigned int> newColors;
    newColors.reserve(originalColor.size());
    for (double c : originalColor){
        unsigned int newColor = lround(c * 255);
        newColors.push_back(newColor);
    }
    return newColors;
}

img::EasyImage createColorRectangle(int &width, int &height){
    img::EasyImage image(width,height);
    for(unsigned int i = 0; i < width; i++)
    {
        for(unsigned int j = 0; j < height; j++)
        {
            image(i,j).red = i;
            image(i,j).green = j;
            image(i,j).blue = (i+j)%width;
        }
    }
    return image;
}

img::EasyImage createBlocks(int &imageWidth, int &imageHeight, int &blocksInX, int &blocksInY, vector<double> &colorWhite , vector<double> &colorBlack, bool &invertColors){
    // Get block width
    int blockWidth = imageWidth / blocksInX;
    // Get block height
    int blockHeight = imageHeight / blocksInY;
    // Make colors set
    enum Color {white,black};
    // Create blockColor
    Color blockColor;
    // Scale colors
    vector<unsigned int>newColorWhite = scaleColors(colorWhite);
    vector<unsigned int>newColorBlack = scaleColors(colorBlack);
    // Create image
    img::EasyImage image(imageWidth,imageHeight);
    // Loop through each pixel
    for(unsigned int i = 0; i < imageWidth; i++)
    {
        for(unsigned int j = 0; j < imageHeight; j++)
        {
            // Calculate block in which the pixel sits
            int coordX = int(floor(i / blockWidth));
            int coordY = int(floor(j / blockHeight));

            if((coordX + coordY) % 2 == 0 && !invertColors){
                blockColor = white;
            }
            else if((coordX + coordY) % 2 == 0 && invertColors){
                blockColor = black;
            }
            else if((coordX + coordY) % 2 != 0 && !invertColors){
                blockColor = black;
            }
            else{
                blockColor = white;
            }

            // Fill block for white
            if (blockColor == white){
                image(i,j).red = newColorWhite.at(0);
                image(i,j).green = newColorWhite.at(1);
                image(i,j).blue = newColorWhite.at(2);
            }
            // Fill block for black
            else {
                image(i,j).red = newColorBlack.at(0);
                image(i,j).green = newColorBlack.at(1);
                image(i,j).blue = newColorBlack.at(2);
            }

        }
    }
    return image;
}

img::EasyImage createQuarterCircle(int &imageWidth, int &imageHeight, int &linesNumber, vector<double> &backgroundColor, vector<double> &lineColor){
    // Scale height
    int heightScale = imageHeight / (linesNumber - 1);
    // Create image
    img::EasyImage image(imageWidth,imageHeight);
    // Scale colors
    vector<unsigned int> newBGColor = scaleColors(backgroundColor);
    vector<unsigned int> newLineColor = scaleColors(lineColor);
    // Set background color
    image.clear(img::Color(newBGColor.at(0),newBGColor.at(1),newBGColor.at(2)));
    // Loop through image height
    for(unsigned int i = 0; i <= imageHeight; i+=heightScale)
    {
        if(i == imageHeight){
            image.draw_line(0,imageHeight-1,imageWidth-1,imageHeight-1,img::Color(newLineColor.at(0),newLineColor.at(1),newLineColor.at(2)));
        }
        else{
            image.draw_line(0, i, i, imageHeight-1, img::Color(newLineColor.at(0),newLineColor.at(1),newLineColor.at(2)));
        }
    }
    return image;
}

img::EasyImage createEye(int &imageWidth, int &imageHeight, int &linesNumber, vector<double> &backgroundColor, vector<double> &lineColor){
    // Scale width
    int widthScale = imageWidth / (linesNumber - 1);
    // Scale height
    int heightScale = imageHeight / (linesNumber - 1);
    // Create image
    img::EasyImage image(imageWidth,imageHeight);
    // Scale colors
    vector<unsigned int> newBGColor = scaleColors(backgroundColor);
    vector<unsigned int> newLineColor = scaleColors(lineColor);
    // Set background color
    image.clear(img::Color(newBGColor.at(0),newBGColor.at(1),newBGColor.at(2)));

    for(unsigned int i = 0; i <= imageHeight; i+=heightScale)
    {
        if(i == imageHeight){
            image.draw_line(0,i-1,i-1,i-1,img::Color(newLineColor.at(0),newLineColor.at(1),newLineColor.at(2)));
        }
        if ( i == imageWidth){
            image.draw_line(i-1,0,i-1,imageHeight-1,img::Color(newLineColor.at(0),newLineColor.at(1),newLineColor.at(2)));
        }
        else{
            image.draw_line(0, i, i, imageWidth-1, img::Color(newLineColor.at(0),newLineColor.at(1),newLineColor.at(2)));
            image.draw_line(i, 0, imageWidth-1, i, img::Color(newLineColor.at(0),newLineColor.at(1),newLineColor.at(2)));
        }
    }

    return image;
}

img::EasyImage createDiamond(int &imageWidth, int &imageHeight, int &linesNumber, vector<double> &backgroundColor, vector<double> &lineColor){
    // Scale width
    int widthScale = imageWidth / (linesNumber - 1);
    // Scale height
    int heightScale = imageHeight / (linesNumber - 1);
    // Find middle point
    int midX = imageWidth / 2;
    int midY = imageHeight / 2;
    // Create image
    img::EasyImage image(imageWidth,imageHeight);
    // Scale colors
    vector<unsigned int> newBGColor = scaleColors(backgroundColor);
    vector<unsigned int> newLineColor = scaleColors(lineColor);
    // Set background color
    image.clear(img::Color(newBGColor.at(0),newBGColor.at(1),newBGColor.at(2)));
    // Loop through image height
    for(unsigned int i = 0; i <= midY; i+=heightScale/2)
    {
        image.draw_line(midX, i, midX-1+i, midX-1, img::Color(newLineColor.at(0),newLineColor.at(1),newLineColor.at(2)));
        image.draw_line(midX, i, midX-i, midX-1, img::Color(newLineColor.at(0),newLineColor.at(1),newLineColor.at(2)));
    }
    // Loop through image width
    for(unsigned int i = 0; i <= midX; i+=widthScale/2)
    {
        image.draw_line(i, midY,  midX-1, midY-1+i, img::Color(newLineColor.at(0),newLineColor.at(1),newLineColor.at(2)));
        image.draw_line(imageWidth-1-i, midY,  midX-1, midY-1+i, img::Color(newLineColor.at(0),newLineColor.at(1),newLineColor.at(2)));
    }

    return image;
}

// Session 1: L-Systems

img::EasyImage draw2DLines (const Lines2D &lines , const int size , vector<double> &backgroundColor){
    // Declare colors vector
    vector<double> originalColors;
    vector<unsigned int> newColors;
    // Scale colors
    vector<unsigned int> bgColor = scaleColors(backgroundColor);
    // Determine x_min , y_min , x_max , y_max
    double x_min = 0, y_min = 0, x_max = 0, y_max = 0;
    // Loop through all the lines
    for( Line2D l : lines) {
        // Assign the points to temp variables
        double p1_x = l.p1.x;
        double p1_y = l.p1.y;
        double p2_x = l.p2.x;
        double p2_y = l.p2.y;
        // Check coordinates of the first point
        if (x_max <= p1_x){
            x_max = p1_x;
        }
        if (x_min >= p1_x){
            x_min = p1_x;
        }
        if (y_max <= p1_y){
            y_max = p1_y;
        }
        if (y_min >= p1_y){
            y_min = p1_y;
        }
        // Check coordinates of the second point
        if (x_max <= p2_x){
            x_max = p2_x;
        }
        if (x_min >= p2_x){
            x_min = p2_x;
        }
        if (y_max <= p2_y){
            y_max = p2_y;
        }
        if (y_min >= p2_y){
            y_min = p2_y;
        }
    }

    // Calculate range along the x-axis and y-axis
    double x_range = x_max - x_min;
    double y_range = y_max - y_min;
    // Calculate image dimensions
    double imageWidth   = size * (x_range / max(x_range,y_range));
    double imageHeight  = size * (y_range / max(x_range,y_range));
    // Create the image file
    img::EasyImage image(lround(imageWidth) , lround(imageHeight));
    image.clear(img::Color(bgColor.at(0) , bgColor.at(1) , bgColor.at(2) ) );
    // Determine the scaling factor
    double d = 0.95 * (imageWidth / x_range);
    // Calculate for x and y
    double DC_x = d * ( (x_min + x_max) / 2 );
    double DC_y = d * ( (y_min + y_max) / 2 );
    double dx = imageWidth / 2 - DC_x;
    double dy = imageHeight / 2 - DC_y;
    // Declare temporary variables
    //int x1 , x2 , y1 , y2;
    // Loop through the lines again
    for (Line2D l : lines) {
        // Multiply all the points by the scaling factor
        l.p1.x  *= d;
        l.p1.y  *= d;
        l.p2.x  *= d;
        l.p2.y  *= d;
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
        originalColors = {l.color.red , l.color.green , l.color.blue};
        newColors = scaleColors(originalColors);
        // Draw the lines
        image.draw_line(x1 , y1 , x2 , y2 , img::Color(newColors.at(0) , newColors.at(1) , newColors.at(2) ) );
    }
    return image;
}

string getEndString(const LParser::LSystem2D &l_system , string &startingString , string &endingString){
    // Replace symbols
    for(char c : startingString){
        // Add the operators
        if(c == '-' || c == '+' || c == '(' || c == ')'){
            endingString += c;
        }
            // Replace the string
        else{
            endingString += l_system.get_replacement(c);
        }
    }
    // Update startingString
    startingString = endingString;
    endingString = "";
    return startingString;
}

Lines2D createSystemLines (const LParser::LSystem2D &l_system , Lines2D &lines , string &startingString , string &endingString ,
        double &startingAngle , double &angle , vector<double> &lineColor, Point2D &currentPoint , int current_c){
    // Create variables to store the current x and y coordinate
    vector<double> current_x;
    vector<double> current_y;
    vector<double> current_angle;
    // Loop through characters in initiating string
    for( int i = current_c; i < startingString.length(); i++) {
        // If angle must increase
        if (startingString[i] == '+') {
            startingAngle += angle;
        }
            // If angle must decrease
        else if (startingString[i] == '-') {
            startingAngle -= angle;
        }
        else if (startingString[i] == '('){
            // Save coordinates, draw everything within the bracket, return to saved coordinates
            current_x.push_back(currentPoint.x);
            current_y.push_back(currentPoint.y);
            current_angle.push_back(startingAngle);
        }
        else if (startingString[i] == ')'){
            currentPoint.x = current_x[current_x.size()-1];
            current_x.pop_back();
            currentPoint.y = current_y[current_y.size()-1];
            current_y.pop_back();
            startingAngle = current_angle[current_angle.size()-1];
            current_angle.pop_back();
        }
            // If we must draw
        else if(l_system.draw(startingString[i])){
            Line2D line(    Point2D(currentPoint.x , currentPoint.y) ,
                            Point2D(currentPoint.x +  cos(startingAngle) , currentPoint.y + sin(startingAngle)) ,
                            Color(lineColor[0] , lineColor[1] , lineColor[2])
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

Lines2D drawSystem2D(const LParser::LSystem2D &l_system, const int &size, vector<double> &backgroundColor,
                     vector<double> &lineColor) {
    // Create the list of lines
    Lines2D lines;

    // Get all the components of the LSystem
    const auto &Alphabet = l_system.get_alphabet();
    const string &initiator = l_system.get_initiator();
    double angle = l_system.get_angle();
    double startingAngle = l_system.get_starting_angle();

    // Convert angles to radians
    angle = ( angle * M_PI ) / 180;
    startingAngle = ( startingAngle * M_PI ) / 180;

    // Get inierations
    unsigned int nr_iterations = l_system.get_nr_iterations();

    // Get initiating string
    string startingString = l_system.get_initiator();
    string endingString;
    // Replace symbols
    for( int i = 0; i < nr_iterations; i++){
        startingString = getEndString(l_system, startingString, endingString);
    }
    Point2D currentPoint(0,0);
    return createSystemLines(l_system,lines,startingString,endingString,startingAngle,angle,lineColor,currentPoint,0);
}

// Session 2 : 3D Lines
// Transformation functions

Matrix scaleFigure(const double scale){
    Matrix S;
    S(1,1) = scale;
    S(2,2) = scale;
    S(3,3) = scale;
    S(4,4) = 1;
    return S;
}

Matrix rotateX(const double angle){
    Matrix M_x;
    M_x(2,2) = cos(angle);
    M_x(2,3) = sin(angle);
    M_x(3,2) = -sin(angle);
    M_x(3,3) = cos(angle);
    return M_x;
}

Matrix rotateY(const double angle){
    Matrix M_y;
    M_y(1,1) = cos(angle);
    M_y(1,3) = -sin(angle);
    M_y(3,1) = sin(angle);
    M_y(3,3) = cos(angle);
    return M_y;
}

Matrix rotateZ(const double angle){
    Matrix M_z;
    M_z(1,1) = cos(angle);
    M_z(1,2) = sin(angle);
    M_z(2,1) = -sin(angle);
    M_z(2,2) = cos(angle);
    return M_z;
}

Matrix translate(const Vector3D &vector){
    Matrix T;
    T(4,1) = vector.x;
    T(4,2) = vector.y;
    T(4,3) = vector.z;
    return T;
}

void toPolar(const Vector3D &point, double &theta , double &phi , double &r){
    // Get the points
    double x = point.x;
    double y = point.y;
    double z = point.z;
    // Calculate r
    r = sqrt(x * x + y * y + z * z);
    // Calculate theta
    theta = atan2(y , x);
    // Calculate phi
    if(r != 0){
        phi = acos(z/r);
    }
    else{
        phi = 0;
    }
}

Matrix eyePointTrans(const Vector3D &eyePoint , double &theta , double &phi , double &r){
    // Make vector
    Vector3D v = Vector3D::vector(0,0,-r);
    toPolar(eyePoint , theta , phi , r);

    Matrix m;
    m(1,1) = -sin(theta);
    m(1,2) = -cos(theta) * cos(phi);
    m(1,3) = cos(theta) * sin(phi);
    m(2,1) = cos(theta);
    m(2,2) = -sin(theta) * cos(phi);
    m(2,3) = sin(theta) * sin(phi);
    m(3,2) = sin(phi);
    m(3,3) = cos(phi);
    m(4,3) = -r;
    // Return matrix for 2 rotations and a translation
    return m;
}

Point2D doProjection(const Vector3D &point , const double d){
    double x_1 = (d * point.x) / -point.z;
    double y_1 = (d * point.y) / -point.z;
    return Point2D(x_1 , y_1);
}

void getLinePointIndex(Face &face, Figure* &f, Lines2D &lines){
    // Create new variables
    int b_index = 0;
    int e_index = 0;
    for (unsigned long i = 0; i < face.point_indexes.size(); ++i) {
        // Get begin and end index
        b_index = face.point_indexes.at(i);
        e_index = face.point_indexes.at(i == face.point_indexes.size() - 1 ? 0 : i + 1);
        // Get points
        Vector3D beginP = f->points[b_index];
        Vector3D endP = f->points[e_index];
        // Convert Vector3D to Point2D
        Point2D newBeginP = doProjection(beginP, 1.0);
        Point2D newEndP = doProjection(endP, 1.0);
        // Create new line
        Line2D newLine(newBeginP, newEndP, f->color);
        lines.push_back(newLine);
    }
}

Lines2D doProjection(Figures3D &figs){
    Lines2D lines;
    for(Figure* f : figs){
        for(Face face : f->faces){
            // Get points index - loop through
            getLinePointIndex(face, f, lines);
        }
    }
    return lines;
}

// Session 3: Figures - Platonic bodies

Figure* createCube(vector<double>&lineColor){
    // Points array
    double Points_T [3][8] = {
            { 1 , -1 , 1 , -1 , 1 , -1 , 1 , -1 },
            { -1 , 1 , 1 , -1 , 1 , -1 , -1 , 1 },
            {-1 , -1 , 1 , 1 , -1 , -1 , 1 , 1  }
    };
    // Faces array
    int Faces_T [4][6] = {
            { 0 , 4 , 1 , 5 , 6 , 0 },
            { 4 , 1 , 5 , 0 , 2 , 5 },
            { 2 , 7 , 3 , 6 , 7 , 1 },
            { 6 , 2 , 7 , 3 , 3 , 4 }
    };
    // Create all the points
    vector<Vector3D> points;
    for (int i = 0; i < 8; ++i) {
        points.push_back(Vector3D::point(Points_T[0][i],Points_T[1][i],Points_T[2][i]));
    }
    // Create all the faces
    vector<Face> faces;
    for (int j = 0; j < 6; ++j) {
        faces.push_back( Face( {Faces_T[0][j] , Faces_T[1][j] , Faces_T[2][j] , Faces_T[3][j]} ) );
    }
    // Create new figure
    Figure* newCube;
    newCube = new Figure( points , faces , Color( lineColor.at(0) , lineColor.at(1) , lineColor.at(2) ) );
    return newCube;
}

Figure* createTetrahedron(vector<double>&lineColor){
    // Points array
    double Points_T [3][4] = {
            { 1 , -1 , 1 , -1 },
            { -1 , 1 , 1 , -1 },
            { -1 , -1 , 1 , 1 }
    };
    // Faces array
    int Faces_T [3][4] = {
            { 0 , 2 , 0 , 0 },
            { 1 , 3 , 3 , 2 },
            { 2 , 2 , 1 , 3 }
    };
    // Create all the points
    vector<Vector3D> points;
    for (int i = 0; i < 4; ++i) {
        points.push_back(Vector3D::point(Points_T[0][i] , Points_T[1][i] , Points_T[2][i] ) );
    }
    // Create faces
    vector<Face> faces;
    for (int j = 0; j < 4; ++j) {
        faces.push_back( Face( { Faces_T[0][j] , Faces_T[1][j] , Faces_T[2][j] } ) );
    }
    // Create new figure
    Figure* newTetra;
    newTetra = new Figure( points , faces , Color( lineColor.at(0) , lineColor.at(1) , lineColor.at(2) ) );
    return newTetra;
}

Figure* createOctahedron(vector<double>&lineColor){
    // Points array
    double Points_T [3][6] = {
            { 1 , 0 ,-1 , 0 , 0 , 0 },
            { 0 , 1 , 0 , -1 , 0 , 0 },
            { 0 , 0 , 0 , 0 , -1 , 1 }
    };
    // Faces array
    int Faces_T [3][8] = {
            { 0 , 1 , 2 , 3 , 1 , 2 , 3 , 0},
            { 1 , 2 , 3 , 0 , 0 , 1 , 2 , 3},
            { 5 , 5 , 5 , 5 , 4 , 4 , 4 , 4}
    };
    // Create all the points
    vector<Vector3D> points;
    for (int i = 0; i < 6; ++i) {
        points.push_back(Vector3D::point(Points_T[0][i] , Points_T[1][i] , Points_T[2][i] ) );
    }
    // Create faces
    vector<Face> faces;
    for (int j = 0; j < 8; ++j) {
        faces.push_back( Face( { Faces_T[0][j] , Faces_T[1][j] , Faces_T[2][j] } ) );
    }
    // Create new figure
    Figure* newOcta;
    newOcta = new Figure( points , faces , Color( lineColor.at(0) , lineColor.at(1) , lineColor.at(2) ) );
    return newOcta;
}

Figure* createIcosahedron(vector<double>&lineColor){
    // Points array
    double Points_T [3][12];
    Points_T[0][0] = 0;
    Points_T[1][0] = 0;
    Points_T[2][0] = sqrt(5) / 2;

    for (int k = 1; k < 6; ++k) {
        Points_T[0][k] = cos(2 * M_PI * (k - 1) / 5 );
        Points_T[1][k] = sin(2 * M_PI * (k - 1) / 5 );
        Points_T[2][k] = 0.5;
    }

    for (int l = 6; l < 11; ++l) {
        Points_T[0][l] = cos(M_PI / 5 + (l - 6) * (2 * M_PI) / 5);
        Points_T[1][l] = sin(M_PI / 5 + (l - 6) * (2 * M_PI) / 5) ;
        Points_T[2][l] = -0.5;
    }

    Points_T[0][11] = 0;
    Points_T[1][11] = 0;
    Points_T[2][11] = -sqrt(5) / 2;

    // Faces array
    int Faces_T [3][20] = {
            { 0 , 0 , 0 , 0 , 0 , 1 , 2 , 2 , 3 , 3 , 4 , 4 , 5 , 5 , 1 , 11 , 11 , 11 , 11 , 11 },
            { 1 , 2 , 3 , 4 , 5 , 6 , 6 , 7 , 7 , 8 , 8 , 9 , 9 , 10 , 10 , 7 , 8 , 9 , 10 , 6 },
            { 2 , 3 , 4 , 5 , 1 , 2 , 7 , 3 , 8 , 4 , 9 , 5 , 10 , 1 , 6 , 6 , 7 , 8 , 9 , 10 }
    };

    // Create all the points
    vector<Vector3D> points;
    for (int i = 0; i < 12; ++i) {
        points.push_back(Vector3D::point(Points_T[0][i] , Points_T[1][i] , Points_T[2][i] ) );
    }
    // Create faces
    vector<Face> faces;
    for (int j = 0; j < 20; ++j) {
        faces.push_back( Face( { Faces_T[0][j] , Faces_T[1][j] , Faces_T[2][j] } ) );
    }
    // Create new figure
    Figure* newIso;
    newIso = new Figure( points , faces , Color( lineColor.at(0) , lineColor.at(1) , lineColor.at(2) ) );
    return newIso;
}

Figure* createDodecahedron(vector<double>&lineColor){
    // Points array
    double Points_T [3][20];
    // Faces array
    int Faces_T [5][12] = {
            { 0 , 0 , 1 , 2 , 3 , 4 , 19 , 19 , 18 , 17 , 16 , 15 },
            { 1 , 5 , 7 , 9 , 11 , 13 , 18 , 14 , 12 , 10 , 8 , 6 },
            { 2 , 6 , 8 , 10 , 12 , 14 , 17 , 13 , 11 , 9 , 7 , 5 },
            { 3 , 7 , 9 , 11 , 13 , 5 , 16 , 12 , 10 , 8 , 6 , 14 },
            { 4 , 1 , 2 , 3 , 4 , 0 , 15 , 18 , 17 , 16 , 15 , 19 }
    };
    // Create dodecahedron
    Figure* newICO = createIcosahedron(lineColor);
    int count = 0;
    for (Face f : newICO->faces){
        Vector3D p4 = Vector3D::point(0 , 0 , 0);
        p4 += newICO->points.at(f.point_indexes.at(0));
        p4 += newICO->points.at(f.point_indexes.at(1));
        p4 += newICO->points.at(f.point_indexes.at(2));
        p4 /= 3;
        Points_T[0][count] = p4.x;
        Points_T[1][count] = p4.y;
        Points_T[2][count] = p4.z;
        count++;
    }
    // Create all the points
    vector<Vector3D> points;
    for (int i = 0; i < 20; ++i) {
        points.push_back(Vector3D::point(Points_T[0][i] , Points_T[1][i] , Points_T[2][i] ) );
    }
    // Create faces
    vector<Face> faces;
    for (int j = 0; j < 12; ++j) {
        faces.push_back( Face( { Faces_T[0][j] , Faces_T[1][j] , Faces_T[2][j] , Faces_T[3][j] , Faces_T[4][j] } ) );
    }
    // Create new figure
    Figure* newDodeca;
    newDodeca = new Figure( points , faces , Color( lineColor.at(0) , lineColor.at(1) , lineColor.at(2) ) );
    return newDodeca;
}

void splitTriangle(Face &originalTriangle , Figure* &originalFigure , vector<Face> &newFaces){
    // Get the points of the triangle
    int index_A = originalTriangle.point_indexes.at(0);
    int index_B = originalTriangle.point_indexes.at(1);
    int index_C = originalTriangle.point_indexes.at(2);
    Vector3D A = originalFigure->points.at(index_A);
    Vector3D B = originalFigure->points.at(index_B);
    Vector3D C = originalFigure->points.at(index_C);
    // Calculate middle points
    Vector3D D = (A + B) / 2;
    Vector3D E = (A + C) / 2;
    Vector3D F = (B + C) / 2;
    // Push points to figure
    originalFigure->points.push_back(D);
    int index_D = static_cast<int>(originalFigure->points.size() - 1);
    originalFigure->points.push_back(E);
    int index_E = static_cast<int>(originalFigure->points.size() - 1);
    originalFigure->points.push_back(F);
    int index_F = static_cast<int>(originalFigure->points.size() - 1);

    // Replace triangles and add the new triangles
    //originalTriangle.point_indexes = { index_A , index_D , index_E};
    newFaces.push_back( Face( {index_A , index_D , index_E} ) );
    newFaces.push_back( Face( {index_B , index_F , index_D} ) );
    newFaces.push_back( Face( {index_C , index_E , index_F} ) );
    newFaces.push_back( Face( {index_D , index_F , index_E} ) );

}

Figure* createSphere(const double &radius , const int &n , vector<double>&lineColor){
    // Create an icosahedron
    Figure* newIcoSphere = createIcosahedron(lineColor);
    // Create temporary faces vector
    vector<Face>newFaces;

    for (int i = 0; i < n; ++i) {

        // Split each triangle into 4
        for (Face f : newIcoSphere->faces){
            splitTriangle( f , newIcoSphere , newFaces);
        }
        // Update the faces vector
        newIcoSphere->faces = newFaces;
        // Reset the temporary vector
        newFaces = {};
    }
    // Rescale all the points
    for (Vector3D &p : newIcoSphere->points){
        p.normalise();
    }
    return newIcoSphere;
}

Figure* createCone(vector<double> &lineColor , const int &n , const double &height){
    // Create points
    vector<Vector3D> points;
    for (int i = 0; i < n + 1; ++i) {
        if (i == n){
            points.push_back( Vector3D::point(0 , 0 , height) );
        }
        else{
            points.push_back ( Vector3D::point(cos( (2 * i * M_PI) / n ) , sin( (2 * i * M_PI) / n ) , 0) );
        }
    }
    // Get indexes vector
    vector<int> indexes;
    for (int i = n - 1; i >= 0; --i) {
        indexes.push_back(i);
    }
    // Create faces
    vector<Face> faces;
    for (int j = 0; j < n + 1; ++j) {
        if (j == n){
            faces.emplace_back(Face(indexes));
        }
        else{
            faces.push_back( Face( { j , (j+1) % n , n } ) );
        }
    }
    // Create new figure
    Figure* newCone;
    newCone = new Figure(points , faces , Color( lineColor.at(0) , lineColor.at(1) , lineColor.at(2) ) );
    return newCone;
}

Figure* createCylinder(vector<double> &lineColor , const int &n , const double &height){
    // Create points
    vector<Vector3D> points;
    // Add points of bottom face
    for (int i = 0; i < n; ++i) {
        points.push_back ( Vector3D::point(cos( (2 * i * M_PI) / n ) , sin( (2 * i * M_PI) / n ) , 0) );
    }
    // Add points of top face
    for (int i = n; i < 2 * n; ++i) {
        points.push_back ( Vector3D::point(cos( (2 * i * M_PI) / n ) , sin( (2 * i * M_PI) / n ) , height) );
    }
    // Create faces
    vector<Face> faces;
    for (int j = 0; j < n; ++j) {
        // Last face loops back to first one
        if (j == n-1){
            faces.push_back( Face( { j , 0 , n , j + n  } ) );
        }
        // Normal face
        else{
            faces.push_back( Face( { j , j+1 , j+n+1 , j + n  } ) );
        }
    }
    // Make new figure
    Figure* newFigure;
    newFigure = new Figure(points , faces , Color(lineColor.at(0) , lineColor.at(1) , lineColor.at(2) ) );
    return newFigure;
}

Figure* createTorus(const double &r , const double &R , const int &n , const int &m , vector<double> &lineColor){
    // Create points vector
    vector<Vector3D> points;
    // Create all the points
    for (unsigned int i = 0; i < n; ++i) {
        for (unsigned int j = 0; j < m; ++j) {
            // Calculate modifiers u and v
            double u = (2 * i * M_PI) / n;
            double v = (2 * j * M_PI) / m;
            // Calculate the coordinates of the new point
            double x_uv = ( R + r * cos(v) ) * cos(u);
            double y_uv = ( R + r * cos(v) ) * sin(u);
            double z_uv = r * sin(v);

            points.push_back( Vector3D::point( x_uv , y_uv , z_uv ) );
        }
    }
    // Create faces vector
    vector<Face> faces;
    for (unsigned int i = 0; i < n; ++i) {
        for (unsigned int j = 0; j < m; ++j) {
            int ind1 = i * n + j;
            int ind2 = ( (i + 1) % n ) * n + j;
            int ind3 = ( (i + 1) % n ) * n + (j + 1) % m;
            int ind4 = i * n + (j + 1) % m;
            faces.push_back( Face( { ind1 , ind2 , ind3 , ind4 } ) );
        }
    }
    // Create new figure
    Figure* newFigure;
    newFigure = new Figure( points , faces , Color( lineColor.at(0) , lineColor.at(1) , lineColor.at(2) ) );
    return newFigure;
}

string getEndString3D(const LParser::LSystem3D &l_system , string &startingString , string &endingString){
    // Replace symbols
    for(char c : startingString){
        // Add the operators
        if(c == '-' || c == '+' || c == '(' || c == ')' || c == '^' || c == '&' || c == '/' || c == '\\' || c == '|'){
            endingString += c;
        }
            // Replace the string
        else{
            endingString += l_system.get_replacement(c);
        }
    }
    // Update startingString
    startingString = endingString;
    endingString = "";
    return startingString;
}

Figure* createSystemLines3D (const LParser::LSystem3D &l_system , string &startingString,
        double &angle , vector<double> &lineColor, Vector3D &currentPoint , int current_c){
    // Create variables to store the current x and y coordinate
    vector<double> current_x;
    vector<double> current_y;
    vector<double> current_z;
    // Create new figure
    Figure* newFigure;
    // Create points and faces vectors
    vector<Vector3D> points;
    vector<Face> faces;
    // Add the origin to the points vector
    points.push_back( Vector3D::point(currentPoint) );
    // Create orientation vectors
    Vector3D H = Vector3D::point( 1 , 0 , 0 );
    Vector3D L = Vector3D::point( 0 , 1 , 0 );
    Vector3D U = Vector3D::point( 0 , 0 , 1 );
    // Create vectors to update the orientations
    Vector3D H_new = H;
    Vector3D L_new = L;
    Vector3D U_new = U;
    // Create containers for orientation vectors
    vector<Vector3D> current_H;
    vector<Vector3D> current_L;
    vector<Vector3D> current_U;
    vector<int> current_index;
    // Loop through characters in initiating string
    int count = 0;
    for( int i = current_c; i < startingString.length(); i++) {
        // Rudder left
        if (startingString[i] == '+') {
            // Calculate the new orientations
            H_new = ( H * cos(angle) + L * sin(angle) );
            L_new = ( L * cos(angle) - H * sin(angle) );
            // Update the orientations
            H = H_new;
            L = L_new;
        }
        // Rudder right
        else if (startingString[i] == '-') {
            // Calculate the new orientations
            H_new = ( H * cos(-angle) + L * sin(-angle) );
            L_new = ( L * cos(-angle) - H * sin(-angle) );
            // Update the orientations
            H = H_new;
            L = L_new;
        }
        // Pitch up
        else if (startingString[i] == '^') {
            // Calculate the new orientations
            H_new = ( H * cos(angle) + U * sin(angle) );
            U_new = ( U * cos(angle) - H * sin(angle) );
            // Update the orientations
            H = H_new;
            U = U_new;
        }
        // Pitch down
        else if (startingString[i] == '&') {
            // Calculate the new orientations
            H_new = ( H * cos(-angle) + U * sin(-angle) );
            U_new = ( U * cos(-angle) - H * sin(-angle) );
            // Update the orientations
            H = H_new;
            U = U_new;
        }
        // Roll left
        else if (startingString[i] == '\\') {
            // Calculate the new orientations
            L_new = ( L * cos(angle) - U * sin(angle) );
            U_new = ( L * sin(angle) + U * cos(angle) );
            // Update the orientations
            L = L_new;
            U = U_new;
        }
        // Roll right
        else if (startingString[i] == '/') {
            // Calculate the new orientations
            L_new = ( L * cos(-angle) - U * sin(-angle) );
            U_new = ( L * sin(-angle) + U * cos(-angle) );
            // Update the orientations
            L = L_new;
            U = U_new;
        }
        // Turn around
        else if (startingString[i] == '|') {
            H = -H;
            L = -L;
        }
        else if (startingString[i] == '('){
            // Save coordinates, draw everything within the bracket, return to saved coordinates
            current_x.push_back(currentPoint.x);
            current_y.push_back(currentPoint.y);
            current_z.push_back(currentPoint.z);
            current_H.push_back(H);
            current_L.push_back(L);
            current_U.push_back(U);
            current_index.push_back(count);
        }
        else if (startingString[i] == ')'){

            currentPoint.x = current_x.back();
            current_x.pop_back();

            currentPoint.y = current_y.back();
            current_y.pop_back();

            currentPoint.z = current_z.back();
            current_z.pop_back();

            H = current_H.back();
            current_H.pop_back();

            L = current_L.back();
            current_L.pop_back();

            U = current_U.back();
            current_U.pop_back();

            count = current_index.back();
            current_index.pop_back();
        }
            // If we must draw
        else if(l_system.draw(startingString[i])){
            // Create new points
            currentPoint += H;
            Vector3D newPoint = Vector3D::point(currentPoint);
            points.push_back(newPoint);
            // Create new faces (the lines)
            faces.push_back( Face( {count , (int)points.size() - 1} ) );
            count = (int)points.size() - 1;
        }
            // If we mustn't draw
        else {
            currentPoint += H;
            Vector3D newPoint = Vector3D::point(currentPoint);
            points.push_back(currentPoint);
        }

    }
    // Create the new figure
    newFigure = new Figure(points , faces , Color( lineColor.at(0), lineColor.at(1), lineColor.at(2) ));
    return newFigure;
}

Figure* drawSystem3D(const LParser::LSystem3D &l_system, const int &size, vector<double> &backgroundColor,
                     vector<double> &lineColor) {
    // Get all the components of the LSystem3D
    const auto &Alphabet = l_system.get_alphabet();
    const string &initiator = l_system.get_initiator();
    double angle = l_system.get_angle();

    // Convert angles to radians
    angle = ( angle * M_PI ) / 180;

    // Get iterations
    unsigned int nr_iterations = l_system.get_nr_iterations();

    // Get initiating string
    string startingString = l_system.get_initiator();
    string endingString;
    // Replace symbols
    for( int i = 0; i < nr_iterations; i++){
        startingString = getEndString3D(l_system, startingString, endingString);
    }
    Vector3D currentPoint = Vector3D::point(0,0,0);
    return createSystemLines3D(l_system,startingString,angle,lineColor,currentPoint,0);
}


Figure* drawLineDrawing(double &scale , double &rotX , double &rotY , double &rotZ , int &nrPoints , int &nrLines ,
                       const ini::Configuration &configuration , vector<double> &lineColor , vector<double> &center ,
                       Matrix &m_eye , int &i){
    // Make new Figure class
    vector<Vector3D> points;
    vector<Face> faces;

    // Make temporary vector;
    vector<double> point_to_add;
    vector<int> line_to_add;
    // Get points
    for (int j = 0; j < nrPoints; j++){
        point_to_add = configuration["Figure"+to_string(i)]["point"+to_string(j)].as_double_tuple_or_die();
        points.push_back(Vector3D::point(point_to_add.at(0) , point_to_add.at(1) , point_to_add.at(2)));
    }
    // Get lines
    for (int j = 0; j < nrLines; j++){
        line_to_add = configuration["Figure"+to_string(i)]["line"+to_string(j)].as_int_tuple_or_die();
        faces.push_back(Face( {line_to_add.at(0) , line_to_add.at(1)} ));
    }
    // Create new figure
    Figure* newFigure;
    newFigure = new Figure(points , faces , Color( lineColor.at(0), lineColor.at(1), lineColor.at(2) ));
    Matrix m = scaleFigure(scale);
    m *= rotateX((rotX * M_PI) / 180);
    m *= rotateY((rotY * M_PI) / 180);
    m *= rotateZ((rotZ * M_PI) / 180);
    m *= translate(Vector3D::point(center.at(0) , center.at(1) , center.at(2)));
    m *= m_eye;
    newFigure->applyTransformation(m);
    return newFigure;
}

Figures3D drawWireframe(int &size , vector<double> &eye , vector<double> &backgroundcolor , int &nrFigures ,
        const ini::Configuration &configuration ){
    // Make the variables
    double theta;
    double phi;
    double r;
    Vector3D eyePoint = Vector3D::point(eye.at(0) , eye.at(1) , eye.at(2));
    Matrix m_eye = eyePointTrans(eyePoint , theta , phi , r);
    // Create figures vector
    Figures3D figures;
    // Get figures
    for (int i = 0; i < nrFigures; i++){
        // Get all attributes
        // Get figure type
        string figure_type = configuration["Figure"+to_string(i)]["type"].as_string_or_die();
        // Get common attributes
        double scale = configuration["Figure"+to_string(i)]["scale"].as_double_or_default(1.0);
        double rotX = configuration["Figure"+to_string(i)]["rotateX"].as_double_or_default(0);
        double rotY = configuration["Figure"+to_string(i)]["rotateY"].as_double_or_default(0);
        double rotZ = configuration["Figure"+to_string(i)]["rotateZ"].as_double_or_default(0);
        vector<double> center = configuration["Figure"+to_string(i)]["center"].as_double_tuple_or_default({0,0,0});
        vector<double> lineColor = configuration["Figure"+to_string(i)]["color"].as_double_tuple_or_default({0,0,0});

        Matrix m = scaleFigure(scale);
        m *= rotateX((rotX * M_PI) / 180);
        m *= rotateY((rotY * M_PI) / 180);
        m *= rotateZ((rotZ * M_PI) / 180);
        m *= translate( Vector3D::point( center.at(0) , center.at(1) , center.at(2) ) );
        m *= m_eye;

        // Figure_type : "LineDrawing"
        if (figure_type == "LineDrawing"){
            // Get attributes specific for LineDrawing
            int nrPoints = configuration["Figure"+to_string(i)]["nrPoints"].as_int_or_die();
            int nrLines = configuration["Figure"+to_string(i)]["nrLines"].as_int_or_die();
            // Create new LineDrawing figure
            Figure* newFigure = drawLineDrawing(scale , rotX , rotY , rotZ , nrPoints , nrLines , configuration ,
                                               lineColor , center , m_eye , i);
            // Add figure to vector of figures
            figures.push_back(newFigure);
        }
        // Figure_type : "Cube"
        else if (figure_type == "Cube"){
            // Create new figure
            Figure* newFigure;
            newFigure = createCube(lineColor);
            newFigure->applyTransformation(m);
            figures.push_back(newFigure);
        }
        // Figure_type == "Tetrahedron"
        else if (figure_type == "Tetrahedron"){
            // Create new figure
            Figure* newFigure;
            newFigure = createTetrahedron(lineColor);
            newFigure->applyTransformation(m);
            figures.push_back(newFigure);
        }
        // Figure_type == "Octahedron"
        else if (figure_type == "Octahedron"){
            // Create new figure
            Figure* newFigure;
            newFigure = createOctahedron(lineColor);
            newFigure->applyTransformation(m);
            figures.push_back(newFigure);
        }
        // Figure_type == "Icosahedron"
        else if (figure_type == "Icosahedron"){
            // Create new figure
            Figure* newFigure;
            newFigure = createIcosahedron(lineColor);
            newFigure->applyTransformation(m);
            figures.push_back(newFigure);
        }
        // Figure_type == "Dodecahedron"
        else if (figure_type == "Dodecahedron"){
            // Create new figure
            Figure* newFigure;
            newFigure = createDodecahedron(lineColor);
            newFigure->applyTransformation(m);
            figures.push_back(newFigure);
        }
        // Figure_type == "Cone"
        else if (figure_type == "Cone"){
            int n = configuration["Figure"+to_string(i)]["n"].as_int_or_die();
            double height = configuration["Figure"+to_string(i)]["height"].as_double_or_die();
            Figure* newFigure;
            newFigure = createCone(lineColor , n , height);
            newFigure->applyTransformation(m);
            figures.push_back(newFigure);
        }
        // Figure_type == "Cylinder"
        else if (figure_type == "Cylinder"){
            const int n = configuration["Figure"+to_string(i)]["n"].as_int_or_die();
            const double height = configuration["Figure" + to_string(i)]["height"].as_double_or_die();
            Figure* newFigure;
            newFigure = createCylinder(lineColor , n , height);
            newFigure->applyTransformation(m);
            figures.push_back(newFigure);
        }
        // Figure_type == "Sphere"
        else if (figure_type == "Sphere"){
            const int n = configuration["Figure"+to_string(i)]["n"].as_int_or_die();
            Figure* newFigure;
            newFigure = createSphere(1.0 , n , lineColor);
            newFigure->applyTransformation(m);
            figures.push_back(newFigure);
        }
        // Figure_type == "Torus"
        else if (figure_type == "Torus"){
            double _R = configuration["Figure"+to_string(i)]["R"].as_double_or_die();
            double _r = configuration["Figure"+to_string(i)]["r"].as_double_or_die();
            int _n = configuration["Figure"+to_string(i)]["n"].as_int_or_die();
            int _m = configuration["Figure"+to_string(i)]["m"].as_int_or_die();
            Figure* newFigure;
            newFigure = createTorus( _r , _R , _n , _m , lineColor);
            newFigure->applyTransformation(m);
            figures.push_back(newFigure);
        }
        // Figure_type == "3DLSystem"
        else if (figure_type == "3DLSystem"){
            string inputFile = configuration["Figure"+to_string(i)]["inputfile"].as_string_or_die();
            // Initialize the parser
            LParser::LSystem3D l_system;
            ifstream input_stream(inputFile);
            input_stream >> l_system;
            input_stream.close();
            Figure* newFigure;
            newFigure = drawSystem3D(l_system , size , backgroundcolor , lineColor);
            newFigure->applyTransformation(m);
            figures.push_back(newFigure);
        }

    }
    return figures;
}

img::EasyImage generate_image(const ini::Configuration &configuration)
{
    // Get type
    string type = configuration["General"]["type"].as_string_or_die();
    // Case : type == "2DLSystem"
    if (type == "2DLSystem"){
        int size = configuration["General"]["size"].as_int_or_die();
        vector<double> backgroundColor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
        string input_filename = configuration["2DLSystem"]["inputfile"].as_string_or_die();
        vector<double> lineColor = configuration["2DLSystem"]["color"].as_double_tuple_or_die();
        // Initialize the parser
        LParser::LSystem2D l_system;
        ifstream input_stream(input_filename);
        input_stream >> l_system;
        input_stream.close();
        return draw2DLines(drawSystem2D(l_system, size, backgroundColor, lineColor) , size , backgroundColor );
    }
    // Case: type == "Wireframe"
    else if (type == "Wireframe"){
        // Get general properties
        int size = configuration["General"]["size"].as_int_or_die();
        vector<double> eye = configuration["General"]["eye"].as_double_tuple_or_die();
        vector<double> backgroundcolor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
        int nrFigures = configuration["General"]["nrFigures"].as_int_or_die();
        // Draw the wireframe
        Figures3D figures = drawWireframe(size , eye , backgroundcolor , nrFigures , configuration);
        Lines2D lines = doProjection(figures);
        return draw2DLines(lines ,
                size , backgroundcolor);

    }

    /*
     * [ImageProperties]
     * width = width_of_image
     * height = height_of_image
     */
    int width = configuration["ImageProperties"]["width"].as_int_or_die(); // Get width
    int height = configuration["ImageProperties"]["height"].as_int_or_die(); // Get height

    // Case : type = "IntroColorRectangle"
    if (type == "IntroColorRectangle"){
        return createColorRectangle(width, height);
    }
    // Case : type == "IntroBlocks"
    else if (type == "IntroBlocks"){
        int blocksInX = configuration["BlockProperties"]["nrXBlocks"].as_int_or_die();
        int blocksInY = configuration["BlockProperties"]["nrYBlocks"].as_int_or_die();
        vector<double> colorWhite = configuration["BlockProperties"]["colorWhite"].as_double_tuple_or_die();
        vector<double> colorBlack = configuration["BlockProperties"]["colorBlack"].as_double_tuple_or_die();
        bool invertColors = configuration["BlockProperties"]["invertColors"].as_bool_or_default(false);
        return createBlocks(width, height, blocksInX, blocksInY, colorWhite, colorBlack, invertColors);
    }
    // Case : type == "IntroLines"
    else if (type == "IntroLines"){
        string figure = configuration["LineProperties"]["figure"].as_string_or_die();
        vector<double> backgroundColor = configuration["LineProperties"]["backgroundcolor"].as_double_tuple_or_die();
        vector<double> lineColor = configuration["LineProperties"]["lineColor"].as_double_tuple_or_die();
        int linesNumber = configuration["LineProperties"]["nrLines"].as_int_or_die();
        // Case: figure == "QuarterCircle"
        if(figure == "QuarterCircle"){
            return createQuarterCircle(width,height,linesNumber,backgroundColor,lineColor);
        }
        // Case: figure == "Eye"
        if(figure == "Eye"){
            return createEye(width,height,linesNumber,backgroundColor,lineColor);
        }
        // Case: figure == "Diamond"
        if(figure == "Diamond"){
            return createDiamond(width,height,linesNumber,backgroundColor,lineColor);
        }

    }

    return img::EasyImage();
}



int main(int argc, char const* argv[])
{
        int retVal = 0;
        try
        {
                std::vector<std::string> args = std::vector<std::string>(argv+1, argv+argc);
                if (args.empty()) {
                        std::ifstream fileIn("filelist");
                        std::string filelistName;
                        while (std::getline(fileIn, filelistName)) {
                                args.push_back(filelistName);
                        }
                }
                for(std::string fileName : args)
                {
                        ini::Configuration conf;
                        try
                        {
                                std::ifstream fin(fileName);
                                fin >> conf;
                                fin.close();
                        }
                        catch(ini::ParseException& ex)
                        {
                                std::cerr << "Error parsing file: " << fileName << ": " << ex.what() << std::endl;
                                retVal = 1;
                                continue;
                        }

                        img::EasyImage image = generate_image(conf);
                        if(image.get_height() > 0 && image.get_width() > 0)
                        {
                                std::string::size_type pos = fileName.rfind('.');
                                if(pos == std::string::npos)
                                {
                                        //filename does not contain a '.' --> append a '.bmp' suffix
                                        fileName += ".bmp";
                                }
                                else
                                {
                                        fileName = fileName.substr(0,pos) + ".bmp";
                                }
                                try
                                {
                                        std::ofstream f_out(fileName.c_str(),std::ios::trunc | std::ios::out | std::ios::binary);
                                        f_out << image;

                                }
                                catch(std::exception& ex)
                                {
                                        std::cerr << "Failed to write image to file: " << ex.what() << std::endl;
                                        retVal = 1;
                                }
                        }
                        else
                        {
                                std::cout << "Could not generate image for " << fileName << std::endl;
                        }
                }
        }
        catch(const std::bad_alloc &exception)
        {
    		//When you run out of memory this exception is thrown. When this happens the return value of the program MUST be '100'.
    		//Basically this return value tells our automated test scripts to run your engine on a pc with more memory.
    		//(Unless of course you are already consuming the maximum allowed amount of memory)
    		//If your engine does NOT adhere to this requirement you risk losing points because then our scripts will
    		//mark the test as failed while in reality it just needed a bit more memory
                std::cerr << "Error: insufficient memory" << std::endl;
                retVal = 100;
        }
        return retVal;
}

#pragma clang diagnostic pop