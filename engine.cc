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
#include <utility>
#include <cassert>
#include <limits>
#include <map>
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

    /**
     * @brief Constructs a Color type object
     * @param red Value of red component
     * @param green Value of green component
     * @param blue Value of blue component
     */
    Color(double red, double green, double blue) : red(red), green(green), blue(blue) {}

    virtual ~Color() {

    }
};

class Point2D {
public:
    double x;
    double y;

    /**
     * @brief Constructs a Point2D type object
     * @param x x-coordinate
     * @param y y-coordinate
     */
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

    double z1;
    double z2;

    /**
     *
     * @param p1 Begin point : Point2D type object
     * @param p2 End point : Point2D type object
     * @param color Line color : Color type object
     */
    Line2D(const Point2D &p1, const Point2D &p2, const Color &color) : p1(p1), p2(p2), color(color) {}

    /**
     *
     * @param p1 Begin point : Point2D type object
     * @param p2 End point : Point2D type object
     * @param color Line color : Color type object
     * @param z1 z-coordinate of begin point
     * @param z2 z-coordinate of end point
     */
    Line2D(const Point2D &p1, const Point2D &p2, const Color &color, double &z1, double &z2) : p1(p1), p2(p2),
                                                                                             color(color),
                                                                                             z1(z1), z2(z2) {}
    /**
     * @brief Constructs a Line2D from a reference object
     * @param refLine A pointer to a Line2D type object
     */
    Line2D(Line2D* refLine) :   p1(refLine->p1) , p2(refLine->p2) ,
                                color(refLine->color) , z1(refLine->z1) ,
                                z2(refLine->z2) {}

    virtual ~Line2D() {

    }
};

class Face {

public:
    // These indexes refer to points in the 'points' vector of the Figure-class
    vector<int> point_indexes;

    Face(const vector<int> &point_indexes) : point_indexes(point_indexes) {}

    Face(Face* refFace) : point_indexes(refFace->point_indexes) {}

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

    Figure(Figure* refFig) : points(refFig->points) , faces(refFig->faces) , color(refFig->color) {}

    void applyTransformation(const Matrix &m){
        // Multiply each vector with the matrix
        for(Vector3D &v : points){
            v *= m;
        }
    }

    virtual ~Figure() {

    }
};

const double posInf = numeric_limits<double>::infinity();
const double negInf = -numeric_limits<double>::infinity();



class ZBuffer : public vector< vector<double> >{
public:

    ZBuffer( const int width , const int height ) {
        for (int i = 0; i < width; ++i) {
            this->push_back(vector<double>(height , posInf) );
        }
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

void getLinePointIndex(Face &face, Figure* &f, Lines2D &lines , const double d){
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
        Point2D newBeginP = doProjection(beginP, d);
        Point2D newEndP = doProjection(endP, d);
        // Create new line
        Line2D newLine(newBeginP, newEndP, f->color , beginP.z , endP.z);
        //newLine.z1 = beginP.z;
        //newLine.z2 = endP.z;
        lines.push_back(newLine);
    }
}

Lines2D doProjection(Figures3D &figs , const double d){
    Lines2D lines;
    for(Figure* f : figs){
        for(Face face : f->faces){
            // Get points index - loop through
            getLinePointIndex(face, f, lines , d);
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
            { 0 , 1 , 0 , 0 },
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
    vector<int> bottom_indexes;
    vector<int> top_indexes;
    // Add points of bottom face
    for (int i = 0; i < n; ++i) {
        points.push_back ( Vector3D::point(cos( (2 * i * M_PI) / n ) , sin( (2 * i * M_PI) / n ) , 0) );
        bottom_indexes.push_back(i);
    }
    // Add points of top face
    for (int i = n; i < 2 * n; ++i) {
        points.push_back ( Vector3D::point(cos( (2 * i * M_PI) / n ) , sin( (2 * i * M_PI) / n ) , height) );
        top_indexes.push_back(i);
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
    // Add bottom face
    faces.emplace_back(bottom_indexes);
    // Add top face
    faces.emplace_back(top_indexes);
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

void insertToPentagon(int index_A , int index_D , int index_I , Vector3D &D , Vector3D &I , map< int , vector<int> > &pentagons , vector<Vector3D> &originalPoints){
    // First pair of points
    if ( pentagons[index_A].empty() ){
        pentagons[index_A].push_back(index_D);
        pentagons[index_A].push_back(index_I);
        return;
    }
    else if( pentagons[index_A].size() == 5 ){
        return;
    }
    // Last pair of points
    else if ( pentagons[index_A].front() == index_D && pentagons[index_A].at(1) != index_I){
        pentagons[index_A].insert(pentagons[index_A].begin() , index_I);
        return;
    }
    else if ( pentagons[index_A].front() == index_I && pentagons[index_A].at(1) != index_D ){
        pentagons[index_A].insert(pentagons[index_A].begin() , index_D);
        return;
    }
        // Second pair of points
    else if (pentagons[index_A].back() == index_D && pentagons[index_A].at(pentagons[index_A].size() - 1) == index_I){
        pentagons[index_A].push_back(index_I);
        return;
    }
    else if (pentagons[index_A].back() == index_I && pentagons[index_A].at(pentagons[index_A].size() - 1) == index_D){
        pentagons[index_A].push_back(index_D);
        return;
    }
        // Pair in the middle
    else{
        for(int j= 0; j < pentagons[index_A].size(); j++){
            if( pentagons[index_A].at(j) == index_D){
                // Insert after this point
                pentagons[index_A].insert(pentagons[index_A].begin() + j + 1 , index_I);
                return;
            }
            else if ( pentagons[index_A].at(j) == index_I){
                // Insert before this point
                pentagons[index_A].insert(pentagons[index_A].begin() + j , index_D);
                return;
            }
        }
        // If line cannot be connected to previously present points
        pentagons[index_A].push_back(index_D);
        return;
    }
}

int isElementOf(vector<Vector3D> &points, Vector3D &p){
    for (unsigned int i = 0; i < points.size(); ++i) {
        if(points.at(i) == p){
            return i;
        }
    }
    return static_cast<int>(points.size());
}

void splitTriangleHexagon(Face &originalTriangle , Figure* &originalFigure , vector<Face> &newFaces , map< int , vector<int> > &pentagons){
    // Get the points of the triangle
    int index_A = originalTriangle.point_indexes.at(0);
    int index_B = originalTriangle.point_indexes.at(1);
    int index_C = originalTriangle.point_indexes.at(2);
    Vector3D A = originalFigure->points.at(index_A);
    Vector3D B = originalFigure->points.at(index_B);
    Vector3D C = originalFigure->points.at(index_C);
    // Calculate new points
    Vector3D D = 2 * (A / 3) + B / 3;
    Vector3D E = 2 * (B / 3) + A / 3;
    Vector3D F = 2 * (B / 3) + C / 3;
    Vector3D G = 2 * (C / 3) + B / 3;
    Vector3D I = 2 * (A / 3) + C / 3;
    Vector3D H = 2 * (C / 3) + A / 3;
    // Push points to figure

    int index_D = isElementOf(originalFigure->points, D);
    if(index_D == originalFigure->points.size()){
        originalFigure->points.push_back(D);
        index_D = static_cast<int>(originalFigure->points.size() - 1);
    }

    int index_E = isElementOf(originalFigure->points, E);
    if(index_E == originalFigure->points.size()){
        originalFigure->points.push_back(E);
        index_E = static_cast<int>(originalFigure->points.size() - 1);
    }

    int index_F = isElementOf(originalFigure->points, F);
    if(index_F == originalFigure->points.size()){
        originalFigure->points.push_back(F);
        index_F = static_cast<int>(originalFigure->points.size() - 1);
    }
    int index_G = isElementOf(originalFigure->points, G);
    if(index_G == originalFigure->points.size()){
        originalFigure->points.push_back(G);
        index_G = static_cast<int>(originalFigure->points.size() - 1);
    }
    int index_H = isElementOf(originalFigure->points, H);
    if(index_H == originalFigure->points.size()){
        originalFigure->points.push_back(H);
        index_H = static_cast<int>(originalFigure->points.size() - 1);
    }
    int index_I = isElementOf(originalFigure->points, I);
    if(index_I == originalFigure->points.size()){
        originalFigure->points.push_back(I);
        index_I = static_cast<int>(originalFigure->points.size() - 1);
    }

    // Create hexagon
    newFaces.push_back( Face( {index_D , index_E , index_F , index_G , index_H , index_I} ) );
    // Insert new points connected to A
    insertToPentagon(index_A , index_D , index_I , D , I , pentagons , originalFigure->points);
    // Insert new points connected to B
    insertToPentagon(index_B , index_F , index_E , F , E , pentagons , originalFigure->points);
    // Insert new points connected to C
    insertToPentagon(index_C , index_H , index_G , H , G , pentagons , originalFigure->points);
}

Figure* createBuckyBall(vector<double> &lineColor){
    /*
     * Een voetbal die bestaat uit 20 zeshoeken en 12 vijfhoeken.
     * Elk van de 20 driehoeken van de icosahedron op te delen in eengelijkzijdige zeshoek en drie driehoeken.
     */
    // Create an icosahedron
    Figure* newIcoSphere = createIcosahedron(lineColor);
    // Create temporary faces vector
    vector<Face>newFaces;
    // Create temporary points vector
    vector<Vector3D> points;
    // Create pentagons
    map< int , vector<int> >pentagons;
    for(int i = 0; i < newIcoSphere->points.size(); i++){
        pentagons[i] = {};
    }
    // Split each triangle into 4
    for (Face &f : newIcoSphere->faces){
        splitTriangleHexagon( f , newIcoSphere , newFaces , pentagons);
    }
    // Create new pentagon faces
    for( auto &f : pentagons ){
        newFaces.push_back(new Face(f.second));
    }
    newIcoSphere->faces = newFaces;
    // Rescale all the points
    for (Vector3D &p : newIcoSphere->points){
        p.normalise();
    }
    return newIcoSphere;

}

Figures3D createMengerSponge(vector<double> &lineColor , int nr_Iterations , Matrix &m ){
    // Start with cube
    Figures3D newSponge = {};
    Figure* newFig = createCube(lineColor);
    newFig->applyTransformation(m);
    newSponge.push_back(newFig);
    // Divide every face into nine squares
    Matrix m_s = scaleFigure(1/3);
    // Create empty translating matrix
    Matrix m_t;
    // Copy begin point
    Figures3D newFractal;
    // For each iteration
    for (int i = 0; i < nr_Iterations; ++i) {
        // For each figure
        for (auto &it : newSponge) {
            // For each point
            for (int j = 0; j < it->points.size(); ++j) {
                // Create corner cube
                // Copy the original figure
                newFig = new Figure(it);
                // Scale the figure
                newFig->applyTransformation(m_s);
                // Get translation matrix
                m_t = translate( it->points.at(j) - newFig->points.at(j) );
                // Translate the points
                newFig->applyTransformation(m_t);
                // Add figure to list of fractals
                newFractal.push_back(newFig);
            }
            // For each new figure made

        }
        newSponge = newFractal;
        newFractal = {};
    }
    return newSponge;
}

void generateFractal(Figures3D &fractal, const int nr_iterations, const double scale){

    Figure* newFig = nullptr;
    // Create scaling matrix
    Matrix m_s = scaleFigure(1/scale);
    // Create empty translating matrix
    Matrix m_t;
    // Copy begin point
    Figures3D newFractal;
    // For each iteration
    for (int i = 0; i < nr_iterations; ++i) {
        // For each figure
        for (auto &it : fractal) {
            // For each point
            for (int j = 0; j < it->points.size(); ++j) {
                // Copy the original figure
                newFig = new Figure(it);
                // Scale the figure
                newFig->applyTransformation(m_s);
                // Get translation matrix
                m_t = translate( it->points.at(j) - newFig->points.at(j) );
                // Translate the points
                newFig->applyTransformation(m_t);
                // Add figure to list of fractals
                newFractal.push_back(newFig);
            }
        }
        fractal = newFractal;
        newFractal = {};
    }
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
        // Figure_type == "BuckyBall"
        else if (figure_type == "BuckyBall"){
            Figure* newFigure;
            newFigure = createBuckyBall(lineColor);
            newFigure->applyTransformation(m);
            figures.push_back(newFigure);
        }
        // Figure_type = "MengerSponge"
        else if (figure_type == "MengerSponge"){
            int nr_Iterations = configuration["Figure"+to_string(i)]["nrIterations"].as_int_or_default(0);
            Figures3D newSponge = createMengerSponge(lineColor , nr_Iterations , m);
            figures.insert(figures.end() , newSponge.begin() , newSponge.end());
        }
        // Figure_type == "FractalCube"
        else if(figure_type == "FractalCube"){
            int nr_Iterations = configuration["Figure"+to_string(i)]["nrIterations"].as_int_or_default(0);
            double fractal_scale = configuration["Figure"+to_string(i)]["fractalScale"].as_double_or_default(1);
            Figure* newFigure = createCube(lineColor);
            newFigure->applyTransformation(m);
            Figures3D newFractal = {newFigure};
            generateFractal(newFractal , nr_Iterations , fractal_scale);
            figures.insert(figures.end() , newFractal.begin() , newFractal.end());
        }
        // Figure_type == "FractalTetrahedron"
        else if(figure_type == "FractalTetrahedron"){
            int nr_Iterations = configuration["Figure"+to_string(i)]["nrIterations"].as_int_or_default(0);
            double fractal_scale = configuration["Figure"+to_string(i)]["fractalScale"].as_double_or_default(1);
            Figure* newFigure = createTetrahedron(lineColor);
            newFigure->applyTransformation(m);
            Figures3D newFractal = {newFigure};
            generateFractal(newFractal , nr_Iterations , fractal_scale);
            figures.insert(figures.end() , newFractal.begin() , newFractal.end());
        }
        // Figure_type == "FractalIcosahedron"
        else if(figure_type == "FractalIcosahedron"){
            int nr_Iterations = configuration["Figure"+to_string(i)]["nrIterations"].as_int_or_default(0);
            double fractal_scale = configuration["Figure"+to_string(i)]["fractalScale"].as_double_or_default(1);
            Figure* newFigure = createIcosahedron(lineColor);
            newFigure->applyTransformation(m);
            Figures3D newFractal = {newFigure};
            generateFractal(newFractal , nr_Iterations , fractal_scale);
            figures.insert(figures.end() , newFractal.begin() , newFractal.end());
        }
        // Figure_type == "FractalOctahedron"
        else if(figure_type == "FractalOctahedron"){
            int nr_Iterations = configuration["Figure"+to_string(i)]["nrIterations"].as_int_or_default(0);
            double fractal_scale = configuration["Figure"+to_string(i)]["fractalScale"].as_double_or_default(1);
            Figure* newFigure = createOctahedron(lineColor);
            newFigure->applyTransformation(m);
            Figures3D newFractal = {newFigure};
            generateFractal(newFractal , nr_Iterations , fractal_scale);
            figures.insert(figures.end() , newFractal.begin() , newFractal.end());
        }
        // Figure_type == "FractalDodecahedron"
        else if(figure_type == "FractalDodecahedron"){
            int nr_Iterations = configuration["Figure"+to_string(i)]["nrIterations"].as_int_or_default(0);
            double fractal_scale = configuration["Figure"+to_string(i)]["fractalScale"].as_double_or_default(1);
            Figure* newFigure = createDodecahedron(lineColor);
            newFigure->applyTransformation(m);
            Figures3D newFractal = {newFigure};
            generateFractal(newFractal , nr_Iterations , fractal_scale);
            figures.insert(figures.end() , newFractal.begin() , newFractal.end());
        }
            // Figure_type == "FractalBuckyBall"
        else if(figure_type == "FractalBuckyBall"){
            int nr_Iterations = configuration["Figure"+to_string(i)]["nrIterations"].as_int_or_default(0);
            double fractal_scale = configuration["Figure"+to_string(i)]["fractalScale"].as_double_or_default(1);
            Figure* newFigure = createBuckyBall(lineColor);
            newFigure->applyTransformation(m);
            Figures3D newFractal = {newFigure};
            generateFractal(newFractal , nr_Iterations , fractal_scale);
            figures.insert(figures.end() , newFractal.begin() , newFractal.end());
        }
    }
    return figures;
}

// Session 4 : Z-Buffering

double calculateInvZ(unsigned int i , unsigned int i_min , unsigned int i_max , const double &z0 , const double &z1){
    // Calculate current p
    double p = (double)i / (double)(i_max - i_min);
    return p / z0 + (1 - p) / z1;
}

void draw_zbuf_line( ZBuffer &zBuffer, img::EasyImage &image,
                     unsigned int x0, unsigned int y0,  double &z0,
                     unsigned int x1, unsigned int y1,  double &z1,
                     const Color &lineColor){

    assert( x0 < image.get_width() && y0 < image.get_height() );
    assert( x1 < image.get_width() && y1 < image.get_height() );

    img::Color newColor = img::Color(lineColor.red , lineColor.green , lineColor.blue);
    // Vertical line
    if (x0 == x1)
    {
        //special case for x0 == x1
        for (unsigned int i = std::min(y0, y1); i <= std::max(y0, y1); i++)
        {
            double inv_z;

            inv_z = y0 < y1 ? calculateInvZ(i - y0 , y0 , y1 , z0 , z1) : calculateInvZ(i - y1 , y1 , y0 , z1 , z0);
            // Check if we can draw
            if (inv_z < zBuffer[x0][i]) {
                // Update z-buffer
                zBuffer[x0][i] = inv_z;
                (image)(x0, i) = newColor;
            }
        }
    }
    // Horizontal line
    else if (y0 == y1)
    {
        //special case for y0 == y1
        for (unsigned int i = std::min(x0, x1); i <= std::max(x0, x1); i++)
        {
            double inv_z;

            inv_z = x0 < x1 ? calculateInvZ(i - x0 , x0 , x1 , z0 , z1) : calculateInvZ(i - x1 , x1 , x0 , z1 , z0);
            // Check if we can draw
            if (inv_z < zBuffer[i][y0]) {
                // Update z-buffer
                zBuffer[i][y0] = inv_z;
                (image)(i, y0) = newColor;
            }
        }
    }
    else
    {
        if (x0 > x1)
        {
            //flip points if x1>x0: we want x0 to have the lowest value
            swap(x0, x1);
            swap(y0, y1);
            swap(z0, z1); //right corners fixed with swap enabled
        }
        // Calculate the coefficient
        double m = ((double) y1 - (double) y0) / ((double) x1 - (double) x0);

        if (-1.0 <= m && m <= 1.0)
        {
            for (unsigned int i = 0; i <= (x1 - x0); i++)
            {
                double inv_z;

                inv_z = calculateInvZ(i , 0 , (x1 - x0) , z1 , z0);
                // Check if we can draw
                if (inv_z < zBuffer[x0 + i][(unsigned int) round(y0 + m * i)]) {
                    // Update z-buffer
                    zBuffer[x0 + i][(unsigned int) round(y0 + m * i)] = inv_z;
                    (image)(x0 + i, (unsigned int) round(y0 + m * i)) = newColor;
                }
            }
        }
        else if (m > 1.0)
        {
            for (unsigned int i = 0; i <= (y1 - y0); i++)
            {
                double inv_z;

                inv_z = calculateInvZ(i , 0 , (y1 - y0) , z1 , z0);
                if (inv_z < zBuffer[(unsigned int) round(x0 + (i / m))][y0 + i]) {
                    // Update z-buffer
                    zBuffer[(unsigned int) round(x0 + (i / m))][y0 + i] = inv_z;
                    (image)((unsigned int) round(x0 + (i / m)), y0 + i) = newColor;
                }
            }
        }
        else if (m < -1.0)
        {
            for (unsigned int i = 0; i <= (y0 - y1); i++)
            {
                double inv_z;

                inv_z = calculateInvZ(i , 0 , (y0 - y1) , z1 , z0);
                if (inv_z < zBuffer[(unsigned int) round(x0 - (i / m))][y0 - i]) {
                    // Update z-buffer
                    zBuffer[(unsigned int) round(x0 - (i / m))][y0 - i] = inv_z;
                    (image)((unsigned int) round(x0 - (i / m)), y0 - i) = newColor;
                }
            }
        }
    }
}

img::EasyImage draw2DZbuffLines (const Lines2D &lines , const int size , vector<double> &backgroundColor){
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
    // Create Z-Buffer
    ZBuffer zBuffer(image.get_width() , image.get_height());
    // Determine the scaling factor
    double d = 0.95 * (imageWidth / x_range);
    // Calculate for x and y
    double DC_x = d * ( (x_min + x_max) / 2 );
    double DC_y = d * ( (y_min + y_max) / 2 );
    double dx = imageWidth / 2 - DC_x;
    double dy = imageHeight / 2 - DC_y;
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
        unsigned int x1 = lround(l.p1.x);
        unsigned int x2 = lround(l.p2.x);
        unsigned int y1 = lround(l.p1.y);
        unsigned int y2 = lround(l.p2.y);
        double z1 = l.z1;
        double z2 = l.z2;
        // Fetch and rescale the colors
        originalColors = {l.color.red , l.color.green , l.color.blue};
        newColors = scaleColors(originalColors);
        // Draw the lines
        draw_zbuf_line( zBuffer , image , x1 , y1 , z1 , x2 , y2 , z2 , Color( newColors.at(0) , newColors.at(1) , newColors.at(2) ) );
    }
    return image;
}

// Session 5 : Z-Buffering met driehoeken

void getProjectedPoints(Figures3D &figs , double &d , double &x_range , double &y_range ,
                        double &x_min, double &y_min, double &x_max, double &y_max){
    for(auto f : figs){
        for(const auto &p : f->points){
            Point2D newPoint = doProjection(p, d);
            // Update min and max
            if(x_max <= newPoint.x){
                x_max = newPoint.x;
            }
            if(x_min >= newPoint.x){
                x_min = newPoint.x;
            }
            if(y_max <= newPoint.y){
                y_max = newPoint.y;
            }
            if(y_min >= newPoint.y){
                y_min = newPoint.y;
            }
        }
    }
    // Calculate range along the x-axis and y-axis
    x_range = x_max - x_min;
    y_range = y_max - y_min;
}

void triangulate(Face &originalFace , vector<Face> &newFaces){
    for (unsigned int i = 1; i < originalFace.point_indexes.size() - 1; ++i) {
        // Create new face
        newFaces.push_back( Face( { originalFace.point_indexes.at(0) , originalFace.point_indexes.at(i) , originalFace.point_indexes.at(i+1) } ) );
    }
}

double calculateIntersection( const int& y_i ,  Point2D const &P , Point2D const &Q){
    return Q.x + ( P.x - Q.x ) * ( ( y_i - Q.y ) / ( P.y - Q.y ) );
}

void draw_zbuf_triag(ZBuffer &zbuf , img::EasyImage &image ,
        Vector3D const &A, Vector3D const &B, Vector3D const &C, double &d, double &dx, double &dy, Color &color){
    // Convert color
    img::Color newColor = img::Color(lround(color.red * 255) , lround(color.green * 255) , lround(color.blue * 255));
    // Projection of the triangle
    Point2D newA = doProjection(A , d);
    newA.x += dx;
    newA.y += dy;
    Point2D newB = doProjection(B , d);
    newB.x += dx;
    newB.y += dy;
    Point2D newC = doProjection(C , d);
    newC.x += dx;
    newC.y += dy;
    // Determine y_min and y_max
    double y_min = min({newA.y , newB.y , newC.y} );
    double y_max = max({newA.y , newB.y , newC.y} );
    // Determine gravity center
    Vector3D G = Vector3D::point( (newA.x + newB.x + newC.x) / 3 , (newA.y + newB.y + newC.y) / 3 , 0);
    double inv_z_G = 1/(3 * A.z) + 1/(3 * B.z) + 1/(3 * C.z);
    // Determine dzdx and dzdy
    Vector3D u = B - A;
    Vector3D v = C - A;
    Vector3D w = Vector3D::cross(u , v);
    double k = w.x * A.x + w.y * A.y + w.z * A.z;
    double dzdx = w.x / ( -d * k );
    double dzdy = w.y / ( -d * k );
    double inv_z;
    // Create variables for x_l and x_r
    double x_l_AB = posInf;
    double x_l_AC = posInf;
    double x_l_BC = posInf;
    double x_r_AB = negInf;
    double x_r_AC = negInf;
    double x_r_BC = negInf;
    int x_l;
    int x_r;
    // Determine which pixels belong in the triangle
    for (int i = (int)round(y_min + 0.5); i <= (int)round(y_max - 0.5); ++i) {
        // Determine x_l and x_r for each line
        if( (i - newA.y) * (i - newB.y) <= 0 && ( newA.y != newB.y ) ){
            x_l_AB = x_r_AB = calculateIntersection( i , newA , newB );
        }
        if( (i - newB.y) * (i - newC.y) <= 0 && ( newB.y != newC.y ) ){
            x_l_BC = x_r_BC = calculateIntersection( i , newB , newC );
        }
        if( (i - newA.y) * (i - newC.y) <= 0 && ( newA.y != newC.y ) ){
            x_l_AC = x_r_AC = calculateIntersection( i , newA , newC );
        }
        // Determine x_l and x_r for y = y_i line
        x_l = static_cast<int>(round(min({x_l_AB, x_l_AC, x_l_BC} ) + 0.5 ));
        x_r = static_cast<int>(round(max({x_r_AB, x_r_AC, x_r_BC} ) - 0.5 ));
        for (int j = x_l; j <= x_r; ++j) {
            // Determine the inv_z value
            inv_z = 1.0001 * inv_z_G + (j - G.x) * dzdx + (i - G.y) * dzdy;
            if(inv_z < zbuf[j][i]){
                zbuf[j][i] = inv_z;
                image(j,i) = newColor;
            }
        }
        // Reset left and right limits
        x_l_AB = posInf;
        x_l_AC = posInf;
        x_l_BC = posInf;
        x_r_AB = negInf;
        x_r_AC = negInf;
        x_r_BC = negInf;
    }

}

img::EasyImage draw2DZbuffTriag (const int &size , vector<double> &backgroundColor , Figures3D &figures){
    // Scale colors
    vector<unsigned int> bgColor = scaleColors(backgroundColor);
    // Triangulate : get the new faces
    vector<Face> newFaces;
    for(auto &fig : figures){
        for(auto &f : fig->faces){
            triangulate(f , newFaces);
        }
        // Replace faces
        fig->faces = newFaces;
        newFaces.clear();
    }
    // Calculate range along the x-axis and y-axis
    double x_min = 0 , y_min = 0 , x_max = 0 , y_max = 0;
    double d = 1.0;
    double x_range;
    double y_range;
    getProjectedPoints(figures , d , x_range , y_range , x_min , y_min , x_max , y_max);
    // Calculate image dimensions
    double imageWidth   = size * ( x_range / max( x_range,y_range ) );
    double imageHeight  = size * ( y_range / max( x_range,y_range ) );
    // Create the image file
    img::EasyImage image(lround(imageWidth) , lround(imageHeight) , img::Color(bgColor.at(0) , bgColor.at(1) , bgColor.at(2) ));
    // Create Z-Buffer
    ZBuffer zBuffer(image.get_width() , image.get_height());
    // Determine the scaling factor
    d = 0.95 * (imageWidth / x_range);
    // Calculate for x and y
    double DC_x = d * ( (x_min + x_max) / 2 );
    double DC_y = d * ( (y_min + y_max) / 2 );
    double dx = imageWidth / 2 - DC_x;
    double dy = imageHeight / 2 - DC_y;
    // Apply the z-buffering algorithm
    for (auto &fig : figures) {
        for(auto &f : fig->faces){
            // Get the point indexes
            int ind_A = f.point_indexes.at(0);
            int ind_B = f.point_indexes.at(1);
            int ind_C = f.point_indexes.at(2);
            // Apply z-buffering algorithm
            draw_zbuf_triag(zBuffer , image , fig->points.at(ind_A) , fig->points.at(ind_B) , fig->points.at(ind_C) ,
                    d , dx , dy , fig->color);
        }
    }
    return image;
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
        Lines2D lines = doProjection(figures , 1.0);
        return draw2DLines(lines ,
                size , backgroundcolor);
    }
    // Case: type == "ZBufferedWireframe"
    else if (type == "ZBufferedWireframe"){
        // Get general properties
        int size = configuration["General"]["size"].as_int_or_die();
        vector<double> eye = configuration["General"]["eye"].as_double_tuple_or_die();
        vector<double> backgroundcolor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
        int nrFigures = configuration["General"]["nrFigures"].as_int_or_die();
        // Draw the wireframe
        Figures3D figures = drawWireframe(size , eye , backgroundcolor , nrFigures , configuration);
        Lines2D lines = doProjection(figures , 1.0);
        // Draw the lines
        return draw2DZbuffLines( lines , size , backgroundcolor);
    }
    // Case: type == "ZBuffering"
    else if(type == "ZBuffering"){
        // Get general properties
        int size = configuration["General"]["size"].as_int_or_die();
        vector<double> eye = configuration["General"]["eye"].as_double_tuple_or_die();
        vector<double> backgroundcolor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
        int nrFigures = configuration["General"]["nrFigures"].as_int_or_die();
        // Draw the wireframe
        Figures3D figures = drawWireframe(size , eye , backgroundcolor , nrFigures , configuration);
        return draw2DZbuffTriag(size , backgroundcolor , figures);
    }


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