#pragma clang diagnostic push
#pragma ide diagnostic ignored "bugprone-narrowing-conversions"
#include "easy_image.h"
#include "ini_configuration.h"
#include "l_parser/l_parser.h"
#include "vector/vector3d.h"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <cmath>
#include <list>

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
    vector<int> point_indexes; // 2 for this exercise , 3+ later

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
            //image(i,j).red = i;
            //image(i,j).green = j;
            //image(i,j).blue = (i+j)%imageHeight;
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
    // Calculate image dimentions
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
            //createSystemLines(l_system,lines,startingString,endingString,startingAngle,angle,lineColor,currentPoint,i+1);
            //lines = createSystemLines(l_system,lines,startingString,endingString,startingAngle,angle,lineColor,currentPoint.x,currentPoint.y);
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

Lines2D drawSystem (const LParser::LSystem2D &l_system , const int &size , vector<double> &backgroundColor , vector<double> &lineColor) {
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

// Session 3 : 3D Lines
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

Lines2D doProjection(const Figures3D &figs){
    Lines2D lines;
    for(Figure *f : figs){
        for(Face face : f->faces){
            // Get points index
            int b_index = face.point_indexes.at(0);
            int e_index = face.point_indexes.at(1);
            // Get points
            Vector3D beginP = f->points[b_index];
            Vector3D endP = f->points[e_index];
            /*
            if(beginP.z == 0){
                continue;
            }
            if(endP.z == 0){
                continue;
            }
             */
            // Convert Vector3D to Point2D
            Point2D newBeginP = doProjection(beginP,1.0);
            Point2D newEndP = doProjection(endP,1.0);
            // Create new line
            Line2D newLine(newBeginP , newEndP , f->color);
            // Push line back to vector
            lines.push_back(newLine);
        }
    }
    return lines;
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
        return draw2DLines( drawSystem(l_system , size , backgroundColor , lineColor ) , size , backgroundColor );
    }

    else if (type == "Wireframe"){
        // Get general properties
        int size = configuration["General"]["size"].as_int_or_die();
        vector<double> eye = configuration["General"]["eye"].as_double_tuple_or_die();
        vector<double> backgroundcolor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
        int nrFigures = configuration["General"]["nrFigures"].as_int_or_die();
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
            string figure_type = configuration["Figure"+to_string(i)]["type"].as_string_or_default("LineDrawing");
            double scale = configuration["Figure"+to_string(i)]["scale"].as_double_or_default(1.0);
            double rotX = configuration["Figure"+to_string(i)]["rotateX"].as_double_or_default(0);
            double rotY = configuration["Figure"+to_string(i)]["rotateY"].as_double_or_default(0);
            double rotZ = configuration["Figure"+to_string(i)]["rotateZ"].as_double_or_default(0);
            vector<double> center = configuration["Figure"+to_string(i)]["center"].as_double_tuple_or_default({0,0,0});
            vector<double> lineColor = configuration["Figure"+to_string(i)]["color"].as_double_tuple_or_default({0,0,0});
            int nrPoints = configuration["Figure"+to_string(i)]["nrPoints"].as_int_or_die();
            int nrLines = configuration["Figure"+to_string(i)]["nrLines"].as_int_or_die();
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
                faces.push_back(Face({line_to_add.at(0) , line_to_add.at(1)}));
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
            // Add figure to vector of figures
            figures.push_back(newFigure);
        }
        return draw2DLines(doProjection(figures) , size , backgroundcolor);
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