#pragma clang diagnostic push
#pragma ide diagnostic ignored "bugprone-narrowing-conversions"

#include <fstream>

#include "ini_configuration.h"
#include "geometry/shapes.h"
#include "lighting/light.h"
#include "transformations.h"
#include "generators/platonicBodies.h"
#include "generators/lSystems3D.h"
#include "generators/draw2d.h"
#include "generators/basic2d.h"
#include "generators/zBuffering.h"
#include "generators/wireframes.h"
#include <corecrt_math_defines.h>

using namespace std;

void createInfinityLight(lights3D &lights, const color &ambientColor,
                         const color &specularColor, int i,
                         const ini::Configuration &configuration,
                         vector<double> &dirPoint, vector<double> &newDiffuse,
                         color &diffuseColor) {
    // Diffuse on infinity
    dirPoint = configuration["light" + to_string(
            i)]["direction"].as_double_tuple_or_die();
    newDiffuse = configuration["light" + to_string(
            i)]["diffuseLight"].as_double_tuple_or_default(
            {0, 0, 0});
    diffuseColor = color(newDiffuse);
    lights.push_back(
            new infLight(ambientColor, diffuseColor, specularColor,
                         dirPoint));
}

void createPointLight(lights3D &lights, const color &ambientColor, double spotAngle,
                      const color &specularColor, int i,
                      const ini::Configuration &configuration, vector<double> &dirPoint,
                      vector<double> &newDiffuse,
                      color &diffuseColor) {// Diffuse with spotlight
    newDiffuse = configuration["light" + to_string(
            i)]["diffuseLight"].as_double_tuple_or_default(
            {0, 0, 0});
    diffuseColor = color(newDiffuse);
    // Position of spotlight in real coordinates
    dirPoint = configuration["light" + to_string(
            i)]["location"].as_double_tuple_or_default({0, 0, 0});
    // Angle of the spot
    spotAngle = configuration["light" + to_string(
            i)]["spotAngle"].as_double_or_default(90.0);
    spotAngle *= (M_PI / 180);
    lights.push_back(
            new pointLight(ambientColor, diffuseColor, specularColor,
                           dirPoint, spotAngle));
}

void parseLights(const ini::Configuration &configuration, lights3D &lights) {
    // Ambient light
    vector<double> newAmbient;
    color ambientColor;
    // Diffuse light
    bool isDiffuse;
    bool infty;
    double spotAngle = 0;
    vector<double> dirPoint;
    vector<double> newDiffuse;
    color diffuseColor;
    // Specular light
    vector<double> newSpecular;
    color specularColor;
    // Number of lights
    int nrLights = configuration["General"]["nrLights"].as_int_or_default(0);

    for (int i = 0; i < nrLights; ++i) {
        // Ambient components
        newAmbient = configuration["light" + to_string(
                i)]["ambientLight"].as_double_tuple_or_default({1, 1, 1});
        ambientColor = color(newAmbient);
        // Specular components
        newSpecular = configuration["light" + to_string(
                i)]["specularLight"].as_double_tuple_or_default({0, 0, 0});
        specularColor = color(newSpecular);
        // Diffuse components
        infty = configuration["light" + to_string(i)]["infinity"].as_bool_if_exists(
                isDiffuse);
        if (infty) {
            if (isDiffuse) {
                createInfinityLight(lights, ambientColor,
                                    specularColor, i, configuration, dirPoint,
                                    newDiffuse, diffuseColor);
                continue;
            } else {
                createPointLight(lights, ambientColor, spotAngle, specularColor, i,
                                 configuration, dirPoint, newDiffuse,
                                 diffuseColor);
            }
        } else {
            lights.push_back(new light(ambientColor, diffuseColor, specularColor));
        }
    }
}

figures3D drawFigures(const ini::Configuration &configuration, vector<double> &eye,
                    int nrFigures, lights3D &lights) {// Draw the wireframe
    auto *newLight = new light(
            light(color(1.0, 1.0, 1.0), color(0, 0, 0), color(0, 0, 0)));
    lights = {newLight};
    figures3D figures = drawWireframe(eye, nrFigures,
                                      configuration, lights);
    return figures;
}

img::EasyImage generateImage(const ini::Configuration &configuration) {
    // Get type
    string type = configuration["General"]["type"].as_string_or_die();
    vector<double> backgroundColor;
    vector<double> eye;
    vector<double> lineColor;
    lights3D lights;
    int nrFigures;
    int size;

    // Case : type == "2DLSystem"
    if (type == "2DLSystem") {
        size = configuration["General"]["size"].as_int_or_die();
        backgroundColor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
        string inputFilename = configuration["2DLSystem"]["inputfile"].as_string_or_die();
        lineColor = configuration["2DLSystem"]["color"].as_double_tuple_or_die();
        // Initialize the parser
        LParser::LSystem2D lSystem;
        ifstream inputStream(inputFilename);
        inputStream >> lSystem;
        inputStream.close();
        return draw2DLines(drawSystem2D(lSystem, size, backgroundColor, lineColor),
                           size, backgroundColor);
    }
        // Case: type == "Wireframe"
    else if (type == "Wireframe") {
        // Get general properties
        size = configuration["General"]["size"].as_int_or_die();
        eye = configuration["General"]["eye"].as_double_tuple_or_die();
        backgroundColor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
        nrFigures = configuration["General"]["nrFigures"].as_int_or_die();
        figures3D figures = drawFigures(configuration, eye, nrFigures, lights);
        lines2D lines = doProjection(figures, 1.0);
        return draw2DLines(lines,
                           size, backgroundColor);
    }
        // Case: type == "ZBufferedWireframe"
    else if (type == "ZBufferedWireframe") {
        // Get general properties
        size = configuration["General"]["size"].as_int_or_die();
        eye = configuration["General"]["eye"].as_double_tuple_or_die();
        backgroundColor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
        nrFigures = configuration["General"]["nrFigures"].as_int_or_die();
        figures3D figures = drawFigures(configuration, eye, nrFigures, lights);
        lines2D lines = doProjection(figures, 1.0);
        return draw2DzBufferLines(lines, size, backgroundColor);

    }
        // Case: type == "ZBuffering"
    else if (type == "ZBuffering") {
        // Get general properties
        size = configuration["General"]["size"].as_int_or_die();
        eye = configuration["General"]["eye"].as_double_tuple_or_die();
        backgroundColor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
        nrFigures = configuration["General"]["nrFigures"].as_int_or_die();
        figures3D figures = drawFigures(configuration, eye, nrFigures, lights);
        return draw2DzBufferTriangle(size, backgroundColor, figures, lights);
    }
        // Case: type == "LightedZBuffering"
    else if (type == "LightedZBuffering") {
        // Get general properties
        size = configuration["General"]["size"].as_int_or_die();
        eye = configuration["General"]["eye"].as_double_tuple_or_die();
        backgroundColor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
        nrFigures = configuration["General"]["nrFigures"].as_int_or_die();
        // Get all the lights
        parseLights(configuration, lights);
        figures3D figures = drawWireframe(eye, nrFigures,
                                          configuration, lights);
        return draw2DzBufferTriangle(size, backgroundColor, figures, lights);
    }


    int width = configuration["ImageProperties"]["width"].as_int_or_die(); // Get width
    int height = configuration["ImageProperties"]["height"].as_int_or_die(); // Get height

    // Case : type = "IntroColorRectangle"
    if (type == "IntroColorRectangle") {
        return createColorRectangle(width, height);
    }
        // Case : type == "IntroBlocks"
    else if (type == "IntroBlocks") {
        int blocksInX = configuration["BlockProperties"]["nrXBlocks"].as_int_or_die();
        int blocksInY = configuration["BlockProperties"]["nrYBlocks"].as_int_or_die();
        vector<double> colorWhite = configuration["BlockProperties"]["colorWhite"].as_double_tuple_or_die();
        vector<double> colorBlack = configuration["BlockProperties"]["colorBlack"].as_double_tuple_or_die();
        bool invertColors = configuration["BlockProperties"]["invertColors"].as_bool_or_default(
                false);
        return createBlocks(width, height, blocksInX, blocksInY, colorWhite, colorBlack,
                            invertColors);
    }
        // Case : type == "IntroLines"
    else if (type == "IntroLines") {
        string figure = configuration["LineProperties"]["figure"].as_string_or_die();
        backgroundColor = configuration["LineProperties"]["backgroundcolor"].as_double_tuple_or_die();
        lineColor = configuration["LineProperties"]["lineColor"].as_double_tuple_or_die();
        int linesNumber = configuration["LineProperties"]["nrLines"].as_int_or_die();
        // Case: figure == "QuarterCircle"
        if (figure == "QuarterCircle") {
            return createQuarterCircle(width, height, linesNumber, backgroundColor,
                                       lineColor);
        }
        // Case: figure == "Eye"
        if (figure == "Eye") {
            return createEye(width, height, linesNumber, backgroundColor, lineColor);
        }
        // Case: figure == "Diamond"
        if (figure == "Diamond") {
            return createDiamond(width, height, linesNumber, backgroundColor,
                                 lineColor);
        }

    }

    return {};
}


int main(int argc, char const *argv[]) {
    int retVal = 0;
    try {
        std::vector<std::string> args = std::vector<std::string>(argv + 1, argv + argc);
        if (args.empty()) {
            std::ifstream fileIn("filelist");
            std::string filelistName;
            while (std::getline(fileIn, filelistName)) {
                args.push_back(filelistName);
            }
        }
        for (std::string fileName: args) {
            ini::Configuration conf;
            try {
                std::ifstream fin(fileName);
                fin >> conf;
                fin.close();
            }
            catch (ini::ParseException &ex) {
                std::cerr << "Error parsing file: " << fileName << ": " << ex.what()
                          << std::endl;
                retVal = 1;
                continue;
            }

            img::EasyImage image = generateImage(conf);
            if (image.get_height() > 0 && image.get_width() > 0) {
                std::string::size_type pos = fileName.rfind('.');
                if (pos == std::string::npos) {
                    //filename does not contain a '.' --> append a '.bmp' suffix
                    fileName += ".bmp";
                } else {
                    fileName = fileName.substr(0, pos) + ".bmp";
                }
                try {
                    std::ofstream fOut(fileName.c_str(),
                                       std::ios::trunc | std::ios::out |
                                       std::ios::binary);
                    fOut << image;

                }
                catch (std::exception &ex) {
                    std::cerr << "Failed to write image to file: " << ex.what()
                              << std::endl;
                    retVal = 1;
                }
            } else {
                std::cout << "Could not generate image for " << fileName << std::endl;
            }
        }
    }
    catch (const std::bad_alloc &exception) {
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