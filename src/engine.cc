#pragma clang diagnostic push
#pragma ide diagnostic ignored "bugprone-narrowing-conversions"

#include <fstream>

#include "ini_configuration.h"
#include "geometry/shapes.h"
#include "lighting/light.h"
#include "transformations.h"
#include "clipping.h"
#include "generators/platonicBodies.h"
#include "generators/lSystems3D.h"
#include "generators/fractals.h"
#include "generators/draw3d.h"
#include "generators/draw2d.h"
#include "generators/basic2d.h"
#include "generators/zBuffering.h"

using namespace std;

figures3D drawWireframe(int &size, vector<double> &eye, vector<double> &backgroundColor,
                        int &nrFigures,
                        const ini::Configuration &configuration, lights3D &lights) {
    // Check if lighting
    string type = configuration["General"]["type"];
    vector<double> ambientCoefficient;
    vector<double> diffuseCoefficient;
    vector<double> specularCoefficient;
    bool lighting = (type == "LightedZBuffering");
    // Check for view cone
    bool viewCone = configuration["General"]["clipping"].as_bool_or_default(false);
    Vector3D viewDir;
    double dNear;
    double dFar;
    double hFov;
    double aspectRatio;
    // Make the variables
    double theta;
    double phi;
    double r;
    Vector3D eyePoint = Vector3D::point(eye.at(0), eye.at(1), eye.at(2));

    // Get view viewCone components
    if (viewCone) {
        // Get view viewCone components
        vector<double> viewDirection = configuration["General"]["viewDirection"].as_double_tuple_or_die();
        dNear = configuration["General"]["dNear"].as_double_or_default(1.0);
        dFar = configuration["General"]["dFar"].as_double_or_default(1000.0);
        hFov = configuration["General"]["hFov"].as_double_or_default(90);
        aspectRatio = configuration["General"]["aspectRatio"];

        viewDir = Vector3D::point(viewDirection.at(0), viewDirection.at(1),
                                  viewDirection.at(2));
    } else {
        viewDir = -eyePoint;
    }
    Matrix mEye = eyePointTrans(eyePoint, viewDir, theta, phi, r);

    // Transform the lights
    for (auto *l: lights) {
        auto infL = dynamic_cast<infLight *>(l);
        if (infL != nullptr) {
            infL->ldVector *= mEye;
            infL->ldVector.normalise();
        }
        auto pntL = dynamic_cast<pointLight *>(l);
        if (pntL != nullptr) {
            pntL->location *= mEye;
        }
    }

    // Create figures vector
    figures3D figures;
    // Get figures
    for (int i = 0; i < nrFigures; i++) {
        // Get all attributes
        // Get figure type
        string figureType = configuration["figure" +
                                          to_string(i)]["type"].as_string_or_die();
        // Get common attributes
        double scale = configuration["figure" +
                                     to_string(i)]["scale"].as_double_or_default(1.0);
        double rotX = configuration["figure" +
                                    to_string(i)]["rotateX"].as_double_or_default(0);
        double rotY = configuration["figure" +
                                    to_string(i)]["rotateY"].as_double_or_default(0);
        double rotZ = configuration["figure" +
                                    to_string(i)]["rotateZ"].as_double_or_default(0);
        vector<double> center = configuration["figure" + to_string(
                i)]["center"].as_double_tuple_or_default({0, 0, 0});

        double reflectionCoefficient = configuration["figure" +
                                                     to_string(
                                                             i)]["reflectionCoefficient"].as_double_or_default(
                0);

        if (lighting) {
            ambientCoefficient = configuration["figure" + to_string(
                    i)]["ambientReflection"].as_double_tuple_or_default(
                    {1, 1, 1});

        } else {
            ambientCoefficient = configuration["figure" + to_string(i)]["color"];
        }
        diffuseCoefficient = configuration["figure" + to_string(
                i)]["diffuseReflection"].as_double_tuple_or_default(
                {0, 0, 0});
        specularCoefficient = configuration["figure" + to_string(
                i)]["specularReflection"].as_double_tuple_or_default(
                {0, 0, 0});

        Matrix m = scaleFigure(scale);
        m *= rotateX((rotX * M_PI) / 180);
        m *= rotateY((rotY * M_PI) / 180);
        m *= rotateZ((rotZ * M_PI) / 180);
        m *= translate(Vector3D::point(center.at(0), center.at(1), center.at(2)));
        m *= mEye;

        // Figure_type : "LineDrawing"
        if (figureType == "LineDrawing") {
            // Get attributes specific for LineDrawing
            int nrPoints = configuration["figure" +
                                         to_string(i)]["nrPoints"].as_int_or_die();
            int nrLines = configuration["figure" +
                                        to_string(i)]["nrLines"].as_int_or_die();
            // Create new LineDrawing figure
            figure *newFigure = drawLineDrawing(scale, rotX, rotY, rotZ, nrPoints,
                                                nrLines, configuration,
                                                ambientCoefficient, center, mEye, i);
            // Add figure to vector of figures
            figures.push_back(newFigure);
        }
            // Figure_type : "Cube"
        else if (figureType == "Cube") {
            // Create new figure
            figure *newFigure;
            newFigure = createCube(ambientCoefficient, diffuseCoefficient,
                                   specularCoefficient, reflectionCoefficient);
            newFigure->applyTransformation(m);
            figures.push_back(newFigure);
        }
            // Figure_type == "Tetrahedron"
        else if (figureType == "Tetrahedron") {
            // Create new figure
            figure *newFigure;
            newFigure = createTetrahedron(ambientCoefficient, diffuseCoefficient,
                                          specularCoefficient,
                                          reflectionCoefficient);
            newFigure->applyTransformation(m);

            figures.push_back(newFigure);
        }
            // Figure_type == "Octahedron"
        else if (figureType == "Octahedron") {
            // Create new figure
            figure *newFigure;
            newFigure = createOctahedron(ambientCoefficient, diffuseCoefficient,
                                         specularCoefficient,
                                         reflectionCoefficient);
            newFigure->applyTransformation(m);
            figures.push_back(newFigure);
        }
            // Figure_type == "Icosahedron"
        else if (figureType == "Icosahedron") {
            // Create new figure
            figure *newFigure;
            newFigure = createIcosahedron(ambientCoefficient, diffuseCoefficient,
                                          specularCoefficient,
                                          reflectionCoefficient);
            newFigure->applyTransformation(m);
            figures.push_back(newFigure);
        }
            // Figure_type == "Dodecahedron"
        else if (figureType == "Dodecahedron") {
            // Create new figure
            figure *newFigure;
            newFigure = createDodecahedron(ambientCoefficient, diffuseCoefficient,
                                           specularCoefficient,
                                           reflectionCoefficient);
            newFigure->applyTransformation(m);
            figures.push_back(newFigure);
        }
            // Figure_type == "Cone"
        else if (figureType == "Cone") {
            int n = configuration["figure" + to_string(i)]["n"].as_int_or_die();
            double height = configuration["figure" +
                                          to_string(i)]["height"].as_double_or_die();
            figure *newFigure;
            newFigure = createCone(n, height, ambientCoefficient, diffuseCoefficient,
                                   specularCoefficient,
                                   reflectionCoefficient);
            newFigure->applyTransformation(m);
            figures.push_back(newFigure);
        }
            // Figure_type == "Cylinder"
        else if (figureType == "Cylinder") {
            const int n = configuration["figure" + to_string(i)]["n"].as_int_or_die();
            const double height = configuration["figure" + to_string(
                    i)]["height"].as_double_or_die();
            figure *newFigure;
            newFigure = createCylinder(n, height, ambientCoefficient,
                                       diffuseCoefficient, specularCoefficient,
                                       reflectionCoefficient);
            newFigure->applyTransformation(m);
            figures.push_back(newFigure);
        }
            // Figure_type == "Sphere"
        else if (figureType == "Sphere") {
            const int n = configuration["figure" + to_string(i)]["n"].as_int_or_die();
            figure *newFigure;
            newFigure = createSphere(1.0, n, ambientCoefficient, diffuseCoefficient,
                                     specularCoefficient,
                                     reflectionCoefficient);
            newFigure->applyTransformation(m);
            figures.push_back(newFigure);
        }
            // Figure_type == "Torus"
        else if (figureType == "Torus") {
            double rBig = configuration["figure" + to_string(i)]["R"].as_double_or_die();
            double torusR = configuration["figure" + to_string(i)]["r"].as_double_or_die();
            int torusN = configuration["figure" + to_string(i)]["n"].as_int_or_die();
            int torusM = configuration["figure" + to_string(i)]["m"].as_int_or_die();
            figure *newFigure;
            newFigure = createTorus(torusR, rBig, torusN, torusM, ambientCoefficient,
                                    diffuseCoefficient, specularCoefficient,
                                    reflectionCoefficient);
            newFigure->applyTransformation(m);
            figures.push_back(newFigure);
        }
            // Figure_type == "3DLSystem"
        else if (figureType == "3DLSystem") {
            string inputFile = configuration["figure" + to_string(
                    i)]["inputfile"].as_string_or_die();
            // Initialize the parser
            LParser::LSystem3D lSystem;
            ifstream inputStream(inputFile);
            inputStream >> lSystem;
            inputStream.close();
            figure *newFigure;
            newFigure = drawSystem3D(lSystem, ambientCoefficient);
            newFigure->applyTransformation(m);
            figures.push_back(newFigure);
        }
            // Figure_type == "BuckyBall"
        else if (figureType == "BuckyBall") {
            figure *newFigure;
            newFigure = createBuckyBall(ambientCoefficient, diffuseCoefficient,
                                        specularCoefficient,
                                        reflectionCoefficient);
            newFigure->applyTransformation(m);
            figures.push_back(newFigure);
        }
            // Figure_type = "MengerSponge"
        else if (figureType == "MengerSponge") {
            int nrIterations = configuration["figure" + to_string(
                    i)]["nrIterations"].as_int_or_default(0);
            figures3D newSponge = createMengerSponge(nrIterations, m,
                                                     ambientCoefficient,
                                                     diffuseCoefficient,
                                                     specularCoefficient,
                                                     reflectionCoefficient);
            figures.insert(figures.end(), newSponge.begin(), newSponge.end());
        }
            // Figure_type == "FractalCube"
        else if (figureType == "FractalCube") {
            int nrIterations = configuration["figure" + to_string(
                    i)]["nrIterations"].as_int_or_default(0);
            double fractalScale = configuration["figure" + to_string(
                    i)]["fractalScale"].as_double_or_default(1);
            figure *newFigure = createCube(ambientCoefficient, diffuseCoefficient,
                                           specularCoefficient,
                                           reflectionCoefficient);
            newFigure->applyTransformation(m);
            figures3D newFractal = {newFigure};
            generateFractal(newFractal, nrIterations, fractalScale);
            figures.insert(figures.end(), newFractal.begin(), newFractal.end());
        }
            // Figure_type == "FractalTetrahedron"
        else if (figureType == "FractalTetrahedron") {
            int nrIterations = configuration["figure" + to_string(
                    i)]["nrIterations"].as_int_or_default(0);
            double fractalScale = configuration["figure" + to_string(
                    i)]["fractalScale"].as_double_or_default(1);
            figure *newFigure = createTetrahedron(ambientCoefficient,
                                                  diffuseCoefficient,
                                                  specularCoefficient,
                                                  reflectionCoefficient);
            newFigure->applyTransformation(m);
            figures3D newFractal = {newFigure};
            generateFractal(newFractal, nrIterations, fractalScale);
            figures.insert(figures.end(), newFractal.begin(), newFractal.end());
        }
            // Figure_type == "FractalIcosahedron"
        else if (figureType == "FractalIcosahedron") {
            int nrIterations = configuration["figure" + to_string(
                    i)]["nrIterations"].as_int_or_default(0);
            double fractalScale = configuration["figure" + to_string(
                    i)]["fractalScale"].as_double_or_default(1);
            figure *newFigure = createIcosahedron(ambientCoefficient,
                                                  diffuseCoefficient,
                                                  specularCoefficient,
                                                  reflectionCoefficient);
            newFigure->applyTransformation(m);
            figures3D newFractal = {newFigure};
            generateFractal(newFractal, nrIterations, fractalScale);
            figures.insert(figures.end(), newFractal.begin(), newFractal.end());
        }
            // Figure_type == "FractalOctahedron"
        else if (figureType == "FractalOctahedron") {
            int nrIterations = configuration["figure" + to_string(
                    i)]["nrIterations"].as_int_or_default(0);
            double fractalScale = configuration["figure" + to_string(
                    i)]["fractalScale"].as_double_or_default(1);
            figure *newFigure = createOctahedron(ambientCoefficient, diffuseCoefficient,
                                                 specularCoefficient,
                                                 reflectionCoefficient);
            newFigure->applyTransformation(m);
            figures3D newFractal = {newFigure};
            generateFractal(newFractal, nrIterations, fractalScale);
            figures.insert(figures.end(), newFractal.begin(), newFractal.end());
        }
            // Figure_type == "FractalDodecahedron"
        else if (figureType == "FractalDodecahedron") {
            int nrIterations = configuration["figure" + to_string(
                    i)]["nrIterations"].as_int_or_default(0);
            double fractalScale = configuration["figure" + to_string(
                    i)]["fractalScale"].as_double_or_default(1);
            figure *newFigure = createDodecahedron(ambientCoefficient,
                                                   diffuseCoefficient,
                                                   specularCoefficient,
                                                   reflectionCoefficient);
            newFigure->applyTransformation(m);
            figures3D newFractal = {newFigure};
            generateFractal(newFractal, nrIterations, fractalScale);
            figures.insert(figures.end(), newFractal.begin(), newFractal.end());
        }
            // Figure_type == "FractalBuckyBall"
        else if (figureType == "FractalBuckyBall") {
            int nrIterations = configuration["figure" + to_string(
                    i)]["nrIterations"].as_int_or_default(0);
            double fractalScale = configuration["figure" + to_string(
                    i)]["fractalScale"].as_double_or_default(1);
            figure *newFigure = createBuckyBall(ambientCoefficient, diffuseCoefficient,
                                                specularCoefficient,
                                                reflectionCoefficient);
            newFigure->applyTransformation(m);
            figures3D newFractal = {newFigure};
            generateFractal(newFractal, nrIterations, fractalScale);
            figures.insert(figures.end(), newFractal.begin(), newFractal.end());
        }
    }
    // Clip view
    if (viewCone) {
        clipView(figures, dNear, dFar, hFov, aspectRatio);
    }
    return figures;
}

// Session 4 : Z-Buffering


img::EasyImage generateImage(const ini::Configuration &configuration) {
    // Get type
    string type = configuration["General"]["type"].as_string_or_die();
    // Case : type == "2DLSystem"
    if (type == "2DLSystem") {
        int size = configuration["General"]["size"].as_int_or_die();
        vector<double> backgroundColor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
        string inputFilename = configuration["2DLSystem"]["inputfile"].as_string_or_die();
        vector<double> lineColor = configuration["2DLSystem"]["color"].as_double_tuple_or_die();
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
        int size = configuration["General"]["size"].as_int_or_die();
        vector<double> eye = configuration["General"]["eye"].as_double_tuple_or_die();
        vector<double> backgroundColor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
        int nrFigures = configuration["General"]["nrFigures"].as_int_or_die();
        // Draw the wireframe
        auto *newLight = new light(
                light(color(1.0, 1.0, 1.0), color(0, 0, 0), color(0, 0, 0)));
        lights3D lights = {newLight};
        figures3D figures = drawWireframe(size, eye, backgroundColor, nrFigures,
                                          configuration, lights);
        lines2D lines = doProjection(figures, 1.0);
        return draw2DLines(lines,
                           size, backgroundColor);
    }
        // Case: type == "ZBufferedWireframe"
    else if (type == "ZBufferedWireframe") {
        // Get general properties
        int size = configuration["General"]["size"].as_int_or_die();
        vector<double> eye = configuration["General"]["eye"].as_double_tuple_or_die();
        vector<double> backgroundColor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
        int nrFigures = configuration["General"]["nrFigures"].as_int_or_die();
        // Draw the wireframe
        auto *newLight = new light(
                light(color(1.0, 1.0, 1.0), color(0, 0, 0), color(0, 0, 0)));
        lights3D lights = {newLight};
        figures3D figures = drawWireframe(size, eye, backgroundColor, nrFigures,
                                          configuration, lights);
        lines2D lines = doProjection(figures, 1.0);
        // Draw the lines
        string btype = configuration["General"]["blurType"].as_string_or_default(
                "None");

        return draw2DzBufferLines(lines, size, backgroundColor);

    }
        // Case: type == "ZBuffering"
    else if (type == "ZBuffering") {
        // Get general properties
        int size = configuration["General"]["size"].as_int_or_die();
        vector<double> eye = configuration["General"]["eye"].as_double_tuple_or_die();
        vector<double> backgroundColor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
        int nrFigures = configuration["General"]["nrFigures"].as_int_or_die();
        // Draw the wireframe
        auto *newLight = new light(
                light(color(1.0, 1.0, 1.0), color(0, 0, 0), color(0, 0, 0)));
        lights3D lights = {newLight};
        figures3D figures = drawWireframe(size, eye, backgroundColor, nrFigures,
                                          configuration, lights);
        //
        return draw2DzBufferTriangle(size, backgroundColor, figures, lights);
    }
        // Case: type == "LightedZBuffering"
    else if (type == "LightedZBuffering") {
        // Get general properties
        int size = configuration["General"]["size"].as_int_or_die();
        vector<double> eye = configuration["General"]["eye"].as_double_tuple_or_die();
        vector<double> backgroundColor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
        int nrFigures = configuration["General"]["nrFigures"].as_int_or_die();
        int nrLights = configuration["General"]["nrLights"].as_int_or_default(0);
        // Get all the lights
        lights3D lights;
        light *newLight;
        // Ambient light
        vector<double> newAmbient;
        color ambientColor;
        // Diffuse light
        bool isDiffuse;
        bool infty;
        double spotAngle;
        vector<double> dirPoint;
        vector<double> newDiffuse;
        color diffuseColor;
        // Specular light
        vector<double> newSpecular;
        color specularColor;

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
                    continue;
                } else {
                    // Diffuse with spotlight
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
            } else {
                newDiffuse = {0, 0, 0};
                lights.push_back(new light(ambientColor, diffuseColor, specularColor));
            }
        }
        figures3D figures = drawWireframe(size, eye, backgroundColor, nrFigures,
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
        vector<double> backgroundColor = configuration["LineProperties"]["backgroundcolor"].as_double_tuple_or_die();
        vector<double> lineColor = configuration["LineProperties"]["lineColor"].as_double_tuple_or_die();
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