#include "wireframes.h"

void getViewDir(bool viewCone, const Vector3D &eyePoint,
                const ini::Configuration &configuration, Vector3D &viewDir,
                double &dNear, double &dFar, double &hFov,
                double &aspectRatio)
{
    // Get view viewCone components
    if (viewCone)
    {
        // Get view viewCone components
        std::vector<double> viewDirection = configuration["General"]["viewDirection"].as_double_tuple_or_die();
        dNear = configuration["General"]["dNear"].as_double_or_default(1.0);
        dFar = configuration["General"]["dFar"].as_double_or_default(1000.0);
        hFov = configuration["General"]["hFov"].as_double_or_default(90);
        aspectRatio = configuration["General"]["aspectRatio"];

        viewDir = Vector3D::point(viewDirection.at(0), viewDirection.at(1),
                                  viewDirection.at(2));
    }
    else
    {
        viewDir = -eyePoint;
    }
}

void transformLights(lights3D &lights, const Matrix &mEye)
{
    // Transform the lights
    for (auto *l : lights)
    {
        auto infL = dynamic_cast<infLight *>(l);
        if (infL != nullptr)
        {
            infL->ldVector *= mEye;
            infL->ldVector.normalise();
        }
        auto pntL = dynamic_cast<pointLight *>(l);
        if (pntL != nullptr)
        {
            pntL->location *= mEye;
        }
    }
}

void parseFigureAttributes(bool lighting, int i, const ini::Configuration &configuration,
                           std::vector<double> &ambientCoefficient,
                           std::vector<double> &diffuseCoefficient,
                           std::vector<double> &specularCoefficient, std::string &figureType,
                           double &scale, double &rotX, double &rotY, double &rotZ,
                           std::vector<double> &center, double &reflectionCoefficient)
{
    figureType = configuration["figure" +
                               std::to_string(i)]["type"]
                     .as_string_or_die();
    scale = configuration["figure" +
                          std::to_string(i)]["scale"]
                .as_double_or_default(
                    1.0);
    rotX = configuration["figure" +
                         std::to_string(i)]["rotateX"]
               .as_double_or_default(0);
    rotY = configuration["figure" +
                         std::to_string(i)]["rotateY"]
               .as_double_or_default(0);
    rotZ = configuration["figure" +
                         std::to_string(i)]["rotateZ"]
               .as_double_or_default(0);
    center = configuration["figure" + std::to_string(
                                          i)]["center"]
                 .as_double_tuple_or_default({0, 0, 0});
    reflectionCoefficient = configuration["figure" + std::to_string(i)]["reflectionCoefficient"].as_double_or_default(0); // Get all attributes
                                                                                                                          // Get figure type
                                                                                                                          // Get common attributes
    if (lighting)
    {
        ambientCoefficient = configuration["figure" + std::to_string(i)]["ambientReflection"].as_double_tuple_or_default({1, 1, 1});
    }
    else
    {
        ambientCoefficient = configuration["figure" + std::to_string(i)]["color"];
    }
    diffuseCoefficient = configuration["figure" + std::to_string(
                                                      i)]["diffuseReflection"]
                             .as_double_tuple_or_default(
                                 {0, 0, 0});
    specularCoefficient = configuration["figure" + std::to_string(
                                                       i)]["specularReflection"]
                              .as_double_tuple_or_default(
                                  {0, 0, 0});
}

figures3D drawWireframe(std::vector<double> &eye, int &nrFigures,
                        const ini::Configuration &configuration, lights3D &lights)
{
    // Check if lighting
    std::string type = configuration["General"]["type"];
    std::vector<double> ambientCoefficient;
    std::vector<double> diffuseCoefficient;
    std::vector<double> specularCoefficient;
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

    getViewDir(viewCone, eyePoint, configuration, viewDir, dNear, dFar, hFov,
               aspectRatio);

    Matrix mEye = eyePointTrans(eyePoint, viewDir, theta, phi, r);
    transformLights(lights, mEye);

    // Create figures vector
    figures3D figures;
    // Get figures
    for (int i = 0; i < nrFigures; i++)
    {
        figure *newFigure;
        std::string figureType;
        double scale;
        double rotX;
        double rotY;
        double rotZ;
        std::vector<double> center;
        double reflectionCoefficient;

        parseFigureAttributes(
            lighting, i, configuration, ambientCoefficient, diffuseCoefficient,
            specularCoefficient, figureType, scale, rotX, rotY, rotZ, center,
            reflectionCoefficient);

        Matrix m = scaleFigure(scale);
        m *= rotateX((rotX * M_PI) / 180);
        m *= rotateY((rotY * M_PI) / 180);
        m *= rotateZ((rotZ * M_PI) / 180);
        m *= translate(Vector3D::point(center.at(0), center.at(1), center.at(2)));
        m *= mEye;

        // Figure_type : "LineDrawing"
        if (figureType == "LineDrawing")
        {
            // Get attributes specific for LineDrawing
            int nrPoints = configuration["figure" +
                                         std::to_string(i)]["nrPoints"]
                               .as_int_or_die();
            int nrLines = configuration["figure" +
                                        std::to_string(i)]["nrLines"]
                              .as_int_or_die();
            // Create new LineDrawing figure
            newFigure = drawLineDrawing(scale, rotX, rotY, rotZ, nrPoints,
                                        nrLines, configuration,
                                        ambientCoefficient, center, mEye, i);
            // Add figure to vector of figures
            figures.push_back(newFigure);
        }
        // Figure_type : "Cube"
        else if (figureType == "Cube")
        {
            // Create new figure
            newFigure = createCube(ambientCoefficient, diffuseCoefficient,
                                   specularCoefficient, reflectionCoefficient);
            newFigure->applyTransformation(m);
            figures.push_back(newFigure);
        }
        // Figure_type == "Tetrahedron"
        else if (figureType == "Tetrahedron")
        {
            // Create new figure
            newFigure = createTetrahedron(ambientCoefficient, diffuseCoefficient,
                                          specularCoefficient,
                                          reflectionCoefficient);
            newFigure->applyTransformation(m);

            figures.push_back(newFigure);
        }
        // Figure_type == "Octahedron"
        else if (figureType == "Octahedron")
        {
            // Create new figure
            newFigure = createOctahedron(ambientCoefficient, diffuseCoefficient,
                                         specularCoefficient,
                                         reflectionCoefficient);
            newFigure->applyTransformation(m);
            figures.push_back(newFigure);
        }
        // Figure_type == "Icosahedron"
        else if (figureType == "Icosahedron")
        {
            // Create new figure
            newFigure = createIcosahedron(ambientCoefficient, diffuseCoefficient,
                                          specularCoefficient,
                                          reflectionCoefficient);
            newFigure->applyTransformation(m);
            figures.push_back(newFigure);
        }
        // Figure_type == "Dodecahedron"
        else if (figureType == "Dodecahedron")
        {
            // Create new figure
            newFigure = createDodecahedron(ambientCoefficient, diffuseCoefficient,
                                           specularCoefficient,
                                           reflectionCoefficient);
            newFigure->applyTransformation(m);
            figures.push_back(newFigure);
        }
        // Figure_type == "Cone"
        else if (figureType == "Cone")
        {
            const int n = configuration["figure" +
                                        std::to_string(i)]["n"]
                              .as_int_or_die();
            const double height = configuration["figure" +
                                                std::to_string(
                                                    i)]["height"]
                                      .as_double_or_die();
            newFigure = createCone(n, height, ambientCoefficient, diffuseCoefficient,
                                   specularCoefficient,
                                   reflectionCoefficient);
            newFigure->applyTransformation(m);
            figures.push_back(newFigure);
        }
        // Figure_type == "Cylinder"
        else if (figureType == "Cylinder")
        {
            const int n = configuration["figure" +
                                        std::to_string(i)]["n"]
                              .as_int_or_die();
            const double height = configuration["figure" + std::to_string(
                                                               i)]["height"]
                                      .as_double_or_die();
            newFigure = createCylinder(n, height, ambientCoefficient,
                                       diffuseCoefficient, specularCoefficient,
                                       reflectionCoefficient);
            newFigure->applyTransformation(m);
            figures.push_back(newFigure);
        }
        // Figure_type == "Sphere"
        else if (figureType == "Sphere")
        {
            const int n = configuration["figure" +
                                        std::to_string(i)]["n"]
                              .as_int_or_die();
            newFigure = createSphere(1.0, n, ambientCoefficient, diffuseCoefficient,
                                     specularCoefficient,
                                     reflectionCoefficient);
            newFigure->applyTransformation(m);
            figures.push_back(newFigure);
        }
        // Figure_type == "Torus"
        else if (figureType == "Torus")
        {
            double majorRadius = configuration["figure" +
                                               std::to_string(
                                                   i)]["R"]
                                     .as_double_or_die();
            double minorRadius = configuration["figure" +
                                               std::to_string(
                                                   i)]["r"]
                                     .as_double_or_die();
            int torusRadialSegments = configuration["figure" +
                                                    std::to_string(
                                                        i)]["n"]
                                          .as_int_or_die();
            int torusTubularSegments = configuration["figure" +
                                                     std::to_string(
                                                         i)]["m"]
                                           .as_int_or_die();
            newFigure = createTorus(minorRadius, majorRadius, torusRadialSegments,
                                    torusTubularSegments, ambientCoefficient,
                                    diffuseCoefficient, specularCoefficient,
                                    reflectionCoefficient);
            newFigure->applyTransformation(m);
            figures.push_back(newFigure);
        }
        // Figure_type == "3DLSystem"
        else if (figureType == "3DLSystem")
        {
            std::string inputFile = configuration["figure" + std::to_string(
                                                                 i)]["inputfile"]
                                        .as_string_or_die();
            // Initialize the parser
            LParser::LSystem3D lSystem;
            std::ifstream inputStream(inputFile);
            inputStream >> lSystem;
            inputStream.close();
            newFigure = drawSystem3D(lSystem, ambientCoefficient);
            newFigure->applyTransformation(m);
            figures.push_back(newFigure);
        }
        // Figure_type == "BuckyBall"
        else if (figureType == "BuckyBall")
        {
            newFigure = createBuckyBall(ambientCoefficient, diffuseCoefficient,
                                        specularCoefficient,
                                        reflectionCoefficient);
            newFigure->applyTransformation(m);
            figures.push_back(newFigure);
        }
        // Figure_type = "MengerSponge"
        else if (figureType == "MengerSponge")
        {
            int nrIterations = configuration["figure" + std::to_string(
                                                            i)]["nrIterations"]
                                   .as_int_or_default(0);
            figures3D newSponge = createMengerSponge(nrIterations, m,
                                                     ambientCoefficient,
                                                     diffuseCoefficient,
                                                     specularCoefficient,
                                                     reflectionCoefficient);
            figures.insert(figures.end(), newSponge.begin(), newSponge.end());
        }
        // Figure_type == "FractalCube"
        else if (figureType == "FractalCube")
        {
            int nrIterations = configuration["figure" + std::to_string(
                                                            i)]["nrIterations"]
                                   .as_int_or_default(0);
            double fractalScale = configuration["figure" + std::to_string(
                                                               i)]["fractalScale"]
                                      .as_double_or_default(1);
            newFigure = createCube(ambientCoefficient, diffuseCoefficient,
                                   specularCoefficient,
                                   reflectionCoefficient);
            newFigure->applyTransformation(m);
            figures3D newFractal = {newFigure};
            generateFractal(newFractal, nrIterations, fractalScale);
            figures.insert(figures.end(), newFractal.begin(), newFractal.end());
        }
        // Figure_type == "FractalTetrahedron"
        else if (figureType == "FractalTetrahedron")
        {
            int nrIterations = configuration["figure" + std::to_string(
                                                            i)]["nrIterations"]
                                   .as_int_or_default(0);
            double fractalScale = configuration["figure" + std::to_string(
                                                               i)]["fractalScale"]
                                      .as_double_or_default(1);
            newFigure = createTetrahedron(ambientCoefficient,
                                          diffuseCoefficient,
                                          specularCoefficient,
                                          reflectionCoefficient);
            newFigure->applyTransformation(m);
            figures3D newFractal = {newFigure};
            generateFractal(newFractal, nrIterations, fractalScale);
            figures.insert(figures.end(), newFractal.begin(), newFractal.end());
        }
        // Figure_type == "FractalIcosahedron"
        else if (figureType == "FractalIcosahedron")
        {
            int nrIterations = configuration["figure" + std::to_string(
                                                            i)]["nrIterations"]
                                   .as_int_or_default(0);
            double fractalScale = configuration["figure" + std::to_string(
                                                               i)]["fractalScale"]
                                      .as_double_or_default(1);
            newFigure = createIcosahedron(ambientCoefficient,
                                          diffuseCoefficient,
                                          specularCoefficient,
                                          reflectionCoefficient);
            newFigure->applyTransformation(m);
            figures3D newFractal = {newFigure};
            generateFractal(newFractal, nrIterations, fractalScale);
            figures.insert(figures.end(), newFractal.begin(), newFractal.end());
        }
        // Figure_type == "FractalOctahedron"
        else if (figureType == "FractalOctahedron")
        {
            int nrIterations = configuration["figure" + std::to_string(
                                                            i)]["nrIterations"]
                                   .as_int_or_default(0);
            double fractalScale = configuration["figure" + std::to_string(
                                                               i)]["fractalScale"]
                                      .as_double_or_default(1);
            newFigure = createOctahedron(ambientCoefficient, diffuseCoefficient,
                                         specularCoefficient,
                                         reflectionCoefficient);
            newFigure->applyTransformation(m);
            figures3D newFractal = {newFigure};
            generateFractal(newFractal, nrIterations, fractalScale);
            figures.insert(figures.end(), newFractal.begin(), newFractal.end());
        }
        // Figure_type == "FractalDodecahedron"
        else if (figureType == "FractalDodecahedron")
        {
            int nrIterations = configuration["figure" + std::to_string(
                                                            i)]["nrIterations"]
                                   .as_int_or_default(0);
            double fractalScale = configuration["figure" + std::to_string(
                                                               i)]["fractalScale"]
                                      .as_double_or_default(1);
            newFigure = createDodecahedron(ambientCoefficient,
                                           diffuseCoefficient,
                                           specularCoefficient,
                                           reflectionCoefficient);
            newFigure->applyTransformation(m);
            figures3D newFractal = {newFigure};
            generateFractal(newFractal, nrIterations, fractalScale);
            figures.insert(figures.end(), newFractal.begin(), newFractal.end());
        }
        // Figure_type == "FractalBuckyBall"
        else if (figureType == "FractalBuckyBall")
        {
            int nrIterations = configuration["figure" + std::to_string(
                                                            i)]["nrIterations"]
                                   .as_int_or_default(0);
            double fractalScale = configuration["figure" + std::to_string(
                                                               i)]["fractalScale"]
                                      .as_double_or_default(1);
            newFigure = createBuckyBall(ambientCoefficient, diffuseCoefficient,
                                        specularCoefficient,
                                        reflectionCoefficient);
            newFigure->applyTransformation(m);
            figures3D newFractal = {newFigure};
            generateFractal(newFractal, nrIterations, fractalScale);
            figures.insert(figures.end(), newFractal.begin(), newFractal.end());
        }
    }
    // Clip view
    if (viewCone)
    {
        clipView(figures, dNear, dFar, hFov, aspectRatio);
    }
    return figures;
}