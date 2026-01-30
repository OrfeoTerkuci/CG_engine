#include "fractals.h"

void insertToPentagon(int indexA, int indexD, int indexI,
                      std::map<int, std::vector<int>> &pentagons)
{
    // The first pair of points
    if (pentagons[indexA].empty())
    {
        pentagons[indexA].push_back(indexD);
        pentagons[indexA].push_back(indexI);
        return;
    }
    else if (pentagons[indexA].size() == 5)
    {
        return;
    }
    // Last pair of points
    else if (pentagons[indexA].front() == indexD &&
             pentagons[indexA].at(1) != indexI)
    {
        pentagons[indexA].insert(pentagons[indexA].begin(), indexI);
        return;
    }
    else if (pentagons[indexA].front() == indexI &&
             pentagons[indexA].at(1) != indexD)
    {
        pentagons[indexA].insert(pentagons[indexA].begin(), indexD);
        return;
    }
    // Second pair of points
    else if (pentagons[indexA].back() == indexD &&
             pentagons[indexA].at(pentagons[indexA].size() - 1) == indexI)
    {
        pentagons[indexA].push_back(indexI);
        return;
    }
    else if (pentagons[indexA].back() == indexI &&
             pentagons[indexA].at(pentagons[indexA].size() - 1) == indexD)
    {
        pentagons[indexA].push_back(indexD);
        return;
    }
    // Pair in the middle
    else
    {
        for (int j = 0; j < pentagons[indexA].size(); j++)
        {
            if (pentagons[indexA].at(j) == indexD)
            {
                // Insert after this point
                pentagons[indexA].insert(pentagons[indexA].begin() + j + 1, indexI);
                return;
            }
            else if (pentagons[indexA].at(j) == indexI)
            {
                // Insert before this point
                pentagons[indexA].insert(pentagons[indexA].begin() + j, indexD);
                return;
            }
        }
        // If line cannot be connected to previously present points
        pentagons[indexA].push_back(indexD);
        return;
    }
}

int isElementOf(std::vector<Vector3D> &points, Vector3D &p)
{
    for (int i = 0; i < points.size(); ++i)
    {
        if (points.at(i) == p)
        {
            return i;
        }
    }
    return static_cast<int>(points.size());
}

void splitTriangleHexagon(face &originalTriangle, figure *&originalFigure,
                          std::vector<face> &newFaces,
                          std::map<int, std::vector<int>> &pentagons,
                          std::vector<Vector3D> &newPoints)
{
    // Get the points of the triangle
    int indexA = originalTriangle.pointIndexes.at(0);
    int indexB = originalTriangle.pointIndexes.at(1);
    int indexC = originalTriangle.pointIndexes.at(2);
    Vector3D a = originalFigure->points.at(indexA);
    Vector3D b = originalFigure->points.at(indexB);
    Vector3D c = originalFigure->points.at(indexC);
    // Calculate new points
    Vector3D d = 2 * (a / 3) + b / 3;
    Vector3D e = 2 * (b / 3) + a / 3;
    Vector3D f = 2 * (b / 3) + c / 3;
    Vector3D g = 2 * (c / 3) + b / 3;
    Vector3D i = 2 * (a / 3) + c / 3;
    Vector3D h = 2 * (c / 3) + a / 3;
    // Push points to figure

    int indexD = isElementOf(newPoints, d);
    if (indexD == newPoints.size())
    {
        newPoints.push_back(d);
        indexD = static_cast<int>(newPoints.size() - 1);
    }

    int indexE = isElementOf(newPoints, e);
    if (indexE == newPoints.size())
    {
        newPoints.push_back(e);
        indexE = static_cast<int>(newPoints.size() - 1);
    }

    int indexF = isElementOf(newPoints, f);
    if (indexF == newPoints.size())
    {
        newPoints.push_back(f);
        indexF = static_cast<int>(newPoints.size() - 1);
    }
    int indexG = isElementOf(newPoints, g);
    if (indexG == newPoints.size())
    {
        newPoints.push_back(g);
        indexG = static_cast<int>(newPoints.size() - 1);
    }
    int indexH = isElementOf(newPoints, h);
    if (indexH == newPoints.size())
    {
        newPoints.push_back(h);
        indexH = static_cast<int>(newPoints.size() - 1);
    }
    int indexI = isElementOf(newPoints, i);
    if (indexI == newPoints.size())
    {
        newPoints.push_back(i);
        indexI = static_cast<int>(newPoints.size() - 1);
    }

    // Create hexagon
    newFaces.push_back(face({indexD, indexE, indexF, indexG, indexH, indexI}));
    // Insert new points connected to a
    insertToPentagon(indexA, indexD, indexI, pentagons);
    // Insert new points connected to b
    insertToPentagon(indexB, indexF, indexE, pentagons);
    // Insert new points connected to c
    insertToPentagon(indexC, indexH, indexG, pentagons);
}

figure *
createBuckyBall(std::vector<double> &ambientCoefficient,
                std::vector<double> &diffuseCoefficient,
                std::vector<double> &specularCoefficient,
                double &reflectionCoefficient)
{

    // Create an icosahedron
    figure *newIcoSphere = createIcosahedron(ambientCoefficient, diffuseCoefficient,
                                             specularCoefficient,
                                             reflectionCoefficient);
    // Create temporary faces std::vector
    std::vector<face> newFaces;
    // Create temporary points std::vector
    std::vector<Vector3D> points;
    // Create pentagons
    std::map<int, std::vector<int>> pentagons;
    for (int i = 0; i < newIcoSphere->points.size(); i++)
    {
        pentagons[i] = {};
    }
    // Split each triangle into 4
    for (face &f : newIcoSphere->faces)
    {
        splitTriangleHexagon(f, newIcoSphere, newFaces, pentagons, points);
    }
    // Create new pentagon faces
    for (auto &f : pentagons)
    {
        newFaces.push_back(new face(f.second));
    }
    newIcoSphere->faces = newFaces;
    newIcoSphere->points = points;
    // Rescale all the points
    for (Vector3D &p : newIcoSphere->points)
    {
        p.normalise();
    }
    return newIcoSphere;
}

figures3D
createMengerSponge(int nrIterations, Matrix &m,
                   std::vector<double> &ambientCoefficient,
                   std::vector<double> &diffuseCoefficient,
                   std::vector<double> &specularCoefficient,
                   double &reflectionCoefficient)
{
    // Start with cube
    figures3D newSponge;
    figure *newFig = createCube(ambientCoefficient, diffuseCoefficient,
                                specularCoefficient,
                                reflectionCoefficient);

    newFig->applyTransformation(m);
    newSponge.push_back(newFig);
    // Divide every face into nine squares
    Matrix mS = scaleFigure((double)1 / 3);
    // Create empty translating matrix
    Matrix mT;
    // Copy begin point
    figures3D newFractal;
    // For each iteration
    for (int i = 0; i < nrIterations; ++i)
    {
        // For each figure
        for (auto &it : newSponge)
        {
            // For each point
            for (int j = 0; j < it->points.size(); ++j)
            {
                // Create corner cube
                // Copy the original figure
                newFig = new figure(it);
                // Scale the figure
                newFig->applyTransformation(mS);
                // Get translation matrix
                mT = translate(it->points.at(j) - newFig->points.at(j));
                // Translate the points
                newFig->applyTransformation(mT);
                // Add figure to list of fractals
                newFractal.push_back(newFig);
            }
        }
        newSponge = newFractal;
        newFractal = {};
    }
    return newSponge;
}

void generateFractal(figures3D &fractal, const int nrIterations, const double scale)
{

    figure *newFig;
    // Create scaling matrix
    Matrix mS = scaleFigure(1 / scale);
    // Create empty translating matrix
    Matrix mT;
    // Copy begin point
    figures3D newFractal;
    // For each iteration
    for (int i = 0; i < nrIterations; ++i)
    {
        // For each figure
        for (auto &it : fractal)
        {
            // For each point
            for (int j = 0; j < it->points.size(); ++j)
            {
                // Copy the original figure
                newFig = new figure(it);
                // Scale the figure
                newFig->applyTransformation(mS);
                // Get translation matrix
                mT = translate(it->points.at(j) - newFig->points.at(j));
                // Translate the points
                newFig->applyTransformation(mT);
                // Add figure to list of fractals
                newFractal.push_back(newFig);
            }
        }
        fractal = newFractal;
        newFractal = {};
    }
}
