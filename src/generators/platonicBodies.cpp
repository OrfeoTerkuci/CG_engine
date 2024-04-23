#include "platonicBodies.h"

figure *createCube(std::vector<double> &ambientCoefficient,
                   std::vector<double> &diffuseCoefficient,
                   std::vector<double> &specularCoefficient,
                   double &reflectionCoefficient)
{
    // Points array
    double pointsT[3][8] = {
        {1, -1, 1, -1, 1, -1, 1, -1},
        {-1, 1, 1, -1, 1, -1, -1, 1},
        {-1, -1, 1, 1, -1, -1, 1, 1}};
    // Faces array
    int facesT[4][6] = {
        {0, 4, 1, 5, 6, 0},
        {4, 1, 5, 0, 2, 5},
        {2, 7, 3, 6, 7, 1},
        {6, 2, 7, 3, 3, 4}};
    // Create all the points
    std::vector<Vector3D> points;
    points.reserve(8);
    for (int i = 0; i < 8; ++i)
    {
        points.push_back(Vector3D::point(pointsT[0][i], pointsT[1][i], pointsT[2][i]));
    }
    // Create all the faces
    std::vector<face> faces;
    faces.reserve(6);
    for (int j = 0; j < 6; ++j)
    {
        faces.push_back(face({facesT[0][j], facesT[1][j], facesT[2][j], facesT[3][j]}));
    }
    // Create new figure
    figure *newCube;
    newCube = new figure(points, faces, ambientCoefficient, diffuseCoefficient,
                         specularCoefficient, reflectionCoefficient);
    return newCube;
}

figure *
createTetrahedron(std::vector<double> &ambientCoefficient,
                  std::vector<double> &diffuseCoefficient,
                  std::vector<double> &specularCoefficient,
                  double &reflectionCoefficient)
{
    // Points array
    double pointsT[3][4] = {
        {1, -1, 1, -1},
        {-1, 1, 1, -1},
        {-1, -1, 1, 1}};
    // Faces array
    int facesT[3][4] = {
        {0, 1, 0, 0},
        {1, 3, 3, 2},
        {2, 2, 1, 3}};
    // Create all the points
    std::vector<Vector3D> points;
    points.reserve(4);
    for (int i = 0; i < 4; ++i)
    {
        points.push_back(
            Vector3D::point(pointsT[0][i], pointsT[1][i], pointsT[2][i]));
    }
    // Create faces
    std::vector<face> faces;
    faces.reserve(4);
    for (int j = 0; j < 4; ++j)
    {
        faces.push_back(face({facesT[0][j], facesT[1][j], facesT[2][j]}));
    }
    // Create new figure
    figure *newTetra;
    newTetra = new figure(points, faces, ambientCoefficient, diffuseCoefficient,
                          specularCoefficient,
                          reflectionCoefficient);
    return newTetra;
}

figure *
createOctahedron(std::vector<double> &ambientCoefficient,
                 std::vector<double> &diffuseCoefficient,
                 std::vector<double> &specularCoefficient,
                 double &reflectionCoefficient)
{
    // Points array
    double pointsT[3][6] = {
        {1, 0, -1, 0, 0, 0},
        {0, 1, 0, -1, 0, 0},
        {0, 0, 0, 0, -1, 1}};
    // Faces array
    int facesT[3][8] = {
        {0, 1, 2, 3, 1, 2, 3, 0},
        {1, 2, 3, 0, 0, 1, 2, 3},
        {5, 5, 5, 5, 4, 4, 4, 4}};
    // Create all the points
    std::vector<Vector3D> points;
    points.reserve(6);
    for (int i = 0; i < 6; ++i)
    {
        points.push_back(
            Vector3D::point(pointsT[0][i], pointsT[1][i], pointsT[2][i]));
    }
    // Create faces
    std::vector<face> faces;
    faces.reserve(8);
    for (int j = 0; j < 8; ++j)
    {
        faces.push_back(face({facesT[0][j], facesT[1][j], facesT[2][j]}));
    }
    // Create new figure
    figure *newOcta;
    newOcta = new figure(points, faces, ambientCoefficient, diffuseCoefficient,
                         specularCoefficient,
                         reflectionCoefficient);
    return newOcta;
}

figure *
createIcosahedron(std::vector<double> &ambientCoefficient,
                  std::vector<double> &diffuseCoefficient,
                  std::vector<double> &specularCoefficient,
                  double &reflectionCoefficient)
{
    // Points array
    double pointsT[3][12];
    pointsT[0][0] = 0;
    pointsT[1][0] = 0;
    pointsT[2][0] = sqrt(5) / 2;

    for (int k = 1; k < 6; ++k)
    {
        pointsT[0][k] = cos(2 * M_PI * (k - 1) / 5);
        pointsT[1][k] = sin(2 * M_PI * (k - 1) / 5);
        pointsT[2][k] = 0.5;
    }

    for (int l = 6; l < 11; ++l)
    {
        pointsT[0][l] = cos(M_PI / 5 + (l - 6) * (2 * M_PI) / 5);
        pointsT[1][l] = sin(M_PI / 5 + (l - 6) * (2 * M_PI) / 5);
        pointsT[2][l] = -0.5;
    }

    pointsT[0][11] = 0;
    pointsT[1][11] = 0;
    pointsT[2][11] = -sqrt(5) / 2;

    // Faces array
    int facesT[3][20] = {
        {0, 0, 0, 0, 0, 1, 2, 2, 3, 3, 4, 4, 5, 5, 1, 11, 11, 11, 11, 11},
        {1, 2, 3, 4, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 7, 8, 9, 10, 6},
        {2, 3, 4, 5, 1, 2, 7, 3, 8, 4, 9, 5, 10, 1, 6, 6, 7, 8, 9, 10}};

    // Create all the points
    std::vector<Vector3D> points;
    points.reserve(12);
    for (int i = 0; i < 12; ++i)
    {
        points.push_back(
            Vector3D::point(pointsT[0][i], pointsT[1][i], pointsT[2][i]));
    }
    // Create faces
    std::vector<face> faces;
    faces.reserve(20);
    for (int j = 0; j < 20; ++j)
    {
        faces.push_back(face({facesT[0][j], facesT[1][j], facesT[2][j]}));
    }
    // Create new figure
    figure *newIso;
    newIso = new figure(points, faces, ambientCoefficient, diffuseCoefficient,
                        specularCoefficient,
                        reflectionCoefficient);
    return newIso;
}

figure *
createDodecahedron(std::vector<double> &ambientCoefficient,
                   std::vector<double> &diffuseCoefficient,
                   std::vector<double> &specularCoefficient,
                   double &reflectionCoefficient)
{
    // Points array
    double pointsT[3][20];
    // Faces array
    int facesT[5][12] = {
        {0, 0, 1, 2, 3, 4, 19, 19, 18, 17, 16, 15},
        {1, 5, 7, 9, 11, 13, 18, 14, 12, 10, 8, 6},
        {2, 6, 8, 10, 12, 14, 17, 13, 11, 9, 7, 5},
        {3, 7, 9, 11, 13, 5, 16, 12, 10, 8, 6, 14},
        {4, 1, 2, 3, 4, 0, 15, 18, 17, 16, 15, 19}};
    // Create dodecahedron
    figure *icosahedron = createIcosahedron(ambientCoefficient, diffuseCoefficient,
                                            specularCoefficient,
                                            reflectionCoefficient);
    int count = 0;
    for (face f : icosahedron->faces)
    {
        Vector3D p4 = Vector3D::point(0, 0, 0);
        p4 += icosahedron->points.at(f.pointIndexes.at(0));
        p4 += icosahedron->points.at(f.pointIndexes.at(1));
        p4 += icosahedron->points.at(f.pointIndexes.at(2));
        p4 /= 3;
        pointsT[0][count] = p4.x;
        pointsT[1][count] = p4.y;
        pointsT[2][count] = p4.z;
        count++;
    }
    // Create all the points
    std::vector<Vector3D> points;
    points.reserve(20);
    for (int i = 0; i < 20; ++i)
    {
        points.push_back(Vector3D::point(pointsT[0][i], pointsT[1][i], pointsT[2][i]));
    }
    // Create faces
    std::vector<face> faces;
    faces.reserve(12);
    for (int j = 0; j < 12; ++j)
    {
        faces.push_back(face({facesT[0][j], facesT[1][j], facesT[2][j], facesT[3][j],
                              facesT[4][j]}));
    }
    // Create new figure
    figure *dodecahedron;
    dodecahedron = new figure(points, faces, ambientCoefficient, diffuseCoefficient,
                              specularCoefficient,
                              reflectionCoefficient);
    return dodecahedron;
}

void splitTriangle(face &originalTriangle, figure *&originalFigure,
                   std::vector<face> &newFaces)
{
    // Get the points of the triangle
    int indexA = originalTriangle.pointIndexes.at(0);
    int indexB = originalTriangle.pointIndexes.at(1);
    int indexC = originalTriangle.pointIndexes.at(2);
    Vector3D a = originalFigure->points.at(indexA);
    Vector3D b = originalFigure->points.at(indexB);
    Vector3D c = originalFigure->points.at(indexC);
    // Calculate middle points
    Vector3D d = (a + b) / 2;
    Vector3D e = (a + c) / 2;
    Vector3D f = (b + c) / 2;
    // Push points to figure
    originalFigure->points.push_back(d);
    int indexD = static_cast<int>(originalFigure->points.size() - 1);
    originalFigure->points.push_back(e);
    int indexE = static_cast<int>(originalFigure->points.size() - 1);
    originalFigure->points.push_back(f);
    int indexF = static_cast<int>(originalFigure->points.size() - 1);

    // Replace triangles and add the new triangles
    // originalTriangle.pointIndexes = { indexA , indexD , indexE};
    newFaces.push_back(face({indexA, indexD, indexE}));
    newFaces.push_back(face({indexB, indexF, indexD}));
    newFaces.push_back(face({indexC, indexE, indexF}));
    newFaces.push_back(face({indexD, indexF, indexE}));
}

figure *
createSphere(const double &radius, const int &n,
             std::vector<double> &ambientCoefficient,
             std::vector<double> &diffuseCoefficient,
             std::vector<double> &specularCoefficient, double &reflectionCoefficient)
{
    // Create an icosahedron
    figure *newIcoSphere = createIcosahedron(ambientCoefficient, diffuseCoefficient,
                                             specularCoefficient,
                                             reflectionCoefficient);
    // Create temporary faces std::vector
    std::vector<face> newFaces;

    for (int i = 0; i < n; ++i)
    {

        // Split each triangle into 4
        for (face f : newIcoSphere->faces)
        {
            splitTriangle(f, newIcoSphere, newFaces);
        }
        // Update the faces std::vector
        newIcoSphere->faces = newFaces;
        // Reset the temporary std::vector
        newFaces = {};
    }
    // Rescale all the points
    for (Vector3D &p : newIcoSphere->points)
    {
        p.normalise();
    }
    return newIcoSphere;
}

figure *createCone(const int &n, const double &height,
                   std::vector<double> &ambientCoefficient,
                   std::vector<double> &diffuseCoefficient,
                   std::vector<double> &specularCoefficient,
                   double &reflectionCoefficient)
{
    // Create points
    std::vector<Vector3D> points;
    for (int i = 0; i < n + 1; ++i)
    {
        if (i == n)
        {
            points.push_back(Vector3D::point(0, 0, height));
        }
        else
        {
            points.push_back(
                Vector3D::point(cos((2 * i * M_PI) / n), sin((2 * i * M_PI) / n),
                                0));
        }
    }
    // Get indexes std::vector
    std::vector<int> indexes;
    for (int i = n - 1; i >= 0; --i)
    {
        indexes.push_back(i);
    }
    // Create faces
    std::vector<face> faces;
    for (int j = 0; j < n + 1; ++j)
    {
        if (j == n)
        {
            faces.emplace_back(indexes);
        }
        else
        {
            faces.push_back(face({j, (j + 1) % n, n}));
        }
    }
    // Create new figure
    figure *newCone;
    newCone = new figure(points, faces, ambientCoefficient, diffuseCoefficient,
                         specularCoefficient,
                         reflectionCoefficient);
    return newCone;
}

figure *createCylinder(const int &n, const double &height,
                       std::vector<double> &ambientCoefficient,
                       std::vector<double> &diffuseCoefficient,
                       std::vector<double> &specularCoefficient,
                       double &reflectionCoefficient)
{
    // Create points
    std::vector<Vector3D> points;
    std::vector<int> bottomIndexes;
    std::vector<int> topIndexes;
    // Add points of bottom face
    for (int i = 0; i < n; ++i)
    {
        points.push_back(
            Vector3D::point(cos((2 * i * M_PI) / n), sin((2 * i * M_PI) / n), 0));
        bottomIndexes.push_back(i);
    }
    // Add points of top face
    for (int i = n; i < 2 * n; ++i)
    {
        points.push_back(
            Vector3D::point(cos((2 * i * M_PI) / n), sin((2 * i * M_PI) / n),
                            height));
        topIndexes.push_back(i);
    }
    // Create faces
    std::vector<face> faces;
    for (int j = 0; j < n; ++j)
    {
        // Last face loops back to first one
        if (j == n - 1)
        {
            faces.push_back(face({j, 0, n, j + n}));
        }
        // Normal face
        else
        {
            faces.push_back(face({j, j + 1, j + n + 1, j + n}));
        }
    }
    // Add bottom face
    faces.emplace_back(bottomIndexes);
    // Add top face
    faces.emplace_back(topIndexes);
    // Make new figure
    figure *newFigure;
    newFigure = new figure(points, faces, ambientCoefficient, diffuseCoefficient,
                           specularCoefficient,
                           reflectionCoefficient);
    return newFigure;
}

figure *createTorus(const double &r, const double &rBig, const int &n, const int &m,
                    std::vector<double> &ambientCoefficient,
                    std::vector<double> &diffuseCoefficient,
                    std::vector<double> &specularCoefficient,
                    double &reflectionCoefficient)
{
    // Create points std::vector
    std::vector<Vector3D> points;
    // Create all the points
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            // Calculate modifiers u and v
            double u = (2 * i * M_PI) / n;
            double v = (2 * j * M_PI) / m;
            // Calculate the coordinates of the new point
            double xUv = (rBig + r * cos(v)) * cos(u);
            double yUv = (rBig + r * cos(v)) * sin(u);
            double zUv = r * sin(v);

            points.push_back(Vector3D::point(xUv, yUv, zUv));
        }
    }
    // Create faces std::vector
    std::vector<face> faces;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            int ind1 = (i * m + j);
            int ind2 = ((i + 1) % n) * m + j;
            int ind3 = ((i + 1) % n) * m + (j + 1) % m;
            int ind4 = i * m + (j + 1) % m;
            faces.push_back(face({ind1, ind2, ind3, ind4}));
        }
    }
    // Create new figure
    figure *newFigure;
    newFigure = new figure(points, faces, ambientCoefficient, diffuseCoefficient,
                           specularCoefficient,
                           reflectionCoefficient);
    return newFigure;
}