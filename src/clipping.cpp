#include "clipping.h"
#include "geometry/shapes.h"
#include "vector/vector3d.h"
#include <vector>

void clipTriangleOneOutNear(Vector3D &a, Vector3D &b, Vector3D &c, int &indB, int &indC,
                            double &dVal, std::vector<Vector3D> &newPoints,
                            std::vector<face> &newFaces)
{
    // a out, b and c in
    Vector3D d;
    Vector3D e;
    int indD;
    int indE;
    // Get new point in AB line
    double p = (dVal - b.z) / (a.z - b.z);
    e = p * a + (1 - p) * b;
    // Get new point in AC line
    p = (dVal - c.z) / (a.z - c.z);
    d = p * a + (1 - p) * c;
    // Make the two new triangles
    // Check if point already exists
    bool bExists = false;
    bool cExists = false;
    bool dExists = false;
    bool eExists = false;
    for (int i = 0; i < newPoints.size(); i++)
    {
        if (b == newPoints.at(i))
        {
            bExists = true;
            indB = i;
        }
        if (c == newPoints.at(i))
        {
            cExists = true;
            indC = i;
        }
        if (d == newPoints.at(i))
        {
            dExists = true;
            indD = i;
        }
        if (e == newPoints.at(i))
        {
            eExists = true;
            indE = i;
        }
    }
    if (!bExists)
    {
        newPoints.push_back(b);
        indB = static_cast<int>(newPoints.size() - 1);
    }
    if (!cExists)
    {
        newPoints.push_back(c);
        indC = static_cast<int>(newPoints.size() - 1);
    }
    if (!dExists)
    {
        newPoints.push_back(d);
        indD = static_cast<int>(newPoints.size() - 1);
    }
    if (!eExists)
    {
        newPoints.push_back(e);
        indE = static_cast<int>(newPoints.size() - 1);
    }
    newFaces.push_back(face({indC, indD, indE}));
    newFaces.push_back(face({indE, indB, indC}));
}

void clipTriangleTwoOutNear(Vector3D &a, Vector3D &b, Vector3D &c, int &indA, double &dVal,
                            std::vector<Vector3D> &newPoints, std::vector<face> &newFaces)
{
    // a in, b and c out
    // Get new point in AB line
    double p = (dVal - b.z) / (a.z - b.z);
    Vector3D d = p * a + (1 - p) * b;
    // Get new point in AC line
    p = (dVal - c.z) / (a.z - c.z);
    Vector3D e = p * a + (1 - p) * c;
    // Make the new triangle
    // Check if point already exists
    int indD;
    int indE;
    bool aExists = false;
    bool dExists = false;
    bool eExists = false;
    for (int i = 0; i < newPoints.size(); i++)
    {
        if (a == newPoints.at(i))
        {
            aExists = true;
            indA = i;
        }
        if (d == newPoints.at(i))
        {
            dExists = true;
            indD = i;
        }
        if (e == newPoints.at(i))
        {
            eExists = true;
            indE = i;
        }
    }
    if (!aExists)
    {
        newPoints.push_back(a);
        indA = static_cast<int>(newPoints.size() - 1);
    }
    if (!dExists)
    {
        newPoints.push_back(d);
        indD = static_cast<int>(newPoints.size() - 1);
    }
    if (!eExists)
    {
        newPoints.push_back(e);
        indE = static_cast<int>(newPoints.size() - 1);
    }
    newFaces.push_back(face({indA, indD, indE}));
}

void clipTriangleOneOutFar(Vector3D &a, Vector3D &b, Vector3D &c, int &indB, int &indC,
                           double &dVal, std::vector<Vector3D> &newPoints,
                           std::vector<face> &newFaces)
{
    // a out, b and c in
    Vector3D d;
    Vector3D e;
    // Get new point in AB line
    double p = (dVal - b.z) / (a.z - b.z);
    d = p * a + (1 - p) * b;
    // Get new point in AC line
    p = (dVal - c.z) / (a.z - c.z);
    e = p * a + (1 - p) * c;
    // Make the two new triangles
    // Check if point already exists
    int indD;
    int indE;
    bool bExists = false;
    bool cExists = false;
    bool dExists = false;
    bool eExists = false;
    for (int i = 0; i < newPoints.size(); i++)
    {
        if (b == newPoints.at(i))
        {
            bExists = true;
        }
        if (c == newPoints.at(i))
        {
            cExists = true;
            indC = i;
        }
        if (d == newPoints.at(i))
        {
            dExists = true;
            indD = i;
        }
        if (e == newPoints.at(i))
        {
            eExists = true;
            indE = i;
        }
    }
    if (!bExists)
    {
        newPoints.push_back(b);
    }
    if (!cExists)
    {
        newPoints.push_back(c);
        indC = static_cast<int>(newPoints.size() - 1);
    }
    if (!dExists)
    {
        newPoints.push_back(d);
        indD = static_cast<int>(newPoints.size() - 1);
    }
    if (!eExists)
    {
        newPoints.push_back(e);
        indE = static_cast<int>(newPoints.size() - 1);
    }
    newFaces.push_back(face({indB, indE, indD}));
    newFaces.push_back(face({indB, indC, indE}));
}

void clipTriangleTwoOutFar(Vector3D &a, Vector3D &b, Vector3D &c, int &indA, double &dVal,
                           std::vector<Vector3D> &newPoints, std::vector<face> &newFaces)
{
    // a in, b and c out
    // Get new point in AB line
    double p = (dVal - b.z) / (a.z - b.z);
    Vector3D e = p * a + (1 - p) * b;
    // Get new point in AC line
    p = (dVal - c.z) / (a.z - c.z);
    Vector3D d = p * a + (1 - p) * c;
    // Make the new triangle
    // Check if point already exists
    int indD;
    int indE;
    bool aExists = false;
    bool dExists = false;
    bool eExists = false;
    for (int i = 0; i < newPoints.size(); i++)
    {
        if (a == newPoints.at(i))
        {
            aExists = true;
            indA = i;
        }
        if (d == newPoints.at(i))
        {
            dExists = true;
            indD = i;
        }
        if (e == newPoints.at(i))
        {
            eExists = true;
            indE = i;
        }
    }
    if (!aExists)
    {
        newPoints.push_back(a);
        indA = static_cast<int>(newPoints.size() - 1);
    }
    if (!dExists)
    {
        newPoints.push_back(d);
        indD = static_cast<int>(newPoints.size() - 1);
    }
    if (!eExists)
    {
        newPoints.push_back(e);
        indE = static_cast<int>(newPoints.size() - 1);
    }
    newFaces.push_back(face({indA, indE, indD}));
}

void clipNear(figure *&originalFigure, std::vector<face> &newFaces, double &dNear)
{

    double dVal = -dNear;
    Vector3D a;
    int indA;
    Vector3D b;
    int indB;
    Vector3D c;
    int indC;

    std::vector<Vector3D> newPoints = {};

    for (auto &figureFace : originalFigure->faces)
    {

        // Get points of the face
        indA = figureFace.pointIndexes.at(0);
        indB = figureFace.pointIndexes.at(1);
        indC = figureFace.pointIndexes.at(2);
        a = originalFigure->points.at(indA);
        b = originalFigure->points.at(indB);
        c = originalFigure->points.at(indC);
        // If all points within range
        bool aExists = false;
        bool bExists = false;
        bool cExists = false;
        if (a.z <= dVal && b.z <= dVal && c.z <= dVal)
        {
            // Check if point already exists
            for (int i = 0; i < newPoints.size(); i++)
            {
                if (a == newPoints.at(i))
                {
                    aExists = true;
                    indA = i;
                }
                if (b == newPoints.at(i))
                {
                    bExists = true;
                    indB = i;
                }
                if (c == newPoints.at(i))
                {
                    cExists = true;
                    indC = i;
                }
            }
            if (!aExists)
            {
                newPoints.push_back(a);
                indA = static_cast<int>(newPoints.size() - 1);
            }
            if (!bExists)
            {
                newPoints.push_back(b);
                indB = static_cast<int>(newPoints.size() - 1);
            }
            if (!cExists)
            {
                newPoints.push_back(c);
                indC = static_cast<int>(newPoints.size() - 1);
            }
            newFaces.push_back(face({indA, indB, indC}));
            continue;
        }
        // All points are out of range
        if (a.z >= dVal && b.z >= dVal && c.z >= dVal)
        {
            continue;
        }
        // If 1 point is out of range
        // If a is out of range
        if (a.z > dVal && b.z <= dVal && c.z <= dVal)
        {
            clipTriangleOneOutNear(a, b, c, indB, indC, dVal, newPoints,
                                   newFaces);
        }
        // If b is out of range
        if (b.z >= dVal && a.z <= dVal && c.z < dVal)
        {
            clipTriangleOneOutNear(b, c, a, indC, indA, dVal, newPoints,
                                   newFaces);
        }
        // If c is out of range
        if (c.z >= dVal && a.z <= dVal && b.z <= dVal)
        {
            clipTriangleOneOutNear(c, a, b, indA, indB, dVal, newPoints,
                                   newFaces);
        }
        // If 2 points are out of range
        // If b and c are out of range
        if (b.z >= dVal && c.z >= dVal && a.z <= dVal)
        {
            clipTriangleTwoOutNear(a, b, c, indA, dVal, newPoints,
                                   newFaces);
        }
        // If c and a are out of range
        if (a.z >= dVal && c.z >= dVal && b.z <= dVal)
        {
            clipTriangleTwoOutNear(b, c, a, indB, dVal, newPoints,
                                   newFaces);
        }
        // If a and b are out of range
        if (a.z >= dVal && b.z >= dVal && c.z <= dVal)
        {
            clipTriangleTwoOutNear(c, a, b, indC, dVal, newPoints,
                                   newFaces);
        }
    }
    // Set the new points
    originalFigure->points = newPoints;
    // Set the new faces
    originalFigure->faces = newFaces;
    newFaces = {};
}

void clipFar(figure *&originalFigure, std::vector<face> &newFaces, double &dFar)
{
    double dVal = -dFar;
    Vector3D a;
    int indA;
    Vector3D b;
    int indB;
    Vector3D c;
    int indC;

    std::vector<Vector3D> newPoints = {};

    for (auto &figureFace : originalFigure->faces)
    {

        // Get points of the face
        indA = figureFace.pointIndexes.at(0);
        indB = figureFace.pointIndexes.at(1);
        indC = figureFace.pointIndexes.at(2);
        a = originalFigure->points.at(indA);
        b = originalFigure->points.at(indB);
        c = originalFigure->points.at(indC);
        // If all points within range
        bool aExists = false;
        bool bExists = false;
        bool cExists = false;
        if (a.z > dVal && b.z > dVal && c.z > dVal)
        {
            // Check if point already exists
            for (int i = 0; i < newPoints.size(); i++)
            {
                if (a == newPoints.at(i))
                {
                    aExists = true;
                    indA = i;
                }
                if (b == newPoints.at(i))
                {
                    bExists = true;
                    indB = i;
                }
                if (c == newPoints.at(i))
                {
                    cExists = true;
                    indC = i;
                }
            }
            if (!aExists)
            {
                newPoints.push_back(a);
                indA = static_cast<int>(newPoints.size() - 1);
            }
            if (!bExists)
            {
                newPoints.push_back(b);
                indB = static_cast<int>(newPoints.size() - 1);
            }
            if (!cExists)
            {
                newPoints.push_back(c);
                indC = static_cast<int>(newPoints.size() - 1);
            }
            newFaces.push_back(face({indA, indB, indC}));
            continue;
        }
        // If all out of range
        if (a.z <= dVal && b.z <= dVal && c.z <= dVal)
        {
            continue;
        }
        // If 1 point is out of range
        // If a is out of range
        if (a.z <= dVal && b.z >= dVal && c.z >= dVal)
        {
            clipTriangleOneOutFar(a, b, c, indB, indC, dVal, newPoints, newFaces);
        }
        // If b is out of range
        if (b.z <= dVal && a.z >= dVal && c.z >= dVal)
        {
            clipTriangleOneOutFar(b, c, a, indC, indA, dVal, newPoints, newFaces);
        }
        // If c is out of range
        if (c.z <= dVal && a.z >= dVal && b.z >= dVal)
        {
            clipTriangleOneOutFar(c, a, b, indA, indB, dVal, newPoints, newFaces);
        }
        // If 2 points are out of range
        // If b and c are out of range
        if (b.z <= dVal && c.z <= dVal && a.z >= dVal)
        {
            clipTriangleTwoOutFar(a, b, c, indA, dVal, newPoints, newFaces);
        }
        // If c and a are out of range
        if (c.z <= dVal && a.z <= dVal && b.z >= dVal)
        {
            clipTriangleTwoOutFar(b, c, a, indB, dVal, newPoints, newFaces);
        }
        // If a and b are out of range
        if (a.z <= dVal && b.z <= dVal && c.z >= dVal)
        {
            clipTriangleTwoOutFar(c, a, b, indC, dVal, newPoints, newFaces);
        }
    }
    // Set the new points
    originalFigure->points = newPoints;
    // Set the new faces
    originalFigure->faces = newFaces;
    newFaces = {};
}

void clipTriangleOneOutLeft(Vector3D &a, Vector3D &b, Vector3D &c, int &indB, int &indC,
                            double &dVal, double &dNear,
                            std::vector<Vector3D> &newPoints,
                            std::vector<face> &newFaces)
{
    // a out , b and c in
    Vector3D d;
    Vector3D e;
    // Get new point in AB line
    double p =
        (b.x * dNear + b.z * dVal) / ((b.x - a.x) * dNear + (b.z - a.z) * dVal);
    if (dVal <= 0)
    {
        e = p * a + (1 - p) * b;
    }
    else
    {
        d = p * a + (1 - p) * b;
    }
    // Get new point in AC line
    p = (c.x * dNear + c.z * dVal) / ((c.x - a.x) * dNear + (c.z - a.z) * dVal);
    if (dVal <= 0)
    {
        d = p * a + (1 - p) * c;
    }
    else
    {
        e = p * a + (1 - p) * c;
    }
    // Make the two new triangles
    bool bExists = false;
    bool cExists = false;
    // Check if point already exists
    int indD;
    int indE;
    bool dExists = false;
    bool eExists = false;
    for (int i = 0; i < newPoints.size(); i++)
    {
        if (b == newPoints.at(i))
        {
            bExists = true;
            indB = i;
        }
        if (c == newPoints.at(i))
        {
            cExists = true;
            indC = i;
        }
        if (d == newPoints.at(i))
        {
            dExists = true;
            indD = i;
        }
        if (e == newPoints.at(i))
        {
            eExists = true;
            indE = i;
        }
    }
    if (!bExists)
    {
        newPoints.push_back(b);
        indB = static_cast<int>(newPoints.size() - 1);
    }
    if (!cExists)
    {
        newPoints.push_back(c);
        indC = static_cast<int>(newPoints.size() - 1);
    }
    if (!dExists)
    {
        newPoints.push_back(d);
        indD = static_cast<int>(newPoints.size() - 1);
    }
    if (!eExists)
    {
        newPoints.push_back(e);
        indE = static_cast<int>(newPoints.size() - 1);
    }

    newFaces.push_back(face({indD, indB, indC}));
    if (dVal <= 0)
    {
        newFaces.push_back(face({indD, indE, indB}));
    }
    else
    {
        newFaces.push_back(face({indD, indC, indE}));
    }
}

void clipTriangleTwoOutLeft(Vector3D &a, Vector3D &b, Vector3D &c, int &indA, double &dVal,
                            double &dNear, std::vector<Vector3D> &newPoints,
                            std::vector<face> &newFaces)
{
    // a in, b and c out
    Vector3D d;
    Vector3D e;
    // Get new point in AB line
    double p =
        (b.x * dNear + b.z * dVal) / ((b.x - a.x) * dNear + (b.z - a.z) * dVal);
    if (dVal <= 0)
    {
        d = p * a + (1 - p) * b;
    }
    else
    {
        e = p * a + (1 - p) * b;
    }
    // Get new point in AC line
    p = (c.x * dNear + c.z * dVal) / ((c.x - a.x) * dNear + (c.z - a.z) * dVal);
    if (dVal <= 0)
    {
        e = p * a + (1 - p) * c;
    }
    else
    {
        d = p * a + (1 - p) * c;
    }
    // Make the new triangle
    // Check if point already exists
    int indD;
    int indE;
    bool aExists = false;
    bool dExists = false;
    bool eExists = false;
    for (int i = 0; i < newPoints.size(); i++)
    {
        if (a == newPoints.at(i))
        {
            aExists = true;
            indA = i;
        }
        if (d == newPoints.at(i))
        {
            dExists = true;
            indD = i;
        }
        if (e == newPoints.at(i))
        {
            eExists = true;
            indE = i;
        }
    }
    if (!aExists)
    {
        newPoints.push_back(a);
        indA = static_cast<int>(newPoints.size() - 1);
    }
    if (!dExists)
    {
        newPoints.push_back(d);
        indD = static_cast<int>(newPoints.size() - 1);
    }
    if (!eExists)
    {
        newPoints.push_back(e);
        indE = static_cast<int>(newPoints.size() - 1);
    }
    if (dVal <= 0)
    {
        newFaces.push_back(face({indA, indD, indE}));
    }
    else
    {
        newFaces.push_back(face({indA, indE, indD}));
    }
}

void clipLeft(figure *&originalFigure, std::vector<face> &newFaces, double &left,
              double &dNear)
{

    double dVal;
    double xA;
    double xB;
    double xC;

    Vector3D a;
    int indA;
    Vector3D b;
    int indB;
    Vector3D c;
    int indC;

    std::vector<Vector3D> newPoints = {};

    for (auto &figureFace : originalFigure->faces)
    {

        // Get points of the face
        indA = figureFace.pointIndexes.at(0);
        indB = figureFace.pointIndexes.at(1);
        indC = figureFace.pointIndexes.at(2);
        a = originalFigure->points.at(indA);
        b = originalFigure->points.at(indB);
        c = originalFigure->points.at(indC);

        xA = -a.x * dNear / a.z;
        xB = -b.x * dNear / b.z;
        xC = -c.x * dNear / c.z;
        // If all points within range
        bool aExists = false;
        bool bExists = false;
        bool cExists = false;
        if (xA >= left && xB >= left && xC >= left)
        {
            // Check if point already exists
            for (int i = 0; i < newPoints.size(); i++)
            {
                if (a == newPoints.at(i))
                {
                    aExists = true;
                    indA = i;
                }
                if (b == newPoints.at(i))
                {
                    bExists = true;
                    indB = i;
                }
                if (c == newPoints.at(i))
                {
                    cExists = true;
                    indC = i;
                }
            }
            if (!aExists)
            {
                newPoints.push_back(a);
                indA = static_cast<int>(newPoints.size() - 1);
            }
            if (!bExists)
            {
                newPoints.push_back(b);
                indB = static_cast<int>(newPoints.size() - 1);
            }
            if (!cExists)
            {
                newPoints.push_back(c);
                indC = static_cast<int>(newPoints.size() - 1);
            }
            newFaces.push_back(face({indA, indB, indC}));
            continue;
        }
        dVal = left;
        if (xA < dVal && xB < dVal && xC < dVal)
        {
            continue;
        }
        // Clip t.o.v left
        dVal = left;
        // If 1 point is out of range
        // If a is out of range
        if (xA < dVal && xB >= dVal && xC >= dVal)
        {
            clipTriangleOneOutLeft(a, b, c, indB, indC, dVal, dNear, newPoints,
                                   newFaces);
        }
        // If b is out of range
        if (xB < dVal && xA >= dVal && xC >= dVal)
        {
            clipTriangleOneOutLeft(b, c, a, indC, indA, dVal, dNear, newPoints,
                                   newFaces);
        }
        // If c is out of range
        if (xC < dVal && xA >= dVal && xB >= dVal)
        {
            clipTriangleOneOutLeft(c, a, b, indA, indB, dVal, dNear, newPoints,
                                   newFaces);
        }
        // If 2 points are out of range
        // If b and c are out of range
        if (xB < dVal && xC < dVal && xA >= dVal)
        {
            clipTriangleTwoOutLeft(a, b, c, indA, dVal, dNear, newPoints,
                                   newFaces);
        }
        // If c and a are out of range
        if (xC < dVal && xA < dVal && xB >= dVal)
        {
            clipTriangleTwoOutLeft(b, c, a, indB, dVal, dNear, newPoints,
                                   newFaces);
        }
        // If a and b are out of range
        if (xA < dVal && xB < dVal && xC >= dVal)
        {
            clipTriangleTwoOutLeft(c, a, b, indC, dVal, dNear, newPoints,
                                   newFaces);
        }
    }
    // Set the new points
    originalFigure->points = newPoints;
    // Set the new faces
    originalFigure->faces = newFaces;
    newFaces = {};
}

void clipRight(figure *&originalFigure, std::vector<face> &newFaces, double &right,
               double &dNear)
{

    double dVal;
    double xA;
    double xB;
    double xC;

    Vector3D a;
    int indA;
    Vector3D b;
    int indB;
    Vector3D c;
    int indC;

    std::vector<Vector3D> newPoints = {};

    for (auto &figureFace : originalFigure->faces)
    {

        // Get points of the face
        indA = figureFace.pointIndexes.at(0);
        indB = figureFace.pointIndexes.at(1);
        indC = figureFace.pointIndexes.at(2);
        a = originalFigure->points.at(indA);
        b = originalFigure->points.at(indB);
        c = originalFigure->points.at(indC);

        xA = -a.x * dNear / a.z;
        xB = -b.x * dNear / b.z;
        xC = -c.x * dNear / c.z;
        bool aExists = false;
        bool bExists = false;
        bool cExists = false;
        // If all points within range
        if (xA <= right && xB <= right && xC <= right)
        {
            // Check if point already exists
            for (int i = 0; i < newPoints.size(); i++)
            {
                if (a == newPoints.at(i))
                {
                    aExists = true;
                    indA = i;
                }
                if (b == newPoints.at(i))
                {
                    bExists = true;
                    indB = i;
                }
                if (c == newPoints.at(i))
                {
                    cExists = true;
                    indC = i;
                }
            }
            if (!aExists)
            {
                newPoints.push_back(a);
                indA = static_cast<int>(newPoints.size() - 1);
            }
            if (!bExists)
            {
                newPoints.push_back(b);
                indB = static_cast<int>(newPoints.size() - 1);
            }
            if (!cExists)
            {
                newPoints.push_back(c);
                indC = static_cast<int>(newPoints.size() - 1);
            }
            newFaces.push_back(face({indA, indB, indC}));
            continue;
        }

        // Clip t.o.v right
        dVal = right;
        // If all points are out of range
        if (xA > dVal && xB > dVal && xC > dVal)
        {
            continue;
        }
        // If 1 point is out of range
        // If a is out of range
        if (xA > dVal && xB <= dVal && xC <= dVal)
        {
            clipTriangleOneOutLeft(a, b, c, indB, indC, dVal, dNear, newPoints,
                                   newFaces);
        }
        // If b is out of range
        if (xB > dVal && xA <= dVal && xC <= dVal)
        {
            clipTriangleOneOutLeft(b, c, a, indC, indA, dVal, dNear, newPoints,
                                   newFaces);
        }
        // If c is out of range
        if (xC > dVal && xA <= dVal && xB <= dVal)
        {
            clipTriangleOneOutLeft(c, a, b, indA, indB, dVal, dNear, newPoints,
                                   newFaces);
        }
        // If 2 points are out of range
        // If b and c are out of range
        if (xB > dVal && xC > dVal && xA <= dVal)
        {
            clipTriangleTwoOutLeft(a, b, c, indA, dVal, dNear, newPoints,
                                   newFaces);
        }
        // If c and a are out of range
        if (xC > dVal && xA > dVal && xB <= dVal)
        {
            clipTriangleTwoOutLeft(b, c, a, indB, dVal, dNear, newPoints,
                                   newFaces);
        }
        // If a and b are out of range
        if (xA > dVal && xB > dVal && xC <= dVal)
        {
            clipTriangleTwoOutLeft(c, a, b, indC, dVal, dNear, newPoints,
                                   newFaces);
        }
    }
    // Set the new points
    originalFigure->points = newPoints;
    // Set the new faces
    originalFigure->faces = newFaces;
    newFaces = {};
}

void clipTriangleOneOutTop(Vector3D &a, Vector3D &b, Vector3D &c, int &indB, int &indC,
                           double &dVal, double &dNear,
                           std::vector<Vector3D> &newPoints,
                           std::vector<face> &newFaces)
{
    // a out, b and c in
    Vector3D d;
    Vector3D e;
    // Get new point in AB line
    double p =
        (b.y * dNear + b.z * dVal) / ((b.y - a.y) * dNear + (b.z - a.z) * dVal);
    if (dVal >= 0)
    {
        e = p * a + (1 - p) * b;
    }
    else
    {
        d = p * a + (1 - p) * b;
    }
    // Get new point in AC line
    p = (c.y * dNear + c.z * dVal) / ((c.y - a.y) * dNear + (c.z - a.z) * dVal);
    if (dVal >= 0)
    {
        d = p * a + (1 - p) * c;
    }
    else
    {
        e = p * a + (1 - p) * c;
    }
    // Make the two new triangles
    bool bExists = false;
    bool cExists = false;
    // Check if point already exists
    int indD;
    int indE;
    bool dExists = false;
    bool eExists = false;
    for (int i = 0; i < newPoints.size(); i++)
    {
        if (b == newPoints.at(i))
        {
            bExists = true;
            indB = i;
        }
        if (c == newPoints.at(i))
        {
            cExists = true;
            indC = i;
        }
        if (d == newPoints.at(i))
        {
            dExists = true;
            indD = i;
        }
        if (e == newPoints.at(i))
        {
            eExists = true;
            indE = i;
        }
    }
    if (!bExists)
    {
        newPoints.push_back(b);
        indB = static_cast<int>(newPoints.size() - 1);
    }
    if (!cExists)
    {
        newPoints.push_back(c);
        indC = static_cast<int>(newPoints.size() - 1);
    }
    if (!dExists)
    {
        newPoints.push_back(d);
        indD = static_cast<int>(newPoints.size() - 1);
    }
    if (!eExists)
    {
        newPoints.push_back(e);
        indE = static_cast<int>(newPoints.size() - 1);
    }

    newFaces.push_back(face({indD, indB, indC}));
    if (dVal >= 0)
    {
        newFaces.push_back(face({indD, indE, indB}));
    }
    else
    {
        newFaces.push_back(face({indD, indC, indE}));
    }
}

void clipTriangleTwoOutTop(Vector3D &a, Vector3D &b, Vector3D &c, int &indA, double &dVal,
                           double &dNear, std::vector<Vector3D> &newPoints,
                           std::vector<face> &newFaces)
{
    Vector3D d;
    Vector3D e;
    // Get new point in AB line
    double p =
        (b.y * dNear + b.z * dVal) / ((b.y - a.y) * dNear + (b.z - a.z) * dVal);
    if (dVal >= 0)
    {
        d = p * a + (1 - p) * b;
    }
    else
    {
        e = p * a + (1 - p) * b;
    }
    // Get new point in AC line
    p = (c.y * dNear + c.z * dVal) / ((c.y - a.y) * dNear + (c.z - a.z) * dVal);
    if (dVal >= 0)
    {
        e = p * a + (1 - p) * c;
    }
    else
    {
        d = p * a + (1 - p) * c;
    }
    // Make the new triangle
    // Check if point already exists
    int indD;
    int indE;
    bool aExists = false;
    bool dExists = false;
    bool eExists = false;
    for (int i = 0; i < newPoints.size(); i++)
    {
        if (a == newPoints.at(i))
        {
            aExists = true;
            indA = i;
        }
        if (d == newPoints.at(i))
        {
            dExists = true;
            indD = i;
        }
        if (e == newPoints.at(i))
        {
            eExists = true;
            indE = i;
        }
    }
    if (!aExists)
    {
        newPoints.push_back(a);
        indA = static_cast<int>(newPoints.size() - 1);
    }
    if (!dExists)
    {
        newPoints.push_back(d);
        indD = static_cast<int>(newPoints.size() - 1);
    }
    if (!eExists)
    {
        newPoints.push_back(e);
        indE = static_cast<int>(newPoints.size() - 1);
    }
    if (dVal >= 0)
    {
        newFaces.push_back(face({indA, indD, indE}));
    }
    else
    {
        newFaces.push_back(face({indA, indE, indD}));
    }
}

void clipTop(figure *&originalFigure, std::vector<face> &newFaces, double &top, double &dNear)
{
    double dVal;
    double yA;
    double yB;
    double yC;

    Vector3D a;
    int indA;
    Vector3D b;
    int indB;
    Vector3D c;
    int indC;

    std::vector<Vector3D> newPoints = {};

    for (auto &figureFace : originalFigure->faces)
    {

        // Get points of the face
        indA = figureFace.pointIndexes.at(0);
        indB = figureFace.pointIndexes.at(1);
        indC = figureFace.pointIndexes.at(2);

        a = originalFigure->points.at(indA);
        b = originalFigure->points.at(indB);
        c = originalFigure->points.at(indC);

        yA = -a.y * dNear / a.z;
        yB = -b.y * dNear / b.z;
        yC = -c.y * dNear / c.z;

        // If all points within range
        bool aExists = false;
        bool bExists = false;
        bool cExists = false;

        if (yA <= top && yB <= top && yC <= top)
        {
            // Check if point already exists
            for (int i = 0; i < newPoints.size(); i++)
            {
                if (a == newPoints.at(i))
                {
                    aExists = true;
                    indA = i;
                }
                if (b == newPoints.at(i))
                {
                    bExists = true;
                    indB = i;
                }
                if (c == newPoints.at(i))
                {
                    cExists = true;
                    indC = i;
                }
            }
            if (!aExists)
            {
                newPoints.push_back(a);
                indA = static_cast<int>(newPoints.size() - 1);
            }
            if (!bExists)
            {
                newPoints.push_back(b);
                indB = static_cast<int>(newPoints.size() - 1);
            }
            if (!cExists)
            {
                newPoints.push_back(c);
                indC = static_cast<int>(newPoints.size() - 1);
            }
            newFaces.push_back(face({indA, indB, indC}));
            continue;
        }
        // All points are out of range
        dVal = top;
        if (yA > dVal && yB > dVal && yC > dVal)
        {
            continue;
        }
        // Clip t.o.v top
        dVal = top;
        // If 1 point is out of range
        // If a is out of range
        if (yA > dVal && yB <= dVal && yC <= dVal)
        {
            clipTriangleOneOutTop(a, b, c, indB, indC, dVal, dNear, newPoints,
                                  newFaces);
        }
        // If b is out of range
        if (yB > dVal && yA <= dVal && yC <= dVal)
        {
            clipTriangleOneOutTop(b, c, a, indC, indA, dVal, dNear, newPoints,
                                  newFaces);
        }
        // If c is out of range
        if (yC > dVal && yA <= dVal && yB <= dVal)
        {
            clipTriangleOneOutTop(c, a, b, indA, indB, dVal, dNear, newPoints,
                                  newFaces);
        }
        // If 2 points are out of range
        // If b and c are out of range
        if (yB > dVal && yC > dVal && yA <= dVal)
        {
            clipTriangleTwoOutTop(a, b, c, indA, dVal, dNear, newPoints,
                                  newFaces);
        }
        // If c and a are out of range
        if (yC > dVal && yA > dVal && yB <= dVal)
        {
            clipTriangleTwoOutTop(b, c, a, indB, dVal, dNear, newPoints,
                                  newFaces);
        }
        // If a and b are out of range
        if (yA > dVal && yB > dVal && yC <= dVal)
        {
            clipTriangleTwoOutTop(c, a, b, indC, dVal, dNear, newPoints,
                                  newFaces);
        }
    }
    // Set the new points
    originalFigure->points = newPoints;
    // Set the new faces
    originalFigure->faces = newFaces;
    newFaces = {};
}

void clipBottom(figure *&originalFigure, std::vector<face> &newFaces, double &bottom,
                double &dNear)
{
    double d;
    double yA;
    double yB;
    double yC;

    Vector3D a;
    int indA;
    Vector3D b;
    int indB;
    Vector3D c;
    int indC;

    std::vector<Vector3D> newPoints = {};

    for (auto &figureFace : originalFigure->faces)
    {

        // Get points of the face
        indA = figureFace.pointIndexes.at(0);
        indB = figureFace.pointIndexes.at(1);
        indC = figureFace.pointIndexes.at(2);

        a = originalFigure->points.at(indA);
        b = originalFigure->points.at(indB);
        c = originalFigure->points.at(indC);

        yA = -a.y * dNear / a.z;
        yB = -b.y * dNear / b.z;
        yC = -c.y * dNear / c.z;

        // If all points within range
        bool aExists = false;
        bool bExists = false;
        bool cExists = false;
        if (yA > bottom && yB > bottom && yC > bottom)
        {
            // Check if point already exists
            for (int i = 0; i < newPoints.size(); i++)
            {
                if (a == newPoints.at(i))
                {
                    aExists = true;
                    indA = i;
                }
                if (b == newPoints.at(i))
                {
                    bExists = true;
                    indB = i;
                }
                if (c == newPoints.at(i))
                {
                    cExists = true;
                    indC = i;
                }
            }
            if (!aExists)
            {
                newPoints.push_back(a);
                indA = static_cast<int>(newPoints.size() - 1);
            }
            if (!bExists)
            {
                newPoints.push_back(b);
                indB = static_cast<int>(newPoints.size() - 1);
            }
            if (!cExists)
            {
                newPoints.push_back(c);
                indC = static_cast<int>(newPoints.size() - 1);
            }
            newFaces.push_back(face({indA, indB, indC}));
            continue;
        }
        // All points are out of range
        d = bottom;
        if (yA <= d && yB <= d && yC <= d)
        {
            continue;
        }
        // Clip t.o.v bottom
        d = bottom;
        // If 1 point is out of range
        // If a is out of range
        if (yA <= d && yB > d && yC > d)
        {
            clipTriangleOneOutTop(a, b, c, indB, indC, d, dNear, newPoints,
                                  newFaces);
        }
        // If b is out of range
        if (yB <= d && yA > d && yC > d)
        {
            clipTriangleOneOutTop(b, c, a, indC, indA, d, dNear, newPoints,
                                  newFaces);
        }
        // If c is out of range
        if (yC <= d && yA > d && yB > d)
        {
            clipTriangleOneOutTop(c, a, b, indA, indB, d, dNear, newPoints,
                                  newFaces);
        }
        // If 2 points are out of range
        // If b and c are out of range
        if (yB <= d && yC <= d && yA > d)
        {
            clipTriangleTwoOutTop(a, b, c, indA, d, dNear, newPoints,
                                  newFaces);
        }
        // If c and a are out of range
        if (yC <= d && yA <= d && yB > d)
        {
            clipTriangleTwoOutTop(b, c, a, indB, d, dNear, newPoints,
                                  newFaces);
        }
        // If a and b are out of range
        if (yA <= d && yB <= d && yC > d)
        {
            clipTriangleTwoOutTop(c, a, b, indC, d, dNear, newPoints,
                                  newFaces);
        }
    }
    // Set the new points
    originalFigure->points = newPoints;
    // Set the new faces
    originalFigure->faces = newFaces;
    newFaces = {};
}

void triangulate(face &originalFace, std::vector<face> &newFaces)
{
    for (unsigned int i = 1; i < originalFace.pointIndexes.size() - 1; ++i)
    {
        // Create new face
        newFaces.push_back(
            face({originalFace.pointIndexes.at(0), originalFace.pointIndexes.at(i),
                  originalFace.pointIndexes.at(i + 1)}));
    }
}

void clipView(figures3D &originalFigures, double &dNear, double &dFar, double &hFov,
              double &aspectRatio)
{
    std::vector<face> newFace = {};
    double right;
    double left;
    double top;
    double bottom;
    // Get right
    hFov *= M_PI / 180;
    right = dNear * tan(hFov / 2);
    left = -right;
    // Get top
    top = right / aspectRatio;
    bottom = -top;
    // Check each figure
    for (auto &f : originalFigures)
    {
        // Clip near/far
        for (auto face : f->faces)
        {
            triangulate(face, newFace);
        }
        f->faces = newFace;
        newFace = {};
        clipNear(f, newFace, dNear);
        // Clip far
        clipFar(f, newFace, dFar);
        // Clip left
        clipLeft(f, newFace, left, dNear);
        // Clip right
        clipRight(f, newFace, right, dNear);
        // Clip top
        clipTop(f, newFace, top, dNear);
        // Clip bottom
        clipBottom(f, newFace, bottom, dNear);
        newFace = {};
    }
}
