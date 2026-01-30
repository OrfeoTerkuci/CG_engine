#include "zBuffering.h"

zBuffer::zBuffer(const int width, const int height)
{
    for (int i = 0; i < width; i++)
    {
        this->push_back(vector<double>(height, posInf));
    }
}

void getProjectedPoints(figures3D &figs, double &d, double &xRange, double &yRange,
                        double &xMin, double &yMin, double &xMax, double &yMax)
{
    for (auto f : figs)
    {
        for (const auto &p : f->points)
        {
            point2D newPoint = doProjection(p, d);
            // Update min and max
            if (xMax <= newPoint.x)
            {
                xMax = newPoint.x;
            }
            if (xMin >= newPoint.x)
            {
                xMin = newPoint.x;
            }
            if (yMax <= newPoint.y)
            {
                yMax = newPoint.y;
            }
            if (yMin >= newPoint.y)
            {
                yMin = newPoint.y;
            }
        }
    }
    // Calculate range along the x-axis and y-axis
    xRange = xMax - xMin;
    yRange = yMax - yMin;
}

double
calculateInvZ(unsigned int i, unsigned int iMin, unsigned int iMax, const double &z0,
              const double &z1)
{
    // Calculate current p
    double p = (double)i / (double)(iMax - iMin);
    return p / z0 + (1 - p) / z1;
}

void drawZBufferLine(zBuffer &zBuffer, img::EasyImage &image,
                     unsigned int x0, unsigned int y0, double &z0,
                     unsigned int x1, unsigned int y1, double &z1,
                     const color &lineColor)
{

    assert(x0 < image.get_width() && y0 < image.get_height());
    assert(x1 < image.get_width() && y1 < image.get_height());

    img::Color newColor = img::Color((int)lineColor.red, (int)lineColor.green,
                                     (int)lineColor.blue);
    // Vertical line
    if (x0 == x1)
    {
        // special case for x0 == x1
        for (unsigned int i = std::min(y0, y1); i <= std::max(y0, y1); i++)
        {
            double invZ;

            invZ = y0 < y1 ? calculateInvZ(i - y0, y0, y1, z0, z1) : calculateInvZ(i - y1, y1, y0, z1, z0);
            // Check if we can draw
            if (invZ < zBuffer[x0][i])
            {
                // Update z-buffer
                zBuffer[x0][i] = invZ;
                (image)(x0, i) = newColor;
            }
        }
    }
    // Horizontal line
    else if (y0 == y1)
    {
        // special case for y0 == y1
        for (unsigned int i = std::min(x0, x1); i <= std::max(x0, x1); i++)
        {
            double invZ;

            invZ = x0 < x1 ? calculateInvZ(i - x0, x0, x1, z0, z1) : calculateInvZ(i - x1, x1, x0, z1, z0);
            // Check if we can draw
            if (invZ < zBuffer[i][y0])
            {
                // Update z-buffer
                zBuffer[i][y0] = invZ;
                (image)(i, y0) = newColor;
            }
        }
    }
    else
    {
        if (x0 > x1)
        {
            // flip points if x1>x0: we want x0 to have the lowest value
            std::swap(x0, x1);
            std::swap(y0, y1);
            std::swap(z0, z1); // right corners fixed with swap enabled
        }
        // Calculate the coefficient
        double m = ((double)y1 - (double)y0) / ((double)x1 - (double)x0);

        if (-1.0 <= m && m <= 1.0)
        {
            for (unsigned int i = 0; i <= (x1 - x0); i++)
            {
                double invZ;

                invZ = calculateInvZ(i, 0, (x1 - x0), z1, z0);
                // Check if we can draw
                if (invZ < zBuffer[x0 + i][(unsigned int)round(y0 + m * i)])
                {
                    // Update z-buffer
                    zBuffer[x0 + i][(unsigned int)round(y0 + m * i)] = invZ;
                    (image)(x0 + i, (unsigned int)round(y0 + m * i)) = newColor;
                }
            }
        }
        else if (m > 1.0)
        {
            for (unsigned int i = 0; i <= (y1 - y0); i++)
            {
                double invZ;

                invZ = calculateInvZ(i, 0, (y1 - y0), z1, z0);
                if (invZ < zBuffer[(unsigned int)round(x0 + (i / m))][y0 + i])
                {
                    // Update z-buffer
                    zBuffer[(unsigned int)round(x0 + (i / m))][y0 + i] = invZ;
                    (image)((unsigned int)round(x0 + (i / m)), y0 + i) = newColor;
                }
            }
        }
        else if (m < -1.0)
        {
            for (unsigned int i = 0; i <= (y0 - y1); i++)
            {
                double invZ;

                invZ = calculateInvZ(i, 0, (y0 - y1), z1, z0);
                if (invZ < zBuffer[(unsigned int)round(x0 - (i / m))][y0 - i])
                {
                    // Update z-buffer
                    zBuffer[(unsigned int)round(x0 - (i / m))][y0 - i] = invZ;
                    (image)((unsigned int)round(x0 - (i / m)), y0 - i) = newColor;
                }
            }
        }
    }
}

img::EasyImage draw2DzBufferLines(const lines2D &lines, const int size,
                                  std::vector<double> &backgroundColor)
{
    // Declare colors vector
    std::vector<double> originalColors;
    std::vector<unsigned int> newColors;
    // Scale colors
    std::vector<unsigned int> bgColor = scaleColors(backgroundColor);
    // Determine xMin , yMin , xMax , yMax
    double xMin = 0, yMin = 0, xMax = 0, yMax = 0;
    // Loop through all the lines
    for (const line2D &l : lines)
    {
        // Assign the points to temp variables
        double p1X = l.p1.x;
        double p1Y = l.p1.y;
        double p2X = l.p2.x;
        double p2Y = l.p2.y;
        // Check coordinates of the first point
        if (xMax <= p1X)
        {
            xMax = p1X;
        }
        if (xMin >= p1X)
        {
            xMin = p1X;
        }
        if (yMax <= p1Y)
        {
            yMax = p1Y;
        }
        if (yMin >= p1Y)
        {
            yMin = p1Y;
        }
        // Check coordinates of the second point
        if (xMax <= p2X)
        {
            xMax = p2X;
        }
        if (xMin >= p2X)
        {
            xMin = p2X;
        }
        if (yMax <= p2Y)
        {
            yMax = p2Y;
        }
        if (yMin >= p2Y)
        {
            yMin = p2Y;
        }
    }

    // Calculate range along the x-axis and y-axis
    double xRange = xMax - xMin;
    double yRange = yMax - yMin;
    // Calculate image dimensions
    double imageWidth = size * (xRange / std::max(xRange, yRange));
    double imageHeight = size * (yRange / std::max(xRange, yRange));
    // Create the image file
    img::EasyImage image(lround(imageWidth), lround(imageHeight));
    image.clear(img::Color(bgColor.at(0), bgColor.at(1), bgColor.at(2)));
    // Create Z-Buffer
    zBuffer zBuffer(static_cast<int>(image.get_width()),
                    static_cast<int>(image.get_height()));
    // Determine the scaling factor
    double d = 0.95 * (imageWidth / xRange);
    // Calculate for x and y
    double dcX = d * ((xMin + xMax) / 2);
    double dcY = d * ((yMin + yMax) / 2);
    double dx = imageWidth / 2 - dcX;
    double dy = imageHeight / 2 - dcY;
    // Loop through the lines again
    for (line2D l : lines)
    {
        // Multiply all the points by the scaling factor
        l.p1.x *= d;
        l.p1.y *= d;
        l.p2.x *= d;
        l.p2.y *= d;
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
        originalColors = {l.lineColor.red, l.lineColor.green, l.lineColor.blue};
        newColors = scaleColors(originalColors);
        // Draw the lines
        drawZBufferLine(zBuffer, image, x1, y1, z1, x2, y2, z2,
                        color(newColors.at(0), newColors.at(1), newColors.at(2)));
    }
    return image;
}

double calculateIntersection(const int &yI, point2D const &p, point2D const &q)
{
    return q.x + (p.x - q.x) * ((yI - q.y) / (p.y - q.y));
}

double getAngle(const Vector3D &a, const Vector3D &b, const Vector3D &c,
                const Vector3D &direction)
{
    Vector3D n = Vector3D::cross(b - a, c - a);
    n.normalise();
    Vector3D l = -direction;
    return Vector3D::dot(l, n);
}

void drawZBufferTriangle(zBuffer &buffer, img::EasyImage &image,
                         Vector3D const &a, Vector3D const &b, Vector3D const &c,
                         double &d, double &dx, double &dy,
                         color &ambientReflection, color &diffuseReflection,
                         color &specularReflection,
                         double &reflectionCoefficient,
                         lights3D &lights)
{
    // Convert color
    color temp;
    color newCol;
    color finalCol;
    double infAngle;
    for (auto *&l : lights)
    {
        // Get ambient light
        temp = l->ambientLight * ambientReflection;
        newCol = newCol + temp;
        // Get diffuse light
        auto infL = dynamic_cast<infLight *>(l);
        if (infL != nullptr)
        {
            infAngle = getAngle(a, b, c, infL->ldVector);
            if (infAngle > 0)
            {
                temp = l->diffuseLight * diffuseReflection * infAngle;
                newCol = newCol + temp;
            }
        }
    }
    img::Color newColor;
    // Projection of the triangle
    point2D newA = doProjection(a, d);
    newA.x += dx;
    newA.y += dy;
    point2D newB = doProjection(b, d);
    newB.x += dx;
    newB.y += dy;
    point2D newC = doProjection(c, d);
    newC.x += dx;
    newC.y += dy;
    // Determine yMin and yMax
    double yMin = std::min({newA.y, newB.y, newC.y});
    double yMax = std::max({newA.y, newB.y, newC.y});
    // Determine gravity center
    Vector3D g = Vector3D::point((newA.x + newB.x + newC.x) / 3,
                                 (newA.y + newB.y + newC.y) / 3, 0);
    double invZg = 1 / (3 * a.z) + 1 / (3 * b.z) + 1 / (3 * c.z);
    // Determine dz_dx and dz_dy
    Vector3D u = b - a;
    Vector3D v = c - a;
    Vector3D w = Vector3D::cross(u, v);
    double k = w.x * a.x + w.y * a.y + w.z * a.z;
    double dzDx = w.x / (-d * k);
    double dzDy = w.y / (-d * k);
    w.normalise();
    double invZ;
    // Create variables for xL and xR
    double xLAb = posInf;
    double xLAc = posInf;
    double xLBc = posInf;
    double xRAb = negInf;
    double xRAc = negInf;
    double xRBc = negInf;
    int xL;
    int xR;
    // Determine which pixels belong in the triangle
    for (int i = (int)round(yMin + 0.5); i <= (int)round(yMax - 0.5); ++i)
    {
        // Determine xL and xR for each line
        if ((i - newA.y) * (i - newB.y) <= 0 && (newA.y != newB.y))
        {
            xLAb = xRAb = calculateIntersection(i, newA, newB);
        }
        if ((i - newB.y) * (i - newC.y) <= 0 && (newB.y != newC.y))
        {
            xLBc = xRBc = calculateIntersection(i, newB, newC);
        }
        if ((i - newA.y) * (i - newC.y) <= 0 && (newA.y != newC.y))
        {
            xLAc = xRAc = calculateIntersection(i, newA, newC);
        }
        // Determine xL and xR for y = y_i line
        xL = static_cast<int>(round(std::min({xLAb, xLAc, xLBc}) + 0.5));
        xR = static_cast<int>(round(std::max({xRAb, xRAc, xRBc}) - 0.5));
        for (int j = xL; j <= xR; ++j)
        {
            // Determine the invZ value
            invZ = 1.0001 * invZg + (j - g.x) * dzDx + (i - g.y) * dzDy;
            if (invZ < buffer[j][i])
            {
                // Get point light color
                finalCol = newCol;
                for (auto *l : lights)
                {
                    // Get original point
                    double x = -(j - dx) / (invZ * d);
                    double y = -(i - dy) / (invZ * d);
                    Vector3D p = Vector3D::point(x, y, 1 / invZ);
                    // Point Lights
                    auto pntL = dynamic_cast<pointLight *>(l);
                    if (pntL != nullptr)
                    {
                        // Get l vector : distance between point and pointLight
                        Vector3D lV = Vector3D::vector(pntL->location - p);
                        lV.normalise();
                        // Get the angle
                        double angle = Vector3D::dot(lV, w);
                        double pntAngle = cos(pntL->spotAngle);
                        if (angle > pntAngle)
                        {
                            temp = pntL->diffuseLight * diffuseReflection *
                                   (1 - ((1 - angle) / (1 - pntAngle)));
                            finalCol = finalCol + temp;
                        }
                        // Calculate the specular light
                        lV = Vector3D::vector(pntL->location - p);
                        lV.normalise();
                        Vector3D origin = Vector3D::point(0, 0, 0);
                        origin -= p;
                        origin.normalise();
                        Vector3D r = 2 * w * angle - lV;
                        r.normalise();
                        angle = Vector3D::dot(r, origin);
                        if (angle > 0)
                        {
                            temp = pntL->specularLight * specularReflection *
                                   pow(angle, reflectionCoefficient);
                            finalCol = finalCol + temp;
                        }
                    }
                    // Specular for light source on infinity
                    auto infL = dynamic_cast<infLight *>(l);
                    if (infL != nullptr)
                    {
                        // Calculate the specular light
                        double angle = infAngle;
                        Vector3D origin = Vector3D::point(0, 0, 0);
                        origin -= p;
                        origin.normalise();
                        Vector3D r = 2 * w * angle + infL->ldVector;
                        r.normalise();
                        angle = Vector3D::dot(r, origin);
                        if (angle > 0)
                        {
                            temp = infL->specularLight * specularReflection *
                                   pow(angle, reflectionCoefficient);
                            finalCol = finalCol + temp;
                        }
                    }
                }
                // Check for overshot
                if (finalCol.red > 1)
                {
                    finalCol.red = 1;
                }
                if (finalCol.green > 1)
                {
                    finalCol.green = 1;
                }
                if (finalCol.blue > 1)
                {
                    finalCol.blue = 1;
                }
                newColor = img::Color(lround(finalCol.red * 255),
                                      lround(finalCol.green * 255),
                                      lround(finalCol.blue * 255));
                buffer[j][i] = invZ;
                image(j, i) = newColor;
            }
        }
        // Reset left and right limits
        xLAb = posInf;
        xLBc = posInf;
        xLAc = posInf;
        xRAb = negInf;
        xRBc = negInf;
        xRAc = negInf;
    }
}

img::EasyImage
draw2DzBufferTriangle(const int &size, std::vector<double> &backgroundColor,
                      figures3D &figures,
                      lights3D &lights)
{
    // Scale colors
    std::vector<unsigned int> bgColor = scaleColors(backgroundColor);
    color ambient;
    color diffuse;
    color specular;
    // Triangulate : get the new faces
    std::vector<face> newFaces;
    for (auto &fig : figures)
    {
        for (auto &f : fig->faces)
        {
            triangulate(f, newFaces);
        }
        // Replace faces
        fig->faces = newFaces;
        newFaces.clear();
    }
    // Calculate range along the x-axis and y-axis
    double xMin = 0, yMin = 0, xMax = 0, yMax = 0;
    double d = 1.0;
    double xRange;
    double yRange;
    getProjectedPoints(figures, d, xRange, yRange, xMin, yMin, xMax, yMax);
    // Calculate image dimensions
    double imageWidth = size * (xRange / std::max(xRange, yRange));
    double imageHeight = size * (yRange / std::max(xRange, yRange));
    // Create the image file
    img::EasyImage image(lround(imageWidth), lround(imageHeight),
                         img::Color(bgColor.at(0), bgColor.at(1), bgColor.at(2)));
    // Create Z-Buffer
    zBuffer zBuffer((int)image.get_width(), (int)image.get_height());
    // Determine the scaling factor
    d = 0.95 * (imageWidth / xRange);
    // Calculate for x and y
    double dcX = d * ((xMin + xMax) / 2);
    double dcY = d * ((yMin + yMax) / 2);
    double dx = imageWidth / 2 - dcX;
    double dy = imageHeight / 2 - dcY;
    // Apply the z-buffering algorithm
    for (auto &fig : figures)
    {
        // r_a
        ambient = fig->ambientReflection;
        // rd_d
        diffuse = fig->diffuseReflection;
        // r_s
        specular = fig->specularReflection;
        for (auto &f : fig->faces)
        {
            // Get the point indexes
            int indA = f.pointIndexes.at(0);
            int indB = f.pointIndexes.at(1);
            int indC = f.pointIndexes.at(2);
            // Apply z-buffering algorithm
            drawZBufferTriangle(zBuffer, image, fig->points.at(indA),
                                fig->points.at(indB), fig->points.at(indC),
                                d, dx, dy, ambient, diffuse, specular,
                                fig->reflectionCoefficient, lights);
        }
    }
    return image;
}
