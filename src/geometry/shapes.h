#ifndef ENGINE_SHAPES_H
#define ENGINE_SHAPES_H

#include <iostream>
#include <vector>

#include "../color/color.h"
#include "../vector/vector3d.h"

class point2D
{
public:
        double x;
        double y;

        /**
         * @brief Constructs a point2D type object
         * @param x x-coordinate
         * @param y y-coordinate
         */
        point2D(double x, double y) : x(x), y(y) {}

        virtual ~point2D() = default;
};

class line2D
{
public:
        point2D p1;
        point2D p2;
        color lineColor;

        double z1;
        double z2;

        /**
         *
         * @param p1 Begin point : point2D type object
         * @param p2 End point : point2D type object
         * @param lineColor Line color : color type object
         */
        line2D(const point2D &p1, const point2D &p2, const color &color) : p1(p1),
                                                                           p2(p2),
                                                                           lineColor(color),
                                                                           z1(1),
                                                                           z2(1) {}

        /**
         *
         * @param p1 Begin point : point2D type object
         * @param p2 End point : point2D type object
         * @param color Line color : color type object
         * @param z1 z-coordinate of begin point
         * @param z2 z-coordinate of end point
         */
        line2D(const point2D &p1, const point2D &p2, const color &color, double &z1,
               double &z2) : p1(p1), p2(p2),
                             lineColor(color),
                             z1(z1), z2(z2) {}

        /**
         * @brief Constructs a line2D from a reference object
         * @param refLine A pointer to a line2D type object
         */
        explicit line2D(line2D *refLine) : p1(refLine->p1), p2(refLine->p2),
                                           lineColor(refLine->lineColor), z1(refLine->z1),
                                           z2(refLine->z2) {}

        virtual ~line2D() = default;
};

class face
{

public:
        // These indexes refer to points in the 'points' vector of the figure-class
        std::vector<int> pointIndexes;

        explicit face(const std::vector<int> &pointIndexes) : pointIndexes(pointIndexes) {}

        face(face *refFace) : pointIndexes(refFace->pointIndexes) {}

        virtual ~face() = default;
};

class figure
{

public:
        std::vector<Vector3D> points;
        std::vector<face> faces;

        color ambientReflection;
        color diffuseReflection;
        color specularReflection;

        double reflectionCoefficient;

        figure(const std::vector<Vector3D> &points, const std::vector<face> &faces,
               const color &figureColor) : points(points),
                                           faces(faces),
                                           ambientReflection(figureColor),
                                           diffuseReflection(
                                               color(0, 0, 0)),
                                           specularReflection(
                                               color(0, 0, 0)),
                                           reflectionCoefficient(0) {}

        figure(const std::vector<Vector3D> &points, const std::vector<face> &faces,
               const color &ambient, const color &diffuse, const color &specular,
               double reflectionCoefficient = 0) : points(points), faces(faces),
                                                   ambientReflection(ambient),
                                                   diffuseReflection(diffuse),
                                                   specularReflection(specular),
                                                   reflectionCoefficient(
                                                       reflectionCoefficient) {}

        figure(const std::vector<Vector3D> &points, const std::vector<face> &faces,
               const std::vector<double> &ambient, const std::vector<double> &diffuse,
               const std::vector<double> &specular,
               double reflectionCoefficient = 0) : points(points), faces(faces),
                                                   ambientReflection(ambient), diffuseReflection(diffuse),
                                                   specularReflection(specular),
                                                   reflectionCoefficient(reflectionCoefficient) {}

        explicit figure(figure *refFig) : points(refFig->points), faces(refFig->faces),
                                          ambientReflection(refFig->ambientReflection),
                                          diffuseReflection(refFig->diffuseReflection),
                                          specularReflection(
                                              refFig->specularReflection),
                                          reflectionCoefficient(
                                              refFig->reflectionCoefficient) {}

        void applyTransformation(const Matrix &m);

        virtual ~figure() = default;
};

typedef std::vector<line2D> lines2D;

typedef std::vector<figure *> figures3D;

#endif // ENGINE_SHAPES_H
