#include "draw3d.h"

figure *
drawLineDrawing(double &scale, double &rotX, double &rotY, double &rotZ, int &nrPoints,
                int &nrLines,
                const ini::Configuration &configuration, std::vector<double> &lineColor,
                std::vector<double> &center,
                Matrix &mEye, int &i)
{
        // Make new figure class
        std::vector<Vector3D> points;
        std::vector<face> faces;

        // Make temporary vector;
        std::vector<double> pointToAdd;
        std::vector<int> lineToAdd;
        // Get points
        for (int j = 0; j < nrPoints; j++)
        {
                pointToAdd = configuration["figure" + std::to_string(i)]["point" + std::to_string(
                                                                                       j)]
                                 .as_double_tuple_or_die();
                points.push_back(
                    Vector3D::point(pointToAdd.at(0), pointToAdd.at(1), pointToAdd.at(2)));
        }
        // Get lines
        for (int j = 0; j < nrLines; j++)
        {
                lineToAdd = configuration["figure" + std::to_string(i)]["line" + std::to_string(
                                                                                     j)]
                                .as_int_tuple_or_die();
                faces.push_back(face({lineToAdd.at(0), lineToAdd.at(1)}));
        }
        // Create new figure
        figure *newFigure;
        newFigure = new figure(points, faces,
                               color(lineColor.at(0), lineColor.at(1), lineColor.at(2)));
        Matrix m = scaleFigure(scale);
        m *= rotateX((rotX * M_PI) / 180);
        m *= rotateY((rotY * M_PI) / 180);
        m *= rotateZ((rotZ * M_PI) / 180);
        m *= translate(Vector3D::point(center.at(0), center.at(1), center.at(2)));
        m *= mEye;
        newFigure->applyTransformation(m);
        return newFigure;
}
