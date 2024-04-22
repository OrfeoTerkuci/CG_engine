//
// Created by oterk on 19/03/2024.
//

#ifndef ENGINE_LIGHTS_H
#define ENGINE_LIGHTS_H

#include <iostream>
#include "../color/color.h"
#include "../vector/vector3d.h"
#include <vector>

class light {
public:
    // The ambient light coefficient
    color ambientLight;
    // The diffuse light coefficient
    color diffuseLight;
    // The specular light coefficient
    color specularLight;

    light(const color &ambientLight, const color &diffuseLight,
          const color &specularLight) : ambientLight(
            ambientLight), diffuseLight(diffuseLight), specularLight(specularLight) {}

    virtual ~light() = default;
};

class infLight : public light {
public:
    // Direction of lighting
    Vector3D ldVector;

    infLight(const color &ambientLight,
                                     const color &diffuseLight,
                                     const color &specularLight,
                                     const Vector3D &ldVector)
            : light(ambientLight, diffuseLight, specularLight), ldVector(ldVector) {}

    infLight(const color &ambientLight, const color &diffuseLight,
             const color &specularLight,
             std::vector<double> &ldVec)
            : light(ambientLight, diffuseLight, specularLight),
              ldVector(Vector3D::vector(ldVec.at(0), ldVec.at(1), ldVec.at(2))) {}

};

class pointLight : public light {
public:
    // The light source position
    Vector3D location;
    // The angle of the spotlight
    double spotAngle;

    pointLight(const color &ambientLight,
                                       const color &diffuseLight,
                                       const color &specularLight,
                                       const Vector3D &location, double spotAngle)
            : light(ambientLight, diffuseLight,
                    specularLight),
              location(location),
              spotAngle(spotAngle) {}

    pointLight(const color &ambientLight, const color &diffuseLight,
               const color &specularLight,
               std::vector<double> &location, double spotAngle) : light(ambientLight,
                                                                        diffuseLight,
                                                                        specularLight),
                                                                  location(
                                                                          Vector3D::point(
                                                                                  location.at(
                                                                                          0),
                                                                                  location.at(
                                                                                          1),
                                                                                  location.at(
                                                                                          2))),
                                                                  spotAngle(
                                                                          spotAngle) {}

};

typedef std::vector<light *> lights3D;

#endif //ENGINE_LIGHTS_H
