#include "shapes.h"

void figure::applyTransformation(const Matrix &m) {
    // Multiply each vector with the matrix
    for (Vector3D &v: points) {
        v *= m;
    }
}