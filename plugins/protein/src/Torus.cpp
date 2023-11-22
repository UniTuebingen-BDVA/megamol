
#ifdef _WIN32
#include <windows.h> // include windows.h to avoid thousands of compile errors even though this class is not depending on Windows
#endif

#include "protein/Torus.h"
#include <cmath>
#include <glm/glm.hpp>

#ifndef M_PI
    #define M_PI 3.141592653589793238462643383279502884L
#endif

Torus::Torus(glm::vec3 center, float radius, int numSegments, int numRings)
        : center(center)
        , radius(radius)    //inner radius
        , numSegments(numSegments)
        , numRings(numRings) {}

void Torus::generateTorus(float probeRadius) {
    float pi = M_PI;
    for (int i = 0; i <= numRings; ++i) {
        for (int j = 0; j <= numSegments; ++j) {
            float theta = 2 * pi * static_cast<float>(j) / static_cast<float>(numSegments);
            float phi = 2 * pi * static_cast<float>(i) / static_cast<float>(numRings);

            float x = (radius + probeRadius * cos(phi)) * cos(theta);
            float y = (radius + probeRadius * cos(phi)) * sin(theta);
            float z = probeRadius * sin(phi);

            // Füge den Vertex hinzu
            addVertex(x, y, z);

            glm::vec3 normal = glm::normalize(glm::vec3(x, y, z));
            addNormal(normal.x, normal.y, normal.z);
        }
    }

    // Erzeuge die Indizes für die Dreiecke
    for (int i = 0; i < numRings; ++i) {
        for (int j = 0; j < numSegments; ++j) {
            // Erster Punkt im Dreieck
            addIndices(i * (numSegments + 1) + j,
                      (i + 1) * (numSegments + 1) + j,
                       i * (numSegments + 1) + j + 1);
            

            // Zweiter Punkt im Dreieck
            addIndices((i + 1) * (numSegments + 1) + j,
                       (i + 1) * (numSegments + 1) + j + 1,
                        i * (numSegments + 1) + j + 1);

            unsigned int currentIndex = i * (numSegments + 1) + j;
            faceIndices.push_back(currentIndex);
        }
    }
}

const float Torus::getDistance(glm::vec3 atomPosition1, glm::vec3 atomPosition2) {
    float distance = glm::distance(atomPosition1, atomPosition2);

    return distance;
}

const glm::vec3 Torus::getTorusAxisUnitVec(glm::vec3 atomPosition1, glm::vec3 atomPosition2) {
    float distance = getDistance(atomPosition1, atomPosition2);
    glm::vec3 axisUnitVec = (atomPosition2 - atomPosition2) / distance;

    return axisUnitVec;
}

const glm::vec3 Torus::getTorusCenter(glm::vec3 atomPosition1, glm::vec3 atomPosition2, float atomRadius1, float atomRadius2, float probeRadius) {
    float distance = getDistance(atomPosition1, atomPosition2);
    glm::vec3 torusCenter = 0.5f * (atomPosition1 + atomPosition2) +
                            0.5f * (atomPosition2 - atomPosition1) *
                            ((std::powf((atomRadius1 + probeRadius), 2)) - (std::powf((atomRadius2 + probeRadius), 2))) / powf(distance, 2);

    return torusCenter;
}

const float Torus::getTorusRadius(float atomRadius1, float atomRadius2, float probeRadius, float distance) {
    float torusRadius = 0.5f * powf((powf((atomRadius1 + atomRadius2 + 2 * probeRadius), 0.5f) - powf(distance, 2)), 0.5f) *
                         powf((powf(distance, 2) - powf((atomRadius1 - atomRadius2), 2)), 0.5f) / distance;

    return torusRadius;
}

const float Torus::getBaseTriangleAngle(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3) {
    glm::vec3 axisUnitVector1 = getTorusAxisUnitVec(atomPosition1, atomPosition2);
    glm::vec3 axisUnitVector2 = getTorusAxisUnitVec(atomPosition1, atomPosition3);

    float dotProduct = glm::dot(axisUnitVector1, axisUnitVector2);

    float length1 = glm::length(axisUnitVector1);
    float length2 = glm::length(axisUnitVector2);

    float baseTriangleAngle_rad = std::acos(dotProduct / (length1 * length2));
    float baseTrinagleAngle_deg = glm::degrees(baseTriangleAngle_rad);

    // rad or deg?
    return baseTriangleAngle_rad;
}

const glm::vec3 Torus::getBasePlaneNormal(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3) {
    float angle = getBaseTriangleAngle(atomPosition1, atomPosition2, atomPosition3);
    glm::vec3 axisUnitVector1 = getTorusAxisUnitVec(atomPosition1, atomPosition2);
    glm::vec3 axisUnitVector2 = getTorusAxisUnitVec(atomPosition1, atomPosition3);

    /*
    float cosTheta = std::cos(angle);
    float sinTheta = std::sin(angle);

    glm::vec3 axisUnitVector2_rotated = rotateVector(cosTheta, sinTheta, axisUnitVector2, axisUnitVector1);

    glm::vec3 basePlaneNormal = glm::cross(axisUnitVector1, axisUnitVector2_rotated);
    */
    glm::vec3 basePlaneNormal = glm::cross(axisUnitVector1, axisUnitVector2) / std::sin(angle);
    return glm::normalize(basePlaneNormal);
}

const glm::vec3 Torus::getTorusBasePointUnitVector(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3) {
    glm::vec3 basePlaneNormal = getBasePlaneNormal(atomPosition1, atomPosition2, atomPosition3);
    glm::vec3 torusaxisUnitVector = getTorusAxisUnitVec(atomPosition1, atomPosition2);

    glm::vec3 torusBasePointUnitVector = glm::cross(basePlaneNormal, torusaxisUnitVector);
    return torusBasePointUnitVector;
}

const glm::vec3 Torus::getBasePoint(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3,
                                    float atomRadius1, float atomRadius2, float atomRadius3, float probeRadius) {
    glm::vec3 torusCenter_12 = getTorusCenter(atomPosition1, atomPosition2, atomRadius1, atomRadius2, probeRadius);
    glm::vec3 torusCenter_13 = getTorusCenter(atomPosition1, atomPosition3, atomRadius1, atomRadius3, probeRadius);
    glm::vec3 torusBasePointUnitVector = getTorusBasePointUnitVector(atomPosition1, atomPosition2, atomPosition3);
    glm::vec3 torusAxisUnitVector_13 = getTorusAxisUnitVec(atomPosition1, atomPosition3);
    float baseTriangleAngle = getBaseTriangleAngle(atomPosition1, atomPosition2, atomPosition3);

    glm::vec3 basePoint = torusCenter_12 +
                          torusBasePointUnitVector * (torusAxisUnitVector_13);

    return glm::vec3(0, 0, 0);
}



const glm::vec3 Torus::rotateVector(float cosTheta, float sinTheta, glm::vec3 rotateVector, glm::vec3 rotateAxis) {
    glm::vec3 rotatedVec = cosTheta * rotateVector +
                           sinTheta + glm::cross(rotateAxis, rotateVector) +
                           (1 - cosTheta) * glm::dot(rotateAxis, rotateVector) * rotateAxis;

    return rotatedVec;
}

void Torus::addVertex(float x, float y, float z) {
    vertices.push_back(glm::vec3(x, y, z));
}

void Torus::addNormal(float nx, float ny, float nz) {
    normals.push_back(glm::vec3(nx, ny, nz));
}

void Torus::addIndices(unsigned int i1, unsigned int i2, unsigned int i3) {
    indices.push_back(i1);
    indices.push_back(i2);
    indices.push_back(i3);
}
