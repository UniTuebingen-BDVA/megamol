
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

    // clear memory of prev arrays
    std::vector<glm::vec3>().swap(vertices);
    std::vector<glm::vec3>().swap(normals);
    std::vector<unsigned int>().swap(indices);
    std::vector<unsigned int>().swap(lineIndices);
    std::vector<unsigned int>().swap(faceIndices);
    int currentIndex = 0;

    for (int i = 0; i <= numRings; ++i) {
        for (int j = 0; j <= numSegments; ++j) {
            float theta = 2 * pi * static_cast<float>(j) / static_cast<float>(numSegments);
            float phi = 2 * pi * static_cast<float>(i) / static_cast<float>(numRings);

            float x = center.x + (radius + probeRadius * cos(phi)) * cos(theta);
            float y = center.y + (radius + probeRadius * cos(phi)) * sin(theta);
            float z = center.z + probeRadius * sin(phi);


            addVertex(x, y, z);

            glm::vec3 normal = glm::normalize(glm::vec3(x - center.x, y - center.y, z - center.z));
            addNormal(normal);

            indices.push_back(vertices.size() - 1);

            lineIndices.push_back(indices.size() - 1);
            lineIndices.push_back((indices.size() - 1 + numSegments) % (numRings * (numSegments + 1)));
        }
    }

    for (int i = 0; i < numRings; ++i) {
        for (int j = 0; j < numSegments; ++j) {

            //first
            faceIndices.push_back(i * (numSegments + 1) + j);
            faceIndices.push_back(i * (numSegments + 1) + (j + 1) % numSegments);
            faceIndices.push_back((i + 1) * (numSegments + 1) + j);

            //second
            faceIndices.push_back((i + 1) * (numSegments + 1) + j);
            faceIndices.push_back(i * (numSegments + 1) + (j + 1) % numSegments);
            faceIndices.push_back((i + 1) * (numSegments + 1) + (j + 1) % numSegments);
            
        }
    }
}

const float Torus::getDistance(glm::vec3 atomPosition1, glm::vec3 atomPosition2) {
    float distance = glm::distance(atomPosition1, atomPosition2);

    return distance; //d_ij
}

const glm::vec3 Torus::getTorusAxisUnitVec(glm::vec3 atomPosition1, glm::vec3 atomPosition2) {
    float distance = getDistance(atomPosition1, atomPosition2);
    glm::vec3 axisUnitVec = (atomPosition2 - atomPosition2) / distance;

    return axisUnitVec; //u_ij
}

const glm::vec3 Torus::getTorusCenter(glm::vec3 atomPosition1, glm::vec3 atomPosition2, float atomRadius1, float atomRadius2, float probeRadius) {
    float distance = getDistance(atomPosition1, atomPosition2);
    glm::vec3 torusCenter = 0.5f * (atomPosition1 + atomPosition2) +
                            0.5f * (atomPosition2 - atomPosition1) *
                            ((std::powf((atomRadius1 + probeRadius), 2)) - (std::powf((atomRadius2 + probeRadius), 2))) / powf(distance, 2);

    return torusCenter; //t_ij
}

const float Torus::getTorusRadius(float atomRadius1, float atomRadius2, float probeRadius, float distance) {
    float torusRadius = 0.5f * powf((powf((atomRadius1 + atomRadius2 + 2 * probeRadius), 0.5f) - powf(distance, 2)), 0.5f) *
                        powf((powf(distance, 2) - powf((atomRadius1 - atomRadius2), 2)), 0.5f) / distance;

    return torusRadius; //r_ij
}

const float Torus::getBaseTriangleAngle(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3) {
    glm::vec3 axisUnitVector1 = getTorusAxisUnitVec(atomPosition1, atomPosition2);
    glm::vec3 axisUnitVector2 = getTorusAxisUnitVec(atomPosition1, atomPosition3);

    float dotProduct = glm::dot(axisUnitVector1, axisUnitVector2);

    /*
    float length1 = glm::length(axisUnitVector1);
    float length2 = glm::length(axisUnitVector2);

    float baseTriangleAngle_rad = std::acos(dotProduct / (length1 * length2));
    float baseTrinagleAngle_deg = glm::degrees(baseTriangleAngle_rad);
    */
    float baseTriangleAngle_rad = std::acos(dotProduct);
    // rad or deg?
    return baseTriangleAngle_rad; //w_ijk
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

    return glm::normalize(basePlaneNormal); //u_ijk
}

const glm::vec3 Torus::getTorusBasePointUnitVector(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3) {
    glm::vec3 basePlaneNormal = getBasePlaneNormal(atomPosition1, atomPosition2, atomPosition3);
    glm::vec3 torusAxisUnitVector = getTorusAxisUnitVec(atomPosition1, atomPosition2);

    glm::vec3 torusBasePointUnitVector = glm::cross(basePlaneNormal, torusAxisUnitVector);

    return torusBasePointUnitVector; //u_tb
}

const glm::vec3 Torus::getBasePoint(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3,
                                    float atomRadius1, float atomRadius2, float atomRadius3, float probeRadius) {
    glm::vec3 torusCenter_12 = getTorusCenter(atomPosition1, atomPosition2, atomRadius1, atomRadius2, probeRadius);
    glm::vec3 torusCenter_13 = getTorusCenter(atomPosition1, atomPosition3, atomRadius1, atomRadius3, probeRadius);
    glm::vec3 torusBasePointUnitVector = getTorusBasePointUnitVector(atomPosition1, atomPosition2, atomPosition3);
    glm::vec3 torusAxisUnitVector_13 = getTorusAxisUnitVec(atomPosition1, atomPosition3);
    float baseTriangleAngle = getBaseTriangleAngle(atomPosition1, atomPosition2, atomPosition3);
    float angleFactor = 1 / baseTriangleAngle;

    glm::vec3 basePoint = torusCenter_12 +
                          torusBasePointUnitVector *
                          (glm::dot(torusAxisUnitVector_13,(torusCenter_13 - torusCenter_12))) *
                          angleFactor;

    return basePoint; //b_ijk
}

const float Torus::getProbeHeight(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3,
                                  float atomRadius1, float atomRadius2, float atomRadius3, float probeRadius){
    glm::vec3 basePoint = getBasePoint(atomPosition1, atomPosition2, atomPosition3, atomRadius1, atomRadius2, atomRadius3, probeRadius);

    float radPow = powf(atomRadius1 + probeRadius, 2);
    float pointPow = powf(glm::distance(basePoint, atomPosition1), 2);
    float probeHeight = sqrtf(radPow - pointPow);

    return probeHeight; //h_ijk
};

const glm::vec3 Torus::getProbePosition(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3,
                                        float atomRadius1, float atomRadius2, float atomRadius3, float probeRadius) {
    //b_ijk, h_ijk, u_ijk
    glm::vec3 basePoint = getBasePoint(atomPosition1, atomPosition2, atomPosition3, atomRadius1, atomRadius2, atomRadius3, probeRadius);
    float probeHeight = getProbeHeight(atomPosition1, atomPosition2, atomPosition3, atomRadius1, atomRadius2, atomRadius3, probeRadius);
    glm::vec3 basePlaneNormal = getBasePlaneNormal(atomPosition1, atomPosition2, atomPosition3);
    //TODO
    glm::vec3 probePosition = basePoint;

    return probePosition; //p_ijk
}

const glm::vec3 Torus::getProbeVertex(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3,
                                      float atomRadius1, float atomRadius2, float atomRadius3, float probeRadius) {
    glm::vec3 probePosition = getProbePosition(atomPosition1, atomPosition2, atomPosition3, atomRadius1, atomRadius2, atomRadius3, probeRadius);
    glm::vec3 dividend = atomRadius1 * probePosition + probeRadius * atomPosition1;

    glm::vec3 vertex = dividend / (atomRadius1 + probeRadius);

    return vertex; //v_pi
}

const glm::vec3 Torus::getContactCircleCenter(glm::vec3 atomPosition1, glm::vec3 atomPosition2, float atomRadius1, float atomRadius2, float probeRadius) {
    glm::vec3 torusCenter = getTorusCenter(atomPosition1, atomPosition2, atomRadius1, atomRadius2, probeRadius);
    glm::vec3 dividend = atomRadius1 * torusCenter + probeRadius * atomPosition1;

    glm::vec3 contactCircleCenter = dividend / (atomRadius1 + probeRadius);

    return contactCircleCenter; //c_ij
}

const float Torus::getContactCircleRadius(glm::vec3 atomPosition1, glm::vec3 atomPosition2, float atomRadius1, float atomRadius2, float probeRadius) {
    float distance = getDistance(atomPosition1, atomPosition2);
    float torusRadius = getTorusRadius(atomRadius1, atomRadius2, probeRadius, distance);

    float contactCircleRadius = (torusRadius * atomRadius1) / (atomRadius1 + probeRadius);

    return contactCircleRadius; //r_c
}

const float Torus::getContactCircleDisplacement(glm::vec3 atomPosition1, glm::vec3 atomPosition2, float atomRadius1, float atomRadius2, float probeRadius) {
    glm::vec3 torusAxisUnitVec = getTorusAxisUnitVec(atomPosition1, atomPosition2);
    glm::vec3 contactCircleCenter = getContactCircleCenter(atomPosition1, atomPosition2, atomRadius1, atomRadius2, probeRadius);

    float contactCirlceDisplacement = glm::dot(torusAxisUnitVec, (contactCircleCenter - atomPosition1));

    return contactCirlceDisplacement;
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

void Torus::addNormal(glm::vec3 normal) {
    normals.push_back(normal);
}

void Torus::addIndices(unsigned int i1, unsigned int i2, unsigned int i3) {
    indices.push_back(i1);
    indices.push_back(i2);
    indices.push_back(i3);
}
