
#ifdef _WIN32
#include <windows.h> // include windows.h to avoid thousands of compile errors even though this class is not depending on Windows
#endif

#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtx/rotate_vector.hpp"
#include "glm/gtx/string_cast.hpp"
#include "protein/Torus.h"
#include <cmath>
#include <glm/glm.hpp>
#include <iostream>


#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884L
#endif

Torus::Torus(glm::vec3 center, float radius, int numSegments, int numRings)
        : center(center)
        , radius(radius) //inner radius
        , numSegments(numSegments)
        , numRings(numRings) {}

void Torus::generateTorus(float probeRadius, glm::vec3 axisUnitVec, float rotationAngle, int offset, int foffset) {
    float pi = M_PI;

    // clear memory of prev arrays
    std::vector<glm::vec3>().swap(vertices);
    std::vector<glm::vec3>().swap(normals);
    std::vector<unsigned int>().swap(indices);
    std::vector<unsigned int>().swap(lineIndices);
    std::vector<unsigned int>().swap(faceIndices);


    for (int i = 0; i <= numRings; ++i) {
        for (int j = 0; j <= numSegments; ++j) {
            float theta = 2 * pi * static_cast<float>(j) / static_cast<float>(numSegments);
            float phi = 2 * pi * static_cast<float>(i) / static_cast<float>(numRings);

            float x = (/*2 + */ radius + probeRadius * cos(phi)) * cos(theta);
            float y = (/*2 + */ radius + probeRadius * cos(phi)) * sin(theta);
            float z = probeRadius * sin(phi);

            glm::vec3 vertex(x, y, z);


            addVertex(x, y, z);

            glm::vec3 normal = glm::normalize(vertex);
            addNormal(normal);

            indices.push_back(vertices.size() - 1 + offset);

            lineIndices.push_back(indices.size() - 1);
            lineIndices.push_back((indices.size() - 1 + numSegments) % (numRings * (numSegments + 1)));

            //if (i < numRings && j < numSegments) {
            //    //point j on ring i
            //    int v1 = i * (numSegments + 1) + j + offset;
            //    //next point (j + 1) on ring i
            //    int v2 = i * (numSegments + 1) + (j + 1) % numSegments + offset;
            //    //point j on ring i + 1
            //    int v3 = (i + 1) * (numSegments + 1) + j + offset;

            //    // Fügen Sie die Face-Indices hinzu
            //    faceIndices.push_back(v1);
            //    faceIndices.push_back(v2);
            //    faceIndices.push_back(v3);

            //    vertexFaceIndex[v1] = glm::vec3(
            //        vertices[(v1 - offset)], vertices[(v1 - offset) + 1], vertices[(v1 - offset) + 2]);
            //    //std::cout << "offset: " << offset << std::endl;
            //    std::cout << "value at " << v1 << " : " << vertexFaceIndex[v1].x << "/" << vertexFaceIndex[v1].y << "/"
            //              << vertexFaceIndex[v1].z << std::endl;

            //    //point j n ring i + 1
            //    int v4 = (i + 1) * (numSegments + 1) + j + offset;
            //    //next point on ring i
            //    int v5 = i * (numSegments + 1) + (j + 1) % numSegments + offset;
            //    //next point (j + 1) on ring i + 1
            //    int v6 = (i + 1) * (numSegments + 1) + (j + 1) % numSegments + offset;

            //    // Fügen Sie die Face-Indices hinzu
            //    faceIndices.push_back(v4);
            //    faceIndices.push_back(v5);
            //    faceIndices.push_back(v6);


            //}
        }
    }


    glm::mat4 translationMatrix = glm::translate(glm::mat4(1.0f), center);
    glm::mat4 rotationMatrix = glm::rotate(glm::mat4(1.0f), rotationAngle, axisUnitVec);
    glm::mat4 transformationMatrix = translationMatrix * rotationMatrix;
    for (int i = 0; i < vertices.size(); ++i) {
        glm::vec4 transformedVertex = transformationMatrix * glm::vec4(vertices[i], 1.0f);
        vertices[i] = glm::vec3(transformedVertex);
    }

    std::cout << "vertices_size: " << vertices.size() << std::endl;

    for (int i = 0; i < numRings; ++i) {
        for (int j = 0; j < numSegments; ++j) {
            //point j on ring i
            int v1 = i * (numSegments + 1) + j + offset;
            //next point (j + 1) on ring i
            int v2 = i * (numSegments + 1) + (j + 1) % numSegments + offset;
            //point j on ring i + 1
            int v3 = (i + 1) * (numSegments + 1) + j + offset;

            // Fügen Sie die Face-Indices hinzu
            faceIndices.push_back(v1);
            faceIndices.push_back(v2);
            faceIndices.push_back(v3);

            std::vector<glm::vec3> vertices_v1 = {vertices[(v1 - offset)], vertices[(v1 - offset) + 1], vertices[(v1 - offset) + 2]};
            vertexFaceIndex[v1] = vertices_v1;

            std::vector<glm::vec3> vertices_v2 = {vertices[(v2 - offset)], vertices[(v2 - offset) + 1], vertices[(v2 - offset) + 2]};
            vertexFaceIndex[v2] = vertices_v2;

            std::vector<glm::vec3> vertices_v3;
            if (v3 - offset + 2 == vertices.size()) {
                vertices_v3 = {vertices[(v3 - offset)], vertices[(v3 - offset) + 1], vertices[0]};
            } else {
                vertices_v3 = {vertices[(v3 - offset)], vertices[(v3 - offset) + 1], vertices[(v3 - offset) + 2]};
            }
            vertexFaceIndex[v3] = vertices_v3;
            //std::cout << v3 - offset << std::endl;
            

            //point j n ring i + 1
            int v4 = (i + 1) * (numSegments + 1) + j + offset;
            //next point on ring i
            int v5 = i * (numSegments + 1) + (j + 1) % numSegments + offset;
            //next point (j + 1) on ring i + 1
            int v6 = (i + 1) * (numSegments + 1) + (j + 1) % numSegments + offset;

            //std::cout << v5 - offset << std::endl;

            // Fügen Sie die Face-Indices hinzu
            faceIndices.push_back(v4);
            faceIndices.push_back(v5);
            faceIndices.push_back(v6);

            std::vector<glm::vec3> vertices_v4;
            if (v4 - offset + 2 == vertices.size()) {
                vertices_v4 = {vertices[(v4 - offset)], vertices[(v4 - offset) + 1], vertices[0]};
            } else {
                vertices_v4 = {vertices[(v4 - offset)], vertices[(v4 - offset) + 1], vertices[(v4 - offset) + 2]};
            }
            vertexFaceIndex[v4] = vertices_v4;

            std::vector<glm::vec3> vertices_v5;
            if (v5 - offset + 2 == vertices.size()) {
                vertices_v5 = {vertices[(v5 - offset)], vertices[(v5 - offset) + 1], vertices[0]};
            } else {
                vertices_v5 = {vertices[(v5 - offset)], vertices[(v5 - offset) + 1], vertices[(v5 - offset) + 2]};
            }
            vertexFaceIndex[v5] = vertices_v5;

            std::vector<glm::vec3> vertices_v6;
            if (v6 - offset + 2 == vertices.size()) {
                vertices_v6 = {vertices[(v6 - offset)], vertices[(v6 - offset) + 1], vertices[0]};
            } else {
                vertices_v6 = {vertices[(v6 - offset)], vertices[(v6 - offset) + 1], vertices[(v6 - offset) + 2]};
            }
            vertexFaceIndex[v6] = vertices_v6;
            //std::cout << v6 - offset << std::endl;

            /*std::cout << "value at " << v1 << " first: " << vertexFaceIndex[v1][0].x << "/" << vertexFaceIndex[v1][0].y << "/"
                      << vertexFaceIndex[v1][0].z << std::endl;
            std::cout << "value at " << v1 << " second: " << vertexFaceIndex[v1][1].x << "/" << vertexFaceIndex[v1][1].y
                      << "/" << vertexFaceIndex[v1][1].z << std::endl;
            std::cout << "value at " << v1 << " third: " << vertexFaceIndex[v1][2].x << "/" << vertexFaceIndex[v1][2].y
                      << "/" << vertexFaceIndex[v1][2].z << std::endl;*/
        }
    }
}
    
    //for (int i = 0; i < numRings; ++i) {
    //    for (int j = 0; j < numSegments; ++j) {

    //        int v1 = i * (numSegments + 1) + j + offset;
    //        int v2 = i * (numSegments + 1) + (j + 1) % numSegments + offset;
    //        int v3 = (i + 1) * (numSegments + 1) + j + offset;

    //        // Fügen Sie die Face-Indices hinzu
    //        faceIndices.push_back(v1);
    //        faceIndices.push_back(v2);
    //        faceIndices.push_back(v3);

    //        // Speichern Sie die Face-Indices für jedes Vertex
    //        vertexFaceIndices[v1] = faceIndices.size() - 3;
    //        vertexFaceIndices[v2] = faceIndices.size() - 2;
    //        vertexFaceIndices[v3] = faceIndices.size() - 1;

    //        int v4 = (i + 1) * (numSegments + 1) + j + offset;
    //        int v5 = i * (numSegments + 1) + (j + 1) % numSegments + offset;
    //        int v6 = (i + 1) * (numSegments + 1) + (j + 1) % numSegments + offset;

    //        // Fügen Sie die Face-Indices hinzu
    //        faceIndices.push_back(v4);
    //        faceIndices.push_back(v5);
    //        faceIndices.push_back(v6);

    //        // Speichern Sie die Face-Indices für jedes Vertex
    //        vertexFaceIndices[v4] = faceIndices.size() - 3;
    //        vertexFaceIndices[v5] = faceIndices.size() - 2;
    //        vertexFaceIndices[v6] = faceIndices.size() - 1;
    //    
    //        ////first
    //        //faceIndices.push_back(i * (numSegments + 1) + j + offset);
    //        //faceIndices.push_back(i * (numSegments + 1) + (j + 1) % numSegments + offset);
    //        //faceIndices.push_back((i + 1) * (numSegments + 1) + j + offset);

    //        ////second
    //        //faceIndices.push_back((i + 1) * (numSegments + 1) + j + offset);
    //        //faceIndices.push_back(i * (numSegments + 1) + (j + 1) % numSegments + offset);
    //        //faceIndices.push_back((i + 1) * (numSegments + 1) + (j + 1) % numSegments + offset);
    //    }
    //}
//}

void Torus::removeTrianglesInsideSphere(glm::vec3 sphereCenter, float sphereRadius, std::vector<glm::vec3> vertices, std::vector<glm::vec3> normals, std::vector<unsigned int> faceIndices) {
    std::vector<glm::vec3> newVertices;
    std::vector<glm::vec3> newNormals;
    std::vector<unsigned int> newIndices;

    for (size_t i = 0; i < faceIndices.size(); i += 3) {
        int v1Index = faceIndices[i];
        int v2Index = faceIndices[i + 1];
        int v3Index = faceIndices[i + 2];

        // Check if indices are within valid range
        if (v1Index < 0 || v1Index >= vertices.size() || v2Index < 0 || v2Index >= vertices.size() || v3Index < 0 ||
            v3Index >= vertices.size()) {
            // Invalid indices, mark this triangle as invalid
            faceIndices[i] = faceIndices[i + 1] = faceIndices[i + 2] = -1;
            continue;
        }

        glm::vec3 v1 = vertices[v1Index];
        glm::vec3 v2 = vertices[v2Index];
        glm::vec3 v3 = vertices[v3Index];

        glm::vec3 triangleCenter = (v1 + v2 + v3) / 3.0f;

        if (isTriangleInsideIco(triangleCenter, sphereCenter, sphereRadius)) {
            // Triangle is inside the sphere, mark it as invalid
            faceIndices[i] = faceIndices[i + 1] = faceIndices[i + 2] = -1;
        }
    }

    // Remove invalid triangles from the vectors
    newVertices.reserve(vertices.size());
    newNormals.reserve(normals.size());
    for (size_t i = 0; i < faceIndices.size(); ++i) {
        int index = faceIndices[i];
        if (index != -1) {
            newVertices.push_back(vertices[index]);
            newNormals.push_back(normals[index]);
        }
    }

    // Update the Torus data with the new vectors
    vertices = newVertices;
    normals = newNormals;
}

void Torus::removeTrianglesInsideSphere2(glm::vec3 sphereCenter, float sphereRadius, std::vector<glm::vec3> vertices,
    std::vector<glm::vec3> normals, std::vector<unsigned int> indices, std::vector<unsigned int> faceIndices) {
    std::vector<glm::vec3> newVertices;
    std::vector<glm::vec3> newNormals;
    std::vector<unsigned int> newIndices;

    for (size_t i = 0; i < faceIndices.size(); i += 3) {
        int v1Index = faceIndices[i];
        int v2Index = faceIndices[i + 1];
        int v3Index = faceIndices[i + 2];

        // Check if indices are within valid range
        if (v1Index < 0 || v1Index >= vertices.size() || v2Index < 0 || v2Index >= vertices.size() || v3Index < 0 ||
            v3Index >= vertices.size()) {
            continue; // Skip invalid indices
        }

        glm::vec3 v1 = vertices[v1Index];
        glm::vec3 v2 = vertices[v2Index];
        glm::vec3 v3 = vertices[v3Index];

        glm::vec3 triangleCenter = (v1 + v2 + v3) / 3.0f;

        if (!isTriangleInsideIco(triangleCenter, sphereCenter, sphereRadius)) {
            // Triangle is not inside the sphere, keep it
            newVertices.push_back(v1);
            newVertices.push_back(v2);
            newVertices.push_back(v3);

            newNormals.push_back(normals[v1Index]);
            newNormals.push_back(normals[v2Index]);
            newNormals.push_back(normals[v3Index]);

            // Adjust indices to match newVertices
            int newIndexBase = static_cast<int>(newVertices.size()) - 3;
            newIndices.push_back(newIndexBase);
            newIndices.push_back(newIndexBase + 1);
            newIndices.push_back(newIndexBase + 2);

            // Update faceIndices to match newVertices
            faceIndices[i] = newIndexBase;
            faceIndices[i + 1] = newIndexBase + 1;
            faceIndices[i + 2] = newIndexBase + 2;
        }
    }

    // Update the Torus data with the new vectors
    vertices = newVertices;
    normals = newNormals;
    indices = newIndices;
}

const bool Torus::isTriangleInsideIco(
    glm::vec3 triangleCenter, const glm::vec3& atomPosition, float atomRadius) {
    //glm::vec3 triangleCenter = (v1 + v2 + v3) / 3.0f;

    float distance = glm::length(triangleCenter - atomPosition);

    return distance < atomRadius;
}

const bool Torus::isPointInsideIco(glm::vec3 atomPosition, float atomRadius, glm::vec3 point){
    float distance = glm::distance(atomPosition, point);

    return distance < atomRadius;
}


void Torus::removeTriangles(const std::vector<unsigned int> trianglesToRemove, std::vector<glm::vec3>& vertices,
    std::vector<unsigned int>& faceIndices) {
    std::vector<glm::vec3> newVertices;
    std::vector<unsigned int> newFaceIndices;

    for (size_t i = 0; i < faceIndices.size(); i += 3) {
        unsigned int vertexIndex1 = faceIndices[i];
        unsigned int vertexIndex2 = faceIndices[i + 1];
        unsigned int vertexIndex3 = faceIndices[i + 2];

        if (std::find(trianglesToRemove.begin(), trianglesToRemove.end(), vertexIndex1) == trianglesToRemove.end() &&
            std::find(trianglesToRemove.begin(), trianglesToRemove.end(), vertexIndex2) == trianglesToRemove.end() &&
            std::find(trianglesToRemove.begin(), trianglesToRemove.end(), vertexIndex3) == trianglesToRemove.end()) {

            newVertices.push_back(vertices[vertexIndex1]);
            newVertices.push_back(vertices[vertexIndex2]);
            newVertices.push_back(vertices[vertexIndex3]);

            newFaceIndices.push_back(newFaceIndices.size());
            newFaceIndices.push_back(newFaceIndices.size() + 1);
            newFaceIndices.push_back(newFaceIndices.size() + 2);
        }
    }

    vertices = newVertices;
    faceIndices = newFaceIndices;
}

const float Torus::getDistance(glm::vec3 atomPosition1, glm::vec3 atomPosition2) {
    float distance = glm::distance(atomPosition1, atomPosition2);

    return distance; //d_ij
}

const glm::vec3 Torus::getTorusAxisUnitVec(glm::vec3 atomPosition1, glm::vec3 atomPosition2) {
    glm::vec3 cuttingPlaneNormal = glm::normalize(atomPosition2 - atomPosition1);
    glm::vec3 referenceNormal = glm::vec3(0.0f, 0.0f, 1.0f);
    glm::vec3 axisUnitVec = glm::normalize(glm::cross(referenceNormal, cuttingPlaneNormal));

    return axisUnitVec; //u_ij
}

const float Torus::getRotationAngle(glm::vec3 atomPosition1, glm::vec3 atomPosition2) {
    glm::vec3 cuttingPlaneNormal = glm::normalize(atomPosition2 - atomPosition1);
    glm::vec3 referenceNormal = glm::vec3(0.0f, 0.0f, 1.0f);
    float rotationAngle = glm::acos(glm::dot(referenceNormal, cuttingPlaneNormal));

    return rotationAngle;
}

const glm::vec3 Torus::getTorusCenter(
    glm::vec3 atomPosition1, glm::vec3 atomPosition2, float atomRadius1, float atomRadius2, float probeRadius) {
    float distance = getDistance(atomPosition1, atomPosition2);
    glm::vec3 torusCenter =
        0.5f * (atomPosition1 + atomPosition2) +
        0.5f * (atomPosition2 - atomPosition1) *
            ((std::powf((atomRadius1 + probeRadius), 2)) - (std::powf((atomRadius2 + probeRadius), 2))) /
            powf(distance, 2);

    return torusCenter; //t_ij
}

const float Torus::getTorusRadius(
    glm::vec3 atomPos1, glm::vec3 atomPos2, float atomRadius1, float atomRadius2, float probeRadius) {
    glm::vec3 vec1_2 = atomPos2 - atomPos1;
    glm::vec3 normal = glm::normalize(glm::cross(vec1_2, glm::vec3(0, 0, 1)));
    glm::vec3 intersecCenter = atomPos1 + 0.5f * vec1_2;

    float averageRadius = (atomRadius1 + atomRadius2) * 0.5f;
    float radius = glm::abs(glm::length(intersecCenter - atomPos1) - atomRadius1) + probeRadius * 2;

    /*float part1 = std::powf(atomRadius1 + atomRadius2 + 2 * probeRadius, 2);
    float part2 = glm::abs(glm::length(atomPos2 - atomPos1));
    float part2sqr = std::powf(glm::abs(glm::length(atomPos2 - atomPos1)), 2);
    float part12 = 1 / powf(part1 - part2sqr, 2);
    float part3 = std::powf(atomRadius1 - atomRadius2, 2);
    float part23 = 1 / std::powf(part2sqr - part3, 2);

    float radius = 0.5f * part12 * part23 / part2;*/

    return radius;
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

const glm::vec3 Torus::getTorusBasePointUnitVector(
    glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3) {
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

    glm::vec3 basePoint = torusCenter_12 + torusBasePointUnitVector *
                                               (glm::dot(torusAxisUnitVector_13, (torusCenter_13 - torusCenter_12))) *
                                               angleFactor;

    return basePoint; //b_ijk
}

const float Torus::getProbeHeight(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3,
    float atomRadius1, float atomRadius2, float atomRadius3, float probeRadius) {
    glm::vec3 basePoint =
        getBasePoint(atomPosition1, atomPosition2, atomPosition3, atomRadius1, atomRadius2, atomRadius3, probeRadius);

    float radPow = powf(atomRadius1 + probeRadius, 2);
    float pointPow = powf(glm::distance(basePoint, atomPosition1), 2);
    float probeHeight = sqrtf(radPow - pointPow);

    return probeHeight; //h_ijk
};

const glm::vec3 Torus::getProbePosition(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3,
    float atomRadius1, float atomRadius2, float atomRadius3, float probeRadius) {
    //b_ijk, h_ijk, u_ijk
    glm::vec3 basePoint =
        getBasePoint(atomPosition1, atomPosition2, atomPosition3, atomRadius1, atomRadius2, atomRadius3, probeRadius);
    float probeHeight =
        getProbeHeight(atomPosition1, atomPosition2, atomPosition3, atomRadius1, atomRadius2, atomRadius3, probeRadius);
    glm::vec3 basePlaneNormal = getBasePlaneNormal(atomPosition1, atomPosition2, atomPosition3);
    //TODO
    glm::vec3 probePosition = basePoint;

    return probePosition; //p_ijk
}

const glm::vec3 Torus::getProbeVertex(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3,
    float atomRadius1, float atomRadius2, float atomRadius3, float probeRadius) {
    glm::vec3 probePosition = getProbePosition(
        atomPosition1, atomPosition2, atomPosition3, atomRadius1, atomRadius2, atomRadius3, probeRadius);
    glm::vec3 dividend = atomRadius1 * probePosition + probeRadius * atomPosition1;

    glm::vec3 vertex = dividend / (atomRadius1 + probeRadius);

    return vertex; //v_pi
}

const glm::vec3 Torus::getContactCircleCenter(
    glm::vec3 atomPosition1, glm::vec3 atomPosition2, float atomRadius1, float atomRadius2, float probeRadius) {
    glm::vec3 torusCenter = getTorusCenter(atomPosition1, atomPosition2, atomRadius1, atomRadius2, probeRadius);
    glm::vec3 dividend = atomRadius1 * torusCenter + probeRadius * atomPosition1;

    glm::vec3 contactCircleCenter = dividend / (atomRadius1 + probeRadius);

    return contactCircleCenter; //c_ij
}

const float Torus::getContactCircleRadius(
    glm::vec3 atomPosition1, glm::vec3 atomPosition2, float atomRadius1, float atomRadius2, float probeRadius) {
    float distance = getDistance(atomPosition1, atomPosition2);
    float torusRadius = getTorusRadius(atomPosition1, atomPosition2, atomRadius1, atomRadius2, probeRadius);

    float contactCircleRadius = (torusRadius * atomRadius1) / (atomRadius1 + probeRadius);

    return contactCircleRadius; //r_c
}

const float Torus::getContactCircleDisplacement(
    glm::vec3 atomPosition1, glm::vec3 atomPosition2, float atomRadius1, float atomRadius2, float probeRadius) {
    glm::vec3 torusAxisUnitVec = getTorusAxisUnitVec(atomPosition1, atomPosition2);
    glm::vec3 contactCircleCenter =
        getContactCircleCenter(atomPosition1, atomPosition2, atomRadius1, atomRadius2, probeRadius);

    float contactCirlceDisplacement = glm::dot(torusAxisUnitVec, (contactCircleCenter - atomPosition1));

    return contactCirlceDisplacement;
}


const glm::vec3 Torus::rotateVector(float cosTheta, float sinTheta, glm::vec3 rotateVector, glm::vec3 rotateAxis) {
    glm::vec3 rotatedVec = cosTheta * rotateVector + sinTheta + glm::cross(rotateAxis, rotateVector) +
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
