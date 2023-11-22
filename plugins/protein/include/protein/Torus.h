#ifndef TORUS_H
#define TORUS_H

#pragma once

#include <glm/glm.hpp>
#include <vector>

class Torus {
public:
    explicit Torus(glm::vec3 center, float radius, int numSegments, int numRings);
    //~Torus() = default;

    void generateTorus(float probeRadius);

    void setCenter(const glm::vec3& newCenter) {
        center = newCenter;
    }

    void setRadius(const float newRadius) {
        radius = newRadius;
    }

    [[nodiscard]] unsigned int getVertexSize() const {
        return (unsigned int)vertices.size() * sizeof(float);
    }

    [[nodiscard]] unsigned int getNormalSize() const {
        return (unsigned int)normals.size() * sizeof(float);
    }

    [[nodiscard]] unsigned int getIndexSize() const {
        return (unsigned int)indices.size() * sizeof(unsigned int);
    }

    [[nodiscard]] unsigned int getIndexCount() const {
        return (unsigned int)indices.size();
    }

    [[nodiscard]] const glm::vec3* getVertices() const {
        if (vertices.empty()) {
            return nullptr;
        } else {
            return &vertices[0];
        }
    }

    const glm::vec3* getNormals() const {
        if (normals.empty()) {
            return nullptr;
        } else {
            return &normals[0];
        }
    }

    const unsigned int* getIndices() const {
        return indices.data();
    }

    static const float getDistance(glm::vec3 atomPosition1, glm::vec3 atomPosition2);
    static const glm::vec3 getTorusAxisUnitVec(glm::vec3 atomPosition1, glm::vec3 atomPosition2);
    static const glm::vec3 getTorusCenter(glm::vec3 atomPosition1, glm::vec3 atomPosition2, float atomRadius1, float atomRadius2, float probeRadius);
    static const float getTorusRadius(float atomRadius1, float atomRadius2, float probeRadius, float distance);

    static const float getBaseTriangleAngle(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3);
    static const glm::vec3 getBasePlaneNormal(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3);
    static const glm::vec3 getTorusBasePointUnitVector(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3);
    static const glm::vec3 getBasePoint(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3,
                                        float atomRadius1, float atomRadius2, float atomRadius3, float probeRadius);


    std::vector<int> faceIndices;

private:
    glm::vec3 center;
    float radius;
    float probeRadius;
    int numSegments;
    int numRings;

    std::vector<glm::vec3> vertices;
    std::vector<glm::vec3> normals;
    std::vector<unsigned int> indices;

    static const glm::vec3 rotateVector(float cosTheta, float sinTheta, glm::vec3 rotateVector, glm::vec3 rotateAxis);

    void addVertex(float x, float y, float z);
    void addNormal(float nx, float ny, float nz);
    void addIndices(unsigned int i1, unsigned int i2, unsigned int i3);
};

#endif
