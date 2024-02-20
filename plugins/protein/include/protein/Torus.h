#ifndef TORUS_H
#define TORUS_H

#pragma once

#include <glm/glm.hpp>
#include <vector>
#include <map>

class Torus {
public:
    explicit Torus(glm::vec3 center, float radius, int numSegments, int numRings);
    //~Torus() = default;

    void generateTorus(float probeRadius, glm::vec3 axisUnitVec, float rotationAngle, int offset, glm::vec3 atomPos1, glm::vec3 atomPos2, float atomRadius1, float atomRadius2);

    void setCenter(const glm::vec3& newCenter) {
        center = newCenter;
    }

    void setRadius(const float newRadius) {
        radius = newRadius;
    }

    [[nodiscard]] unsigned int getVertexSize() const {
        return (unsigned int)vertices.size() * sizeof(float);
    }

      [[nodiscard]] unsigned int getVertexCount() const {
        return (unsigned int)vertices.size();
    }

    [[nodiscard]] unsigned int getNormalSize() const {
        return (unsigned int)normals.size() * sizeof(float);
    }

    [[nodiscard]] unsigned int getNormalCount() const {
        return (unsigned int)normals.size();
    }

    [[nodiscard]] unsigned int getIndexSize() const {
        return (unsigned int)indices.size() * sizeof(unsigned int);
    }

    [[nodiscard]] unsigned int getIndexCount() const {
        return (unsigned int)indices.size();
    }

    [[nodiscard]] unsigned int getLineIndicesSize() const {
        return (unsigned int)lineIndices.size() * sizeof(unsigned int);
    }

    [[nodiscard]] unsigned int getLineIndicesCount() const {
        return (unsigned int)lineIndices.size();
    }

    [[nodiscard]] unsigned int getFaceIndicesCount() const {
        return (unsigned int)faceIndices.size();
    }

    [[nodiscard]] const glm::vec3* getVertices() const {
        if (vertices.empty()) {
            return nullptr;
        } else {
            //return &temp[0];
            return &vertices[0];
        }
    }

    [[nodiscard]] const bool getVertices_bool(int i) const {
        if (vertices_vs_check.empty()) {
            return false;
        } else {
            return &vertices_vs_check[i].second;
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

    const unsigned int* getLineIndices() const {
        return lineIndices.data();
    }

    const unsigned int* getFaceIndices() const {
        return faceIndices.data();
    }

    static const glm::vec3 getContactPoint(const glm::vec3& probePos, glm::vec3 atomPos, float atomRad);
    static const glm::vec3 getProbePos(const glm::vec3& torusCenter, float torusRad, float probeRad,
                                       const glm::vec3& torusPoint, float rotationAngle, glm::vec3 axisUnitVec);
    static const std::pair<glm::vec3, float> getVisibilitySphere(glm::vec3 contactPoint, glm::vec3 probePos, glm::vec3 atomPos1, glm::vec3 atomPos2, glm::vec3 torusCenter);


    static const float getDistance(glm::vec3 atomPosition1, glm::vec3 atomPosition2);
    static const glm::vec3 getTorusAxisUnitVec(glm::vec3 atomPosition1, glm::vec3 atomPosition2);
    static const glm::vec3 getTorusCenter(glm::vec3 atomPosition1, glm::vec3 atomPosition2, float atomRadius1, float atomRadius2, float probeRadius);
    static const float getRotationAngle(glm::vec3 atomPosition1, glm::vec3 atomPosition2);
    static const float getTorusRadius(glm::vec3 atomPos1, glm::vec3 atomPos2, float atomRadius1, float atomRadius2, float probeRadius);

    static const float getBaseTriangleAngle(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3);
    static const glm::vec3 getBasePlaneNormal(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3);
    static const glm::vec3 getTorusBasePointUnitVector(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3);
    static const glm::vec3 getBasePoint(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3,
                                        float atomRadius1, float atomRadius2, float atomRadius3, float probeRadius);

    static const float getProbeHeight(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3,
                                      float atomRadius1, float atomRadius2, float atomRadius3, float probeRadius);
    static const glm::vec3 getProbePosition(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3,
                                      float atomRadius1, float atomRadius2, float atomRadius3, float probeRadius);
    static const glm::vec3 getProbeVertex(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3,
                                          float atomRadius1, float atomRadius2, float atomRadius3, float probeRadius);

    static const glm::vec3 getContactCircleCenter(glm::vec3 atomPosition1, glm::vec3 atomPosition2, float atomRadius1, float atomRadius2, float probeRadius);
    static const float getContactCircleRadius(glm::vec3 atomPosition1, glm::vec3 atomPosition2, float atomRadius1, float atomRadius2, float probeRadius);
    static const float getContactCircleDisplacement(glm::vec3 atomPosition1, glm::vec3 atomPosition2, float atomRadius1, float atomRadius2, float probeRadius);
    static const glm::vec3 getConcaveArcPlaneNormal(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3,
                                                    float atomRadius1, float atomRadius2, float atomRadius3, float probeRadius);
    static const float getConcaveTriangleAngle(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3,
                                               float atomRadius1, float atomRadius2, float atomRadius3, float probeRadius);
    static const float getConvexFaceAngle(float concaveTriangleAngle);
    static const float getSaddleWrapAngle(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3, glm::vec3 atomPosition4,
                                          float atomRadius1, float atomRadius2, float atomRadius3, float atomRadius4, float probeRadius);
    static const float getSaddleWithAngle(float contactCircleRadius, float contactCircleDisplacement);
    std::vector<glm::vec3> vertices;
    std::vector<std::pair<glm::vec3, bool>> vertices_vs_check;
    std::vector<glm::vec3> normals;
    std::vector<unsigned int> indices;
    std::vector<unsigned int> faceIndices;
    std::map<int, int> vertexFaceIndices;
    std::map<unsigned int, std::vector<glm::vec3>> vertexFaceIndex;
    std::vector<glm::vec3> edgeVertices;

private:
    glm::vec3 center;
    float radius;
    float probeRadius;
    int numSegments;
    int numRings;

    
    std::vector<unsigned int> lineIndices;

    static const glm::vec3 rotateVector(float cosTheta, float sinTheta, glm::vec3 rotateVector, glm::vec3 rotateAxis);

    void addVertex(float x, float y, float z);
    void addNormal(glm::vec3 normal);
    void addIndices(unsigned int i1, unsigned int i2, unsigned int i3);
};

#endif
