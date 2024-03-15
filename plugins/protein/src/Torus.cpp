
#ifdef _WIN32
#include <windows.h> // include windows.h to avoid thousands of compile errors even though this class is not depending on Windows
#endif

#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtx/rotate_vector.hpp"
#include "glm/gtx/string_cast.hpp"
#include <glm/gtx/norm.hpp>
#include "protein/Torus.h"
#include <cmath>
#include <glm/glm.hpp>
#include <iostream>


#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884L
#endif


Torus::Torus(glm::vec3 center, float radius, int numSegments, int numRings)
        : center(center)
        , radius(radius) 
        , numSegments(numSegments)
        , numRings(numRings) {}

void Torus::generateTorus(float probeRadius, glm::vec3 axisUnitVec, float rotationAngle, int offset, glm::vec3 atomPos1, glm::vec3 atomPos2, float atomRadius1, float atomRadius2) {
    float pi = M_PI;

    // clear memory of prev arrays
    std::vector<glm::vec3>().swap(vertices);
    std::vector<glm::vec3>().swap(normals);
    std::vector<unsigned int>().swap(indices);
    std::vector<unsigned int>().swap(lineIndices);
    std::vector<unsigned int>().swap(faceIndices);

    //dummy point for visibility sphere
    float theta = 2 * pi * static_cast<float>(0) / static_cast<float>(numSegments);
    float phi = 2 * pi * static_cast<float>(0) / static_cast<float>(numRings);

    float x = (radius + probeRadius * cos(phi)) * cos(theta);
    float y = (radius + probeRadius * cos(phi)) * sin(theta);
    float z = probeRadius * sin(phi);

    //calculate visibility sphere
    glm::vec3 probePos = getProbePos(center, radius, probeRadius, glm::vec3(x, y, z), rotationAngle, axisUnitVec);
    glm::vec3 contactPoint(0, 0, 0);
    std::pair<glm::vec3, float> visibilitySphere;
    if (atomRadius1 <= atomRadius2) {
        contactPoint = getContactPoint(probePos, atomPos1, atomRadius1);
        visibilitySphere = getVisibilitySphere(contactPoint, probePos, atomPos1, atomPos2, center);
        std::cout << "used atom1 for contact point" << std::endl;
    } else {
        contactPoint = getContactPoint(probePos, atomPos2, atomRadius2);
        visibilitySphere = getVisibilitySphere(contactPoint, probePos, atomPos2, atomPos1, center);
        std::cout << "used atom2 for contact point" << std::endl;
    }


    for (int i = 0; i <= numRings; ++i) {
        for (int j = 0; j <= numSegments; ++j) {
            //calculate (x, y, z) still need to rotate
            float theta = 2 * pi * static_cast<float>(j) / static_cast<float>(numSegments);
            float phi = 2 * pi * static_cast<float>(i) / static_cast<float>(numRings);

            float x = (radius + probeRadius * cos(phi)) * cos(theta);
            float y = (radius + probeRadius * cos(phi)) * sin(theta);
            float z = probeRadius * sin(phi);

            glm::vec3 vertex(x, y, z);
            addVertex(x, y, z);

            //tried to interpolate the normals to get a smooth transition from stitching triangles to torus
            //not working as expected
            /*
            glm::vec3 probePos = getProbePos(center, radius, probeRadius, vertex, rotationAngle, axisUnitVec);
            glm::vec3 contactPoint(0, 0, 0);
            glm::vec3 contactNormal(0, 0, 0);
            float dist = 0.0f;
            float dist1 = glm::distance(vertex, atomPos1);
            float dist2 = glm::distance(vertex, atomPos2);
            if (dist1 < dist2) {
                contactPoint = getContactPoint(probePos, atomPos1, atomRadius1);
                contactNormal = glm::normalize(contactPoint - atomPos1);
                dist = dist1;
            } else {
                contactPoint = getContactPoint(probePos, atomPos2, atomRadius2);
                contactNormal = glm::normalize(contactPoint - atomPos2);
                dist = dist2;
            }
            float interpolationDist = radius / 10;

            float interpolationFactor = glm::clamp(dist / interpolationDist, 0.0f, 1.0f);

            glm::vec3 normal = glm::normalize(glm::mix(contactNormal, glm::normalize(vertex), interpolationFactor));
            */

            indices.push_back(vertices.size() - 1 + offset);

            glm::vec3 normal(glm::normalize(vertex));
            addNormal(normal);
        }
    }

    //transform torus
    glm::mat4 translationMatrix = glm::translate(glm::mat4(1.0f), center);
    glm::mat4 rotationMatrix = glm::rotate(glm::mat4(1.0f), rotationAngle, axisUnitVec);
    glm::mat4 transformationMatrix = translationMatrix * rotationMatrix;
    for (int i = 0; i < vertices.size(); ++i) {
        glm::vec4 transformedVertex = transformationMatrix * glm::vec4(vertices[i], 1.0f);
        vertices[i] = glm::vec3(transformedVertex);
    }

    //tranform normals
    glm::mat4 inverseTranspose = glm::inverse(glm::transpose(transformationMatrix));
    for (int i = 0; i < normals.size(); ++i) {
        glm::vec4 tranformedNormal = inverseTranspose * glm::vec4(normals[i], 0.0f);
        normals[i] = glm::normalize(glm::vec3(tranformedNormal));
    }

    //mark vertices inside the visibility sphere
    for (int i = 0; i < vertices.size(); ++i) {
        float distanceToVSCenter = glm::distance(visibilitySphere.first, vertices[i]);

        if (distanceToVSCenter <= visibilitySphere.second) {
            vertices_vs_check.push_back(std::pair(vertices[i], true));
        } else {
            vertices_vs_check.push_back(std::pair(vertices[i], false));
        }
    }

    //std::cout << "vertices_size: " << vertices.size() << std::endl;

    //add face indices only if all three points of a triangle are inside the visibility sphere
    for (int i = 0; i < numRings; ++i) {
        for (int j = 0; j < numSegments; ++j) {
            
            //point j on ring i
            int v1 = i * (numSegments + 1) + j + offset;
            //next point (j + 1) on ring i
            int v2 = i * (numSegments + 1) + (j + 1) % numSegments + offset;
            //point j on ring i + 1
            int v3 = (i + 1) * (numSegments + 1) + j + offset;

            if (vertices_vs_check[v1 - offset].second && vertices_vs_check[v2 - offset].second && vertices_vs_check[v3 - offset].second) {
                    faceIndices.push_back(v3);
                    faceIndices.push_back(v2);
                    faceIndices.push_back(v1);
            }

            std::vector<glm::vec3> vertices_v1 = {vertices[(v1 - offset) % vertices.size()],
                                                    vertices[(v1 - offset + 1) % vertices.size()],
                                                    vertices[(v1 - offset + 2) % vertices.size()]};
            vertexFaceIndex[v1] = vertices_v1;

            std::vector<glm::vec3> vertices_v2 = {vertices[(v2 - offset) % vertices.size()],
                                                    vertices[(v2 - offset + 1) % vertices.size()],
                                                    vertices[(v2 - offset + 2) % vertices.size()]};
            vertexFaceIndex[v2] = vertices_v2;

            std::vector<glm::vec3> vertices_v3 = {vertices[(v3 - offset) % vertices.size()],
                                                    vertices[(v3 - offset + 1) % vertices.size()],
                                                    vertices[(v3 - offset + 2) % vertices.size()]};
            vertexFaceIndex[v3] = vertices_v3;

            //point j on ring i + 1
            int v4 = (i + 1) * (numSegments + 1) + j + offset;
            //next point on ring i
            int v5 = i * (numSegments + 1) + (j + 1) % numSegments + offset;
            //next point (j + 1) on ring i + 1
            int v6 = (i + 1) * (numSegments + 1) + (j + 1) % numSegments + offset;

            if (vertices_vs_check[v4 - offset].second && vertices_vs_check[v5 - offset].second && vertices_vs_check[v6 - offset].second) {
                    faceIndices.push_back(v6);
                    faceIndices.push_back(v5);
                    faceIndices.push_back(v4);
            }

            std::vector<glm::vec3> vertices_v4 = {vertices[(v4 - offset) % vertices.size()],
                                                    vertices[(v4 - offset + 1) % vertices.size()],
                                                    vertices[(v4 - offset + 2) % vertices.size()]};
            vertexFaceIndex[v4] = vertices_v4;

            std::vector<glm::vec3> vertices_v5 = {vertices[(v5 - offset) % vertices.size()],
                                                    vertices[(v5 - offset + 1) % vertices.size()],
                                                    vertices[(v5 - offset + 2) % vertices.size()]};
            vertexFaceIndex[v5] = vertices_v5;

            std::vector<glm::vec3> vertices_v6 = {vertices[(v6 - offset) % vertices.size()],
                                                    vertices[(v6 - offset + 1) % vertices.size()],
                                                    vertices[(v6 - offset + 2) % vertices.size()]};
            vertexFaceIndex[v6] = vertices_v6;


            //save edge vertices of torus
            //a vertex is edge vertex if it is inside the visibility sphere and one of its neighbours isn't
            //no offset because we only search in local vertices
            int currentIndex = i * (numSegments + 1) + j; 
            if (vertices_vs_check[currentIndex].second) {     
                std::vector<int> neighbors = {
                    ((i - 1 + numRings) % numRings) * (numSegments + 1) + j, 
                    ((i + 1) % numRings) * (numSegments + 1) + j,            
                    i * (numSegments + 1) + ((j - 1 + numSegments) % numSegments), 
                    i * (numSegments + 1) + ((j + 1) % numSegments) 
                };
                for (auto& neighborIndex : neighbors) {
                    if (!vertices_vs_check[neighborIndex].second) { 
                        edgeVertices.push_back(vertices[currentIndex]); 
                        break;                                         
                    }
                }
            }

        }
    }
    /*for (auto elem : edgeVertices) {
        std::cout << elem.x << "/" << elem.y << "/" << elem.z << std::endl;
    }*/
    //std::cout << "FaceIndices: " << getFaceIndicesCount() << std::endl;
    /*std::cout << "EDGE VERTICES OF TORUS:" << std::endl;
    for (auto elem : edgeVertices) {
        std::cout << elem.x << "/" << elem.y << "/" << elem.z << std::endl;
    }
    std::cout << "END OF EDGE VERTICES" << std::endl;*/
    //std::cout << "edgeVertices.size() (torus) = " << edgeVertices.size() << std::endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////
///// VISIBILITY SPHERE
////////////////////////////////////////////////////////////////////////////////////////////////
///// Interactive Visualization of Molecular Surface Dynamics (M. Krone, K. Bidmon, T. Ertl) (Fig. 7 Visibility Sphere)
////////////////////////////////////////////////////////////////////////////////////////////////
const glm::vec3 Torus::getProbePos(const glm::vec3& torusCenter, float torusRad, float probeRad, const glm::vec3& torusPoint, float rotationAngle, glm::vec3 axisUnitVec) {
    //probe center on base plane, need to transform
    glm::vec3 tubeCenter(0.0f, torusRad, 0.0f);

    //same as torus transformationMatrix
    glm::mat4 translationMatrix = glm::translate(glm::mat4(1.0f), torusCenter);
    glm::mat4 rotationMatrix = glm::rotate(glm::mat4(1.0f), rotationAngle, axisUnitVec);
    glm::mat4 transformationMatrix = translationMatrix * rotationMatrix;

    //tranform probe center
    glm::vec4 transformedTubeCenter = transformationMatrix * glm::vec4(tubeCenter, 1.0f);
    glm::vec3 probePos = glm::vec3(transformedTubeCenter);
    //std::cout << "probePos: " << probePos.x << "/" << probePos.y << "/" << probePos.z << std::endl;
    
    return probePos;
}

const glm::vec3 Torus::getContactPoint(const glm::vec3& probePos, glm::vec3 atomPos, float atomRad) {
    //x = ((p - a_i) / (|p - a_i|)) * r_i
    glm::vec3 temp = probePos - atomPos;
    glm::vec3 tempNorm = glm::normalize(temp);
    //std::cout << "probePos - atomPos = " << temp.x << "/" << temp.y << "/" << temp.z << std::endl;
    float tempLength = glm::distance(probePos, atomPos);
    //std::cout << "dist(probePos, atomPos) = " << tempLength << std::endl;

    //glm::vec3 contactPoint = (temp * atomRad) / tempLength;
    glm::vec3 contactPoint = atomPos + tempNorm * atomRad;
    //std::cout << "contactPoint: " << contactPoint.x << "/" << contactPoint.y << "/" << contactPoint.z << std::endl;

    return contactPoint;
}


const std::pair<glm::vec3, float> Torus::getVisibilitySphere(glm::vec3 contactPoint, glm::vec3 probePos, glm::vec3 atomPos1, glm::vec3 atomPos2, glm::vec3 torusCenter) {

    //p = probePos (vec3)
    glm::vec3 p = probePos;
    //std::cout << "p = " << p.x << "/" << p.y << "/" << p.z << std::endl;

    //a_i = atomPos1 (vec3)
    glm::vec3 a_i = atomPos1;
    //std::cout << "atomPos1 = " << a_i.x << "/" << a_i.y << "/" << a_i.z << std::endl;

    //a_j = atomPos2 (vec3)
    glm::vec3 a_j = atomPos2;
    //std::cout << "atomPos2 = " << a_j.x << "/" << a_j.y << "/" << a_j.z << std::endl;

    //x = contactPoint (vec3)
    glm::vec3 x = contactPoint;
    //std::cout << "x = " << x.x << "/" << x.y << "/" << x.z << std::endl;

    //t_ij = torusCenter (vec3)
    glm::vec3 t_ij = torusCenter;
    //std::cout << "torusCenter = " << t_ij.x << "/" << t_ij.y << "/" << t_ij.z << std::endl;

    //c' = (|p - a_i| / (|p - a_j| + |p - a_i|)) * (a_j - a_i) (vec3)
    //=> added (+ t_ij)
    glm::vec3 c_ = ((glm::distance(p, a_i) / (glm::distance(p, a_j) + glm::distance(p, a_i))) * (a_j - a_i)) + t_ij;
    //std::cout << "c' = " << c_.x << "/" << c_.y << "/" << c_.z << std::endl;

    //c = (c' + a_i) - t_ij (vec3)
    glm::vec3 c = (c_ + a_i) - t_ij;

    //r_vs = |x - c'| (float)
    //more accurate with |x - c|
    float r_vs = glm::distance(x, c /*c_*/) /* - 0.1375f*/;

    //std::cout << "vsCenter = " << c.x << "/" << c.y << "/" << c.z << std::endl;
    //std::cout << "vsRad = " << r_vs << std::endl;

    std::pair<glm::vec3, float> visibilitySphere(c, r_vs);

    return visibilitySphere;
}
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////
///// BASIC TORUS CALCULATIONS
////////////////////////////////////////////////////////////////////////////////////////////////
///// Analytical Molecular Surface Calculation (M. L Connolly) (Surface definition equarions)
////////////////////////////////////////////////////////////////////////////////////////////////
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

    /*
    std::cout << "TorusRadiusCalc:" << std::endl;
    std::cout << "AtomRad1: " << atomRadius1 << std::endl;
    std::cout << "AtomRad2: " << atomRadius2 << std::endl;
    std::cout << "AtomPos1: " << atomPos1.x << "/" << atomPos1.y << "/" << atomPos1.z << std::endl;
    std::cout << "AtomPos2: " << atomPos2.x << "/" << atomPos2.y << "/" << atomPos2.z << std::endl;
    std::cout << "ProbeRad: " << probeRadius << std::endl;
    */

    //(r_i + r_j + 2r_p)^2
    float part1 = std::powf(atomRadius1 + atomRadius2 + (2 * probeRadius), 2);
    //std::cout << "(r_i + r_j + 2r_p)^2 = " << part1 << std::endl;

    //d_ij
    float part2 = glm::abs(glm::length(atomPos2 - atomPos1));
    //std::cout << "d_ij = " << part2 << std::endl;

    //d_ij^2
    float part2sqr = std::powf(/*glm::abs(glm::length(atomPos2 - atomPos1))*/ part2 , 2);
    //std::cout << "d_ij^2 = " << part2sqr << std::endl;

    //[(r_i + r_j + 2r_p)^2 - d_ij^2]^(1/2)
    float part12 = std::sqrtf(part1 - part2sqr);
    //std::cout << "[(r_i + r_j + 2r_p)^2 - d_ij^2]^(1/2) = " << part12 << std::endl;

    //(r_i - r_j)^2
    float part3 = std::powf(std::abs(atomRadius1 - atomRadius2), 2);
    //std::cout << "(r_i - r_j)^2 = " << part3 << std::endl;

    //[d_ij^2 - (r_i - r_j)^2]^(1/2)
    //float part23 = 1 / std::powf(part2sqr - part3, 2);
    float part23 = std::sqrtf(part2sqr - part3);
    //std::cout << "[d_ij^2 - (r_i - r_j)^2]^(1/2) = " << part23 << std::endl;

    float radius = 0.5f * part12 * part23 / part2;
    //std::cout << "radius: " << radius << std::endl;

    return radius;
}


const float Torus::getBaseTriangleAngle(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3) {
    glm::vec3 axisUnitVector1 = getTorusAxisUnitVec(atomPosition1, atomPosition2);
    glm::vec3 axisUnitVector2 = getTorusAxisUnitVec(atomPosition1, atomPosition3);

    float dotProduct = glm::dot(axisUnitVector1, axisUnitVector2);

    float baseTriangleAngle_rad = std::acos(dotProduct);
    // rad or deg?
    return baseTriangleAngle_rad; //w_ijk
}

const glm::vec3 Torus::getBasePlaneNormal(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3) {
    float angle = getBaseTriangleAngle(atomPosition1, atomPosition2, atomPosition3);
    glm::vec3 axisUnitVector1 = getTorusAxisUnitVec(atomPosition1, atomPosition2);
    glm::vec3 axisUnitVector2 = getTorusAxisUnitVec(atomPosition1, atomPosition3);

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
    //t_ij
    glm::vec3 torusCenter_12 = getTorusCenter(atomPosition1, atomPosition2, atomRadius1, atomRadius2, probeRadius);

    //t_ik
    glm::vec3 torusCenter_13 = getTorusCenter(atomPosition1, atomPosition3, atomRadius1, atomRadius3, probeRadius);

    //u_tb
    glm::vec3 torusBasePointUnitVector = getTorusBasePointUnitVector(atomPosition1, atomPosition2, atomPosition3);

    //u_ik
    glm::vec3 torusAxisUnitVector_13 = getTorusAxisUnitVec(atomPosition1, atomPosition3);

    //w_ijk
    float baseTriangleAngle = getBaseTriangleAngle(atomPosition1, atomPosition2, atomPosition3);

    //(sin(w_ijk))^-1
    float angleFactor = 1 / sin(baseTriangleAngle);


    //b_ijk = t_ij + u_tb[u_ik • (t_ik - t_ij)] x (sin(w_ijk))^-1
    glm::vec3 basePoint = torusCenter_12 + torusBasePointUnitVector *
                                           (glm::dot(torusAxisUnitVector_13, (torusCenter_13 - torusCenter_12))) *
                                           angleFactor;

    return basePoint; //b_ijk
}

const float Torus::getProbeHeight(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3,
                                  float atomRadius1, float atomRadius2, float atomRadius3, float probeRadius) {
    //b_ijk
    glm::vec3 basePoint = getBasePoint(atomPosition1, atomPosition2, atomPosition3, atomRadius1, atomRadius2, atomRadius3, probeRadius);

    //(r_i + r_p)^2
    float radPow = powf(atomRadius1 + probeRadius, 2);

    //|b_ijk - a_i|^2
    float pointPow = powf(glm::distance(basePoint, atomPosition1), 2);

    //h_ijk = [(r_i + r_p)^2 - |b_ijk - a_i|^2]^(1/2)
    float probeHeight = sqrtf(radPow - pointPow);

    return probeHeight; //h_ijk
};

const glm::vec3 Torus::getProbePosition(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3,
                                        float atomRadius1, float atomRadius2, float atomRadius3, float probeRadius) {

    /*
    glm::vec3 v1 = atomPosition2 - atomPosition1;
    glm::vec3 v2 = atomPosition3 - atomPosition1;
    glm::vec3 normal = glm::normalize(glm::cross(v1, v2));

    glm::vec3 triangleCenter = (atomPosition1 + atomPosition2 + atomPosition3) / 3.0f;

    float radiusSum = (atomRadius1 + atomRadius2 + atomRadius3) / 3.0f;
    float distanceToPlane = probeRadius + radiusSum;

    glm::vec3 probeCenter = triangleCenter + normal * distanceToPlane;

    bool isCorrect = glm::distance(probeCenter, atomPosition1) == probeRadius + atomRadius1 &&
                     glm::distance(probeCenter, atomPosition2) == probeRadius + atomRadius2 &&
                     glm::distance(probeCenter, atomPosition3) == probeRadius + atomRadius3;

    if (isCorrect) {
        return probeCenter;
    } else {
        std::cout << "WRONG PROBE CENTER" << std::endl;
        return probeCenter;
    }
    */

    
    //b_ijk, h_ijk, u_ijk
    glm::vec3 basePoint =
        getBasePoint(atomPosition1, atomPosition2, atomPosition3, atomRadius1, atomRadius2, atomRadius3, probeRadius);
    float probeHeight =
        getProbeHeight(atomPosition1, atomPosition2, atomPosition3, atomRadius1, atomRadius2, atomRadius3, probeRadius);
    glm::vec3 basePlaneNormal = getBasePlaneNormal(atomPosition1, atomPosition2, atomPosition3);
    //TODO
    glm::vec3 probePosition = basePoint - (probeHeight * basePlaneNormal);

    return probePosition; //p_ijk
    
}

const glm::vec3 Torus::getProbeVertex(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3,
                                      float atomRadius1, float atomRadius2, float atomRadius3, float probeRadius) {
    //p_ijk
    glm::vec3 probePosition = getProbePosition(atomPosition1, atomPosition2, atomPosition3, atomRadius1, atomRadius2, atomRadius3, probeRadius);

    //(r_i * p_ijk + r_p * a_i)
    glm::vec3 dividend = atomRadius1 * probePosition + probeRadius * atomPosition1;

    //v_pi = (r_i * p_ijk + r_p * a_i) / (r_i + r_p)
    glm::vec3 vertex = dividend / (atomRadius1 + probeRadius);

    return vertex; //v_pi
}

const glm::vec3 Torus::getContactCircleCenter(glm::vec3 atomPosition1, glm::vec3 atomPosition2, float atomRadius1, float atomRadius2, float probeRadius) {
    //t_ij
    glm::vec3 torusCenter = getTorusCenter(atomPosition1, atomPosition2, atomRadius1, atomRadius2, probeRadius);

    //(r_i * t_ij + r_p * a_i)
    glm::vec3 dividend = atomRadius1 * torusCenter + probeRadius * atomPosition1;


    //c_ij = (r_i * t_ij + r_p * a_i) / (r_i + r_p)
    glm::vec3 contactCircleCenter = dividend / (atomRadius1 + probeRadius);

    return contactCircleCenter; //c_ij
}

const float Torus::getContactCircleRadius(glm::vec3 atomPosition1, glm::vec3 atomPosition2, float atomRadius1, float atomRadius2, float probeRadius) {
    float distance = getDistance(atomPosition1, atomPosition2);
    float torusRadius = getTorusRadius(atomPosition1, atomPosition2, atomRadius1, atomRadius2, probeRadius);

    //r_c = r_ij * r_i / (r_i + r_p)
    float contactCircleRadius = (torusRadius * atomRadius1) / (atomRadius1 + probeRadius);

    return contactCircleRadius; //r_c
}

const float Torus::getContactCircleDisplacement(glm::vec3 atomPosition1, glm::vec3 atomPosition2, float atomRadius1, float atomRadius2, float probeRadius) {
    glm::vec3 torusAxisUnitVec = getTorusAxisUnitVec(atomPosition1, atomPosition2);
    glm::vec3 contactCircleCenter =getContactCircleCenter(atomPosition1, atomPosition2, atomRadius1, atomRadius2, probeRadius);

    //d_c = u_ij • (c_ij - a_i)
    float contactCirlceDisplacement = glm::dot(torusAxisUnitVec, (contactCircleCenter - atomPosition1));

    return contactCirlceDisplacement; //d_c
}

const glm::vec3 Torus::getConcaveArcPlaneNormal(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3,
                                                float atomRadius1, float atomRadius2, float atomRadius3, float probeRadius) {
    glm::vec3 probePosition = getProbePosition(atomPosition1, atomPosition2, atomPosition3, atomRadius1, atomRadius2, atomRadius3, probeRadius);
    glm::vec3 torusCenter_12 = getTorusCenter(atomPosition1, atomPosition2, atomRadius1, atomRadius2, probeRadius);
    glm::vec3 torusAxisUnitVec_12 = getTorusAxisUnitVec(atomPosition1, atomPosition2);
    float torusRadius_12 = getTorusRadius(atomPosition1, atomPosition2, atomRadius1, atomRadius2, probeRadius);

    //(p_ijk - t_ij) x u_ij
    glm::vec3 dividend = glm::cross((probePosition - torusCenter_12), torusAxisUnitVec_12);

    //n_ijk = ((p_ijk - t_ij) x u_ij) / r_ij
    glm::vec3 concaveArcPlaneNormal = dividend / torusRadius_12;

    return concaveArcPlaneNormal;
}

const float Torus::getConcaveTriangleAngle(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3,
                                    float atomRadius1, float atomRadius2, float atomRadius3, float probeRadius) {
    //atom1 := i; atom2 := j; atom3 := k
    glm::vec3 n_ijk = getConcaveArcPlaneNormal(atomPosition1, atomPosition2, atomPosition3, atomRadius1, atomRadius2, atomRadius3, probeRadius);
    glm::vec3 n_ikj = getConcaveArcPlaneNormal(atomPosition1, atomPosition3, atomPosition2, atomRadius1, atomRadius3, atomRadius2, probeRadius);

    //beta = arccos(n_ijk • n_ikj);
    float beta = acos(glm::dot(n_ijk, n_ikj));

    return beta;
}

const float Torus::getConvexFaceAngle(float concaveTriangleAngle) {
    float convexFaceAngle = M_PI - concaveTriangleAngle;

    return convexFaceAngle;
}

const float Torus::getSaddleWrapAngle(glm::vec3 atomPosition1, glm::vec3 atomPosition2, glm::vec3 atomPosition3, glm::vec3 atomPosition4,
                                      float atomRadius1, float atomRadius2, float atomRadius3, float atomRadius4, float probeRadius) {
    //atom1 := i; atom2 := j; atom3 := k; atom4 := l
    glm::vec3 n_ijk = getConcaveArcPlaneNormal(atomPosition1, atomPosition2, atomPosition3, atomRadius1, atomRadius2, atomRadius3, probeRadius);
    glm::vec3 n_ijl = getConcaveArcPlaneNormal(atomPosition1, atomPosition2, atomPosition4, atomRadius1, atomRadius2, atomRadius4, probeRadius);
    glm::vec3 u_ij = getTorusAxisUnitVec(atomPosition1, atomPosition2);

    //n_ijk x n_ijl • u_ij
    float test = glm::dot(glm::cross(n_ijk, n_ijl), u_ij);

    float phi;
    if (test >= 0) {
        phi = acos(glm::dot(n_ijk, n_ijl));
    } else {
        phi = acos(glm::dot(n_ijk, n_ijl)) + 2 * M_PI;
    }

    return phi;
}

const float Torus::getSaddleWithAngle(float contactCircleRadius, float contactCircleDisplacement) {
    float theta = atan(contactCircleDisplacement / contactCircleRadius);

    return theta;
}
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////
///// HELPER CALCULATIONS
////////////////////////////////////////////////////////////////////////////////////////////////
const glm::vec3 Torus::rotateVector(float cosTheta, float sinTheta, glm::vec3 rotateVector, glm::vec3 rotateAxis) {
    glm::vec3 rotatedVec = cosTheta * rotateVector + sinTheta + glm::cross(rotateAxis, rotateVector) +
                           (1 - cosTheta) * glm::dot(rotateAxis, rotateVector) * rotateAxis;

    return rotatedVec;
}
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////
///// ADD FUNCTIONS
////////////////////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
