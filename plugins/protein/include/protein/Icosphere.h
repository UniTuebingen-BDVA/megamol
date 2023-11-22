///////////////////////////////////////////////////////////////////////////////
// Icosphere.h
// ===========
// Polyhedron subdividing icosahedron (20 tris) by N-times iteration
// The icosphere with N=1 (default) has 80 triangles by subdividing a triangle
// of icosahedron into 4 triangles. If N=0, it is identical to icosahedron.
//
//  AUTHOR: Song Ho Ahn (song.ahn@gmail.com)
// CREATED: 2018-07-23
// UPDATED: 2019-12-28
///////////////////////////////////////////////////////////////////////////////

#ifndef GEOMETRY_ICOSPHERE_H
#define GEOMETRY_ICOSPHERE_H

#include <map>
#include <vector>
#include <glm/glm.hpp>

class Icosphere {
public:
    // ctor/dtor
    explicit Icosphere(float radius = 1.0f, int subdivision = 1, bool smooth = false);
    ~Icosphere() = default;

    // getters/setters
    [[nodiscard]] float getRadius() const {
        return radius;
    }
    void setRadius(float radius);
    [[nodiscard]] int getSubdivision() const {
        return subdivision;
    }
    void setSubdivision(int subdivision);
    [[nodiscard]] bool getSmooth() const {
        return smooth;
    }
    void setSmooth(bool smooth);

    // for vertex data
    [[nodiscard]] unsigned int getVertexCount() const {
        return (unsigned int)vertices.size();
    }
    [[nodiscard]] unsigned int getNormalCount() const {
        return (unsigned int)normals.size();
    }
    [[nodiscard]] unsigned int getIndexCount() const {
        return (unsigned int)indices.size();
    }
    [[nodiscard]] unsigned int getLineIndexCount() const {
        return (unsigned int)lineIndices.size();
    }
    [[nodiscard]] unsigned int getTriangleCount() const {
        return getIndexCount() / 3;
    }

    [[nodiscard]] unsigned int getVertexSize() const {
        return (unsigned int)vertices.size() * sizeof(float);
    } // # of bytes
    [[nodiscard]] unsigned int getNormalSize() const {
        return (unsigned int)normals.size() * sizeof(float);
    }
    [[nodiscard]] unsigned int getIndexSize() const {
        return (unsigned int)indices.size() * sizeof(unsigned int);
    }
    [[nodiscard]] unsigned int getLineIndexSize() const {
        return (unsigned int)lineIndices.size() * sizeof(unsigned int);
    }

    [[nodiscard]] const glm::vec3* getVertices() const {
        if (vertices.empty())
        {
            return nullptr;
        }
        else
        {
            return &vertices[0];
        }
        
        //return vertices.data();
    }
    const glm::vec3* getNormals() const {
        if (normals.empty())
        {
            return nullptr;
        }
        else
        {
            return &normals[0];
        }
    }

    const unsigned int* getIndices() const {
        return indices.data();
    }
    const unsigned int* getLineIndices() const {
        return lineIndices.data();
    }

    // for interleaved vertices: V/N/T
    unsigned int getInterleavedVertexCount() const {
        return getVertexCount();
    } // # of vertices
    unsigned int getInterleavedVertexSize() const {
        return (unsigned int)interleavedVertices.size() * sizeof(float);
    } // # of bytes
    int getInterleavedStride() const {
        return interleavedStride;
    } // should be 32 bytes
    const float* getInterleavedVertices() const {
        return interleavedVertices.data();
    }

    // draw in VertexArray mode
    void draw() const;
    void drawLines(const float lineColor[4]) const;
    void drawWithLines(const float lineColor[4]) const;

    // debug
    void printSelf() const;


protected:
private:
    // static functions
    static void computeFaceNormal(const glm::vec3 v1, const glm::vec3 v2, const glm::vec3 v3, glm::vec3& normal);
    static void computeVertexNormal(const glm::vec3 v, glm::vec3& normal);
    static float computeScaleForLength(const glm::vec3, float length);
    static void computeHalfVertex(const glm::vec3 v1, const glm::vec3 v2, float length, glm::vec3& newV);
    static bool isOnLineSegment(const float a[2], const float b[2], const glm::vec2 c);

    // member functions
    void updateRadius();
    std::vector<glm::vec3> computeIcosahedronVertices();
    void buildVerticesFlat();
    void buildVerticesSmooth();
    void subdivideVerticesFlat();
    void subdivideVerticesSmooth();
    void buildInterleavedVertices();
    void addVertex(float x, float y, float z);
    void addVertices(const glm::vec3 v1, const glm::vec3 v2, const glm::vec3 v3);
    void addNormal(float nx, float ny, float nz);
    void addNormals(const glm::vec3 n1, const glm::vec3 n2, const glm::vec3 n3);
    void addIndices(unsigned int i1, unsigned int i2, unsigned int i3);
    void addSubLineIndices(
        unsigned int i1, unsigned int i2, unsigned int i3, unsigned int i4, unsigned int i5, unsigned int i6);
    unsigned int addSubVertexAttribs(const glm::vec3 v, const  glm::vec3 n);

    // memeber vars
    float radius; // circumscribed radius
    int subdivision;
    bool smooth;

    std::vector<glm::vec3> vertices;
    std::vector<unsigned int> vertex_index;

    std::vector<glm::vec3> normals;
    std::vector<unsigned int> indices;
    std::vector<unsigned int> lineIndices;
    std::map<std::pair<float, float>, unsigned int> sharedIndices; // indices of shared vertices, key is tex coord (s,t)

    // interleaved
    std::vector<float> interleavedVertices;
    int interleavedStride; // # of bytes to hop to the next vertex (should be 32 bytes)
};

#endif
