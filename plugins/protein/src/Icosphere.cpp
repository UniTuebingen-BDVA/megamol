///////////////////////////////////////////////////////////////////////////////
// Icosphere.cpp
// =============
// Polyhedron subdividing icosahedron (20 tris) by N-times iteration
// The icosphere with N=1 (default) has 80 triangles by subdividing a triangle
// of icosahedron into 4 triangles. If N=0, it is identical to icosahedron.
//
//  AUTHOR: Song Ho Ahn (song.ahn@gmail.com)
// CREATED: 2018-07-23
// UPDATED: 2019-12-18
///////////////////////////////////////////////////////////////////////////////

#ifdef _WIN32
#include <windows.h>    // include windows.h to avoid thousands of compile errors even though this class is not depending on Windows
#endif

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#include <iostream>
#include <iomanip>
#include <cmath>
#include "protein/Icosphere.h"

#include <vector>
#include <tuple>
#include <algorithm>
#include <glm/glm.hpp>


// constants //////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////
// ctor
///////////////////////////////////////////////////////////////////////////////
Icosphere::Icosphere(float radius, int sub, bool smooth) : radius(radius), subdivision(sub), smooth(smooth), interleavedStride(32)
{
    if(smooth)
        buildVerticesSmooth();
    else
        buildVerticesFlat();
}



///////////////////////////////////////////////////////////////////////////////
// setters
///////////////////////////////////////////////////////////////////////////////
void Icosphere::setRadius(float radius)
{
    this->radius = radius;
    updateRadius(); // update vertex positions only
}

void Icosphere::setSubdivision(int iteration)
{
    this->subdivision = iteration;
    // rebuild vertices
    if(smooth)
        buildVerticesSmooth();
    else
        buildVerticesFlat();
}

void Icosphere::setSmooth(bool smooth)
{
    if(this->smooth == smooth)
        return;

    this->smooth = smooth;
    if(smooth)
        buildVerticesSmooth();
    else
        buildVerticesFlat();
}



///////////////////////////////////////////////////////////////////////////////
// print itself
///////////////////////////////////////////////////////////////////////////////
void Icosphere::printSelf() const
{

    std::cout << "===== Icosphere =====\n"
              << "        Radius: " << radius << "\n"
              << "   Subdivision: " << subdivision << "\n"
              << "    Smoothness: " << (smooth ? "true" : "false") << "\n"
              << "Triangle Count: " << getTriangleCount() << "\n"
              << "   Index Count: " << getIndexCount() << "\n"
              << "  Vertex Count: " << getVertexCount() << "\n"
              << "  Normal Count: " << getNormalCount() << "\n" << std::endl;
}



///////////////////////////////////////////////////////////////////////////////
// draw a icosphere in VertexArray mode
// OpenGL RC must be set before calling it
///////////////////////////////////////////////////////////////////////////////
void Icosphere::draw() const
{
    // interleaved array
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glVertexPointer(3, GL_FLOAT, interleavedStride, &interleavedVertices[0]);
    glNormalPointer(GL_FLOAT, interleavedStride, &interleavedVertices[3]);

    glDrawElements(GL_TRIANGLES, (unsigned int)indices.size(), GL_UNSIGNED_INT, indices.data());

    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
}



///////////////////////////////////////////////////////////////////////////////
// draw lines only
// the caller must set the line width before call this
///////////////////////////////////////////////////////////////////////////////
void Icosphere::drawLines(const float lineColor[4]) const
{
    // set line colour
    glColor4fv(lineColor);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,   lineColor);

    // draw lines with VA
    glDisable(GL_LIGHTING);
    glDisable(GL_TEXTURE_2D);
    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(3, GL_FLOAT, 0, vertices.data());

    glDrawElements(GL_LINES, (unsigned int)lineIndices.size(), GL_UNSIGNED_INT, lineIndices.data());

    glDisableClientState(GL_VERTEX_ARRAY);
    glEnable(GL_LIGHTING);
    glEnable(GL_TEXTURE_2D);
}



///////////////////////////////////////////////////////////////////////////////
// draw a icosphere surfaces and lines on top of it
// the caller must set the line width before call this
///////////////////////////////////////////////////////////////////////////////
void Icosphere::drawWithLines(const float lineColor[4]) const
{
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0, 1.0f); // move polygon backward
    this->draw();
    glDisable(GL_POLYGON_OFFSET_FILL);

    // draw lines with VA
    drawLines(lineColor);
}



///////////////////////////////////////////////////////////////////////////////
// update vertex positions only
///////////////////////////////////////////////////////////////////////////////
void Icosphere::updateRadius()
{
    float scale = computeScaleForLength(vertices[0], radius);

    std::size_t i, j;
    std::size_t count = vertices.size();
    for (i = 0, j = 0; i < count; i++, j += 6) {
        /*
        for (int k = 0; k < 3; k++)
        {
            std::get<k>(vertices[i]) *= scale;
        }
        */
        vertices[i] *= scale;
        

        //TODO
        // for interleaved array
        interleavedVertices[j] *= scale;
        interleavedVertices[j + 1] *= scale;
        interleavedVertices[j + 2] *= scale;
    }
}


///////////////////////////////////////////////////////////////////////////////
// compute 12 vertices of icosahedron using spherical coordinates
// The north pole is at (0, 0, r) and the south pole is at (0,0,-r).
// 5 vertices are placed by rotating 72 deg at elevation 26.57 deg (=atan(1/2))
// 5 vertices are placed by rotating 72 deg at elevation -26.57 deg
///////////////////////////////////////////////////////////////////////////////
std::vector<glm::vec3> Icosphere::computeIcosahedronVertices()
{
    const float PI = acos(-1);
    const float H_ANGLE = PI / 180 * 72;   // 72 degree = 360 / 5
    const float V_ANGLE = atanf(1.0f / 2); // elevation = 26.565 degree

    std::vector<glm::vec3> vertices(12);
    int i2;                                // indices
    float z, xy;                           // coords
    float hAngle1 = -PI / 2 - H_ANGLE / 2; // start from -126 deg at 2nd row
    float hAngle2 = -PI / 2;               // start from -90 deg at 3rd row

    vertices[0] = glm::vec3(0.0f, 0.0f, radius);

    // 10 vertices at 2nd and 3rd rows
    for (int i = 1; i <= 5; ++i) {
        i2 = i + 5;

        z = radius * sinf(V_ANGLE); // elevaton
        xy = radius * cosf(V_ANGLE);

        /*
        std::tuple<float, float, float> xyz1(xy * cosf(hAngle1), xy * sinf(hAngle1), z);
        std::tuple<float, float, float> xyz2(xy * cosf(hAngle2), xy * sinf(hAngle2), -z);
        */

        vertices[i] = glm::vec3(xy * cosf(hAngle1), xy * sinf(hAngle1), z);
        vertices[i2] = glm::vec3(xy * cosf(hAngle2), xy * sinf(hAngle2), -z);

        // next horizontal angles
        hAngle1 += H_ANGLE;
        hAngle2 += H_ANGLE;
    }

    // the last bottom vertex (0, 0, -r)
    vertices[11] = glm::vec3(0.0f, 0.0f, -radius);
    

    return vertices; 
}



///////////////////////////////////////////////////////////////////////////////
// generate vertices with flat shading
// each triangle is independent (no shared vertices)
///////////////////////////////////////////////////////////////////////////////
void Icosphere::buildVerticesFlat()
{
    //const float S_STEP = 1 / 11.0f;         // horizontal texture step
    //const float T_STEP = 1 / 3.0f;          // vertical texture step
    const float S_STEP = 186 / 2048.0f;     // horizontal texture step
    const float T_STEP = 322 / 1024.0f;     // vertical texture step

    // compute 12 vertices of icosahedron
    std::vector<glm::vec3> tmpVertices = computeIcosahedronVertices();

    // clear memory of prev arrays
    std::vector<glm::vec3>().swap(vertices);
    std::vector<glm::vec3>().swap(normals);
    std::vector<unsigned int>().swap(indices);
    std::vector<unsigned int>().swap(lineIndices);

    glm::vec3 v0, v1, v2, v3, v4, v11;          // vertex positions
    glm::vec3 n;                                         // face normal
    unsigned int index = 0;

    // compute and add 20 tiangles of icosahedron first
    v0 = tmpVertices[0];       // 1st vertex
    v11 = tmpVertices[11]; // 12th vertex
    for(int i = 1; i <= 5; ++i)
    {
        // 4 vertices in the 2nd row
        v1 = tmpVertices[i];
        //if i is less than 5 use i + 1 as next neighbour, otherwise go back to i
        v2 = tmpVertices[(i < 5) ? (i + 1) : i];
        //go to next row
        v3 = tmpVertices[i + 5];
        //analogous to v2 for the next row
        v4 = tmpVertices[((i + 5) < 10) ? (i + 6) : 6];

        // add a triangle in 1st row
        Icosphere::computeFaceNormal(v0, v1, v2, n);
        addVertices(v0, v1, v2);
        addNormals(n, n, n);
        addIndices(index, index+1, index+2);

        // add 2 triangles in 2nd row
        Icosphere::computeFaceNormal(v1, v3, v2, n);
        addVertices(v1, v3, v2);
        addNormals(n, n, n);
        addIndices(index+3, index+4, index+5);

        Icosphere::computeFaceNormal(v2, v3, v4, n);
        addVertices(v2, v3, v4);
        addNormals(n, n, n);
        addIndices(index+6, index+7, index+8);

        // add a triangle in 3rd row
        Icosphere::computeFaceNormal(v3, v11, v4, n);
        addVertices(v3, v11, v4);
        addNormals(n, n, n);
        addIndices(index+9, index+10, index+11);

        // add 6 edge lines per iteration
        //  i
        //  /   /   /   /   /       : (i, i+1)                              //
        // /__ /__ /__ /__ /__                                              //
        // \  /\  /\  /\  /\  /     : (i+3, i+4), (i+3, i+5), (i+4, i+5)    //
        //  \/__\/__\/__\/__\/__                                            //
        //   \   \   \   \   \      : (i+9,i+10), (i+9, i+11)               //
        //    \   \   \   \   \                                             //
        lineIndices.push_back(index);       // (i, i+1)
        lineIndices.push_back(index+1);       // (i, i+1)
        lineIndices.push_back(index+3);     // (i+3, i+4)
        lineIndices.push_back(index+4);
        lineIndices.push_back(index+3);     // (i+3, i+5)
        lineIndices.push_back(index+5);
        lineIndices.push_back(index+4);     // (i+4, i+5)
        lineIndices.push_back(index+5);
        lineIndices.push_back(index+9);     // (i+9, i+10)
        lineIndices.push_back(index+10);
        lineIndices.push_back(index+9);     // (i+9, i+11)
        lineIndices.push_back(index+11);

        // next index
        index += 12;
    }

    // subdivide icosahedron
    subdivideVerticesFlat();

    // generate interleaved vertex array as well
    buildInterleavedVertices();
}


///////
//EVTL TODO
///////

///////////////////////////////////////////////////////////////////////////////
// generate vertices with smooth shading
// NOTE: The north and south pole vertices cannot be shared for smooth shading
// because they have same position and normal, but different texcoords per face
// And, the first vertex on each row is also not shared.
///////////////////////////////////////////////////////////////////////////////
void Icosphere::buildVerticesSmooth() {
    //const float S_STEP = 1 / 11.0f;         // horizontal texture step
    //const float T_STEP = 1 / 3.0f;          // vertical texture step
    const float S_STEP = 186 / 2048.0f; // horizontal texture step
    const float T_STEP = 322 / 1024.0f; // vertical texture step

    // compute 12 vertices of icosahedron
    // NOTE: v0 (top), v11(bottom), v1, v6(first vert on each row) cannot be
    // shared for smooth shading (they have different texcoords)
    std::vector<glm::vec3> tmpVertices = computeIcosahedronVertices();

    // clear memory of prev arrays
    std::vector<glm::vec3>().swap(vertices);
    std::vector<glm::vec3>().swap(normals);
    std::vector<unsigned int>().swap(indices);
    std::vector<unsigned int>().swap(lineIndices);
    std::map<std::pair<float, float>, unsigned int>().swap(sharedIndices);

    glm::vec3 v;    // vertex
    glm::vec3 n;  // normal
    float scale; // scale factor for normalization

    // smooth icosahedron has 14 non-shared (0 to 13) and
    // 8 shared vertices (14 to 21) (total 22 vertices)
    //        OLD                   //        NEW                  //        WRONG
    //  00  01  02  03  04          //  00  00  00  00  00         //  00  00  00  00  00         //v0                     //tmpv0
    //  /\  /\  /\  /\  /\          //  /\  /\  /\  /\  /\         //  /\  /\  /\  /\  /\          
    // /  \/  \/  \/  \/  \         // /  \/  \/  \/  \/  \        // /  \/  \/  \/  \/  \    
    //10--14--15--16--17--11        //02--04--05--06--07--02       //02--06--07--08--09--03       //v2 v4 v5 v6 v7 v2      //tmpv1 tmpv2 tmpv3 tmpv4 tmpv5
    // \  /\  /\  /\  /\  /\        // \  /\  /\  /\  /\  /\       // \  /\  /\  /\  /\  /\        
    //  \/  \/  \/  \/  \/  \       //  \/  \/  \/  \/  \/  \      //  \/  \/  \/  \/  \/  \       
    //  12--18--19--20--21--13      //  03--08--09--10--11--03     //  04--10--11--12--13--05     //v3 v8 v9 v10 v11 v3    //tmpv6 tmpv7 tmpv8 tmpv9 tmpv10
    //   \  /\  /\  /\  /\  /       //   \  /\  /\  /\  /\  /      //   \  /\  /\  /\  /\  / 
    //    \/  \/  \/  \/  \/        //    \/  \/  \/  \/  \/       //    \/  \/  \/  \/  \/
    //    05  06  07  08  09        //    01  01  01  01  01       //    01  01  01  01  01       //v1                     //tmpv11
    // add 14 non-shared vertices first (index from 0 to 13)
    addVertex(tmpVertices[0][0], tmpVertices[0][1], tmpVertices[0][2]); //v0 (top)
    addNormal(0, 0, 1);

    addVertex(tmpVertices[11][0], tmpVertices[11][1], tmpVertices[11][2]); //v1 (bottom)
    addNormal(0, 0, -1);

    v = tmpVertices[1];  //v2
    Icosphere::computeVertexNormal(v, n);
    addVertex(v[0], v[1], v[2]);
    addNormal(n[0], n[1], n[2]);

    v = tmpVertices[6]; //v3
    Icosphere::computeVertexNormal(v, n);
    addVertex(v[0], v[1], v[2]);
    addNormal(n[0], n[1], n[2]);

    v = tmpVertices[2]; //v4
    Icosphere::computeVertexNormal(v, n);
    addVertex(v[0], v[1], v[2]);
    addNormal(n[0], n[1], n[2]);

    v = tmpVertices[3]; //v5
    Icosphere::computeVertexNormal(v, n);
    addVertex(v[0], v[1], v[2]);
    addNormal(n[0], n[1], n[2]);

    v = tmpVertices[4]; //v6
    //scale = Icosphere::computeScaleForLength(v, 1);
    //n[0] = v[0] * scale;    n[1] = v[1] * scale;    n[2] = v[2] * scale;
    Icosphere::computeVertexNormal(v, n);
    addVertex(v[0], v[1], v[2]);
    addNormal(n[0], n[1], n[2]);

    v = tmpVertices[5]; //v7
    Icosphere::computeVertexNormal(v, n);
    addVertex(v[0], v[1], v[2]);
    addNormal(n[0], n[1], n[2]);

    v = tmpVertices[7]; //v8
    Icosphere::computeVertexNormal(v, n);
    addVertex(v[0], v[1], v[2]);
    addNormal(n[0], n[1], n[2]);

    v = tmpVertices[8]; //v9
    Icosphere::computeVertexNormal(v, n);
    addVertex(v[0], v[1], v[2]);
    addNormal(n[0], n[1], n[2]);

    v = tmpVertices[9]; //v10
    Icosphere::computeVertexNormal(v, n);
    addVertex(v[0], v[1], v[2]);
    addNormal(n[0], n[1], n[2]);

    v = tmpVertices[10]; //v11
    Icosphere::computeVertexNormal(v, n);
    addVertex(v[0], v[1], v[2]);
    addNormal(n[0], n[1], n[2]);

    //build index list for icosahedron (20 triangles)  
    //first row
    addIndices( 0,  2,  4);
    addIndices( 0,  4,  5);
    addIndices( 0,  5,  6);
    addIndices( 0,  6,  7);
    addIndices( 0,  7,  2);
    //second row
    addIndices( 2,  3,  4);
    addIndices( 3,  8,  4);
    addIndices( 4,  8,  5);
    addIndices( 8,  9,  5);
    addIndices( 5,  9,  6);
    addIndices( 9, 10,  6);
    addIndices( 6, 10,  7);
    addIndices(10, 11,  7);
    addIndices( 7, 11,  2);
    addIndices(11,  3,  2);
    //third row
    addIndices( 1,  8,  3);
    addIndices( 1,  9,  8);
    addIndices( 1, 10,  9);
    addIndices( 1, 11, 10);
    addIndices( 1,  3, 11);
    
    //add edge lines of icosahedron   
    //first row
    lineIndices.push_back(0);   lineIndices.push_back(2);
    lineIndices.push_back(0);   lineIndices.push_back(4);
    lineIndices.push_back(0);   lineIndices.push_back(5);
    lineIndices.push_back(0);   lineIndices.push_back(6);
    lineIndices.push_back(0);   lineIndices.push_back(7);

    lineIndices.push_back(2);   lineIndices.push_back(4);
    lineIndices.push_back(4);   lineIndices.push_back(5);
    lineIndices.push_back(5);   lineIndices.push_back(6);
    lineIndices.push_back(6);   lineIndices.push_back(7);
    lineIndices.push_back(7);   lineIndices.push_back(2);
    //second row
    lineIndices.push_back(2);   lineIndices.push_back(3);
    lineIndices.push_back(3);   lineIndices.push_back(4);
    lineIndices.push_back(4);   lineIndices.push_back(8);
    lineIndices.push_back(8);   lineIndices.push_back(5);
    lineIndices.push_back(5);   lineIndices.push_back(9);
    lineIndices.push_back(9);   lineIndices.push_back(6);
    lineIndices.push_back(6);   lineIndices.push_back(10);
    lineIndices.push_back(10);  lineIndices.push_back(7);
    lineIndices.push_back(7);   lineIndices.push_back(11);
    lineIndices.push_back(11);  lineIndices.push_back(2);

    lineIndices.push_back(3);   lineIndices.push_back(8);
    lineIndices.push_back(8);   lineIndices.push_back(9);
    lineIndices.push_back(9);   lineIndices.push_back(10);
    lineIndices.push_back(10);  lineIndices.push_back(11);
    lineIndices.push_back(11);  lineIndices.push_back(3);
    //third row
    lineIndices.push_back(1);   lineIndices.push_back(3);
    lineIndices.push_back(1);   lineIndices.push_back(8);
    lineIndices.push_back(1);   lineIndices.push_back(9);
    lineIndices.push_back(1);   lineIndices.push_back(10);
    lineIndices.push_back(1);   lineIndices.push_back(11);
    
    // subdivide icosahedron
    subdivideVerticesSmooth();

    // generate interleaved vertex array as well
    buildInterleavedVertices();
}



///////////////////////////////////////////////////////////////////////////////
// divide a trinage into 4 sub triangles and repeat N times
// If subdivision=0, do nothing.
///////////////////////////////////////////////////////////////////////////////
void Icosphere::subdivideVerticesFlat()
{
    std::vector<glm::vec3> tmpVertices;
    std::vector<unsigned int> tmpIndices;
    int indexCount;
    glm::vec3 *v1, *v2, *v3;          // ptr to original vertices of a triangle
    glm::vec3 tmpv1, tmpv2, tmpv3;
    glm::vec3 newV1, newV2, newV3; // new vertex positions
    glm::vec3 normal;                    // new face normal
    unsigned int index = 0;             // new index value
    int i, j;

    // iteration
    for(i = 1; i <= subdivision; ++i)
    {
        // copy prev arrays
        tmpVertices = vertices;
        tmpIndices = indices;

        // clear prev arrays
        vertices.clear();
        normals.clear();
        indices.clear();
        lineIndices.clear();

        index = 0;
        indexCount = (int)tmpIndices.size();
        for(j = 0; j < indexCount; j += 3)
        {
            // get 3 vertice and texcoords of a triangle
            v1 = &tmpVertices[tmpIndices[j]];
            v2 = &tmpVertices[tmpIndices[j + 1]];
            v3 = &tmpVertices[tmpIndices[j + 2]];

            tmpv1 = *v1;
            tmpv2 = *v2;
            tmpv3 = *v3;

            // get 3 new vertices by spliting half on each edge
            computeHalfVertex(tmpv1, tmpv2, radius, newV1);
            computeHalfVertex(tmpv2, tmpv3, radius, newV2);
            computeHalfVertex(tmpv1, tmpv3, radius, newV3);

            // add 4 new triangles
            addVertices(tmpv1, newV1, newV3);
            computeFaceNormal(tmpv1, newV1, newV3, normal);
            addNormals(normal, normal, normal);
            addIndices(index, index+1, index+2);

            addVertices(newV1, tmpv2, newV2);
            computeFaceNormal(newV1, tmpv2, newV2, normal);
            addNormals(normal, normal, normal);
            addIndices(index+3, index+4, index+5);

            addVertices(newV1, newV2, newV3);
            computeFaceNormal(newV1, newV2, newV3, normal);
            addNormals(normal, normal, normal);
            addIndices(index+6, index+7, index+8);

            addVertices(newV3, newV2, tmpv3);
            computeFaceNormal(newV3, newV2, tmpv3, normal);
            addNormals(normal, normal, normal);
            addIndices(index+9, index+10, index+11);

            // add new line indices per iteration
            addSubLineIndices(index, index+1, index+4, index+5, index+11, index+9); //CCW

            // next index
            index += 12;
        }
    }
}



///////////////////////////////////////////////////////////////////////////////
// divide a trianlge (v1-v2-v3) into 4 sub triangles by adding middle vertices
// (newV1, newV2, newV3) and repeat N times
// If subdivision=0, do nothing.
//         v1           //
//        / \           //
// newV1 *---* newV3    //
//      / \ / \         //
//    v2---*---v3       //
//        newV2         //
///////////////////////////////////////////////////////////////////////////////
void Icosphere::subdivideVerticesSmooth()
{
    std::vector<unsigned int> tmpIndices;
    int indexCount;
    unsigned int i1, i2, i3;            // indices from original triangle
    glm::vec3 *v1, *v2, *v3; // ptr to original vertices of a triangle
    glm::vec3 tmpv1, tmpv2, tmpv3;
    glm::vec3 newV1, newV2, newV3; // new subdivided vertex positions
    glm::vec3 newN1, newN2, newN3; // new subdivided normals
    unsigned int newI1, newI2, newI3;   // new subdivided indices
    int i, j;

    // iteration for subdivision
    for(i = 1; i <= subdivision; ++i)
    {
        // copy prev indices
        tmpIndices = indices;

        // clear prev arrays
        indices.clear();
        lineIndices.clear();

        indexCount = (int)tmpIndices.size();
        for(j = 0; j < indexCount; j += 3)
        {
            // get 3 indices of each triangle
            i1 = tmpIndices[j];
            i2 = tmpIndices[j+1];
            i3 = tmpIndices[j+2];

            // get 3 vertex attribs from prev triangle
            v1 = &vertices[i1];
            v2 = &vertices[i2];
            v3 = &vertices[i3];

            tmpv1 = *v1;
            tmpv2 = *v2;
            tmpv3 = *v3;

            // get 3 new vertex attribs by spliting half on each edge
            computeHalfVertex(tmpv1, tmpv2, radius, newV1);
            computeHalfVertex(tmpv2, tmpv3, radius, newV2);
            computeHalfVertex(tmpv1, tmpv3, radius, newV3);
            computeVertexNormal(newV1, newN1);
            computeVertexNormal(newV2, newN2);
            computeVertexNormal(newV3, newN3);

            // add new vertices/normals/texcoords to arrays
            // It will check if it is shared/non-shared and return index
            newI1 = addSubVertexAttribs(newV1, newN1);
            newI2 = addSubVertexAttribs(newV2, newN2);
            newI3 = addSubVertexAttribs(newV3, newN3);

            // add 4 new triangle indices
            addIndices(i1, newI1, newI3);
            addIndices(newI1, i2, newI2);
            addIndices(newI1, newI2, newI3);
            addIndices(newI3, newI2, i3);

            // add new line indices
            addSubLineIndices(i1, newI1, i2, newI2, i3, newI3); //CCW
        }
    }
}



///////////////////////////////////////////////////////////////////////////////
// generate interleaved vertices: V/N/T
// stride must be 32 bytes
///////////////////////////////////////////////////////////////////////////////
void Icosphere::buildInterleavedVertices()
{
    std::vector<float>().swap(interleavedVertices);

    std::size_t i, j;
    std::size_t count = vertices.size();
    for(i = 0; i < count; i += 3)
    {
        interleavedVertices.push_back(vertices[i][0]);
        interleavedVertices.push_back(vertices[i][1]);
        interleavedVertices.push_back(vertices[i][2]);

        interleavedVertices.push_back(normals[i][0]);
        interleavedVertices.push_back(normals[i][1]);
        interleavedVertices.push_back(normals[i][2]);

    }
}



///////////////////////////////////////////////////////////////////////////////
// add single vertex to array
///////////////////////////////////////////////////////////////////////////////
void Icosphere::addVertex(float x, float y, float z)
{
    glm::vec3 xyz(x, y, z);
    vertices.emplace_back(xyz);
}


///////////////////////////////////////////////////////////////////////////////
// add 3 vertices of a triangle to array
// only add the index of vertice in vector 'vertices'
///////////////////////////////////////////////////////////////////////////////
void Icosphere::addVertices(const glm::vec3 v1, const glm::vec3 v2, const glm::vec3 v3) {
    auto compare_vec3 = [&](const glm::vec3& vec1, const glm::vec3& vec2)
    { return vec1 == vec2; };

    unsigned int index;
    std::vector<glm::vec3> tmpVec = {v1, v2, v3};

    for (const glm::vec3& input_vec : tmpVec) 
     {
         auto vec3_search = std::find_if(vertices.begin(), vertices.end(),
             [&](const glm::vec3& vec)
         {
             return compare_vec3(vec, input_vec);
         });

         if (vec3_search != vertices.end())
         {
             index = static_cast<unsigned int>(std::distance(vertices.begin(), vec3_search));
         }

         vertex_index.emplace_back(index);
     }
 }




///////////////////////////////////////////////////////////////////////////////
// add single normal to array
///////////////////////////////////////////////////////////////////////////////
void Icosphere::addNormal(float nx, float ny, float nz)
{
    glm::vec3 nxyz = {nx, ny, nz};
    normals.push_back(nxyz);
}



///////////////////////////////////////////////////////////////////////////////
// add 3 normals of a triangle to array
///////////////////////////////////////////////////////////////////////////////
void Icosphere::addNormals(const glm::vec3 n1, const glm::vec3 n2, const glm::vec3 n3)
{
    normals.push_back(n1);
    normals.push_back(n2);
    normals.push_back(n3);
}




///////////////////////////////////////////////////////////////////////////////
// add 3 indices to array
///////////////////////////////////////////////////////////////////////////////
void Icosphere::addIndices(unsigned int i1, unsigned int i2, unsigned int i3)
{
    indices.push_back(i1);
    indices.push_back(i2);
    indices.push_back(i3);
}



///////////////////////////////////////////////////////////////////////////////
// add 7 sub edge lines per triangle to array using 6 indices (CCW)
//     i1                                           //
//     /            : (i1, i2)                      //
//   i2---i6        : (i2, i6)                      //
//   / \  /         : (i2, i3), (i2, i4), (i6, i4)  //
// i3---i4---i5     : (i3, i4), (i4, i5)            //
///////////////////////////////////////////////////////////////////////////////
void Icosphere::addSubLineIndices(unsigned int i1,
                                  unsigned int i2,
                                  unsigned int i3,
                                  unsigned int i4,
                                  unsigned int i5,
                                  unsigned int i6)
{
    lineIndices.push_back(i1);      // i1 - i2
    lineIndices.push_back(i2);
    lineIndices.push_back(i2);      // i2 - i6
    lineIndices.push_back(i6);
    lineIndices.push_back(i2);      // i2 - i3
    lineIndices.push_back(i3);
    lineIndices.push_back(i2);      // i2 - i4
    lineIndices.push_back(i4);
    lineIndices.push_back(i6);      // i6 - i4
    lineIndices.push_back(i4);
    lineIndices.push_back(i3);      // i3 - i4
    lineIndices.push_back(i4);
    lineIndices.push_back(i4);      // i4 - i5
    lineIndices.push_back(i5);
}



///////////////////////////////////////////////////////////////////////////////
// add a subdivided vertex attribs (vertex, normal, texCoord) to arrays, then
// return its index value
// If it is a shared vertex, remember its index, so it can be re-used
///////////////////////////////////////////////////////////////////////////////
unsigned int Icosphere::addSubVertexAttribs(const glm::vec3 v, const glm::vec3 n) {
    unsigned int index;     // return value;

    addVertex(v[0], v[1], v[2]);
    addNormal(n[0], n[1], n[2]);
    index = vertices.size() - 1;
  
    return index;
}




// static functions ===========================================================
///////////////////////////////////////////////////////////////////////////////
// return face normal (4th param) of a triangle v1-v2-v3
// if a triangle has no surface (normal length = 0), then return a zero vector
///////////////////////////////////////////////////////////////////////////////
void Icosphere::computeFaceNormal(const glm::vec3 v1, const glm::vec3 v2, const glm::vec3 v3, glm::vec3& n) {
    const float EPSILON = 0.000001f;

    // default return value (0, 0, 0)
    n = {0, 0, 0};
    //n[0] = n[1] = n[2] = 0;

    // find 2 edge vectors: v1-v2, v1-v3
    float ex1 = v2[0] - v1[0];
    float ey1 = v2[1] - v1[1];
    float ez1 = v2[2] - v1[2];
    float ex2 = v3[0] - v1[0];
    float ey2 = v3[1] - v1[1];
    float ez2 = v3[2] - v1[2];

    // cross product: e1 x e2
    float nx, ny, nz;
    nx = ey1 * ez2 - ez1 * ey2;
    ny = ez1 * ex2 - ex1 * ez2;
    nz = ex1 * ey2 - ey1 * ex2;

    // normalize only if the length is > 0
    float length = sqrtf(nx * nx + ny * ny + nz * nz);
    if(length > EPSILON)
    {
        // normalize
        float lengthInv = 1.0f / length;
        n[0] = nx * lengthInv;
        n[1] = ny * lengthInv;
        n[2] = nz * lengthInv;
    }
}



///////////////////////////////////////////////////////////////////////////////
// return vertex normal (2nd param) by mormalizing the vertex vector
///////////////////////////////////////////////////////////////////////////////
void Icosphere::computeVertexNormal(const glm::vec3 v, glm::vec3& normal) {
    // normalize
    float scale = Icosphere::computeScaleForLength(v, 1);
    normal[0] = v[0] * scale;
    normal[1] = v[1] * scale;
    normal[2] = v[2] * scale;
}



///////////////////////////////////////////////////////////////////////////////
// get the scale factor for vector to resize to the given length of vector
///////////////////////////////////////////////////////////////////////////////
float Icosphere::computeScaleForLength(const glm::vec3 v, float length)
{
    // and normalize the vector then re-scale to new radius
    return length / glm::length(v);
}



///////////////////////////////////////////////////////////////////////////////
// find middle point of 2 vertices
// NOTE: new vertex must be resized, so the length is equal to the given length
///////////////////////////////////////////////////////////////////////////////
void Icosphere::computeHalfVertex(const glm::vec3 v1, const glm::vec3 v2, float length, glm::vec3& newV) {
    newV[0] = v1[0] + v2[0];
    newV[1] = v1[1] + v2[1];
    newV[2] = v1[2] + v2[2];
    float scale = Icosphere::computeScaleForLength(newV, length);
    newV[0] *= scale;
    newV[1] *= scale;
    newV[2] *= scale;
}




///////////////////////////////////////////////////////////////////////////////
// determine a point c is on the line segment a-b
///////////////////////////////////////////////////////////////////////////////
bool Icosphere::isOnLineSegment(const float a[2], const float b[2], const glm::vec2 c) {
    const float EPSILON = 0.0001f;

    // cross product must be 0 if c is on the line
    float cross = ((b[0] - a[0]) * (c[1] - a[1])) - ((b[1] - a[1]) * (c[0] - a[0]));
    if(cross > EPSILON || cross < -EPSILON)
        return false;

    // c must be within a-b
    if((c[0] > a[0] && c[0] > b[0]) || (c[0] < a[0] && c[0] < b[0]))
        return false;
    if((c[1] > a[1] && c[1] > b[1]) || (c[1] < a[1] && c[1] < b[1]))
        return false;

    return true;    // all passed, it is on the line segment
}
