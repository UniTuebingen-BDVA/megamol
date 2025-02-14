//
// StreamlineRenderer.cpp
//
// Copyright (C) 2013 by University of Stuttgart (VISUS).
// All rights reserved.
//
// Created on: Jun 11, 2013
//     Author: scharnkn
//

#include "StreamlineRenderer.h"

#include "CUDAFieldTopology.cuh"
#include "VBODataCall.h"
#include "cuda_error_check.h"
#include "ogl_error_check.h"
#include "protein_calls/Interpol.h"
#include "protein_calls/VTIDataCall.h"

#include "vislib/math/Cuboid.h"
#include "vislib/math/Vector.h"
#include "vislib/math/mathfunctions.h"

#include "mmcore/CoreInstance.h"
#include "mmcore/param/EnumParam.h"
#include "mmcore/param/FloatParam.h"
#include "mmcore/param/IntParam.h"
#include "mmcore_gl/utility/ShaderFactory.h"
#include "mmstd/renderer/CallClipPlane.h"

#include "vislib_gl/graphics/gl/IncludeAllGL.h"
#include <cstdlib>
#include <cuda_gl_interop.h>

using namespace megamol;
using namespace megamol::protein_cuda;
using namespace megamol::core;
using namespace megamol::core::utility::log;

typedef vislib::math::Vector<float, 3> Vec3f;
const Vec3f StreamlineRenderer::uniformColor = Vec3f(0.88f, 0.86f, 0.39f);


/*
 * StreamlineRenderer::StreamlineRenderer
 */
StreamlineRenderer::StreamlineRenderer(void)
        : Renderer3DModuleGL()
        ,
        /* Caller slots */
        fieldDataCallerSlot("getFieldData", "Connects the renderer with the field data")
        , getClipPlaneSlot("getClipPlane", "Provides the clip plane")
        ,
        /* Streamline integration parameters */
        nStreamlinesSlot("nStreamlines", "Set the number of streamlines")
        , streamlineMaxStepsSlot("nSteps", "Set the number of steps for streamline integration")
        , streamlineStepSlot("step", "Set stepsize for the streamline integration")
        , seedClipZSlot("seedClipZ", "Starting z value for the clipping plane")
        , seedIsoSlot("seedIso", "Iso value for the seed point")
        , renderModeSlot("renderMode", "Set rendermode for the streamlines")
        , streamtubesThicknessSlot("tubesScl", "The scale factor for the streamtubes thickness")
        , minColSlot("minCol", "Minimum color value")
        , maxColSlot("maxCol", "Maximum color value")
        , triggerComputeGradientField(true)
        , triggerComputeStreamlines(true) {

    // Data caller for volume data
    this->fieldDataCallerSlot.SetCompatibleCall<protein_calls::VTIDataCallDescription>();
    this->MakeSlotAvailable(&this->fieldDataCallerSlot);

    // Data caller for clipping plane
    view::CallClipPlaneDescription ccpd;
    this->getClipPlaneSlot.SetCallback(ccpd.ClassName(), ccpd.FunctionName(0), &StreamlineRenderer::requestPlane);
    this->MakeSlotAvailable(&this->getClipPlaneSlot);


    /* Streamline integration parameters */

    // Set the number of streamlines
    this->nStreamlines = 10;
    this->nStreamlinesSlot.SetParameter(new core::param::IntParam(this->nStreamlines, 0));
    this->MakeSlotAvailable(&this->nStreamlinesSlot);

    // Set the number of steps for streamline integration
    this->streamlineMaxSteps = 10;
    this->streamlineMaxStepsSlot.SetParameter(new core::param::IntParam(this->streamlineMaxSteps, 0));
    this->MakeSlotAvailable(&this->streamlineMaxStepsSlot);

    // Set the step size for streamline integration
    this->streamlineStep = 1.0f;
    this->streamlineStepSlot.SetParameter(new core::param::FloatParam(this->streamlineStep, 0.1f));
    this->MakeSlotAvailable(&this->streamlineStepSlot);

    // Set the step size for streamline integration
    this->seedClipZ = 0.5f;
    this->seedClipZSlot.SetParameter(new core::param::FloatParam(this->seedClipZ, 0.0f));
    this->MakeSlotAvailable(&this->seedClipZSlot);

    // Set the step size for streamline integration
    this->seedIso = 0.5f;
    this->seedIsoSlot.SetParameter(new core::param::FloatParam(this->seedIso));
    this->MakeSlotAvailable(&this->seedIsoSlot);


    /* Streamline render parameters */

    this->renderMode = NONE;
    param::EnumParam* rm = new param::EnumParam(int(this->renderMode));
    rm->SetTypePair(NONE, "None");
    rm->SetTypePair(LINES, "Lines");
    rm->SetTypePair(ILLUMINATED_LINES, "Illuminated lines");
    rm->SetTypePair(TUBES, "Stream tubes");
    this->renderModeSlot << rm;
    this->MakeSlotAvailable(&this->renderModeSlot);

    // Set the streamtubes thickness slot
    this->streamtubesThickness = 1.0f;
    this->streamtubesThicknessSlot.SetParameter(new core::param::FloatParam(this->streamtubesThickness, 0.0f));
    this->MakeSlotAvailable(&this->streamtubesThicknessSlot);

    // Set the streamtubes min color
    this->minCol = -1.0f;
    this->minColSlot.SetParameter(new core::param::FloatParam(this->minCol));
    this->MakeSlotAvailable(&this->minColSlot);

    // Set the streamtubes max color
    this->maxCol = 1.0f;
    this->maxColSlot.SetParameter(new core::param::FloatParam(this->maxCol));
    this->MakeSlotAvailable(&this->maxColSlot);
}


/*
 * StreamlineRenderer::~StreamlineRenderer
 */
StreamlineRenderer::~StreamlineRenderer(void) {
    this->Release();
}


/*
 * StreamlineRenderer::create
 */
bool StreamlineRenderer::create(void) {

    using namespace vislib_gl::graphics::gl;

    // Init extensions
    /*if (!ogl_IsVersionGEQ(2, 0) ||
        !areExtsAvailable("GL_ARB_vertex_shader GL_ARB_vertex_program GL_ARB_shader_objects GL_EXT_gpu_shader4 "
                          "GL_EXT_geometry_shader4 GL_EXT_bindable_uniform GL_ARB_draw_buffers GL_ARB_copy_buffer "
                          "GL_ARB_vertex_buffer_object")) {
        return false;
    }*/

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE_ARB);
    glEnable(GL_VERTEX_PROGRAM_TWO_SIDE);

    auto const shader_options = msf::ShaderFactoryOptionsOpenGL(GetCoreInstance()->GetShaderPaths());

    try {
        this->tubeShader =
            core::utility::make_glowl_shader("tubeShader", shader_options, "protein_cuda/streamlines/tube.vert.glsl",
                "protein_cuda/streamlines/tube.geom.glsl", "protein_cuda/streamlines/tube.frag.glsl");
        glProgramParameteriEXT(tubeShader->getHandle(), GL_GEOMETRY_INPUT_TYPE_EXT, GL_LINES_ADJACENCY);
        glProgramParameteriEXT(tubeShader->getHandle(), GL_GEOMETRY_OUTPUT_TYPE_EXT, GL_TRIANGLE_STRIP);
        glProgramParameteriEXT(tubeShader->getHandle(), GL_GEOMETRY_VERTICES_OUT_EXT, 200);
        // TODO Note from shader factory migration: No idea if relink is ok or setting parameter must happen before initial link within glowl.
        glLinkProgram(tubeShader->getHandle());

        this->illumShader = core::utility::make_glowl_shader("illumShader", shader_options,
            "protein_cuda/streamlines/illuminated.vert.glsl", "protein_cuda/streamlines/illuminated.geom.glsl",
            "protein_cuda/streamlines/illuminated.frag.glsl");
        glProgramParameteriEXT(illumShader->getHandle(), GL_GEOMETRY_INPUT_TYPE_EXT, GL_LINES_ADJACENCY);
        glProgramParameteriEXT(illumShader->getHandle(), GL_GEOMETRY_OUTPUT_TYPE_EXT, GL_LINE_STRIP);
        glProgramParameteriEXT(illumShader->getHandle(), GL_GEOMETRY_VERTICES_OUT_EXT, 200);
        // TODO Note from shader factory migration: No idea if relink is ok or setting parameter must happen before initial link within glowl.
        glLinkProgram(illumShader->getHandle());

    } catch (std::exception& e) {
        Log::DefaultLog.WriteError(("StreamlineRenderer: " + std::string(e.what())).c_str());
        return false;
    }

    return true;
}


/*
 * StreamlineRenderer::release
 */
void StreamlineRenderer::release(void) {
    this->tubeShader.reset();
    this->illumShader.reset();
}


/*
 * StreamlineRenderer::GetExtents
 */
bool StreamlineRenderer::GetExtents(mmstd_gl::CallRender3DGL& call) {
    // Extent of volume data

    protein_calls::VTIDataCall* vtiCall = this->fieldDataCallerSlot.CallAs<protein_calls::VTIDataCall>();
    if (!(*vtiCall)(protein_calls::VTIDataCall::CallForGetExtent)) {
        return false;
    }

    vtiCall->SetCalltime(call.Time());
    vtiCall->SetFrameID(static_cast<int>(call.Time()));
    if (!(*vtiCall)(protein_calls::VTIDataCall::CallForGetData)) {
        return false;
    }

    //this->bbox.SetObjectSpaceBBox(vtiCall->GetWholeExtent());
    this->bbox.SetObjectSpaceBBox(vtiCall->AccessBoundingBoxes().ObjectSpaceBBox());
    //this->bbox.MakeScaledWorld(scale);
    call.AccessBoundingBoxes() = vtiCall->AccessBoundingBoxes();
    call.SetTimeFramesCount(vtiCall->FrameCount());

    return true;
}


/*
 * StreamlineRenderer::Render
 */
bool StreamlineRenderer::Render(mmstd_gl::CallRender3DGL& call) {

    // Update parameters
    this->updateParams();

    protein_calls::VTIDataCall* vtiCall = this->fieldDataCallerSlot.CallAs<protein_calls::VTIDataCall>();
    if (vtiCall == NULL) {
        return false;
    }

    // Get volume data
    vtiCall->SetCalltime(call.Time());
    if (!(*vtiCall)(protein_calls::VTIDataCall::CallForGetData)) {
        return false;
    }

    // (Re)compute streamlines if necessary
    if (this->triggerComputeStreamlines) {

        float zHeight = (vtiCall->GetGridsize().GetZ() - 1) * vtiCall->GetSpacing().GetZ();
        this->genSeedPoints(vtiCall, zHeight * this->seedClipZ, this->seedIso); // Isovalues

        //printf("height: %f, clip %f %f %f %f\n", zHeight, zHeight*0.2, zHeight*0.4, zHeight*0.6,zHeight*0.8);

        if (!this->strLines.InitStreamlines(
                this->streamlineMaxSteps, this->nStreamlines, CUDAStreamlines::BIDIRECTIONAL)) {
            return false;
        }

        // Integrate streamlines
        if (!this->strLines.IntegrateRK4(this->seedPoints.PeekElements(), this->streamlineStep,
                (float*)vtiCall->GetPointDataByIdx(1, 0), // TODO Do not hardcode array
                make_int3(vtiCall->GetGridsize().GetX(), vtiCall->GetGridsize().GetY(), vtiCall->GetGridsize().GetZ()),
                make_float3(vtiCall->GetOrigin().GetX(), vtiCall->GetOrigin().GetY(), vtiCall->GetOrigin().GetZ()),
                make_float3(
                    vtiCall->GetSpacing().GetX(), vtiCall->GetSpacing().GetY(), vtiCall->GetSpacing().GetZ()))) {
            return false;
        }

        // Sample the density field to the alpha component
        if (!this->strLines.SampleScalarFieldToAlpha(
                (float*)vtiCall->GetPointDataByIdx(0, 0), // TODO do not hardcode array
                make_int3(vtiCall->GetGridsize().GetX(), vtiCall->GetGridsize().GetY(), vtiCall->GetGridsize().GetZ()),
                make_float3(vtiCall->GetOrigin().GetX(), vtiCall->GetOrigin().GetY(), vtiCall->GetOrigin().GetZ()),
                make_float3(
                    vtiCall->GetSpacing().GetX(), vtiCall->GetSpacing().GetY(), vtiCall->GetSpacing().GetZ()))) {
            return false;
        }

        // Set RGB color value
        //        if (!this->strLines.SetUniformRGBColor(make_float3(1.0, 0.0, 1.0))) {
        //            return false;
        //        }

        if (!this->strLines.SampleVecFieldToRGB((float*)vtiCall->GetPointDataByIdx(1, 0), // TODO do not hardcode array
                make_int3(vtiCall->GetGridsize().GetX(), vtiCall->GetGridsize().GetY(), vtiCall->GetGridsize().GetZ()),
                make_float3(vtiCall->GetOrigin().GetX(), vtiCall->GetOrigin().GetY(), vtiCall->GetOrigin().GetZ()),
                make_float3(
                    vtiCall->GetSpacing().GetX(), vtiCall->GetSpacing().GetY(), vtiCall->GetSpacing().GetZ()))) {
            return false;
        }

        this->triggerComputeStreamlines = false;
    }

    glPushMatrix();

    // Render streamlines
    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_LINE_SMOOTH);
    glColor3f(0.0f, 1.0f, 1.0f);

    if (this->renderMode == TUBES) {
        this->tubeShader->use();
        glUniform1fARB(this->tubeShader->getUniformLocation("streamTubeThicknessScl"), this->streamtubesThickness);
        glUniform1fARB(this->tubeShader->getUniformLocation("minColTexValue"), this->minCol); // TODO Parameter?
        glUniform1fARB(this->tubeShader->getUniformLocation("maxColTexValue"), this->maxCol); // TODO Paremeter?
        if (!this->strLines.RenderLineStripWithColor()) {
            return false;
        }
        glUseProgram(0);
    } else if (this->renderMode == LINES) {
        glColor3f(StreamlineRenderer::uniformColor.GetX(), StreamlineRenderer::uniformColor.GetY(),
            StreamlineRenderer::uniformColor.GetZ());
        if (!this->strLines.RenderLineStrip()) {
            return false;
        }
    } else if (this->renderMode == ILLUMINATED_LINES) {
        this->illumShader->use();
        if (!this->strLines.RenderLineStripWithColor()) {
            return false;
        }
        glUseProgram(0);
    }

    glPopMatrix();

    return true;
}


/*
 * StreamlineRenderer::genSeedPoints
 */
void StreamlineRenderer::genSeedPoints(protein_calls::VTIDataCall* vti, float zClip, float isoval) {


    float posZ = vti->GetOrigin().GetZ() + zClip; // Start above the lower boundary

    float xMax = vti->GetOrigin().GetX() + vti->GetSpacing().GetX() * (vti->GetGridsize().GetX() - 1);
    float yMax = vti->GetOrigin().GetY() + vti->GetSpacing().GetY() * (vti->GetGridsize().GetY() - 1);
    //float zMax = vti->GetOrigin().GetZ() + vtiCall->GetSpacing().GetZ()*(vtiCall->GetGridsize().GetZ()-1);
    float xMin = vti->GetOrigin().GetX();
    float yMin = vti->GetOrigin().GetY();

    // Initialize random seed
    srand(static_cast<unsigned int>(time(NULL)));
    this->seedPoints.SetCount(0);
    //for (size_t cnt = 0; cnt < this->nStreamlines; ++cnt) {
    while (this->seedPoints.Count() / 3 < this->nStreamlines) {
        Vec3f pos;
        pos.SetX(vti->GetOrigin().GetX() + (float(rand() % 10000) / 10000.0f) * (xMax - xMin));
        pos.SetY(vti->GetOrigin().GetY() + (float(rand() % 10000) / 10000.0f) * (yMax - yMin));
        pos.SetZ(posZ + (float(rand() % 10000) / 10000.0f) * (10));
        //        pos.SetZ(posZ);
        //printf("Random pos %f %f %f\n", pos.GetX(), pos.GetY(), pos.GetZ());

        float sample = this->sampleFieldAtPosTrilin(
            vti, make_float3(pos.GetX(), pos.GetY(), pos.GetZ()), (float*)vti->GetPointDataByIdx(0, 0));

        // Sample density value
        //if (vislib::math::Abs(sample - isoval) < 0.05) {
        if ((sample - isoval) > 0.00) {
            this->seedPoints.Add(pos.GetX());
            this->seedPoints.Add(pos.GetY());
            this->seedPoints.Add(pos.GetZ());
        }
    }
}


/*
 * StreamlineRenderer::sampleFieldAtPosTrilin
 */
float StreamlineRenderer::sampleFieldAtPosTrilin(protein_calls::VTIDataCall* vtiCall, float3 pos, float* field_D) {

    int3 c;
    float3 f;

    int3 gridSize_D =
        make_int3(vtiCall->GetGridsize().GetX(), vtiCall->GetGridsize().GetY(), vtiCall->GetGridsize().GetZ());
    float3 gridOrg_D =
        make_float3(vtiCall->GetOrigin().GetX(), vtiCall->GetOrigin().GetY(), vtiCall->GetOrigin().GetZ());
    float3 gridDelta_D =
        make_float3(vtiCall->GetSpacing().GetX(), vtiCall->GetSpacing().GetY(), vtiCall->GetSpacing().GetZ());

    // Get id of the cell containing the given position and interpolation
    // coefficients
    f.x = (pos.x - gridOrg_D.x) / gridDelta_D.x;
    f.y = (pos.y - gridOrg_D.y) / gridDelta_D.y;
    f.z = (pos.z - gridOrg_D.z) / gridDelta_D.z;
    c.x = (int)(f.x);
    c.y = (int)(f.y);
    c.z = (int)(f.z);
    f.x = f.x - (float)c.x; // alpha
    f.y = f.y - (float)c.y; // beta
    f.z = f.z - (float)c.z; // gamma

    c.x = vislib::math::Clamp(c.x, int(0), gridSize_D.x - 2);
    c.y = vislib::math::Clamp(c.y, int(0), gridSize_D.y - 2);
    c.z = vislib::math::Clamp(c.z, int(0), gridSize_D.z - 2);

    // Get values at corners of current cell
    float s[8];
    s[0] = field_D[gridSize_D.x * (gridSize_D.y * (c.z + 0) + (c.y + 0)) + c.x + 0];
    s[1] = field_D[gridSize_D.x * (gridSize_D.y * (c.z + 0) + (c.y + 0)) + c.x + 1];
    s[2] = field_D[gridSize_D.x * (gridSize_D.y * (c.z + 0) + (c.y + 1)) + c.x + 0];
    s[3] = field_D[gridSize_D.x * (gridSize_D.y * (c.z + 0) + (c.y + 1)) + c.x + 1];
    s[4] = field_D[gridSize_D.x * (gridSize_D.y * (c.z + 1) + (c.y + 0)) + c.x + 0];
    s[5] = field_D[gridSize_D.x * (gridSize_D.y * (c.z + 1) + (c.y + 0)) + c.x + 1];
    s[6] = field_D[gridSize_D.x * (gridSize_D.y * (c.z + 1) + (c.y + 1)) + c.x + 0];
    s[7] = field_D[gridSize_D.x * (gridSize_D.y * (c.z + 1) + (c.y + 1)) + c.x + 1];

    // Use trilinear interpolation to sample the volume
    return protein_calls::Interpol::Trilin(s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7], f.x, f.y, f.z);
}


/*
 * StreamlineRenderer::requestPlane
 */
bool StreamlineRenderer::requestPlane(Call& call) {
    view::CallClipPlane* ccp = dynamic_cast<view::CallClipPlane*>(&call);

    if (ccp == NULL)
        return false;

    this->clipPlane.Set(vislib::math::Point<float, 3>(0.0, 0.0, this->seedClipZ * this->bbox.ObjectSpaceBBox().Depth()),
        Vec3f(0.0, 0.0, 1.0));
    this->clipPlane.Distance(vislib::math::Point<float, 3>(0.0f, 0.0f, 0.0f));

    ccp->SetColour(0, 0, 0, 0);
    ccp->SetPlane(this->clipPlane);

    return true;
}


/*
 * StreamlineRenderer::updateParams
 */
void StreamlineRenderer::updateParams() {

    /* Streamline integration parameters */

    // Set the number of steps for streamline integration
    if (this->nStreamlinesSlot.IsDirty()) {
        this->nStreamlines = this->nStreamlinesSlot.Param<core::param::IntParam>()->Value();
        this->nStreamlinesSlot.ResetDirty();
        this->triggerComputeStreamlines = true;
    }

    // Set the number of steps for streamline integration
    if (this->streamlineMaxStepsSlot.IsDirty()) {
        this->streamlineMaxSteps = this->streamlineMaxStepsSlot.Param<core::param::IntParam>()->Value();
        this->streamlineMaxStepsSlot.ResetDirty();
        this->triggerComputeStreamlines = true;
    }

    // Set the step size for streamline integration
    if (this->streamlineStepSlot.IsDirty()) {
        this->streamlineStep = this->streamlineStepSlot.Param<core::param::FloatParam>()->Value();
        this->streamlineStepSlot.ResetDirty();
        this->triggerComputeStreamlines = true;
    }

    // Set the epsilon for the streamline termination
    if (this->seedClipZSlot.IsDirty()) {
        this->seedClipZ = this->seedClipZSlot.Param<core::param::FloatParam>()->Value();
        this->seedClipZSlot.ResetDirty();
        this->triggerComputeStreamlines = true;
    }

    // Set the epsilon for the streamline termination
    if (this->seedIsoSlot.IsDirty()) {
        this->seedIso = this->seedIsoSlot.Param<core::param::FloatParam>()->Value();
        this->seedIsoSlot.ResetDirty();
        this->triggerComputeStreamlines = true;
    }


    /* Streamlines render mode */

    // parameter refresh
    if (this->renderModeSlot.IsDirty()) {
        this->renderMode = static_cast<RenderModes>(int(this->renderModeSlot.Param<param::EnumParam>()->Value()));
        this->renderModeSlot.ResetDirty();
    }

    // Set the streamtubes thickness slot
    if (this->streamtubesThicknessSlot.IsDirty()) {
        this->streamtubesThickness = this->streamtubesThicknessSlot.Param<core::param::FloatParam>()->Value();
        this->streamtubesThicknessSlot.ResetDirty();
    }

    // Set minimum color value
    if (this->minColSlot.IsDirty()) {
        this->minCol = this->minColSlot.Param<core::param::FloatParam>()->Value();
        this->minColSlot.ResetDirty();
    }

    // Set maximum color value
    if (this->maxColSlot.IsDirty()) {
        this->maxCol = this->maxColSlot.Param<core::param::FloatParam>()->Value();
        this->maxColSlot.ResetDirty();
    }
}
