/*
 * GLSLVolumeRenderer.cpp
 *
 * Copyright (C) 2008 by Universitaet Stuttgart (VISUS).
 * Alle Rechte vorbehalten.
 */


#define _USE_MATH_DEFINES 1

#include <cmath>
#include <ctime>
#include <iostream>

#include <glm/gtc/matrix_transform.hpp>

#include "GLSLVolumeRenderer.h"
#include "mmcore/CoreInstance.h"
#include "mmcore/param/BoolParam.h"
#include "mmcore/param/ColorParam.h"
#include "mmcore/param/EnumParam.h"
#include "mmcore/param/FloatParam.h"
#include "mmcore/param/IntParam.h"
#include "mmcore/param/StringParam.h"
#include "mmcore/param/Vector3fParam.h"
#include "mmcore/utility/ColourParser.h"
#include "mmcore_gl/utility/ShaderFactory.h"
#include "mmstd/renderer/AbstractCallRender.h"
#include "protein_calls/ProteinColor.h"
#include "vislib/OutOfRangeException.h"
#include "vislib/String.h"
#include "vislib/StringConverter.h"
#include "vislib/Trace.h"
#include "vislib/assert.h"
#include "vislib/math/Point.h"
#include "vislib/math/Quaternion.h"
#include "vislib/sys/ASCIIFileBuffer.h"
#include "vislib/sys/File.h"
#include "vislib_gl/graphics/gl/IncludeAllGL.h"
#include "vislib_gl/graphics/gl/glverify.h"


using namespace megamol;
using namespace megamol::core;
using namespace megamol::protein_gl;
using namespace megamol::protein_calls;
using namespace megamol::core::utility::log;


/*
 * protein::GLSLVolumeRenderer::GLSLVolumeRenderer (CTOR)
 */
protein_gl::GLSLVolumeRenderer::GLSLVolumeRenderer(void)
        : mmstd_gl::Renderer3DModuleGL()
        , protDataCallerSlot("getData", "Connects the volume rendering with data storage")
        , protRendererCallerSlot("renderProtein", "Connects the volume rendering with a protein renderer")
        , coloringModeParam("coloringMode", "Coloring Mode")
        , volIsoValueParam("volIsoValue", "Isovalue for isosurface rendering")
        , volFilterRadiusParam("volFilterRadius", "Filter Radius for volume generation")
        , volDensityScaleParam("volDensityScale", "Density scale factor for volume generation")
        , volIsoOpacityParam("volIsoOpacity", "Opacity of isosurface")
        , volClipPlaneFlagParam("volClipPlane", "Enable volume clipping")
        , volClipPlane0NormParam("clipPlane0Norm", "Volume clipping plane 0 normal")
        , volClipPlane0DistParam("clipPlane0Dist", "Volume clipping plane 0 distance")
        , volClipPlaneOpacityParam("clipPlaneOpacity", "Volume clipping plane opacity")
        , colorTableFileParam("colorTableFilename", "The filename of the color table.")
        , minGradColorParam("minGradColor", "The color for the minimum value for gradient coloring")
        , midGradColorParam("midGradColor", "The color for the middle value for gradient coloring")
        , maxGradColorParam("maxGradColor", "The color for the maximum value for gradient coloring")
        , renderVolumeParam("renderVolumeFlag", "Render Volume")
        , renderProteinParam("renderProteinFlag", "Render Protein")
        , currentFrameId(0)
        , atomCount(0)
        , proteinFBO(nullptr)
        , volumeTex(0)
        , volumeSize(256)
        , volFBO(0)
        , volFilterRadius(1.75f)
        , volDensityScale(1.0f)
        , width(0)
        , height(0)
        , volRayTexWidth(0)
        , volRayTexHeight(0)
        , volRayStartTex(0)
        , volRayLengthTex(0)
        , volRayDistTex(0)
        , renderIsometric(true)
        , isoValue(0.5f)
        , volIsoOpacity(0.4f)
        , volClipPlaneFlag(false)
        , volClipPlaneOpacity(0.4f)
        , forceUpdateVolumeTexture(true)
        , renderUtils() {
    this->protDataCallerSlot.SetCompatibleCall<MolecularDataCallDescription>();
    this->MakeSlotAvailable(&this->protDataCallerSlot);

    // set renderer caller slot
    this->protRendererCallerSlot.SetCompatibleCall<mmstd_gl::CallRender3DGLDescription>();
    this->MakeSlotAvailable(&this->protRendererCallerSlot);

    // --- set the coloring mode ---
    this->SetColoringMode(protein_calls::ProteinColor::ColoringMode::ELEMENT);
    param::EnumParam* cm = new param::EnumParam(int(this->currentColoringMode));
    MolecularDataCall* mol = new MolecularDataCall();
    unsigned int cCnt;
    protein_calls::ProteinColor::ColoringMode cMode;
    for (cCnt = 0; cCnt < static_cast<uint32_t>(ProteinColor::ColoringMode::MODE_COUNT); ++cCnt) {
        cMode = static_cast<ProteinColor::ColoringMode>(cCnt);
        cm->SetTypePair(static_cast<int>(cMode), ProteinColor::GetName(cMode).c_str());
    }
    delete mol;
    this->coloringModeParam << cm;
    this->MakeSlotAvailable(&this->coloringModeParam);

    // --- set up parameters for isovalue ---
    this->volIsoValueParam.SetParameter(new param::FloatParam(this->isoValue));
    this->MakeSlotAvailable(&this->volIsoValueParam);
    // --- set up parameter for volume filter radius ---
    this->volFilterRadiusParam.SetParameter(new param::FloatParam(this->volFilterRadius, 0.0f));
    this->MakeSlotAvailable(&this->volFilterRadiusParam);
    // --- set up parameter for volume density scale ---
    this->volDensityScaleParam.SetParameter(new param::FloatParam(this->volDensityScale, 0.0f));
    this->MakeSlotAvailable(&this->volDensityScaleParam);
    // --- set up parameter for isosurface opacity ---
    this->volIsoOpacityParam.SetParameter(new param::FloatParam(this->volIsoOpacity, 0.0f, 1.0f));
    this->MakeSlotAvailable(&this->volIsoOpacityParam);

    // set default clipping plane
    this->volClipPlane.Clear();
    this->volClipPlane.Add(vislib::math::Vector<double, 4>(0.0, 1.0, 0.0, 0.0));

    // --- set up parameter for volume clipping ---
    this->volClipPlaneFlagParam.SetParameter(new param::BoolParam(this->volClipPlaneFlag));
    this->MakeSlotAvailable(&this->volClipPlaneFlagParam);
    // --- set up parameter for volume clipping plane normal ---
    vislib::math::Vector<float, 3> cp0n(static_cast<float>(this->volClipPlane[0].PeekComponents()[0]),
        static_cast<float>(this->volClipPlane[0].PeekComponents()[1]),
        static_cast<float>(this->volClipPlane[0].PeekComponents()[2]));
    this->volClipPlane0NormParam.SetParameter(new param::Vector3fParam(cp0n));
    this->MakeSlotAvailable(&this->volClipPlane0NormParam);
    // --- set up parameter for volume clipping plane distance ---
    float d = static_cast<float>(this->volClipPlane[0].PeekComponents()[3]);
    this->volClipPlane0DistParam.SetParameter(new param::FloatParam(d, 0.0f, 1.0f));
    this->MakeSlotAvailable(&this->volClipPlane0DistParam);
    // --- set up parameter for clipping plane opacity ---
    this->volClipPlaneOpacityParam.SetParameter(new param::FloatParam(this->volClipPlaneOpacity, 0.0f, 1.0f));
    this->MakeSlotAvailable(&this->volClipPlaneOpacityParam);

    // fill color table with default values and set the filename param
    std::string filename("colors.txt");
    ProteinColor::ReadColorTableFromFile(filename, this->fileColorTable);
    this->colorTableFileParam.SetParameter(new param::StringParam(filename));
    this->MakeSlotAvailable(&this->colorTableFileParam);

    // the color for the minimum value (gradient coloring
    this->minGradColorParam.SetParameter(new param::ColorParam("#146496"));
    this->MakeSlotAvailable(&this->minGradColorParam);

    // the color for the middle value (gradient coloring
    this->midGradColorParam.SetParameter(new param::ColorParam("#f0f0f0"));
    this->MakeSlotAvailable(&this->midGradColorParam);

    // the color for the maximum value (gradient coloring
    this->maxGradColorParam.SetParameter(new param::ColorParam("#ae3b32"));
    this->MakeSlotAvailable(&this->maxGradColorParam);

    // --- set up parameter for volume rendering ---
    this->renderVolumeParam.SetParameter(new param::BoolParam(true));
    this->MakeSlotAvailable(&this->renderVolumeParam);

    // --- set up parameter for protein rendering ---
    this->renderProteinParam.SetParameter(new param::BoolParam(true));
    this->MakeSlotAvailable(&this->renderProteinParam);

    // fill rainbow color table
    ProteinColor::MakeRainbowColorTable(100, this->rainbowColors);

    // initialize vertex and color array for tests
    c = new float[3 * NUM];
    p = new float[4 * NUM];
    for (int cnt = 0; cnt < NUM; cnt++) {
        c[cnt * 3 + 0] = 1.0f - float(cnt) / float(NUM);
        c[cnt * 3 + 1] = 0.0f;
        c[cnt * 3 + 2] = 0.0f;
        p[cnt * 4 + 0] = 0.0f;
        p[cnt * 4 + 1] = 0.0f;
        p[cnt * 4 + 2] = 0.0f;
        p[cnt * 4 + 3] = 1.0f;
    }
}


/*
 * protein::GLSLVolumeRenderer::~GLSLVolumeRenderer (DTOR)
 */
protein_gl::GLSLVolumeRenderer::~GLSLVolumeRenderer(void) {
    this->Release();
}


/*
 * protein::GLSLVolumeRenderer::release
 */
void protein_gl::GLSLVolumeRenderer::release(void) {}


/*
 * protein::GLSLVolumeRenderer::create
 */
bool protein_gl::GLSLVolumeRenderer::create(void) {

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    glEnable(GL_VERTEX_PROGRAM_TWO_SIDE);
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);

    using namespace vislib_gl::graphics::gl;

    auto const shaderOptions = msf::ShaderFactoryOptionsOpenGL(GetCoreInstance()->GetShaderPaths());

    try {
        updateVolumeShader = core::utility::make_glowl_shader("updateVolumeShader", shaderOptions,
            "protein_gl/volume/updateVolume.vert.glsl", "protein_gl/volume/updateVolume.frag.glsl");

        volRayStartShader = core::utility::make_glowl_shader("volRayStartShader", shaderOptions,
            "protein_gl/volume/rayStart.vert.glsl", "protein_gl/volume/rayStart.frag.glsl");

        volRayStartEyeShader = core::utility::make_glowl_shader("volRayStartEyeShader", shaderOptions,
            "protein_gl/volume/rayStartEye.vert.glsl", "protein_gl/volume/rayStartEye.frag.glsl");

        volRayLengthShader = core::utility::make_glowl_shader("volRayLengthShader", shaderOptions,
            "protein_gl/volume/rayStart.vert.glsl", "protein_gl/volume/rayLength.frag.glsl");

        volumeShader = core::utility::make_glowl_shader(
            "volumeShader", shaderOptions, "protein_gl/volume/volume.vert.glsl", "protein_gl/volume/volume.frag.glsl");

        colorWriterShader = core::utility::make_glowl_shader("colorWriterShader", shaderOptions,
            "protein_gl/volume/colorWriter.vert.glsl", "protein_gl/volume/colorWriter.frag.glsl");

    } catch (std::exception& e) {
        Log::DefaultLog.WriteError(("GLSLVolumeRenderer: " + std::string(e.what())).c_str());
        return false;
    }

    // Initialize render utils
    if (!renderUtils.InitPrimitiveRendering(GetCoreInstance()->GetShaderPaths())) {
        megamol::core::utility::log::Log::DefaultLog.WriteError(
            "Couldn't initialize primitive rendering. [%s, %s, line %d]\n", __FILE__, __FUNCTION__, __LINE__);
        return false;
    }

    return true;
}


/**********************************************************************
 * 'render'-functions
 **********************************************************************/


/*
 * protein::ProteinRenderer::GetExtents
 */
bool protein_gl::GLSLVolumeRenderer::GetExtents(mmstd_gl::CallRender3DGL& call) {
    mmstd_gl::CallRender3DGL* cr3d = dynamic_cast<mmstd_gl::CallRender3DGL*>(&call);
    if (cr3d == NULL)
        return false;

    MolecularDataCall* mol = this->protDataCallerSlot.CallAs<MolecularDataCall>();

    float scale, xoff, yoff, zoff;
    vislib::math::Cuboid<float> boundingBox;
    vislib::math::Point<float, 3> bbc;

    // try to get the bounding box from the active data call
    if (mol) {
        // try to call the molecular data
        if (!(*mol)(MolecularDataCall::CallForGetExtent))
            return false;
        // get the bounding box
        boundingBox = mol->AccessBoundingBoxes().ObjectSpaceBBox();
        // set the frame count
        cr3d->SetTimeFramesCount(mol->FrameCount());

    } else {
        return false;
    }

    bbc = boundingBox.CalcCenter();
    xoff = -bbc.X();
    yoff = -bbc.Y();
    zoff = -bbc.Z();
    if (!vislib::math::IsEqual(boundingBox.LongestEdge(), 0.0f)) {
        scale = 2.0f / boundingBox.LongestEdge();
    } else {
        scale = 1.0f;
    }

    BoundingBoxes_2& bbox = cr3d->AccessBoundingBoxes();
    bbox.SetBoundingBox(boundingBox);

    // get the pointer to CallRender3DGL (protein renderer)
    mmstd_gl::CallRender3DGL* protrencr3d = this->protRendererCallerSlot.CallAs<mmstd_gl::CallRender3DGL>();
    vislib::math::Point<float, 3> protrenbbc;
    if (protrencr3d) {
        (*protrencr3d)(core::view::AbstractCallRender::FnGetExtents);
        BoundingBoxes_2& protrenbb = protrencr3d->AccessBoundingBoxes();
        this->protrenScale = protrenbb.BoundingBox().Width() / boundingBox.Width();
        // this->protrenTranslate = ( protrenbb.ObjectSpaceBBox().CalcCenter() - bbc) * scale;
        if (mol) {
            this->protrenTranslate.Set(xoff, yoff, zoff);
            this->protrenTranslate *= scale;
        } else {
            this->protrenTranslate = (protrenbb.BoundingBox().CalcCenter() - bbc) * scale;
        }
    }
    return true;
}

/*
 * protein::GLSLVolumeRenderer::Render
 */
bool protein_gl::GLSLVolumeRenderer::Render(mmstd_gl::CallRender3DGL& call) {
#if 0
    // generate volume, if necessary
    if( !glIsTexture( this->volumeTex) ) {
        // from CellVis: cellVis.cpp, initGL
        glGenTextures( 1, &this->volumeTex);
        glBindTexture( GL_TEXTURE_3D, this->volumeTex);
        glTexImage3D( GL_TEXTURE_3D, 0,
                      //GL_LUMINANCE32F_ARB,
                      GL_RGBA16F,
                      this->volumeSize, this->volumeSize, this->volumeSize, 0,
                      //GL_LUMINANCE, GL_FLOAT, 0);
                      GL_RGBA, GL_FLOAT, 0);
        GLint param = GL_LINEAR;
        glTexParameteri( GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, param);
        glTexParameteri( GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, param);
        GLint mode = GL_CLAMP_TO_EDGE;
        //GLint mode = GL_REPEAT;
        glTexParameteri( GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, mode);
        glTexParameteri( GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, mode);
        glTexParameteri( GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, mode);
        glBindTexture( GL_TEXTURE_3D, 0);
        CHECK_FOR_OGL_ERROR();
    }
    // generate FBO, if necessary
    if( !glIsFramebufferEXT( this->volFBO ) ) {
        glGenFramebuffersEXT( 1, &this->volFBO);
        CHECK_FOR_OGL_ERROR();
    }

    // counter variable
    unsigned int z = 1;

    // store current frame buffer object ID
    GLint prevFBO;
    glGetIntegerv( GL_FRAMEBUFFER_BINDING_EXT, &prevFBO);

    glMatrixMode( GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glMatrixMode( GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    // store old viewport
    GLint viewport[4];
    glGetIntegerv( GL_VIEWPORT, viewport);
    // set viewport
    glViewport( 0, 0, this->volumeSize, this->volumeSize);

    glDisable(GL_DEPTH_TEST);

    glBindFramebufferEXT( GL_FRAMEBUFFER_EXT, this->volFBO);
    glFramebufferTexture2DEXT( GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);

    glDisable( GL_BLEND);
    // clear 3d texture
    this->colorWriterShader.Enable();
    glUniform4f( this->colorWriterShader.ParameterLocation( "color"), 0.1f, 0.1f, 0.1f, 0.0f);
    for( z = 0; z < this->volumeSize; ++z) {
        // attach texture slice to FBO
        glFramebufferTexture3DEXT( GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0, GL_TEXTURE_3D, this->volumeTex, 0, z);
        glRecti(-1, -1, 1, 1);
    }
    this->colorWriterShader.Disable();
    
    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ONE);

    // ========================================================================
    //  TEST RENDERING AND BLENDING PERFORMANCE
    // ========================================================================
    glPointSize( 100.0f);
    glEnable(GL_POINT_SIZE);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
    glVertexPointer( 4, GL_FLOAT, 0, p);
    glColorPointer( 3, GL_FLOAT, 0, c);
    for( z = 0; z < this->volumeSize; z++ ) {
        // attach texture slice to FBO
        glFramebufferTexture3DEXT( GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_3D, this->volumeTex, 0, z);
        // draw
        glDrawArrays( GL_POINTS, 0, NUM);
    }
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
    glDisable(GL_POINT_SIZE);
    
    // restore viewport
    glViewport( viewport[0], viewport[1], viewport[2], viewport[3]);

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();

    glBindFramebufferEXT( GL_FRAMEBUFFER_EXT, prevFBO);

    glDisable( GL_BLEND);
    glEnable( GL_DEPTH_TEST);
    glDepthMask( GL_TRUE);
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    return true;
#else
    // cast the call to Render3D
    mmstd_gl::CallRender3DGL* cr3d = dynamic_cast<mmstd_gl::CallRender3DGL*>(&call);
    if (!cr3d)
        return false;
    // get the pointer to CallRender3DGL (protein renderer)
    mmstd_gl::CallRender3DGL* protrencr3d = this->protRendererCallerSlot.CallAs<mmstd_gl::CallRender3DGL>();
    // get pointer to MolecularDataCall
    MolecularDataCall* mol = this->protDataCallerSlot.CallAs<MolecularDataCall>();

    // get camera information
    cameraInfo = cr3d->GetCamera();

    // =============== Query Camera View Dimensions ===============
    if (static_cast<unsigned int>(cr3d->GetFramebuffer()->getWidth()) != this->width ||
        static_cast<unsigned int>(cr3d->GetFramebuffer()->getHeight()) != this->height) {
        this->width = static_cast<unsigned int>(cr3d->GetFramebuffer()->getWidth());
        this->height = static_cast<unsigned int>(cr3d->GetFramebuffer()->getHeight());
    }

    // (re)create the fbo, if necessary
    if (this->proteinFBO == nullptr ||
        (this->proteinFBO->getWidth() != this->width || this->proteinFBO->getHeight() != this->height)) {

        this->proteinFBO = std::make_shared<glowl::FramebufferObject>(this->width, this->height);
        this->proteinFBO->createColorAttachment(GL_RGBA16F, GL_RGBA, GL_FLOAT);
    }


    if (this->renderProteinParam.Param<param::BoolParam>()->Value()) {
        // =============== Protein Rendering ===============
        // start rendering to the FBO for protein rendering
        this->proteinFBO->bind();
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        if (protrencr3d) {
            // setup and call protein renderer
            glPushMatrix();
            glTranslatef(this->protrenTranslate.X(), this->protrenTranslate.Y(), this->protrenTranslate.Z());
            // glScalef( this->protrenScale, this->protrenScale, this->protrenScale);
            *protrencr3d = *cr3d;
            protrencr3d->SetFramebuffer(this->proteinFBO); // TODO: Handle incoming buffers!
            (*protrencr3d)();
            glPopMatrix();
        }
        // stop rendering to the FBO for protein rendering by re-enabling the output buffer
        cr3d->GetFramebuffer()->bind();
    }

    // =============== Refresh all parameters ===============
    this->ParameterRefresh(cr3d);

    // get the call time
    this->callTime = cr3d->Time();

    if (mol) {
        // set frame ID and call data
        mol->SetFrameID(static_cast<int>(callTime));
        if (!(*mol)(MolecularDataCall::CallForGetData))
            return false;
        // check if atom count is zero
        if (mol->AtomCount() == 0)
            return true;

        this->colorLookupTable = {glm::make_vec3(this->minGradColorParam.Param<param::ColorParam>()->Value().data()),
            glm::make_vec3(this->midGradColorParam.Param<param::ColorParam>()->Value().data()),
            glm::make_vec3(this->maxGradColorParam.Param<param::ColorParam>()->Value().data())};

        ProteinColor::MakeColorTable(*mol, this->currentColoringMode, this->atomColorTable, this->colorLookupTable,
            this->fileColorTable, this->rainbowColors, nullptr, nullptr, true);
    }


    unsigned int cpCnt;
    for (cpCnt = 0; cpCnt < this->volClipPlane.Count(); ++cpCnt) {
        glClipPlane(GL_CLIP_PLANE0 + cpCnt, this->volClipPlane[cpCnt].PeekComponents());
    }

    // =============== Volume Rendering ===============
    bool retval = false;
    // try to start volume rendering using protein data
    if (mol) {
        retval = this->RenderMolecularData(cr3d, mol);
    }

    // unlock the current frame
    if (mol)
        mol->Unlock();

    return retval;
#endif
}


/*
 * Volume rendering using molecular data.
 */
bool protein_gl::GLSLVolumeRenderer::RenderMolecularData(mmstd_gl::CallRender3DGL* call, MolecularDataCall* mol) {

    // check last atom count with current atom count
    if (this->atomCount != mol->AtomCount()) {
        this->atomCount = mol->AtomCount();
        this->forceUpdateVolumeTexture = true;
    }

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);

    glPushMatrix();

    // translate scene for volume ray casting
    if (!vislib::math::IsEqual(mol->AccessBoundingBoxes().ObjectSpaceBBox().LongestEdge(), 0.0f)) {
        this->scale = 2.0f / mol->AccessBoundingBoxes().ObjectSpaceBBox().LongestEdge();
    } else {
        this->scale = 1.0f;
    }
    vislib::math::Vector<float, 3> trans(mol->AccessBoundingBoxes().ObjectSpaceBBox().GetSize().PeekDimension());
    trans *= -this->scale * 0.5f;
    glTranslatef(trans.GetX(), trans.GetY(), trans.GetZ());

    // ------------------------------------------------------------
    // --- Volume Rendering                                     ---
    // --- update & render the volume                           ---
    // ------------------------------------------------------------

    vislib::StringA paramSlotName(call->PeekCallerSlot()->Parent()->FullName());
    paramSlotName += "::anim::play";
    param::ParamSlot* paramSlot = dynamic_cast<param::ParamSlot*>(this->FindNamedObject(paramSlotName, true).get());
    if (paramSlot->Param<param::BoolParam>()->Value() || this->forceUpdateVolumeTexture) {
        this->UpdateVolumeTexture(mol);
        this->forceUpdateVolumeTexture = false;
    }

    if (this->renderProteinParam.Param<param::BoolParam>()->Value()) {
        call->GetFramebuffer()->bind();
        // Create 2D orthographic mvp matrix
        glm::mat4 ortho = glm::ortho(0.0f, 1.0f, 0.0f, 1.0f, -1.0f, 1.0f);
        glm::vec3 pos_bottom_left = {0.0, 0.0, 1.0f};
        glm::vec3 pos_upper_left = {0.0, 1.0, 1.0f};
        glm::vec3 pos_upper_right = {1.0, 1.0, 1.0f};
        glm::vec3 pos_bottom_right = {1.0, 0.0, 1.0f};
        renderUtils.Push2DColorTexture(this->proteinFBO->getColorAttachment(0)->getName(), pos_bottom_left,
            pos_upper_left, pos_upper_right, pos_bottom_right);
        renderUtils.DrawTextures(ortho, {call->GetFramebuffer()->getWidth(), call->GetFramebuffer()->getHeight()});
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
    }

    unsigned int cpCnt;
    if (this->volClipPlaneFlag) {
        for (cpCnt = 0; cpCnt < this->volClipPlane.Count(); ++cpCnt) {
            glEnable(GL_CLIP_PLANE0 + cpCnt);
        }
    }

    if (this->renderVolumeParam.Param<param::BoolParam>()->Value())
        this->RenderVolume(mol->AccessBoundingBoxes().ObjectSpaceBBox());

    if (this->volClipPlaneFlag) {
        for (cpCnt = 0; cpCnt < this->volClipPlane.Count(); ++cpCnt) {
            glDisable(GL_CLIP_PLANE0 + cpCnt);
        }
    }

    glDisable(GL_VERTEX_PROGRAM_POINT_SIZE);

    glDisable(GL_DEPTH_TEST);

    glPopMatrix();

    return true;
}


/*
 * refresh parameters
 */
void protein_gl::GLSLVolumeRenderer::ParameterRefresh(mmstd_gl::CallRender3DGL* call) {

    // parameter refresh
    if (this->coloringModeParam.IsDirty()) {
        this->SetColoringMode(
            static_cast<ProteinColor::ColoringMode>(int(this->coloringModeParam.Param<param::EnumParam>()->Value())));
        this->coloringModeParam.ResetDirty();
        this->forceUpdateVolumeTexture = true;
    }
    // volume parameters
    if (this->volIsoValueParam.IsDirty()) {
        this->isoValue = this->volIsoValueParam.Param<param::FloatParam>()->Value();
        this->volIsoValueParam.ResetDirty();
    }
    if (this->volFilterRadiusParam.IsDirty()) {
        this->volFilterRadius = this->volFilterRadiusParam.Param<param::FloatParam>()->Value();
        this->volFilterRadiusParam.ResetDirty();
        this->forceUpdateVolumeTexture = true;
    }
    if (this->volDensityScaleParam.IsDirty()) {
        this->volDensityScale = this->volDensityScaleParam.Param<param::FloatParam>()->Value();
        this->volDensityScaleParam.ResetDirty();
        this->forceUpdateVolumeTexture = true;
    }
    if (this->volIsoOpacityParam.IsDirty()) {
        this->volIsoOpacity = this->volIsoOpacityParam.Param<param::FloatParam>()->Value();
        this->volIsoOpacityParam.ResetDirty();
    }
    if (this->volClipPlaneFlagParam.IsDirty()) {
        this->volClipPlaneFlag = this->volClipPlaneFlagParam.Param<param::BoolParam>()->Value();
        this->volClipPlaneFlagParam.ResetDirty();
    }

    // get current clip plane normal
    vislib::math::Vector<float, 3> cp0n((float)this->volClipPlane[0].PeekComponents()[0],
        (float)this->volClipPlane[0].PeekComponents()[1], (float)this->volClipPlane[0].PeekComponents()[2]);
    // get current clip plane distance
    float cp0d = (float)this->volClipPlane[0].PeekComponents()[3];

    // check clip plane normal parameter
    if (this->volClipPlane0NormParam.IsDirty()) {
        // overwrite clip plane normal
        cp0n = this->volClipPlane0NormParam.Param<param::Vector3fParam>()->Value();
        // normalize clip plane normal, if necessary and set normalized clip plane normal to parameter
        if (!vislib::math::IsEqual<float>(cp0n.Length(), 1.0f)) {
            cp0n.Normalise();
            this->volClipPlane0NormParam.Param<param::Vector3fParam>()->SetValue(cp0n);
        }
        this->volClipPlane0NormParam.ResetDirty();
    }
    // compute maximum extent
    vislib::math::Cuboid<float> bbox(call->AccessBoundingBoxes().BoundingBox());
    vislib::math::Vector<float, 3> tmpVec;
    float d, maxD, minD;
    // 1
    tmpVec.Set(bbox.GetLeftBottomBack().X(), bbox.GetLeftBottomBack().Y(), bbox.GetLeftBottomBack().Z());
    maxD = minD = cp0n.Dot(tmpVec);
    // 2
    tmpVec.Set(bbox.GetRightBottomBack().X(), bbox.GetRightBottomBack().Y(), bbox.GetRightBottomBack().Z());
    d = cp0n.Dot(tmpVec);
    if (minD > d)
        minD = d;
    if (maxD < d)
        maxD = d;
    // 3
    tmpVec.Set(bbox.GetLeftBottomFront().X(), bbox.GetLeftBottomFront().Y(), bbox.GetLeftBottomFront().Z());
    d = cp0n.Dot(tmpVec);
    if (minD > d)
        minD = d;
    if (maxD < d)
        maxD = d;
    // 4
    tmpVec.Set(bbox.GetRightBottomFront().X(), bbox.GetRightBottomFront().Y(), bbox.GetRightBottomFront().Z());
    d = cp0n.Dot(tmpVec);
    if (minD > d)
        minD = d;
    if (maxD < d)
        maxD = d;
    // 5
    tmpVec.Set(bbox.GetLeftTopBack().X(), bbox.GetLeftTopBack().Y(), bbox.GetLeftTopBack().Z());
    d = cp0n.Dot(tmpVec);
    if (minD > d)
        minD = d;
    if (maxD < d)
        maxD = d;
    // 6
    tmpVec.Set(bbox.GetRightTopBack().X(), bbox.GetRightTopBack().Y(), bbox.GetRightTopBack().Z());
    d = cp0n.Dot(tmpVec);
    if (minD > d)
        minD = d;
    if (maxD < d)
        maxD = d;
    // 7
    tmpVec.Set(bbox.GetLeftTopFront().X(), bbox.GetLeftTopFront().Y(), bbox.GetLeftTopFront().Z());
    d = cp0n.Dot(tmpVec);
    if (minD > d)
        minD = d;
    if (maxD < d)
        maxD = d;
    // 8
    tmpVec.Set(bbox.GetRightTopFront().X(), bbox.GetRightTopFront().Y(), bbox.GetRightTopFront().Z());
    d = cp0n.Dot(tmpVec);
    if (minD > d)
        minD = d;
    if (maxD < d)
        maxD = d;
    // check clip plane distance
    if (this->volClipPlane0DistParam.IsDirty()) {
        cp0d = this->volClipPlane0DistParam.Param<param::FloatParam>()->Value();
        cp0d = minD + (maxD - minD) * cp0d;
        this->volClipPlane0DistParam.ResetDirty();
    }

    // set clip plane normal and distance to current clip plane
    this->volClipPlane[0].Set(cp0n.X(), cp0n.Y(), cp0n.Z(), cp0d);

    // check clip plane opacity parameter
    if (this->volClipPlaneOpacityParam.IsDirty()) {
        this->volClipPlaneOpacity = this->volClipPlaneOpacityParam.Param<param::FloatParam>()->Value();
        this->volClipPlaneOpacityParam.ResetDirty();
    }

    // update color table
    if (this->colorTableFileParam.IsDirty()) {
        ProteinColor::ReadColorTableFromFile(
            this->colorTableFileParam.Param<param::StringParam>()->Value(), this->colorLookupTable);
        this->colorTableFileParam.ResetDirty();
        this->forceUpdateVolumeTexture = true;
    }
}


/*
 * Create a volume containing all molecule atoms
 */
void protein_gl::GLSLVolumeRenderer::UpdateVolumeTexture(MolecularDataCall* mol) {
    // generate volume, if necessary
    if (!glIsTexture(this->volumeTex)) {
        // from CellVis: cellVis.cpp, initGL
        glGenTextures(1, &this->volumeTex);
        glBindTexture(GL_TEXTURE_3D, this->volumeTex);
        glTexImage3D(GL_TEXTURE_3D, 0,
            // GL_LUMINANCE32F_ARB,
            GL_RGBA16F, this->volumeSize, this->volumeSize, this->volumeSize, 0,
            // GL_LUMINANCE, GL_FLOAT, 0);
            GL_RGBA, GL_FLOAT, 0);
        GLint param = GL_LINEAR;
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, param);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, param);
        GLint mode = GL_CLAMP_TO_EDGE;
        // GLint mode = GL_REPEAT;
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, mode);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, mode);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, mode);
        glBindTexture(GL_TEXTURE_3D, 0);
    }
    // generate FBO, if necessary
    if (!glIsFramebufferEXT(this->volFBO)) {
        glGenFramebuffersEXT(1, &this->volFBO);
    }

    // counter variable
    unsigned int z;

    // store current frame buffer object ID
    GLint prevFBO;
    glGetIntegerv(GL_FRAMEBUFFER_BINDING_EXT, &prevFBO);

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    // store old viewport
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    // set viewport
    glViewport(0, 0, this->volumeSize, this->volumeSize);

    glDisable(GL_BLEND);
    glDisable(GL_DEPTH_TEST);

    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, this->volFBO);
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);

    // clear 3d texture
    this->colorWriterShader->use();
    glUniform4f(this->colorWriterShader->getUniformLocation("color"), 0.1f, 0.1f, 0.1f, 0.0f);
    for (z = 0; z < this->volumeSize; ++z) {
        // attach texture slice to FBO
        glFramebufferTexture3DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0, GL_TEXTURE_3D, this->volumeTex, 0, z);
        glRecti(-1, -1, 1, 1);
    }
    glUseProgram(0);

    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ONE);

// TEST BLENDING PERFORMANCE
#if 0
    for( z = 0; z < this->volumeSize; ++z ) {
        // attach texture slice to FBO
        glFramebufferTexture3DEXT( GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_3D, this->volumeTex, 0, z);
        glBegin( GL_POINTS);
        glVertex2f(1, 1);
        glEnd();
    }
#endif

#if 1
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);

    // scale[i] = 1/extent[i] --- extent = size of the bbox
    this->volScale[0] = 1.0f / (mol->AccessBoundingBoxes().ObjectSpaceBBox().Width() * this->scale);
    this->volScale[1] = 1.0f / (mol->AccessBoundingBoxes().ObjectSpaceBBox().Height() * this->scale);
    this->volScale[2] = 1.0f / (mol->AccessBoundingBoxes().ObjectSpaceBBox().Depth() * this->scale);
    // scaleInv = 1 / scale = extend
    this->volScaleInv[0] = 1.0f / this->volScale[0];
    this->volScaleInv[1] = 1.0f / this->volScale[1];
    this->volScaleInv[2] = 1.0f / this->volScale[2];

    this->updateVolumeShader->use();
    vislib::math::Vector<float, 3> orig(
        mol->AccessBoundingBoxes().ObjectSpaceBBox().GetLeftBottomBack().PeekCoordinates());
    orig = (orig + this->translation) * this->scale;
    vislib::math::Vector<float, 3> nullVec(0.0f, 0.0f, 0.0f);

    // set shader params
    glUniform1f(this->updateVolumeShader->getUniformLocation("filterRadius"), this->volFilterRadius);
    glUniform1f(this->updateVolumeShader->getUniformLocation("densityScale"), this->volDensityScale);
    glUniform3fv(this->updateVolumeShader->getUniformLocation("scaleVol"), 1, this->volScale);
    glUniform3fv(this->updateVolumeShader->getUniformLocation("scaleVolInv"), 1, this->volScaleInv);
    glUniform3f(this->updateVolumeShader->getUniformLocation("invVolRes"), 1.0f / float(this->volumeSize),
        1.0f / float(this->volumeSize), 1.0f / float(this->volumeSize));
    glUniform3fv(this->updateVolumeShader->getUniformLocation("translate"), 1, orig.PeekComponents());
    glUniform1f(this->updateVolumeShader->getUniformLocation("volSize"), float(this->volumeSize));

    float* atoms = new float[mol->AtomCount() * 4];
    int atomCnt;
#pragma omp parallel for
    for (atomCnt = 0; atomCnt < static_cast<int>(mol->AtomCount()); ++atomCnt) {
        atoms[atomCnt * 4 + 0] = (mol->AtomPositions()[3 * atomCnt + 0] + this->translation.X()) * this->scale;
        atoms[atomCnt * 4 + 1] = (mol->AtomPositions()[3 * atomCnt + 1] + this->translation.Y()) * this->scale;
        atoms[atomCnt * 4 + 2] = (mol->AtomPositions()[3 * atomCnt + 2] + this->translation.Z()) * this->scale;
        atoms[atomCnt * 4 + 3] = mol->AtomTypes()[mol->AtomTypeIndices()[atomCnt]].Radius() * this->scale;
    }

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
    glVertexPointer(4, GL_FLOAT, 0, atoms);
    glColorPointer(3, GL_FLOAT, 0, this->atomColorTable.data());
    for (z = 0; z < this->volumeSize; ++z) {
        // attach texture slice to FBO
        glFramebufferTexture3DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_3D, this->volumeTex, 0, z);
        glUniform1f(
            this->updateVolumeShader->getUniformLocation("sliceDepth"), (float(z) + 0.5f) / float(this->volumeSize));
        // set vertex and color pointers and draw them
        glDrawArrays(GL_POINTS, 0, mol->AtomCount());
    }
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);

    delete[] atoms;

    glUseProgram(0);
#endif

    // restore viewport
    glViewport(viewport[0], viewport[1], viewport[2], viewport[3]);

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();

    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, prevFBO);

    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    writeVolumeRAW();
}


/*
 * draw the volume
 */
void protein_gl::GLSLVolumeRenderer::RenderVolume(vislib::math::Cuboid<float> boundingbox) {
    const float stepWidth = 1.0f / (2.0f * float(this->volumeSize));
    glDisable(GL_BLEND);

    GLint prevFBO;
    glGetIntegerv(GL_FRAMEBUFFER_BINDING_EXT, &prevFBO);

    this->RayParamTextures(boundingbox);

    glDisable(GL_DEPTH_TEST);
    glDepthMask(GL_FALSE);

    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, prevFBO);

    // glUseProgram(_app->shader->volume->progId);
    this->volumeShader->use();

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // glUniform4fv(_app->shader->paramsCvolume.scaleVol, 1, vd->scale);
    glUniform4fv(this->volumeShader->getUniformLocation("scaleVol"), 1, this->volScale);
    // glUniform4fv(_app->shader->paramsCvolume.scaleVolInv, 1, vd->scaleInv);
    glUniform4fv(this->volumeShader->getUniformLocation("scaleVolInv"), 1, this->volScaleInv);
    // glUniform1f(_app->shader->paramsCvolume.stepSize, stepWidth);
    glUniform1f(this->volumeShader->getUniformLocation("stepSize"), stepWidth);

    // glUniform1f(_app->shader->paramsCvolume.alphaCorrection, _app->volStepSize/512.0f);
    // TODO: what is the correct value for volStepSize??
    glUniform1f(this->volumeShader->getUniformLocation("alphaCorrection"), this->volumeSize / 256.0f);
    // glUniform1i(_app->shader->paramsCvolume.numIterations, 255);
    glUniform1i(this->volumeShader->getUniformLocation("numIterations"), 255);
    // glUniform2f(_app->shader->paramsCvolume.screenResInv, 1.0f/_fboWidth, 1.0f/_fboHeight);
    glUniform2f(
        this->volumeShader->getUniformLocation("screenResInv"), 1.0f / float(this->width), 1.0f / float(this->height));

    // bind depth texture
    glUniform1i(this->volumeShader->getUniformLocation("volumeSampler"), 0);
    glUniform1i(this->volumeShader->getUniformLocation("transferRGBASampler"), 1);
    glUniform1i(this->volumeShader->getUniformLocation("rayStartSampler"), 2);
    glUniform1i(this->volumeShader->getUniformLocation("rayLengthSampler"), 3);

    glUniform1f(this->volumeShader->getUniformLocation("isoValue"), this->isoValue);
    glUniform1f(this->volumeShader->getUniformLocation("isoOpacity"), this->volIsoOpacity);
    glUniform1f(this->volumeShader->getUniformLocation("clipPlaneOpacity"), this->volClipPlaneOpacity);

    // transfer function
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_1D, 0);
    // ray start positions
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, this->volRayStartTex);
    // ray direction and length
    glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_2D, this->volRayLengthTex);

    // volume texture
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_3D, this->volumeTex);

    // draw a screen-filling quad
    glRectf(-1.0f, -1.0f, 1.0f, 1.0f);

    glUseProgram(0);

    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);
    glDisable(GL_BLEND);

    // restore depth buffer
    // glFramebufferTexture2DEXT( GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D,
    // this->proteinFBO.GetDepthTextureID(), 0);
}

/*
 * write the parameters of the ray to the textures
 */
void protein_gl::GLSLVolumeRenderer::RayParamTextures(vislib::math::Cuboid<float> boundingbox) {

    GLint param = GL_NEAREST;
    GLint mode = GL_CLAMP_TO_EDGE;

    // generate / resize ray start texture for volume ray casting
    if (!glIsTexture(this->volRayStartTex)) {
        glGenTextures(1, &this->volRayStartTex);
        glBindTexture(GL_TEXTURE_2D, this->volRayStartTex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, this->width, this->height, 0, GL_RGBA, GL_FLOAT, 0);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, param);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, param);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, mode);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, mode);
    } else if (this->width != this->volRayTexWidth || this->height != this->volRayTexHeight) {
        glBindTexture(GL_TEXTURE_2D, this->volRayStartTex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, this->width, this->height, 0, GL_RGBA, GL_FLOAT, 0);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, param);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, param);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, mode);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, mode);
    }
    // generate / resize ray length texture for volume ray casting
    if (!glIsTexture(this->volRayLengthTex)) {
        glGenTextures(1, &this->volRayLengthTex);
        glBindTexture(GL_TEXTURE_2D, this->volRayLengthTex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, this->width, this->height, 0, GL_RGBA, GL_FLOAT, 0);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, param);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, param);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, mode);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, mode);
    } else if (this->width != this->volRayTexWidth || this->height != this->volRayTexHeight) {
        glBindTexture(GL_TEXTURE_2D, this->volRayLengthTex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, this->width, this->height, 0, GL_RGBA, GL_FLOAT, 0);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, param);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, param);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, mode);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, mode);
    }
    // generate / resize ray distance texture for volume ray casting
    if (!glIsTexture(this->volRayDistTex)) {
        glGenTextures(1, &this->volRayDistTex);
        glBindTexture(GL_TEXTURE_2D, this->volRayDistTex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, this->width, this->height, 0, GL_RGBA, GL_FLOAT, 0);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, param);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, param);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, mode);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, mode);
    } else if (this->width != this->volRayTexWidth || this->height != this->volRayTexHeight) {
        glBindTexture(GL_TEXTURE_2D, this->volRayDistTex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, this->width, this->height, 0, GL_RGBA, GL_FLOAT, 0);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, param);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, param);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, mode);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, mode);
    }
    glBindTexture(GL_TEXTURE_2D, 0);
    // set vol ray dimensions
    this->volRayTexWidth = this->width;
    this->volRayTexHeight = this->height;

    GLuint db[2] = {GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1};
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, this->volFBO);

    // -------- ray start ------------
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, this->volRayStartTex, 0);
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, this->volRayDistTex, 0);

    // glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClear(GL_DEPTH_BUFFER_BIT);
    // draw to two rendertargets (the second does not need to be cleared)
    glDrawBuffers(2, db);
    // CHECK_FRAMEBUFFER_STATUS();

    // draw near clip plane
    glDisable(GL_DEPTH_TEST);
    glDepthMask(GL_FALSE);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    glColor4f(0.0f, 0.0f, 0.0f, 0.0f);

    // the shader transforms camera coords back to object space
    this->volRayStartEyeShader->use();

    float fovy = 0.0;
    float near_plane = 0.0;
    float far_plane = 0.0;
    try {
        auto cam_intrinsics = cameraInfo.get<megamol::core::view::Camera::PerspectiveParameters>();
        near_plane = cam_intrinsics.near_plane;
        far_plane = cam_intrinsics.far_plane;
        fovy = cam_intrinsics.fovy;
    } catch (...) {
        megamol::core::utility::log::Log::DefaultLog.WriteError(
            "GLSLVolumeRenderer - Error when getting perspective camera intrinsics");
    }

    float u = near_plane * tan(fovy);
    float r = (this->width / this->height) * u;

    glBegin(GL_QUADS);
    // glVertex3f(-r, -u, -_nearClip);
    glVertex3f(-r, -u, -near_plane);
    glVertex3f(r, -u, -near_plane);
    glVertex3f(r, u, -near_plane);
    glVertex3f(-r, u, -near_plane);
    glEnd();

    glUseProgram(0);

    glDrawBuffers(1, db);

    // glUseProgram(_app->shader->volRayStart->progId);
    this->volRayStartShader->use();

    // ------------ !useSphere && iso -------------
    vislib::math::Vector<float, 3> trans(boundingbox.GetSize().PeekDimension());
    trans *= this->scale * 0.5f;
    if (this->renderIsometric) {
        glUniform3f(this->volRayStartShader->getUniformLocation("translate"), 0.0f, 0.0f, 0.0f);
    } else {
        glUniform3fv(this->volRayStartShader->getUniformLocation("translate"), 1, trans.PeekComponents());
    }

    glDepthMask(GL_TRUE);
    glEnable(GL_DEPTH_TEST);

    glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_TRUE);
    glColor4f(0.0f, 0.0f, 0.0f, 1.0f);

    glEnable(GL_CULL_FACE);

    // draw nearest backfaces
    glCullFace(GL_FRONT);

    // enableClipPlanesVolume();

    // draw bBox
    this->DrawBoundingBox(boundingbox);

    // draw nearest frontfaces
    glCullFace(GL_BACK);
    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);

    // draw bBox
    this->DrawBoundingBox(boundingbox);

    glUseProgram(0);

    // --------------------------------
    // -------- ray length ------------
    // --------------------------------
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, this->volFBO);
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, this->volRayLengthTex, 0);

    // get clear color
    float clearCol[4];
    glGetFloatv(GL_COLOR_CLEAR_VALUE, clearCol);
    glClearColor(0, 0, 0, 0);
    glClearDepth(0.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearDepth(1.0f);
    glDrawBuffers(2, db);
    glClearColor(clearCol[0], clearCol[1], clearCol[2], clearCol[3]);

    // glUseProgram(_app->shader->volRayLength->progId);
    this->volRayLengthShader->use();

    glUniform1i(this->volRayLengthShader->getUniformLocation("sourceTex"), 0);
    glUniform1i(this->volRayLengthShader->getUniformLocation("depthTex"), 1);
    glUniform2f(this->volRayLengthShader->getUniformLocation("screenResInv"), 1.0f / float(this->width),
        1.0f / float(this->height));
    glUniform2f(this->volRayLengthShader->getUniformLocation("zNearFar"), near_plane, far_plane);

    if (this->renderIsometric) {
        glUniform3f(this->volRayLengthShader->getUniformLocation("translate"), 0.0f, 0.0f, 0.0f);
    } else {
        glUniform3fv(this->volRayLengthShader->getUniformLocation("translate"), 1, trans.PeekComponents());
    }
    glUniform1f(this->volRayLengthShader->getUniformLocation("scale"), this->scale);

    glActiveTexture(GL_TEXTURE1);
    // glBindTexture( GL_TEXTURE_2D, _depthTexId[0]);
    this->proteinFBO->bindDepthbuffer();
    glActiveTexture(GL_TEXTURE0);
    // glBindTexture( GL_TEXTURE_2D, _volRayStartTex);
    glBindTexture(GL_TEXTURE_2D, this->volRayStartTex);

    // draw farthest backfaces
    glCullFace(GL_FRONT);
    glDepthFunc(GL_GREATER);

    // draw bBox
    this->DrawBoundingBox(boundingbox);

    glUseProgram(0);

    glDrawBuffers(1, db);
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, 0, 0);

    glDepthFunc(GL_LESS);
    glCullFace(GL_BACK);
    glDisable(GL_CULL_FACE);

    // disableClipPlanes();

    // DEBUG check texture values
    /*
    float *texdata = new float[this->width*this->height];
    float max = 0.0f;
    memset( texdata, 0, sizeof(float)*(this->width*this->height));
    glBindTexture( GL_TEXTURE_2D, this->volRayLengthTex);
    glGetTexImage( GL_TEXTURE_2D, 0, GL_ALPHA, GL_FLOAT, texdata);
    glBindTexture( GL_TEXTURE_2D, 0);
    for( unsigned int z = 1; z <= this->width*this->height; ++z ) {
        std::cout << texdata[z-1] << " ";
        max = max < texdata[z-1] ? texdata[z-1] : max;
        if( z%this->width == 0 )
            std::cout << std::endl;
    }
    delete[] texdata;
    */
}

/*
 * Draw the bounding box.
 */
void protein_gl::GLSLVolumeRenderer::DrawBoundingBox(vislib::math::Cuboid<float> boundingbox) {

    // vislib::math::Vector<float, 3> position( protein->BoundingBox().GetSize().PeekDimension() );
    vislib::math::Vector<float, 3> position(boundingbox.GetSize().PeekDimension());
    position *= this->scale;

    glBegin(GL_QUADS);
    {
        // back side
        glVertex3f(0.0f, 0.0f, 0.0f);
        glVertex3f(0.0f, position.GetY(), 0.0f);
        glVertex3f(position.GetX(), position.GetY(), 0.0f);
        glVertex3f(position.GetX(), 0.0f, 0.0f);

        // front side
        glVertex3f(0.0f, 0.0f, position.GetZ());
        glVertex3f(position.GetX(), 0.0f, position.GetZ());
        glVertex3f(position.GetX(), position.GetY(), position.GetZ());
        glVertex3f(0.0f, position.GetY(), position.GetZ());

        // top side
        glVertex3f(0.0f, position.GetY(), 0.0f);
        glVertex3f(0.0f, position.GetY(), position.GetZ());
        glVertex3f(position.GetX(), position.GetY(), position.GetZ());
        glVertex3f(position.GetX(), position.GetY(), 0.0f);

        // bottom side
        glVertex3f(0.0f, 0.0f, 0.0f);
        glVertex3f(position.GetX(), 0.0f, 0.0f);
        glVertex3f(position.GetX(), 0.0f, position.GetZ());
        glVertex3f(0.0f, 0.0f, position.GetZ());

        // left side
        glVertex3f(0.0f, 0.0f, 0.0f);
        glVertex3f(0.0f, 0.0f, position.GetZ());
        glVertex3f(0.0f, position.GetY(), position.GetZ());
        glVertex3f(0.0f, position.GetY(), 0.0f);

        // right side
        glVertex3f(position.GetX(), 0.0f, 0.0f);
        glVertex3f(position.GetX(), position.GetY(), 0.0f);
        glVertex3f(position.GetX(), position.GetY(), position.GetZ());
        glVertex3f(position.GetX(), 0.0f, position.GetZ());
    }
    glEnd();


    // draw slice for volume clipping
    if (this->volClipPlaneFlag)
        this->drawClippedPolygon(boundingbox);
}

/*
 * draw the clipped polygon for correct clip plane rendering
 */
void GLSLVolumeRenderer::drawClippedPolygon(vislib::math::Cuboid<float> boundingbox) {
    if (!this->volClipPlaneFlag)
        return;

    // vislib::math::Vector<float, 3> position( protein->BoundingBox().GetSize().PeekDimension() );
    vislib::math::Vector<float, 3> position(boundingbox.GetSize().PeekDimension());
    position *= this->scale;

    // check for each clip plane
    float vcpd;
    for (int i = 0; i < (int)this->volClipPlane.Count(); ++i) {
        slices.setupSingleSlice(this->volClipPlane[i].PeekComponents(), position.PeekComponents());
        float d = 0.0f;
        vcpd = static_cast<float>(this->volClipPlane[i].PeekComponents()[3]);
        glBegin(GL_TRIANGLE_FAN);
        slices.drawSingleSlice(-(-d + vcpd - 0.0001f));
        glEnd();
    }
}

/*
 * write the current volume as a raw file
 */
void GLSLVolumeRenderer::writeVolumeRAW() {
    unsigned int numVoxel = this->volumeSize * this->volumeSize * this->volumeSize;

    float* volume = new float[numVoxel];

    glBindTexture(GL_TEXTURE_3D, this->volumeTex);
    glGetTexImage(GL_TEXTURE_3D, 0, GL_ALPHA, GL_FLOAT, volume);
    glBindTexture(GL_TEXTURE_3D, 0);

    // write array to file
    FILE* foutRaw = fopen("test.raw", "wb");
    if (!foutRaw) {
        std::cout << "could not open file for writing." << std::endl;
    } else {
        fwrite(volume, sizeof(float), numVoxel, foutRaw);
        fclose(foutRaw);
    }

    // get and write colors
    float* color = new float[4 * numVoxel];

    glBindTexture(GL_TEXTURE_3D, this->volumeTex);
    glGetTexImage(GL_TEXTURE_3D, 0, GL_RGBA, GL_FLOAT, color);
    glBindTexture(GL_TEXTURE_3D, 0);

    // write array to file
    foutRaw = fopen("color.raw", "wb");
    if (!foutRaw) {
        std::cout << "could not open file for writing." << std::endl;
    } else {
        fwrite(color, sizeof(float), 4 * numVoxel, foutRaw);
        fclose(foutRaw);
    }

    delete[] color;
}
