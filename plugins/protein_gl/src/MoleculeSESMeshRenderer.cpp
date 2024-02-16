
/*
 * MoleculeSESMeshRenderer.cpp
 *
 * Copyright (C) 2009-2021 by Universitaet Stuttgart (VIS). Alle Rechte vorbehalten.
 */


#define _USE_MATH_DEFINES 1

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884L
#endif

#include "MoleculeSESMeshRenderer.h"
#include "glm/gtx/string_cast.hpp"
#include "mmcore/param/BoolParam.h"
#include "mmcore/param/ColorParam.h"
#include "mmcore/param/EnumParam.h"
#include "mmcore/param/FilePathParam.h"
#include "mmcore/param/FloatParam.h"
#include "mmcore/param/IntParam.h"
#include "mmcore/param/StringParam.h"
#include "mmcore/utility/ColourParser.h"
#include "mmcore_gl/utility/ShaderFactory.h"
#include "mmstd/light/DistantLight.h"
#include "mmstd/light/PointLight.h"
#include "protein_calls/ProteinColor.h"
#include "vislib/OutOfRangeException.h"
#include "vislib/StringConverter.h"
#include "vislib/StringTokeniser.h"
#include "vislib/Trace.h"
#include "vislib/assert.h"
#include "vislib/sys/ASCIIFileBuffer.h"
#include "vislib/sys/File.h"
#include "vislib_gl/graphics/gl/IncludeAllGL.h"
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <cmath>

#include <fstream>

#include "protein/Icosphere.h"
#include "protein/Torus.h"

using namespace megamol;
using namespace megamol::core;
using namespace megamol::geocalls_gl;
using namespace megamol::protein;
using namespace megamol::protein_calls;
using namespace megamol::protein_gl;
using namespace megamol::core::utility::log;

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884L
#endif
float pi = M_PI;

/*
 * MoleculeSESMeshRenderer::MoleculeSESMeshRenderer
 */
MoleculeSESMeshRenderer::MoleculeSESMeshRenderer()
        : Renderer3DModuleGL()
        , molDataCallerSlot("getData", "Connects the protein SES rendering with protein data storage")
        , getLightsSlot("getLights", "Connects the protein SES rendering with light sources")
        , bsDataCallerSlot("getBindingSites", "Connects the molecule rendering with binding site data storage")
        , getTriangleDataSlot("gettriadata", "The slot publishing the generated triangle data")
        , coloringModeParam0("color::coloringMode0", "The first coloring mode.")
        , coloringModeParam1("color::coloringMode1", "The second coloring mode.")
        , cmWeightParam("color::colorWeighting", "The weighting of the two coloring modes.")
        , minGradColorParam("color::minGradColor", "The color for the minimum value for gradient coloring")
        , midGradColorParam("color::midGradColor", "The color for the middle value for gradient coloring")
        , maxGradColorParam("color::maxGradColor", "The color for the maximum value for gradient coloring")
        , colorTableFileParam("color::colorTableFilename", "The filename of the color table.")
        , probeRadiusSlot("probeRadius", "The probe radius for the surface computation")
        , datahash(0)
        , computeSesPerMolecule(false)
        , sphereColorBuffer_(nullptr)
        , sphereVertexBuffer_(nullptr)
        , pointLightBuffer_(nullptr)
        , directionalLightBuffer_(nullptr)
        , atomCount_(0)
        , curMDChash(-1) {
    this->molDataCallerSlot.SetCompatibleCall<MolecularDataCallDescription>();
    this->molDataCallerSlot.SetNecessity(core::AbstractCallSlotPresentation::Necessity::SLOT_REQUIRED);
    this->MakeSlotAvailable(&this->molDataCallerSlot);
    this->getLightsSlot.SetCompatibleCall<core::view::light::CallLightDescription>();
    this->getLightsSlot.SetNecessity(core::AbstractCallSlotPresentation::Necessity::SLOT_REQUIRED);
    this->MakeSlotAvailable(&this->getLightsSlot);
    this->bsDataCallerSlot.SetCompatibleCall<BindingSiteCallDescription>();
    this->MakeSlotAvailable(&this->bsDataCallerSlot);

    // the triangle data out slot
    this->getTriangleDataSlot.SetCallback(
        CallTriMeshDataGL::ClassName(), "GetData", &MoleculeSESMeshRenderer::getTriangleDataCallback);
    this->getTriangleDataSlot.SetCallback(
        CallTriMeshDataGL::ClassName(), "GetExtent", &MoleculeSESMeshRenderer::getExtentCallback);
    this->MakeSlotAvailable(&this->getTriangleDataSlot);

    // set epsilon value for float-comparison
    this->epsilon = vislib::math::FLOAT_EPSILON;
    // set probe radius
    this->probeRadius = 1.4f;


    this->probeRadiusSlot.SetParameter(new param::FloatParam(1.4f, 0.1f));
    this->MakeSlotAvailable(&this->probeRadiusSlot);

    // coloring modes
    this->currentColoringMode0 = ProteinColor::ColoringMode::CHAIN;
    this->currentColoringMode1 = ProteinColor::ColoringMode::ELEMENT;
    param::EnumParam* cm0 = new param::EnumParam(int(this->currentColoringMode0));
    param::EnumParam* cm1 = new param::EnumParam(int(this->currentColoringMode1));
    MolecularDataCall* mol = new MolecularDataCall();
    BindingSiteCall* bs = new BindingSiteCall();
    unsigned int cCnt;
    ProteinColor::ColoringMode cMode;
    for (cCnt = 0; cCnt < static_cast<uint32_t>(ProteinColor::ColoringMode::MODE_COUNT); ++cCnt) {
        cMode = static_cast<ProteinColor::ColoringMode>(cCnt);
        cm0->SetTypePair(static_cast<int>(cMode), ProteinColor::GetName(cMode).c_str());
        cm1->SetTypePair(static_cast<int>(cMode), ProteinColor::GetName(cMode).c_str());
    }
    delete mol;
    delete bs;
    this->coloringModeParam0 << cm0;
    this->coloringModeParam1 << cm1;
    this->MakeSlotAvailable(&this->coloringModeParam0);
    this->MakeSlotAvailable(&this->coloringModeParam1);

    // Color weighting parameter
    this->cmWeightParam.SetParameter(new param::FloatParam(0.5f, 0.0f, 1.0f));
    this->MakeSlotAvailable(&this->cmWeightParam);

    // the color for the minimum value (gradient coloring
    this->minGradColorParam.SetParameter(new param::ColorParam("#146496"));
    this->MakeSlotAvailable(&this->minGradColorParam);

    // the color for the middle value (gradient coloring
    this->midGradColorParam.SetParameter(new param::ColorParam("#f0f0f0"));
    this->MakeSlotAvailable(&this->midGradColorParam);

    // the color for the maximum value (gradient coloring
    this->maxGradColorParam.SetParameter(new param::ColorParam("#ae3b32"));
    this->MakeSlotAvailable(&this->maxGradColorParam);

    // fill color table with default values and set the filename param
    std::string filename("colors.txt");
    ProteinColor::ReadColorTableFromFile(filename, this->fileLookupTable);
    this->colorTableFileParam.SetParameter(
        new param::FilePathParam(filename, core::param::FilePathParam::FilePathFlags_::Flag_File_ToBeCreated));
    this->MakeSlotAvailable(&this->colorTableFileParam);

    // fill rainbow color table
    ProteinColor::MakeRainbowColorTable(100, this->rainbowColors);

    // width and height of the screen
    this->width = 0;
    this->height = 0;

    this->preComputationDone = false;

    auto defparams = deferredProvider_.getUsedParamSlots();
    for (const auto& param : defparams) {
        this->MakeSlotAvailable(param);
    }
}


/*
 * MoleculeSESMeshRenderer::~MoleculeSESMeshRenderer
 */
MoleculeSESMeshRenderer::~MoleculeSESMeshRenderer(void) {
    this->Release();
}


/*
 * protein::MoleculeSESMeshRenderer::release
 */
void MoleculeSESMeshRenderer::release(void) {}


/*
 * MoleculeSESMeshRenderer::create
 */
bool MoleculeSESMeshRenderer::create(void) {

    // glEnable( GL_NORMALIZE);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE_ARB);
    glEnable(GL_VERTEX_PROGRAM_TWO_SIDE);

    /* WIP
    CoreInstance* ci = this->GetCoreInstance();
    if (!ci)
        return false;
    */

    try {
        auto const shdr_options = core::utility::make_path_shader_options(
            frontend_resources.get<megamol::frontend_resources::RuntimeConfig>());

        sphereShader_ = core::utility::make_shared_glowl_shader("sphere", shdr_options,
            std::filesystem::path("protein_gl/moleculeses/mses_sphere.vert.glsl"),
            std::filesystem::path("protein_gl/moleculeses/mses_sphere.frag.glsl"));

    } catch (glowl::GLSLProgramException const& ex) {
        megamol::core::utility::log::Log::DefaultLog.WriteError("[MoleculeSESMeshRenderer] %s", ex.what());
    } catch (std::exception const& ex) {
        megamol::core::utility::log::Log::DefaultLog.WriteError(
            "[MoleculeSESMeshRenderer] Unable to compile shader: Unknown exception: %s", ex.what());
    } catch (...) {
        megamol::core::utility::log::Log::DefaultLog.WriteError(
            "[MoleculeSESMeshRenderer] Unable to compile shader: Unknown exception.");
    }

    // create the buffer objects
    sphereVertexBuffer_ = std::make_unique<glowl::BufferObject>(GL_ARRAY_BUFFER, nullptr, 0, GL_DYNAMIC_DRAW);
    sphereColorBuffer_ = std::make_unique<glowl::BufferObject>(GL_ARRAY_BUFFER, nullptr, 0, GL_DYNAMIC_DRAW);

    pointLightBuffer_ = std::make_unique<glowl::BufferObject>(GL_SHADER_STORAGE_BUFFER, nullptr, 0, GL_DYNAMIC_DRAW);
    directionalLightBuffer_ =
        std::make_unique<glowl::BufferObject>(GL_SHADER_STORAGE_BUFFER, nullptr, 0, GL_DYNAMIC_DRAW);

    glGenVertexArrays(1, &vertexArraySphere_);
    glBindVertexArray(vertexArraySphere_);

    sphereVertexBuffer_->bind();
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, nullptr);

    sphereColorBuffer_->bind();
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, nullptr);

    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);

    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    for (int i = 0; i < 6; ++i) {
        glDisableVertexAttribArray(i);
    }

       // setup all the deferred stuff
    deferredProvider_.setup(frontend_resources.get<megamol::frontend_resources::RuntimeConfig>());

    return true;
}

/*
 * MoleculeSESMeshRenderer::GetExtents
 */
bool MoleculeSESMeshRenderer::GetExtents(mmstd_gl::CallRender3DGL& call) {

    MolecularDataCall* mol = this->molDataCallerSlot.CallAs<MolecularDataCall>();
    if (mol == nullptr)
        return false;
    if (!(*mol)(1))
        return false;

    call.AccessBoundingBoxes() = mol->AccessBoundingBoxes();
    call.SetTimeFramesCount(mol->FrameCount());

    return true;
}

/*
 * MoleculeSESMeshRenderer::Render
 */
bool MoleculeSESMeshRenderer::Render(mmstd_gl::CallRender3DGL& call) {
    // temporary variables
    unsigned int cntRS = 0;

    // get camera information
    this->camera = call.GetCamera();
    view_ = this->camera.getViewMatrix();
    proj_ = this->camera.getProjectionMatrix();
    invview_ = glm::inverse(view_);
    transview_ = glm::transpose(view_);
    invproj_ = glm::inverse(proj_);
    invtransview_ = glm::transpose(invview_);
    mvp_ = proj_ * view_;
    mvpinverse_ = glm::inverse(mvp_);
    mvptranspose_ = glm::transpose(mvp_);

    fbo_ = call.GetFramebuffer();
    deferredProvider_.setFramebufferExtents(fbo_->getWidth(), fbo_->getHeight());

    std::array<int, 2> resolution = {fbo_->getWidth(), fbo_->getHeight()};

    float callTime = call.Time();

    // get pointer to CallProteinData
    MolecularDataCall* mol = this->molDataCallerSlot.CallAs<MolecularDataCall>();
    // if something went wrong --> return
    if (!mol)
        return false;

    // execute the call
    mol->SetFrameID(static_cast<int>(callTime));
    if (!(*mol)(MolecularDataCall::CallForGetData))
        return false;

    // get pointer to BindingSiteCall
    BindingSiteCall* bs = this->bsDataCallerSlot.CallAs<BindingSiteCall>();
    if (bs) {
        (*bs)(BindingSiteCall::CallForGetData);
    }

    // framebuffer object for rendering / deferred shading
    fbo_->bind();
    bool externalfbo = false;
    if (fbo_->getNumColorAttachments() == 3) {
        externalfbo = true;
    } else {
        deferredProvider_.bindDeferredFramebufferToDraw();
    }

    // ==================== check parameters ====================
    this->UpdateParameters(mol, bs);

    // ==================== Precomputations ====================
    this->probeRadius = this->probeRadiusSlot.Param<param::FloatParam>()->Value();

    // init the reduced surfaces
    if (this->reducedSurface.empty()) {
        time_t t = clock();
        // create the reduced surface
        unsigned int molIds;
        if (!this->computeSesPerMolecule) {
            this->reducedSurface.push_back(new ReducedSurface(mol, this->probeRadius));
            this->reducedSurface.back()->ComputeReducedSurface();
        } else {
            for (molIds = 0; molIds < mol->MoleculeCount(); ++molIds) {
                this->reducedSurface.push_back(new ReducedSurface(molIds, mol, this->probeRadius));
                this->reducedSurface.back()->ComputeReducedSurface();
            }
        }
        megamol::core::utility::log::Log::DefaultLog.WriteInfo(
            "%s: RS computed in: %f s\n", this->ClassName(), (double(clock() - t) / double(CLOCKS_PER_SEC)));
    }

    if (!this->preComputationDone) {
        this->colorLookupTable = {glm::make_vec3(this->minGradColorParam.Param<param::ColorParam>()->Value().data()),
            glm::make_vec3(this->midGradColorParam.Param<param::ColorParam>()->Value().data()),
            glm::make_vec3(this->maxGradColorParam.Param<param::ColorParam>()->Value().data())};

        // compute the color table
        ProteinColor::MakeWeightedColorTable(*mol, this->currentColoringMode0, this->currentColoringMode1,
            this->cmWeightParam.Param<param::FloatParam>()->Value(), // weight for the first cm
            1.0f - this->cmWeightParam.Param<param::FloatParam>()->Value(), this->atomColorTable,
            this->colorLookupTable, this->fileLookupTable, this->rainbowColors, nullptr, nullptr, true);

        // set the precomputation of the data as done
        this->preComputationDone = true;
    }

    bool virtualViewportChanged = false;
    if (static_cast<unsigned int>(std::get<0>(resolution)) != this->width ||
        static_cast<unsigned int>(std::get<1>(resolution)) != this->height) {
        this->width = static_cast<unsigned int>(std::get<0>(resolution));
        this->height = static_cast<unsigned int>(std::get<1>(resolution));
        virtualViewportChanged = true;
    }

    // ==================== Start actual rendering ====================

    this->RenderAtoms(mol);
    if (!this->reducedSurface.empty())
        this->RenderReducedSurface(this->reducedSurface[0]);

    if (externalfbo) {
        fbo_->bind();
    } else {
        deferredProvider_.resetToPreviousFramebuffer();
        deferredProvider_.draw(call, this->getLightsSlot.CallAs<core::view::light::CallLight>());
    }

    // unlock the current frame
    mol->Unlock();

    return true;
}

/*
 * update parameters
 */
void MoleculeSESMeshRenderer::UpdateParameters(const MolecularDataCall* mol, const BindingSiteCall* bs) {
    // variables
    bool recomputeColors = false;

    if (atomCount_ != mol->AtomCount()) {
        atomCount_ = mol->AtomCount();
        reducedSurface.clear();
        this->preComputationDone = false;
    }

    // ==================== check parameters ====================
    if (this->coloringModeParam0.IsDirty() || this->coloringModeParam1.IsDirty() || this->cmWeightParam.IsDirty()) {
        this->currentColoringMode0 =
            static_cast<ProteinColor::ColoringMode>(this->coloringModeParam0.Param<param::EnumParam>()->Value());
        this->currentColoringMode1 =
            static_cast<ProteinColor::ColoringMode>(this->coloringModeParam1.Param<param::EnumParam>()->Value());

        this->colorLookupTable = {glm::make_vec3(this->minGradColorParam.Param<param::ColorParam>()->Value().data()),
            glm::make_vec3(this->midGradColorParam.Param<param::ColorParam>()->Value().data()),
            glm::make_vec3(this->maxGradColorParam.Param<param::ColorParam>()->Value().data())};

        ProteinColor::MakeWeightedColorTable(*mol, this->currentColoringMode0, this->currentColoringMode1,
            this->cmWeightParam.Param<param::FloatParam>()->Value(), // weight for the first cm
            1.0f - this->cmWeightParam.Param<param::FloatParam>()->Value(), this->atomColorTable,
            this->colorLookupTable, this->fileLookupTable, this->rainbowColors, nullptr, nullptr, true);

        this->preComputationDone = false;
        this->coloringModeParam0.ResetDirty();
        this->coloringModeParam1.ResetDirty();
        this->cmWeightParam.ResetDirty();
    }
    // color table param
    if (this->colorTableFileParam.IsDirty()) {
        ProteinColor::ReadColorTableFromFile(
            this->colorTableFileParam.Param<param::FilePathParam>()->Value(), this->fileLookupTable);
        this->colorTableFileParam.ResetDirty();
        recomputeColors = true;
    }
    // check probe radius parameter slot
    if (this->probeRadiusSlot.IsDirty()) {
        this->probeRadius = this->probeRadiusSlot.Param<param::FloatParam>()->Value();
        this->reducedSurface.clear();
        this->preComputationDone = false;
        this->probeRadiusSlot.ResetDirty();
    }

    if (recomputeColors) {
        this->preComputationDone = false;
    }
}

/*
 * Render the molecular surface using GPU raycasting
 */
void MoleculeSESMeshRenderer::RenderAtoms(const MolecularDataCall* mol) {

    bool virtualViewportChanged = false;
    if (static_cast<unsigned int>(fbo_->getWidth()) != this->width ||
        static_cast<unsigned int>(fbo_->getHeight()) != this->height) {
        this->width = static_cast<unsigned int>(fbo_->getWidth());
        this->height = static_cast<unsigned int>(fbo_->getHeight());
        virtualViewportChanged = true;
    }

    // set viewport
    glm::vec4 viewportStuff;
    viewportStuff[0] = 0.0f;
    viewportStuff[1] = 0.0f;
    viewportStuff[2] = static_cast<float>(fbo_->getWidth());
    viewportStuff[3] = static_cast<float>(fbo_->getHeight());
    if (viewportStuff[2] < 1.0f)
        viewportStuff[2] = 1.0f;
    if (viewportStuff[3] < 1.0f)
        viewportStuff[3] = 1.0f;
    viewportStuff[2] = 2.0f / viewportStuff[2];
    viewportStuff[3] = 2.0f / viewportStuff[3];

    glm::vec3 camdir = camera.get<core::view::Camera::Pose>().direction;
    glm::vec3 right = camera.get<core::view::Camera::Pose>().right;
    glm::vec3 up = camera.get<core::view::Camera::Pose>().up;
    float nearplane = camera.get<core::view::Camera::NearPlane>();
    float farplane = camera.get<core::view::Camera::FarPlane>();

    unsigned int cntRS;

    /////////////////////////////////////
    // ray cast the spheres on the GPU //
    /////////////////////////////////////
    std::vector<float> tmpAtoms;
    tmpAtoms.resize(mol->AtomCount() * 4);
    for (auto i = 0; i < mol->AtomCount(); i++) {
        tmpAtoms[i * 4 + 0] = mol->AtomPositions()[i * 3 + 0];
        tmpAtoms[i * 4 + 1] = mol->AtomPositions()[i * 3 + 1];
        tmpAtoms[i * 4 + 2] = mol->AtomPositions()[i * 3 + 2];
        //tmpAtoms[i * 4 + 3] = mol->AtomTypes()[mol->AtomTypeIndices()[i]].Radius();
        tmpAtoms[i * 4 + 3] = 0.2f;
    }
    this->sphereVertexBuffer_->rebuffer(tmpAtoms.data(), tmpAtoms.size() * sizeof(float));
    this->sphereColorBuffer_->rebuffer(this->atomColorTable.data(), this->atomColorTable.size() * 3 * sizeof(float));
    glBindVertexArray(vertexArraySphere_);

    sphereShader_->use();

    // set shader variables
    sphereShader_->setUniform("viewAttr", viewportStuff);
    sphereShader_->setUniform("camIn", camdir);
    sphereShader_->setUniform("camRight", right);
    sphereShader_->setUniform("camUp", up);
    sphereShader_->setUniform("zValues", nearplane, farplane);
    sphereShader_->setUniform("view", view_);
    sphereShader_->setUniform("proj", proj_);
    sphereShader_->setUniform("viewInverse", invview_);
    sphereShader_->setUniform("mvp", mvp_);
    sphereShader_->setUniform("mvpinverse", mvpinverse_);
    sphereShader_->setUniform("mvptransposed", mvptranspose_);

    //glDrawArrays(GL_POINTS, 0, (unsigned int)mol->AtomCount());

    // disable sphere shader
    glUseProgram(0);
    glBindVertexArray(0);
}

/*
 * Render the reduced surface for debugging/testing
 */
void MoleculeSESMeshRenderer::RenderReducedSurface(protein::ReducedSurface* rs) {
    glBegin(GL_TRIANGLES);
    for (auto faceIdx = 0; faceIdx < rs->GetRSFaceCount(); faceIdx++) {
        glVertex3fv(rs->GetRSFace(faceIdx)->GetVertex1()->GetPosition().PeekComponents());
        glVertex3fv(rs->GetRSFace(faceIdx)->GetVertex2()->GetPosition().PeekComponents());
        glVertex3fv(rs->GetRSFace(faceIdx)->GetVertex3()->GetPosition().PeekComponents());
    }
    glEnd();
}


/*
 * MoleculeSESMeshRenderer::deinitialise
 */
void MoleculeSESMeshRenderer::deinitialise() {}

    /*
 * MoleculeSESMeshRenderer::getDataCallback
 */
bool MoleculeSESMeshRenderer::getTriangleDataCallback(core::Call& caller) {

    // controls
    isFlatShading = false;
    isStitching = true;
    isDebug = true;


    auto* ctmd = dynamic_cast<CallTriMeshDataGL*>(&caller);
    if (ctmd == nullptr)
        return false;

    auto* mol = this->molDataCallerSlot.CallAs<MolecularDataCall>();
    if (mol == nullptr)
        return false;
    if (!(*mol)(1))
        return false;

    int atomCnt = mol->AtomCount();
    if (isDebug) {
        atomCnt = 2;
    }


    ctmd->SetFrameCount(1);
    ctmd->AccessBoundingBoxes() = mol->AccessBoundingBoxes();

    // if no mesh is available, generate a new mesh object
    if (this->triaMesh.empty()) {
        this->triaMesh.push_back(new CallTriMeshDataGL::Mesh());
    }

    // reset vectors if data has changed (PDB data)
    if (curMDChash != static_cast<int>(mol->DataHash())) {
        curMDChash = static_cast<int>(mol->DataHash());
        vertex.resize(0);
        normal.resize(0);
        color.resize(0);
        face.resize(0);
    }

    if (this->vertex.size() == 0) {
        // fill the triangle mesh with data
        if (this->reducedSurface.empty())
            return false;
        int subdivision_level = 2;
        auto ico = new Icosphere(1, subdivision_level, !isFlatShading);
        std::vector<std::vector<unsigned int>> multipleVertexPerIco = getMultipleVertices(ico);
        int vertex_counter = (int)(ico->getVertexCount() * atomCnt);
        std::vector<bool> muss_raus;
        //std::vector<unsigned int> facesTempStore;

        //std::vector<std::vector<unsigned int>> edgeVerticesPerAtom;
        //Hier befinden sich alle Vertices die Teile einer Kante eines Atoms sind


        for (auto i = 0; i < atomCnt; i++) {
            for (int j = 0; j < ico->getVertexCount(); j++) {
                //Vertex-Koordinaten
                vertex.push_back(ico->getVertices()[j][0] * mol->AtomTypes()[mol->AtomTypeIndices()[i]].Radius() +
                                 mol->AtomPositions()[3 * i + 0]); // x
                vertex.push_back(ico->getVertices()[j][1] * mol->AtomTypes()[mol->AtomTypeIndices()[i]].Radius() +
                                 mol->AtomPositions()[3 * i + 1]); // y
                vertex.push_back(ico->getVertices()[j][2] * mol->AtomTypes()[mol->AtomTypeIndices()[i]].Radius() +
                                 mol->AtomPositions()[3 * i + 2]); // z

                //normal[ico->getVertexCount() * 3 * i + 3 * j + 0] = ico->getNormals()[3 * j + 0];
                normal.push_back(ico->getNormals()[j][0]);
                normal.push_back(ico->getNormals()[j][1]);
                normal.push_back(ico->getNormals()[j][2]);


                /*color.push_back(0.6f);
                color.push_back(0.6f);
                color.push_back(0.6f);*/
                if (i == 0) {
                    color.push_back(1.0f);
                    color.push_back(0.0f);
                    color.push_back(0.0f);
                } else if (i == 1) {
                    color.push_back(0.0f);
                    color.push_back(1.0f);
                    color.push_back(0.0f);
                } else {
                    color.push_back(0.0f);
                    color.push_back(0.0f);
                    color.push_back(1.0f);
                }

                /*color.push_back(0.0f);
                color.push_back(0.0f);
                color.push_back(0.0f);*/

                muss_raus.push_back(false);
            }
        }

        for (auto i = 0; i < atomCnt; i++) {
            //Grenzen für Interatkion
            unsigned int lower_bound = vertex_counter / atomCnt * i;
            unsigned int upper_bound = vertex_counter / atomCnt * (i + 1);
            //Koord. aktuelles Atom
            float atom_x = mol->AtomPositions()[(3 * i) + 0];
            float atom_y = mol->AtomPositions()[(3 * i) + 1];
            float atom_z = mol->AtomPositions()[(3 * i) + 2];
            for (int j = 0; j < vertex_counter; ++j) {
                if (j >= lower_bound && j < upper_bound) {
                    continue;
                }
                //prüfe Kollision
                if (std::sqrt((atom_x - vertex[(j * 3) + 0]) * (atom_x - vertex[(j * 3) + 0]) +
                              (atom_y - vertex[(j * 3) + 1]) * (atom_y - vertex[(j * 3) + 1]) +
                              (atom_z - vertex[(j * 3) + 2]) * (atom_z - vertex[(j * 3) + 2])) <
                    mol->AtomTypes()[mol->AtomTypeIndices()[i]].Radius()) {
                    //color.at(j*3) = 1.0f;
                    muss_raus.at(j) = true;

                    CutTri cutTriangle = {//wandle lokale Indizes zu globalen Indizes um
                        ico->getVertexCount() * i + ico->getIndices()[j + 0],
                        ico->getVertexCount() * i + ico->getIndices()[j + 1],
                        ico->getVertexCount() * i + ico->getIndices()[j + 2], i};

                    int atomIndex1 = i;
                    int atomIndex2 = j / ico->getVertexCount();

                    atomCollisions.insert({atomIndex1, atomIndex2});
                }
            }
        }

        std::tuple<unsigned int, int> bla;
        bla = {1, 2};
        int xx = std::get<1>(bla);

        //TEST KUGEL ZUM DEBUGGEN in grün
        /*
    std::vector<int> zusatz = {79, 80, 337, 248, 86, 87, 89, 93, 95, 32, 97, 17, 110, 111, 116, 120}; //84 == 80 247 oder 333

        for (int i : zusatz) {


            for (int j = 0; j < ico->getVertexCount(); j++) {
                vertex.push_back(ico->getVertices()[j][0] * 0.01 + vertex.at(i * 3 + 0)); // x
                vertex.push_back(ico->getVertices()[j][1] * 0.01 + vertex.at(i * 3 + 1));
                vertex.push_back(ico->getVertices()[j][2] * 0.01 + vertex.at(i * 3 + 2));

                //normal[ico->getVertexCount() * 3 * i + 3 * j + 0] = ico->getNormals()[3 * j + 0];
                normal.push_back(ico->getNormals()[j][0]);
                normal.push_back(ico->getNormals()[j][1]);
                normal.push_back(ico->getNormals()[j][2]);

                if (i == 80) {
                    color.push_back(1);
                    color.push_back(0);
                    color.push_back(0);
                } else if (i == 337) {
                    color.push_back(0);
                    color.push_back(0);
                    color.push_back(1);
                } else if (i == 17) {
                    color.push_back(0);
                    color.push_back(1);
                    color.push_back(1);
                } else {
                    color.push_back(0);
                    color.push_back(1);
                    color.push_back(0);
                }
            }
        }
        */
        int face_counter = (int)(ico->getIndexCount() * atomCnt);
        std::vector<std::vector<std::tuple<unsigned int, int>>> referenceToOtherVertice(
            ico->getVertexCount() * atomCnt);
        /* referenceToOtherVertice: vector of all vertices with a "tuple" in it.
         * 1. field is the index of the connected vertex
         * 2. field is the number of faces of the edge
         */
        std::vector<bool> besucht;
        std::vector<unsigned int> zuAtom;
        std::vector<std::vector<unsigned int>> edgesInfo;
        /* edgesInfo:
         * erstes element, referenz zu erstem vertice
         * zwrites element, referenz zu zweitem vertice
         * drittes element, Anzahl faces der Kante
         */

        // Gehe jedes Atom durch           | zusatz ist debug
        // for (auto i = 0; i < atomCnt + zusatz.size(); i++) {
        for (auto i = 0; i < atomCnt; i++) {

            std::vector<unsigned int> atomEdgeVerticeVector;

            int redCounter = 0;
            // Gehe in jedem atom alle Ecken durch
            for (int j = 0; j < ico->getIndexCount(); j += 3) {

                // for DEBUG: set atomCnt to number of atoms that are going to be cut (e.g., reduce to 1 or 2)
                if (i >= atomCnt) {                                                       //DEBUG CODE FÜR ZUSATZ KUGEL
                    face.push_back(ico->getVertexCount() * i + ico->getIndices()[j + 0]); //Vertice 1 von face
                    face.push_back(ico->getVertexCount() * i + ico->getIndices()[j + 1]); //Vertice 2 von face
                    face.push_back(ico->getVertexCount() * i + ico->getIndices()[j + 2]); //Vertice 3 von face
                }

                // checks for each Vertex of a face, weather it is in an atom or not
                //Alle Vertices sind außerhalb anderen Kugeln, also werden sie gezeichnet!
                else if (!(muss_raus.at(ico->getVertexCount() * i + ico->getIndices()[j + 0]) ||
                             muss_raus.at(ico->getVertexCount() * i + ico->getIndices()[j + 1]) ||
                             muss_raus.at(ico->getVertexCount() * i + ico->getIndices()[j + 2]))) {
                    face.push_back(ico->getVertexCount() * i + ico->getIndices()[j + 0]); //Vertice 1 von face
                    face.push_back(ico->getVertexCount() * i + ico->getIndices()[j + 1]); //Vertice 2 von face
                    face.push_back(ico->getVertexCount() * i + ico->getIndices()[j + 2]); //Vertice 3 von face
                } // Nicht alle Vertices sind außerhalb. Schaue, dass nur Vertices die 2/3 drin sind ausgewählt werden.
                // Wenn alle vertices drin liegen: ignore
                // Wenn 2 Vertices drin liegen:Ignore, da
                else {
                    //                              |Richtige Vertices pro Atom| Finde Richtige Vertices aus liste
                    bool staysIn0 = !muss_raus.at(ico->getVertexCount() * i + ico->getIndices()[j + 0]);
                    bool staysIn1 = !muss_raus.at(ico->getVertexCount() * i + ico->getIndices()[j + 1]);
                    bool staysIn2 = !muss_raus.at(ico->getVertexCount() * i + ico->getIndices()[j + 2]);

                    //2 out of 3
                    // Sind IMMER miteinander verbunden UND folgen aufeinander, d.h. ich kann nix weglassen
                    if (staysIn0 ? (staysIn1 || staysIn2) : (staysIn1 && staysIn2)) {
                        if (staysIn0) {
                            zuAtom.push_back(ico->getVertexCount() * i + ico->getIndices()[j + 0]);
                            //color.at((ico->getVertexCount() * i + ico->getIndices()[j + 0])*3) = 1;

                            atomEdgeVerticeVector.push_back(ico->getVertexCount() * i + ico->getIndices()[j + 0]);
                            if (staysIn1) {
                                referenceToOtherVertice.at(ico->getVertexCount() * i + ico->getIndices()[j + 0]);

                                referenceToOtherVertice.at(ico->getVertexCount() * i + ico->getIndices()[j + 0])
                                    .emplace_back(ico->getVertexCount() * i + ico->getIndices()[j + 1], 1);
                            }
                            if (staysIn2) {
                                referenceToOtherVertice.at(ico->getVertexCount() * i + ico->getIndices()[j + 0])
                                    .emplace_back(ico->getVertexCount() * i + ico->getIndices()[j + 2], 1);
                            }
                        }
                        if (staysIn1) {
                            zuAtom.push_back(ico->getVertexCount() * i + ico->getIndices()[j + 1]);
                            //color.at((ico->getVertexCount() * i + ico->getIndices()[j + 1])*3) = 1;

                            atomEdgeVerticeVector.push_back(ico->getVertexCount() * i + ico->getIndices()[j + 1]);
                            if (staysIn0) {
                                referenceToOtherVertice.at(ico->getVertexCount() * i + ico->getIndices()[j + 1])
                                    .emplace_back(ico->getVertexCount() * i + ico->getIndices()[j + 0], 1);
                            }
                            if (staysIn2) {
                                referenceToOtherVertice.at(ico->getVertexCount() * i + ico->getIndices()[j + 1])
                                    .emplace_back(ico->getVertexCount() * i + ico->getIndices()[j + 2], 1);
                            }
                        }
                        if (staysIn2) {
                            zuAtom.push_back(ico->getVertexCount() * i + ico->getIndices()[j + 2]);
                            //color.at((ico->getVertexCount() * i + ico->getIndices()[j + 2])*3) = 1;

                            atomEdgeVerticeVector.push_back(ico->getVertexCount() * i + ico->getIndices()[j + 2]);
                            if (staysIn0) {
                                referenceToOtherVertice.at(ico->getVertexCount() * i + ico->getIndices()[j + 2])
                                    .emplace_back(ico->getVertexCount() * i + ico->getIndices()[j + 0], 1);
                            }
                            if (staysIn1) {
                                referenceToOtherVertice.at(ico->getVertexCount() * i + ico->getIndices()[j + 2])
                                    .emplace_back(ico->getVertexCount() * i + ico->getIndices()[j + 1], 1);
                            }
                        }
                    }

                    //TODO: Hier wohl die erweiterte Datenstruktur rein bringen.
                    // Was kann ich hier rein bringen
                    // Nehme 3 Listen : Done
                    //  Fülle die Listen mit Daten
                    /*
                     *  +------------------+
                        |     Vertice      |
                        +------------------+
                        | - besucht: bool  |
                        | - index: int[4]  |
                        | - index_atom: int|
                        +------------------+

                     */
                }
            }
            edgeVerticesPerAtom.push_back(atomEdgeVerticeVector);
        }

        // Hier habe ich jetzt alle Vertices der Atome.
        // Gehe diese durch und finde nächsten Vertice auf anderen Atomen, nach Verticeliste

        /* TODO:
         *  Hier muss der neue Algorithmus rein.
         *  Komme an, schaue nächste Vertices an.
         */
        std::cout << "vertex.size() / 3 = " << vertex.size() / 3 << std::endl;
        
        std::vector<float> torusVertices;
        std::vector<float> torusNormals;
        std::vector<unsigned int> torusIndices;
        std::vector<unsigned int> torusLineIndices;
        std::vector<unsigned int> torusFaceIndices;

        int highestIcosphereIndex = -1;

        //get highest ico index for torus offset
        for (int i = 0; i < ico->getIndexCount(); ++i) {
            highestIcosphereIndex = std::max(highestIcosphereIndex, static_cast<int>(ico->getIndices()[i]));
        }

        //offset for torus vertices
        int offset = vertex.size() / 3;

        if (offset <= highestIcosphereIndex) {
            offset = highestIcosphereIndex + 1;
        }

        //set subdivision for torus
        int numSegments = 60;
        int numRings = 40;

        bool faceRemoval = false;
        bool withStitching = false;

        //generate torus for every colliding atom pair
        for (auto elem : atomCollisions) {

            //values of colliding atoms
            unsigned int atomIndex1 = elem.first;
            unsigned int atomIndex2 = elem.second;

            std::cout << "Atompaar: " << atomIndex1 << "/" << atomIndex2 << std::endl;

            float atomRadius1 = mol->AtomTypes()[mol->AtomTypeIndices()[atomIndex1]].Radius();
            float atomRadius2 = mol->AtomTypes()[mol->AtomTypeIndices()[atomIndex2]].Radius();

            std::cout << "atomRadius1: " << atomRadius1 << std::endl;
            std::cout << "atomRadius2: " << atomRadius2 << std::endl;

            glm::vec3 atomPos1(mol->AtomPositions()[3 * atomIndex1 + 0], mol->AtomPositions()[3 * atomIndex1 + 1],
                mol->AtomPositions()[3 * atomIndex1 + 2]);

            glm::vec3 atomPos2(mol->AtomPositions()[3 * atomIndex2 + 0], mol->AtomPositions()[3 * atomIndex2 + 1],
                mol->AtomPositions()[3 * atomIndex2 + 2]);


            //values for torus generation
            float distance = Torus::getDistance(atomPos1, atomPos2);
            glm::vec3 center = Torus::getTorusCenter(atomPos1, atomPos2, atomRadius1, atomRadius2, probeRadius);
            float radius = Torus::getTorusRadius(atomPos1, atomPos2, atomRadius1, atomRadius2, probeRadius);
            float rotationAngle = Torus::getRotationAngle(atomPos1, atomPos2);
            glm::vec3 axisUnitVec = Torus::getTorusAxisUnitVec(atomPos1, atomPos2);

            //torus generation
            Torus torus(center, radius, numSegments, numRings);
            torus.generateTorus(probeRadius, axisUnitVec, rotationAngle, offset, atomPos1, atomPos2,
                atomRadius1, atomRadius2);

            //Visibility Sphere: used to remove vertices of atoms inside of visibility sphere
            //dummy point for visibility sphere
            float theta = 2 * pi * static_cast<float>(0) / static_cast<float>(numSegments);
            float phi = 2 * pi * static_cast<float>(0) / static_cast<float>(numRings);

            float x = (radius + probeRadius * cos(phi)) * cos(theta);
            float y = (radius + probeRadius * cos(phi)) * sin(theta);
            float z = probeRadius * sin(phi);

            //calculate visibility sphere
            glm::vec3 probePos =
                torus.getProbePos(center, radius, probeRadius, glm::vec3(x, y, z), rotationAngle, axisUnitVec);
            glm::vec3 contactPoint(0, 0, 0);
            std::pair<glm::vec3, float> visibilitySphere;

            //use smaller atom to calculate contact point, maybe unnecessary(?)
            if (atomRadius1 <= atomRadius2) {
                contactPoint = torus.getContactPoint(probePos, atomPos1, atomRadius1);
                visibilitySphere = torus.getVisibilitySphere(contactPoint, probePos, atomPos1, atomPos2, center);
            } else {
                contactPoint = torus.getContactPoint(probePos, atomPos2, atomRadius2);
                visibilitySphere = torus.getVisibilitySphere(contactPoint, probePos, atomPos2, atomPos1, center);
            }
            glm::vec3 vsCenter = visibilitySphere.first;
            float vsRadius = visibilitySphere.second;
            std::cout << "modified vsRadius: " << vsRadius << std::endl;
            std::cout << "vsCenter: " << vsCenter.x << "/" << vsCenter.y << "/" << vsCenter.z << std::endl;
            std::vector<float> stitchingVertex;

            //remove all vertices of atoms that are inside of the visibility sphere
            std::vector<float> edgeVerticesAtom;
            if (faceRemoval) {
                //std::cout << "faceSize before removal: " << face.size() << std::endl;
                //std::cout << "vertex.size() = " << vertex.size() << std::endl;
                std::vector<unsigned int> indicesToRemove;
                for (int i = 0; i < vertex.size(); i += 3) {
                    float x = vertex[i];
                    float y = vertex[i + 1];
                    float z = vertex[i + 2];

                    glm::vec3 currentVertex(x, y, z);
                    float dist = glm::distance(currentVertex, vsCenter);

                    //if vertex is inside of vs, add its index in face to indicesToRemove
                    if (dist <= vsRadius) {
                        indicesToRemove.push_back(i / 3);
                    }
                }
                std::sort(indicesToRemove.rbegin(), indicesToRemove.rend());


                std::vector<bool> removeMarker(face.size(), false);
                for (auto index : indicesToRemove) {
                    for (size_t i = 0; i < face.size(); ++i) {
                        if (face[i] == index) {
                            //if at least one vertex of a triangle is inside of vs => remove that triangle
                            //three values in face are one triangle
                            //need to make sure to remove the right values in face
                            //possibilities:
                            //indexToRemove, +1          , +2
                            //-1           ,indexToRemove, +1
                            //-2           ,-1           ,indexToRemove

                            int modResult = i % 3;
                            removeMarker[i - modResult] = true;
                            removeMarker[i - modResult + 1] = true;
                            removeMarker[i - modResult + 2] = true;
                        }
                    }
                }

                std::vector<unsigned int> newFace;
                for (size_t i = 0; i < face.size(); ++i) {
                    if (!removeMarker[i]) {
                        //only keep all completely visible triangles
                        newFace.push_back(face[i]);
                    }
                }
                face = std::move(newFace);

                //look for edge vertices for each atom, needed for stitching
                std::vector<EdgeVerticesOfAtoms> atomEdges;
                for (int i = 0; i < face.size(); i += 3) {
                    //add all edges of all visible triangles

                    //std::cout << "triangle " << i / 3 << ": " << face[i] << "/" << face[i + 1] << "/" << face[i + 2] << std::endl;
                    EdgeVerticesOfAtoms edge1 = EdgeVerticesOfAtoms(face[i], face[i + 1]);
                    EdgeVerticesOfAtoms edge2 = EdgeVerticesOfAtoms(face[i + 1], face[i + 2]);
                    EdgeVerticesOfAtoms edge3 = EdgeVerticesOfAtoms(face[i + 2], face[i]);

                    atomEdges.push_back(edge1);
                    atomEdges.push_back(edge2);
                    atomEdges.push_back(edge3);
                }

                //count occurences of each vertex index in each edge
                std::unordered_map<unsigned int, int> vertexCount;
                for (const auto& edge : atomEdges) {
                    vertexCount[edge.index1]++;
                    vertexCount[edge.index2]++;
                }

                //if a vertex has less than six occurences it's an edge vertex
                std::vector<unsigned int> edgeIndices;
                for (auto elem : vertexCount) {
                    //std::cout << elem.first << " existiert " << elem.second << " mal" << std::endl;
                    if (elem.second < 6) {
                        edgeIndices.push_back(elem.first);
                    }
                }

                //convert indices to vertex values
                for (auto elem : edgeIndices) {
                    //std::cout << elem << std::endl;
                    edgeVerticesAtom.push_back(vertex[elem * 3]);
                    edgeVerticesAtom.push_back(vertex[(elem * 3) + 1]);
                    edgeVerticesAtom.push_back(vertex[(elem * 3) + 2]);
                }
            }
            
            //add generated vertices to global vertex structure
            int acceptedVertices = 0;
            for (int i = 0; i < torus.getVertexCount(); i++) {
                //need to add all vertices, remove torus triangles outside of vs by only adding those to face, ignore next comment
                //only add vertex, color and normal if vertex is inside visibility sphere
                if (torus.getVertices_bool(i) == true) {
                    vertex.push_back(torus.getVertices()[i][0]);
                    vertex.push_back(torus.getVertices()[i][1]);
                    vertex.push_back(torus.getVertices()[i][2]);
                    torusVertices.push_back(torus.getVertices()[i][0]);
                    torusVertices.push_back(torus.getVertices()[i][1]);
                    torusVertices.push_back(torus.getVertices()[i][2]);
                    /*color.push_back(1.0f);
                    color.push_back(0.0f);
                    color.push_back(0.0f);*/
                    color.push_back(0.6f);
                    color.push_back(0.6f);
                    color.push_back(0.6f);
                    normal.push_back(torus.getNormals()[i][0]);
                    normal.push_back(torus.getNormals()[i][1]);
                    normal.push_back(torus.getNormals()[i][2]);
                    acceptedVertices++;
                } else {
                    vertex.push_back(torus.getVertices()[i][0]);
                    vertex.push_back(torus.getVertices()[i][1]);
                    vertex.push_back(torus.getVertices()[i][2]);
                    color.push_back(0.0f);
                    color.push_back(0.0f);
                    color.push_back(0.0f);
                    normal.push_back(torus.getNormals()[i][0]);
                    normal.push_back(torus.getNormals()[i][1]);
                    normal.push_back(torus.getNormals()[i][2]);
                }
            }
            std::cout << "Accepted vertices: " << acceptedVertices << std::endl;

            //not necessary
            for (int i = 0; i < torus.getIndexCount(); i++) {
                torusIndices.push_back(torus.getIndices()[i]);
            }

            //not necessary
            for (int i = 0; i < torus.getLineIndicesCount(); i++) {
                torusLineIndices.push_back(torus.getLineIndices()[i]);
            }

            //add generated face indices to global face index structure
            for (int i = 0; i < torus.getFaceIndicesCount(); i++) {
                face.push_back(torus.getFaceIndices()[i]);
            }

            //get all edge vertices of a torus => less possible targets for stitching
            for (int i = 0; i < torus.edgeVertices.size(); i++) {
                torusEdgeVertices.push_back(torus.edgeVertices[i]);
                //std::cout << "torusEdgeVertices: " << torus.edgeVertices[i].x << "/" << torus.edgeVertices[i].y << "/" << torus.edgeVertices[i].z << std::endl;
            }


            //stitching
            if (withStitching) {
                for (int i = 0; i < edgeVerticesAtom.size(); i += 3) {
                    float x = edgeVerticesAtom[i];
                    float y = edgeVerticesAtom[i + 1];
                    float z = edgeVerticesAtom[i + 2];

                    //take current and next vertex of edge vertices of atoms, calculate the closest vertex of torus edge vertices
                    //form two triangles out of those four vertices

                    //problem when atom changes

                    glm::vec3 currentVertex(x, y, z);
                    //std::cout << "currentVertex: " << currentVertex.x << "/" << currentVertex.y << "/" << currentVertex.z << std::endl;

                    glm::vec3 closestVertexToCurrent = findClosestVertex(currentVertex, torusEdgeVertices);
                    //std::cout << "closestVertexToCurrent: " << closestVertexToCurrent.x << "/" << closestVertexToCurrent.y << "/" << closestVertexToCurrent.z << std::endl;

                    glm::vec3 nextVertex;
                    if (i + 3 >= edgeVerticesAtom.size()) {
                        nextVertex = glm::vec3(edgeVerticesAtom[0], edgeVerticesAtom[1], edgeVerticesAtom[2]);
                    } else {
                        nextVertex =
                            glm::vec3(edgeVerticesAtom[i + 3], edgeVerticesAtom[i + 4], edgeVerticesAtom[i + 5]);
                    }
                    //std::cout << "nextVertex: " << nextVertex.x << "/" << nextVertex.y << "/" << nextVertex.z << std::endl;

                    glm::vec3 closestVertexToNext = findClosestVertex(nextVertex, torusEdgeVertices);
                    //std::cout << "closestVertexToNext: " << closestVertexToNext.x << "/" << closestVertexToNext.y << "/" << closestVertexToNext.z << std::endl;

                    
                    std::vector<unsigned int> tempIndices;
                    for (const auto& point : {currentVertex, closestVertexToCurrent, closestVertexToNext, nextVertex}) {
                        vertex.push_back(point.x);
                        vertex.push_back(point.y);
                        vertex.push_back(point.z);

                        color.push_back(0.6f);
                        color.push_back(0.6f);
                        color.push_back(0.6f);

                        glm::vec3 vertexNormal = glm::normalize(point);
                        normal.push_back(vertexNormal.x);
                        normal.push_back(vertexNormal.y);
                        normal.push_back(vertexNormal.z);

                        unsigned int tempIndex = vertex.size() / 3 - 1;
                        tempIndices.push_back(tempIndex);
                        //std::cout << "tempIndex: " << tempIndex << std::endl;
                    }

                    face.push_back(tempIndices[1]);
                    face.push_back(tempIndices[2]);
                    face.push_back(tempIndices[0]);

                    face.push_back(tempIndices[3]);
                    face.push_back(tempIndices[0]);
                    face.push_back(tempIndices[2]);
                    
                }
                
            }
            //update offset
            offset += torus.getVertexCount();
        }

        
        delete ico;
    }


    //TODO: Kann rausoptimiert werden.
    auto* vertex2 = new float[vertex.size()]{};
    for (int i = 0; i < vertex.size(); ++i) {
        vertex2[i] = vertex.at(i);
    }
    auto* normal2 = new float[normal.size()]{};
    for (int i = 0; i < normal.size(); ++i) {
        normal2[i] = normal.at(i);
    }
    auto* color2 = new float[color.size()]{};
    for (int i = 0; i < color.size(); ++i) {
        color2[i] = color.at(i);
    }
    auto* face2 = new unsigned int[face.size()];
    for (int i = 0; i < face.size(); ++i) {
        face2[i] = face.at(i);
    }


    this->triaMesh[0]->SetVertexData(vertex.size() / 3, vertex2, normal2, color2, nullptr, true);
    this->triaMesh[0]->SetTriangleData(face.size() / 3, face2, true);
    this->triaMesh[0]->SetMaterial(nullptr);



    // set triangle mesh to caller
    if (this->triaMesh[0]->GetVertexCount() > 0) {
        ctmd->SetDataHash(this->datahash);
        ctmd->SetObjects(1, this->triaMesh[0]);
        ctmd->SetUnlocker(nullptr);
        return true;
    } else {
        return false;
    }
}




/*
 * MoleculeSESMeshRenderer::getExtentCallback
 */
bool MoleculeSESMeshRenderer::getExtentCallback(core::Call& caller) {
    auto* ctmd = dynamic_cast<CallTriMeshDataGL*>(&caller);
    if (ctmd == nullptr)
        return false;

    auto* mol = this->molDataCallerSlot.CallAs<MolecularDataCall>();
    if (mol == nullptr)
        return false;
    if (!(*mol)(1))
        return false;

    ctmd->SetFrameCount(1);
    ctmd->AccessBoundingBoxes() = mol->AccessBoundingBoxes();

    return true;
}
/*
 * MoleculeSESMeshRenderer::findNearestVertice
 */
std::vector<unsigned int> MoleculeSESMeshRenderer::findNearestVertice(const std::vector<std::vector<unsigned int>>& edgelord,
    unsigned int& referenceIndex0, const std::vector<float>& vertex, int index) {
    /* Int referenceIndex ist vetice index, edgelord ist
     * edgelord ist kanten vector
     *
     * gl point
     */
    //TODO: Baue Winkel ein
    float winkelFaktor = 1;
    float referenceIndexX = (vertex.at(referenceIndex0 * 3 + 0) + vertex.at(referenceIndex0 * 3 + 0)) / 2;
    float referenceIndexY = (vertex.at(referenceIndex0 * 3 + 1) + vertex.at(referenceIndex0 * 3 + 1)) / 2;
    float referenceIndexZ = (vertex.at(referenceIndex0 * 3 + 2) + vertex.at(referenceIndex0 * 3 + 2)) / 2;
    unsigned int indexOfNearestVertex = 0;
    unsigned int indexOfSecondNearestVertex = 0;
    unsigned int indexOfThirdNearestVertex = 0;
    float nearestDistance = FLT_MAX;
    float secondNearestDistance = FLT_MAX;
    float thirdNearestDistance = FLT_MAX;
    for (int j = 0; j < edgelord.size(); j++) {
        //Gehe alle vectoren drch
        if (index != j) {
            for (unsigned int i : edgelord[j]) {
                float dist =
                    std::sqrt((referenceIndexX - vertex.at(i * 3 + 0)) * (referenceIndexX - vertex.at(i * 3 + 0)) +
                              (referenceIndexY - vertex.at(i * 3 + 1)) * (referenceIndexY - vertex.at(i * 3 + 1)) +
                              (referenceIndexZ - vertex.at(i * 3 + 2)) * (referenceIndexZ - vertex.at(i * 3 + 2)));

                if (dist > 0 && dist * winkelFaktor < nearestDistance) {
                    if (i != indexOfNearestVertex){
                        thirdNearestDistance = secondNearestDistance;
                        secondNearestDistance = nearestDistance;
                        nearestDistance = dist;

                        indexOfThirdNearestVertex = indexOfSecondNearestVertex;
                        indexOfSecondNearestVertex = indexOfNearestVertex;
                        indexOfNearestVertex = i;
                    }}
                //TODO: hier auch nochmal nach winkel schauen
                else if (dist < secondNearestDistance){
                    if (i != indexOfNearestVertex) {
                        thirdNearestDistance = secondNearestDistance;
                        secondNearestDistance = dist;

                        indexOfThirdNearestVertex = indexOfSecondNearestVertex;
                        indexOfSecondNearestVertex = i;
                    }
                }
                else if (dist < thirdNearestDistance) {
                    if (i != indexOfSecondNearestVertex) {

                        thirdNearestDistance = dist;

                        indexOfThirdNearestVertex = i;
                    }
                }

            }
        }
    }
    return {indexOfNearestVertex, indexOfSecondNearestVertex, indexOfThirdNearestVertex};
}

glm::vec3 MoleculeSESMeshRenderer::findClosestVertex(const glm::vec3& inputVertex, const std::vector<glm::vec3>& vertexList) {
    float minimalDist = std::numeric_limits<float>::max();
    glm::vec3 closestVertex(0, 0, 0);

    //with input std::vector<float>&
    /*for (int i = 0; i < vertexList.size(); i+=3) {
        glm::vec3 vertex(vertexList[i], vertexList[i + 1], vertexList[i + 2]);

        float dist = glm::distance(inputVertex, vertex);
        if (dist < minimalDist) {
            minimalDist = dist;
            closestVertex = vertex;
        }
    }*/

    //with input std::vector<glm::vec3>&
    for (const auto& elem : vertexList) {
        float dist = glm::distance(inputVertex, elem);

        if (dist < minimalDist) {
            minimalDist = dist;
            closestVertex = elem;
        }
    }

    return closestVertex;
}

void MoleculeSESMeshRenderer::cleanupFace(std::vector<unsigned int>& face, int indexToRemove) {
    for (size_t i = 0; i < face.size(); ++i) {
        //std::cout << face[i] << std::endl;
        if (face[i] == indexToRemove) {
            //std::cout << "found i: " << i << std::endl;
            //determine the range of indices to remove based on the current index and modulo 3
            size_t startIndex = (i % 3 == 0) ? i : ((i % 3 == 1) ? i - 1 : i - 2);
            size_t endIndex = (i % 3 == 2) ? i : ((i % 3 == 1) ? i + 1 : i + 2);

            //check bounds and adjust if necessary
            startIndex = std::max(startIndex, static_cast<size_t>(0));
            endIndex = std::min(endIndex, face.size() - 1);

            // Erase the elements from startIndex to endIndex inclusive
            //std::cout << "start removing..." << std::endl;
            face.erase(face.begin() + startIndex, face.begin() + endIndex + 1);
            //std::cout << "removed" << std::endl;
            // Adjust the index after erase
            if (startIndex == 0) {
                i = -1; 
            } else {
                i = startIndex - 1;
            }
        } else {
            continue;
        }
    }
}

void MoleculeSESMeshRenderer::cleanupFace2(std::vector<unsigned int>& face, std::vector<unsigned int>& indicesToRemove) {
    // Erstellt eine Kopie von `face`, die aktualisiert wird
    std::vector<unsigned int> updatedFace;

    // Erstellt ein Set für schnellen Zugriff, um zu überprüfen, ob ein Index entfernt wurde
    std::set<int> indicesToRemoveSet(indicesToRemove.begin(), indicesToRemove.end());

    // Geht durch jedes Dreieck (angenommen, dass `face` in Dreiergruppen von Indizes organisiert ist)
    for (size_t i = 0; i < face.size(); i += 3) {
        // Überprüft, ob eines der Dreiecksgesichter einen zu entfernenden Vertex verwendet
        if (indicesToRemoveSet.find(face[i] / 3) == indicesToRemoveSet.end() &&
            indicesToRemoveSet.find(face[i + 1] / 3) == indicesToRemoveSet.end() &&
            indicesToRemoveSet.find(face[i + 2] / 3) == indicesToRemoveSet.end()) {
            // Wenn keiner der Vertizes des Gesichts entfernt wurde, fügt es dem aktualisierten `face` hinzu
            updatedFace.push_back(face[i]);
            updatedFace.push_back(face[i + 1]);
            updatedFace.push_back(face[i + 2]);
        }
    }

    // Aktualisiert die Indizes in `updatedFace` basierend auf den gelöschten Vertizes
    for (auto& index : updatedFace) {
        int decrement = 0;
        for (int removedIndex : indicesToRemove) {
            // Verringert den Index, wenn er größer als der entfernte Index ist
            if (index > removedIndex) {
                decrement++;
            }
        }
        index -= decrement;
    }

    // Ersetzt `face` durch `updatedFace`
    face = std::move(updatedFace);
}

std::vector<std::vector<unsigned int>> MoleculeSESMeshRenderer::getMultipleVertices(Icosphere* pIcosphere) {
    std::vector<std::vector<unsigned int>> multipleVertexVector(pIcosphere->getVertexCount());
    //1953 torusVertexCount 
    for (int i = 0; i < pIcosphere->getVertexCount(); ++i) {
        // doppelte for schleife, gehe jedes vertice durch und schaue die vertices an
        for (int j = i+1; j < pIcosphere->getVertexCount(); ++j) {
            if (pIcosphere->getVertices()[i*3] == pIcosphere->getVertices()[j*3] &&
                pIcosphere->getVertices()[i*3+1] == pIcosphere->getVertices()[j*3+1] &&
                pIcosphere->getVertices()[i*3+2] == pIcosphere->getVertices()[j*3+2]){
                multipleVertexVector.at(i).push_back(j);
                multipleVertexVector.at(j).push_back(i);
            }
        }
    }
    return multipleVertexVector;
}
std::vector<unsigned int> MoleculeSESMeshRenderer::findVector(const std::vector<unsigned int>& edgelord) {
    // Finde Winkel zu anderem Punkt
    /*  Gegeben: Punkt A, Punkt B; A auf Atom A, B auf Atom B
     * gegeben Schnittkreis
     * Berechne nächsten Punkt auf Kreis zu Punkt A,
     * Strecke Kreisebene Punkt A
     * Strecke Punkt A und Punkt B
     * Berechne Winkel zwischen den zwei Strecken
     * gebe winkel zurück
     */
    return {};
}
// stichworte rein schreiben, was muss rein, aus den stichworten sätze machen und dann sätze verbessern
/* anfangshürde klein machen und stichpunkte, 2,3 worte für jeden Absatz etc
 * eltech language tool
 *
 */
