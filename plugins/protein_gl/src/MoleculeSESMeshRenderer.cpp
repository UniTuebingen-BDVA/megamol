
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


#include "protein/Icosphere.h"
#include "protein/Torus.h"

using namespace megamol;
using namespace megamol::core;
using namespace megamol::geocalls_gl;
using namespace megamol::protein;
using namespace megamol::protein_calls;
using namespace megamol::protein_gl;
using namespace megamol::core::utility::log;

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
    if(isDebug) {
        atomCnt = 3;
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

        std::vector<std::vector<unsigned int>>
            edgeVerticesPerAtom; //Hier befinden sich alle Vertices die Teile einer Kante eines Atoms sind


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


                color.push_back(0.6f);
                color.push_back(0.6f);
                color.push_back(0.6f);

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
                    
                    CutTri cutTriangle = {
                        //wandle lokale Indizes zu globalen Indizes um
                        ico->getVertexCount() * i + ico->getIndices()[j + 0],
                        ico->getVertexCount() * i + ico->getIndices()[j + 1],
                        ico->getVertexCount() * i + ico->getIndices()[j + 2],
                        i
                    };

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
                if (i >= atomCnt) {                                               //DEBUG CODE FÜR ZUSATZ KUGEL
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
        std::vector<float> torusVertices;
        std::vector<float> torusNormals;
        std::vector<unsigned int> torusIndices;
        std::vector<unsigned int> torusLineIndices;
        std::vector<unsigned int> torusFaceIndices;

        int highestIcosphereIndex = -1;

        // Iteriere durch die Indizes der Icosphere-Vertizes und finde den höchsten Index
        for (int i = 0; i < ico->getIndexCount(); ++i) {
            highestIcosphereIndex = std::max(highestIcosphereIndex, static_cast<int>(ico->getIndices()[i]));
        }

        int offset = vertex.size() / 3;
        int torusFaceOffset = face.size();
        //int torusVertexOffset = highestIcosphereIndex + 1;

        if (offset <= highestIcosphereIndex) {
            offset = highestIcosphereIndex + /* torusVertexOffset +*/ 1;
        }


        int numSegments = 30;
        int numRings = 20;
        for (auto elem : atomCollisions) {

            

            //if (torusVertexOffset <= highestIcosphereIndex) {
            //    torusVertexOffset = highestIcosphereIndex + /* torusVertexOffset +*/ 1;
            //}

            unsigned int atomIndex1 = elem.first;
            unsigned int atomIndex2 = elem.second;

            std::cout << atomIndex1 << "/" << atomIndex2 << std::endl;

            probeRadius = 0.6f;

            float atomRadius1 = mol->AtomTypes()[mol->AtomTypeIndices()[atomIndex1]].Radius();
            float atomRadius2 = mol->AtomTypes()[mol->AtomTypeIndices()[atomIndex2]].Radius();

            glm::vec3 atomPos1(mol->AtomPositions()[3 * atomIndex1 + 0], mol->AtomPositions()[3 * atomIndex1 + 1],
                mol->AtomPositions()[3 * atomIndex1 + 2]);

            glm::vec3 atomPos2(mol->AtomPositions()[3 * atomIndex2 + 0], mol->AtomPositions()[3 * atomIndex2 + 1],
                mol->AtomPositions()[3 * atomIndex2 + 2]);

            float distance = Torus::getDistance(atomPos1, atomPos2);
            glm::vec3 center = Torus::getTorusCenter(atomPos1, atomPos2, atomRadius1, atomRadius2, probeRadius);
            float radius = Torus::getTorusRadius(atomPos1, atomPos2, atomRadius1, atomRadius2, probeRadius);
            float rotationAngle = Torus::getRotationAngle(atomPos1, atomPos2);
            //std::cout << "center: " << center.x << "/" << center.y << "/" << center.z << std::endl;

            glm::vec3 axisUnitVec = Torus::getTorusAxisUnitVec(atomPos1, atomPos2);
            std::cout << "axisUnitVec: " << axisUnitVec.x << "/" << axisUnitVec.y << "/" << axisUnitVec.z << std::endl;

            Torus torus(center, radius, numSegments, numRings);
            torus.generateTorus(probeRadius, axisUnitVec, rotationAngle, offset, torusFaceOffset);

            std::cout << "radius: " << radius << std::endl;

            for (int i = 0; i < atomCnt; ++i) {
                glm::vec3 atomPosition(mol->AtomPositions()[3 * i + 0], mol->AtomPositions()[3 * i + 1], mol->AtomPositions()[3 * i + 2]);
                float atomRadius = mol->AtomTypes()[mol->AtomTypeIndices()[i]].Radius();

                for (auto elem : torus.vertexFaceIndex) {
                    bool first = torus.isPointInsideIco(atomPosition, atomRadius, elem.second[0]);
                    std::cout << first << std::endl;
                }

                //for (size_t i = 0; i < torus.faceIndices.size(); ++i) {
                //    unsigned int elem = torus.faceIndices[i] - offset;

                //    bool first_point_inside = torus.isPointInsideIco(atomPosition, atomRadius, torus.vertexFaceIndex[elem][0]);
                //    bool second_point_inside = torus.isPointInsideIco(atomPosition, atomRadius, torus.vertexFaceIndex[elem][1]);
                //    bool third_point_inside = torus.isPointInsideIco(atomPosition, atomRadius, torus.vertexFaceIndex[elem][2]);

                //    if (first_point_inside || second_point_inside || third_point_inside) {
                //        torus.faceIndices[i] = -1;
                //    }
                //    //std::cout << elem - offset << std::endl;
                //}
            }
            /*for (auto elem : torus.vertexFaceIndex) {
                std::cout << elem.first << std::endl;
            }*/

            //for (int i = 0; i < atomCnt; ++i) {
            //    glm::vec3 atomPosition(mol->AtomPositions()[3 * i + 0],
            //        mol->AtomPositions()[3 * i + 1], mol->AtomPositions()[3 * i + 2]);
            //    float atomRadius = mol->AtomTypes()[mol->AtomTypeIndices()[i]].Radius();
            //    //torus.removeTrianglesInsideSphere(atomPosition, atomRadius, torus.vertices, torus.normals, torus.faceIndices);
            //    torus.removeTrianglesInsideSphere2(atomPosition, atomRadius, torus.vertices, torus.normals, torus.indices, torus.faceIndices);

            //}


            /*std::vector<glm::vec3> updatedVertices;
            std::vector<unsigned int> updatedIndices;*/
            /*
                Prüfen, ob 1/2/3 vertex innerhalb Atom, wenn ja => raus
                ODER
                Prüfen, ob Zentrum von Dreieck innnerhalb Atom, wenn ja => raus
            */



            //for (int i = 0; i < torus.faceIndices.size(); i += 3) {
            //    if (i + 2 < torus.faceIndices.size()) {
            //        int v1Index = torus.faceIndices[i + 0];
            //        int v2Index = torus.faceIndices[i + 1];
            //        int v3Index = torus.faceIndices[i + 2];

            //        // Extrahieren Sie die Face-Indices für jedes Vertex aus der Map
            //        const int& v1FaceIndices = torus.vertexFaceIndices[v1Index];
            //        const int& v2FaceIndices = torus.vertexFaceIndices[v2Index];
            //        const int& v3FaceIndices = torus.vertexFaceIndices[v3Index];

            //        glm::vec3 v1 = torus.vertices[v1FaceIndices];
            //        glm::vec3 v2 = torus.vertices[v2FaceIndices];
            //        glm::vec3 v3 = torus.vertices[v3FaceIndices];

            //        // Überprüfen Sie, ob das Dreieck innerhalb eines Atoms liegt
            //        if (!torus.isTriangleInsideIco(v1, v2, v3, atomPos1, atomRadius1)) {
            //            // Fügen Sie die aktualisierten Vertices und Indices hinzu, wenn das Dreieck nicht entfernt werden soll
            //            updatedVertices.push_back(torus.vertices[v1FaceIndices]);
            //            updatedVertices.push_back(torus.vertices[v2FaceIndices]);
            //            updatedVertices.push_back(torus.vertices[v3FaceIndices]);

            //            updatedIndices.push_back(updatedVertices.size() - 3);
            //            updatedIndices.push_back(updatedVertices.size() - 2);
            //            updatedIndices.push_back(updatedVertices.size() - 1);
            //        }
            //    }
            //}
            /*std::cout << "verticesSize = " << torus.vertices.size() << std::endl;

            std::cout << "faceIndexSize = " << torus.faceIndices.size() << std::endl;*/
            // Aktualisieren Sie den ursprünglichen Vektor mit den aktualisierten Werten
            /*torus.vertices = updatedVertices;
            torus.indices = updatedIndices;*/
            
            for (int i = 0; i < torus.getVertexCount(); i++) {
                /*torusVertices.push_back(torus.getVertices()[i][0]);
                torusVertices.push_back(torus.getVertices()[i][1]);
                torusVertices.push_back(torus.getVertices()[i][2]);*/
                vertex.push_back(torus.getVertices()[i][0]);
                vertex.push_back(torus.getVertices()[i][1]);
                vertex.push_back(torus.getVertices()[i][2]);

                color.push_back(1.0f);
                color.push_back(0.0f);
                color.push_back(0.0f);

                /*torusNormals.push_back(torus.getNormals()[i][0]);
                torusNormals.push_back(torus.getNormals()[i][1]);
                torusNormals.push_back(torus.getNormals()[i][2]);*/
                normal.push_back(torus.getNormals()[i][0]);
                normal.push_back(torus.getNormals()[i][1]);
                normal.push_back(torus.getNormals()[i][2]);
                
            }

            for (int i = 0; i < torus.getIndexCount(); i++) {
                torusIndices.push_back(torus.getIndices()[i]);
            }

            for (int i = 0; i < torus.getLineIndicesCount(); i++) {
                torusLineIndices.push_back(torus.getLineIndices()[i]);
            }

            for (int i = 0; i < torus.getFaceIndicesCount(); i++) {
                /*torusFaceIndices.push_back(torus.getFaceIndices()[i]);*/
                face.push_back(torus.getFaceIndices()[i]);
            }

            offset += torus.getVertexCount();
        }

        
    for (int atom = 0; atom < edgeVerticesPerAtom.size(); ++atom) {
        for (int vertices = 0; vertices < edgeVerticesPerAtom[atom].size(); vertices = vertices + 2) {
            unsigned int a = edgeVerticesPerAtom[atom].at(vertices + 0); // 79
            unsigned int b = edgeVerticesPerAtom[atom].at(vertices + 1); // 80
            std::vector<unsigned int> finde_79 = findNearestVertice(edgeVerticesPerAtom, edgeVerticesPerAtom[atom].at(vertices + 0),
                vertex, atom);
            std::vector<unsigned int> finde_80 = findNearestVertice(edgeVerticesPerAtom, edgeVerticesPerAtom[atom].at(vertices + 1),
                vertex, atom);
            // schau, ob 337 79 verbindung schon da ist
            bool verbindung_Vorhanden = false;
            for (int i = 0; i < 3; ++i) { //liste aus 3 nächsten Vertices
                for (auto & j : referenceToOtherVertice.at(a)) { // Liste aus Vertices mit Verbindung zum aktuellen Vertice
                    if (finde_79.at(i) == std::get<0>(j)) {
                        verbindung_Vorhanden = true;
                    }
                }
            }
            if (!verbindung_Vorhanden){ //Baue Verbindung zum nächsten Vertice, beschreibe das in die Faces liste, mache dreieck fertig
                std::tuple first = {finde_79.at(0), 1};
                std::tuple second = {a, 1};
                referenceToOtherVertice.at(a).emplace_back(first);
                referenceToOtherVertice.at(finde_79.at(0)).emplace_back(second);
                face.push_back(a); // 79
                face.push_back(b);
                face.push_back(finde_79.at(0)); //337
            }

            referenceToOtherVertice.at(79);
            //337, 248
            std::vector<unsigned int> finde_248 = findNearestVertice(edgeVerticesPerAtom, edgeVerticesPerAtom[1].at(4),
                vertex, 1);
            std::vector<unsigned int> finde_341 = findNearestVertice(edgeVerticesPerAtom, edgeVerticesPerAtom[1].at(9),
                vertex, 1);

                bool found = false;
                int index = INT_MIN;

                for (std::tuple<unsigned int, int> i : referenceToOtherVertice.at(79)) {
                    if (std::get<0>(i) == 337) {}
                }

                // Marco: this part was disabled and not able to compile in last commit of Julian Wolf (may still has mistakes!)
                // sets the additional triangles for sewing the Ico spheres together
                if (isStitching) {                    
                    int firstNearestVertex =
                        findNearestVertice(edgeVerticesPerAtom, edgeVerticesPerAtom[atom].at(vertices + 0), vertex, atom)[0];
                    face.push_back(edgeVerticesPerAtom[atom].at(vertices + 0)); // Vertice 1
                    face.push_back(edgeVerticesPerAtom[atom].at(vertices + 1)); // Vertice 2
                    face.push_back(firstNearestVertex);
                }

            }
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
