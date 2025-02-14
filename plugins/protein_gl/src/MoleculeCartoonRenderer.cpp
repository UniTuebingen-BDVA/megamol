/*
 * MoleculeCartoonRenderer.cpp
 *
 * Copyright (C) 2008 by Universitaet Stuttgart (VISUS).
 * Alle Rechte vorbehalten.
 */


#define _USE_MATH_DEFINES 1

#include "MoleculeCartoonRenderer.h"

#include <math.h>
#include <time.h>

#include "mmcore/CoreInstance.h"
#include "mmcore/param/BoolParam.h"
#include "mmcore/param/ColorParam.h"
#include "mmcore/param/EnumParam.h"
#include "mmcore/param/FilePathParam.h"
#include "mmcore/param/FloatParam.h"
#include "mmcore/param/StringParam.h"
#include "mmcore/utility/ColourParser.h"
#include "mmcore/view/Camera.h"
#include "mmcore_gl/utility/ShaderFactory.h"
#include "mmstd/light/DistantLight.h"
#include "protein_calls/ProteinColor.h"
#include "vislib/OutOfRangeException.h"
#include "vislib/String.h"
#include "vislib/StringConverter.h"
#include "vislib/Trace.h"
#include "vislib/assert.h"
#include "vislib/math/Quaternion.h"
#include "vislib/sys/ASCIIFileBuffer.h"
#include "vislib/sys/File.h"
#include "vislib_gl/graphics/gl/IncludeAllGL.h"

using namespace megamol;
using namespace megamol::core;
using namespace megamol::protein;
using namespace megamol::protein_gl;
using namespace megamol::protein_calls;

/*
 * MoleculeCartoonRenderer::MoleculeCartoonRenderer (CTOR)
 */
MoleculeCartoonRenderer::MoleculeCartoonRenderer(void)
        : mmstd_gl::Renderer3DModuleGL()
        , molDataCallerSlot("getdata", "Connects the protein rendering with protein data storage")
        , getLightsSlot("getlights", "Connects the protein rendering with light sources")
        , molRendererCallerSlot("renderMolecule", "Connects the cartoon rendering with another molecule renderer")
        , bsDataCallerSlot("getBindingSites", "Connects the molecule rendering with binding site data storage")
        , renderingModeParam("renderingMode", "Rendering Mode")
        , coloringModeParam0("color::coloringMode0", "The first coloring mode.")
        , coloringModeParam1("color::coloringMode1", "The second coloring mode.")
        , cmWeightParam("color::colorWeighting", "The weighting of the two coloring modes.")
        , stickColoringModeParam("color::stickColoringMode", "Stick Coloring Mode")
        , smoothCartoonColoringParam("color::smoothCartoonColoring", "Use smooth coloring with Cartoon representation")
        , colorTableFileParam("color::colorTableFilename", "The filename of the color table.")
        , minGradColorParam("color::minGradColor", "The color for the minimum value for gradient coloring")
        , midGradColorParam("color::midGradColor", "The color for the middle value for gradient coloring")
        , maxGradColorParam("color::maxGradColor", "The color for the maximum value for gradient coloring")
        , stickRadiusParam("stickRadius", "The radius for stick rendering")
        , interpolParam("posInterpolation", "Enable positional interpolation between frames")
        , compareParam("comparison::compare", "Enable comparing between two different molecules")
        , proteinOnlyParam("proteinOnly", "Render only the protein")
        , tubeRadiusParam("tubeScl", "Scale factor for the tubes when rendering in tubes only mode.")
        , recomputeAlwaysParam("alwaysRecompute", "Shall the positions be recomputed every frame?")
        , currentFrameId(0)
        , atomCount(0) {
    this->molDataCallerSlot.SetCompatibleCall<MolecularDataCallDescription>();
    this->molDataCallerSlot.SetNecessity(core::AbstractCallSlotPresentation::Necessity::SLOT_REQUIRED);
    this->MakeSlotAvailable(&this->molDataCallerSlot);

    this->getLightsSlot.SetCompatibleCall<core::view::light::CallLightDescription>();
    this->getLightsSlot.SetNecessity(core::AbstractCallSlotPresentation::Necessity::SLOT_REQUIRED);
    this->MakeSlotAvailable(&this->getLightsSlot);

    this->molRendererCallerSlot.SetCompatibleCall<mmstd_gl::CallRender3DGLDescription>();
    this->MakeSlotAvailable(&this->molRendererCallerSlot);


    this->bsDataCallerSlot.SetCompatibleCall<BindingSiteCallDescription>();
    this->MakeSlotAvailable(&this->bsDataCallerSlot);

    // fill color table with default values and set the filename param
    std::string filename("colors.txt");
    ProteinColor::ReadColorTableFromFile(filename, this->fileLookupTable);
    this->colorTableFileParam.SetParameter(
        new param::FilePathParam(filename, core::param::FilePathParam::FilePathFlags_::Flag_File_ToBeCreated));
    this->MakeSlotAvailable(&this->colorTableFileParam);

    // coloring modes
    this->currentColoringMode0 = ProteinColor::ColoringMode::CHAIN;
    this->currentColoringMode1 = ProteinColor::ColoringMode::SECONDARY_STRUCTURE;
    param::EnumParam* cm0 = new param::EnumParam(int(this->currentColoringMode0));
    param::EnumParam* cm1 = new param::EnumParam(int(this->currentColoringMode1));
    param::EnumParam* scm = new param::EnumParam(int(ProteinColor::ColoringMode::ELEMENT));
    MolecularDataCall* mol = new MolecularDataCall();
    BindingSiteCall* bs = new BindingSiteCall();
    unsigned int cCnt;
    ProteinColor::ColoringMode cMode;
    for (cCnt = 0; cCnt < static_cast<int>(ProteinColor::ColoringMode::MODE_COUNT); ++cCnt) {
        cMode = static_cast<ProteinColor::ColoringMode>(cCnt);
        cm0->SetTypePair(static_cast<int>(cMode), ProteinColor::GetName(cMode).c_str());
        cm1->SetTypePair(static_cast<int>(cMode), ProteinColor::GetName(cMode).c_str());
        scm->SetTypePair(static_cast<int>(cMode), ProteinColor::GetName(cMode).c_str());
    }
    delete mol;
    delete bs;
    this->coloringModeParam0 << cm0;
    this->coloringModeParam1 << cm1;
    this->stickColoringModeParam << scm;
    this->MakeSlotAvailable(&this->coloringModeParam0);
    this->MakeSlotAvailable(&this->coloringModeParam1);
    this->MakeSlotAvailable(&this->stickColoringModeParam);

    // Color weighting parameter
    this->cmWeightParam.SetParameter(new param::FloatParam(0.5f, 0.0f, 1.0f));
    this->MakeSlotAvailable(&this->cmWeightParam);

    // --- set the render mode ---
    // SetRenderMode(CARTOON);
    // SetRenderMode(CARTOON_SIMPLE);
    // SetRenderMode(CartoonRenderMode::CARTOON_LINE);
    this->currentRenderMode = CartoonRenderMode::CARTOON_GPU;
    param::EnumParam* rm = new param::EnumParam(int(this->currentRenderMode));
    rm->SetTypePair(static_cast<int>(CartoonRenderMode::CARTOON), "Cartoon Hybrid");
    rm->SetTypePair(static_cast<int>(CartoonRenderMode::CARTOON_SIMPLE), "Cartoon Hybrid (simple)");
    rm->SetTypePair(static_cast<int>(CartoonRenderMode::CARTOON_GPU), "Cartoon GPU");
    rm->SetTypePair(static_cast<int>(CartoonRenderMode::CARTOON_CPU), "Cartoon CPU");
    rm->SetTypePair(static_cast<int>(CartoonRenderMode::CARTOON_LINE), "Cartoon Lines");
    this->renderingModeParam << rm;
    this->MakeSlotAvailable(&this->renderingModeParam);

    // --- set smooth coloring for cartoon rendering ---
    this->smoothCartoonColoringMode = false;
    this->smoothCartoonColoringParam.SetParameter(new param::BoolParam(this->smoothCartoonColoringMode));
    this->MakeSlotAvailable(&this->smoothCartoonColoringParam);

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
    this->stickRadiusParam.SetParameter(new param::FloatParam(0.3f, 0.0f));
    this->MakeSlotAvailable(&this->stickRadiusParam);

    this->compareParam.SetParameter(new param::BoolParam(false));
    this->MakeSlotAvailable(&this->compareParam);

    this->proteinOnlyParam.SetParameter(new param::BoolParam(false));
    this->MakeSlotAvailable(&this->proteinOnlyParam);

    // --- set the radius for the cartoon rednering mode ---
    this->radiusCartoon = 0.2f;

    // --- initialize all pointers and variables for cartoon ---
    this->vertTube = new float[1];
    this->colorsParamsTube = new float[1];
    this->vertArrow = new float[1];
    this->colorsParamsArrow = new float[1];
    this->vertHelix = new float[1];
    this->colorsParamsHelix = new float[1];
    this->normalArrow = new float[1];
    this->normalHelix = new float[1];
    this->normalTube = new float[1];

    // hybrid CARTOON render mode was not prepared yet
    this->prepareCartoonHybrid = true;
    // CPU CARTOON render mode was not prepared yet
    this->prepareCartoonCPU = true;
    // CARTOON LINE render mode was not prepared yet
    this->prepareCartoonLine = true;

    // set default value for spline segments per amino acid
    this->numberOfSplineSeg = 6;
    // set default value for tube segments
    this->numberOfTubeSeg = 6;

    // fill rainbow color table
    ProteinColor::MakeRainbowColorTable(100, this->rainbowColors);

    // en-/disable positional interpolation
    this->interpolParam.SetParameter(new param::BoolParam(true));
    this->MakeSlotAvailable(&this->interpolParam);

    // Tubes scale parameter
    this->tubeRadiusParam.SetParameter(new param::FloatParam(1.0f, 0.0f));
    this->MakeSlotAvailable(&this->tubeRadiusParam);

    this->recomputeAlwaysParam.SetParameter(new param::BoolParam(false));
    this->MakeSlotAvailable(&this->recomputeAlwaysParam);
}


/*
 * MoleculeCartoonRenderer::~MoleculeCartoonRenderer (DTOR)
 */
MoleculeCartoonRenderer::~MoleculeCartoonRenderer(void) {
    this->Release();
}


/*
 * MoleculeCartoonRenderer::release
 */
void MoleculeCartoonRenderer::release(void) {}


/*
 * MoleculeCartoonRenderer::create
 */
bool MoleculeCartoonRenderer::create(void) {
    using megamol::core::utility::log::Log;

    {
        // TODO Note from shader factory migration:
        // The old shaders did not set a #version, which is required with the shader factory. The problem is, geometry
        // shaders require at least version 150, but the code within the geometry shaders uses a lot of deprecated
        // functions only available before version 150. So some refactoring is needed within this module.
        // Please do this, whoever wants to use this stuff!
        Log::DefaultLog.WriteError("MoleculeCartoonRenderer is currently broken and requires some refactoring after "
                                   "the shader factory migration. If you are a maintainer or user of this module "
                                   "please do this. For more info see comment above this error message in code.");
        return false;
    }


    glEnable(GL_DEPTH_TEST);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE_ARB);
    glEnable(GL_VERTEX_PROGRAM_TWO_SIDE);

    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE_ARB);

    using namespace vislib_gl::graphics::gl;

    auto const shader_options = msf::ShaderFactoryOptionsOpenGL(GetCoreInstance()->GetShaderPaths());

    try {
        ////////////////////////////////////////////////////
        // load the shader sources for the cartoon shader //
        ////////////////////////////////////////////////////

        this->cartoonShader = core::utility::make_glowl_shader("cartoonShader", shader_options,
            "protein_gl/protein/cartoon_cartoon.vert.glsl", "protein_gl/protein/cartoon_cartoon.geom.glsl",
            "protein_gl/protein/cartoon_cartoon.frag.glsl");
        // set GL_TRIANGLES_ADJACENCY_EXT primitives as INPUT
        glProgramParameteriEXT(cartoonShader->getHandle(), GL_GEOMETRY_INPUT_TYPE_EXT, GL_TRIANGLES_ADJACENCY_EXT);
        // set TRIANGLE_STRIP as OUTPUT
        glProgramParameteriEXT(cartoonShader->getHandle(), GL_GEOMETRY_OUTPUT_TYPE_EXT, GL_TRIANGLE_STRIP);
        // set maximum number of vertices to be generated by geometry shader to
        glProgramParameteriEXT(cartoonShader->getHandle(), GL_GEOMETRY_VERTICES_OUT_EXT, 200);

        // TODO Note from shader factory migration: No idea if relink is ok or setting parameter must happen before initial link within glowl.
        glLinkProgram(cartoonShader->getHandle());

        /////////////////////////////////////////////////
        // load the shader sources for the tube shader //
        /////////////////////////////////////////////////

        this->tubeShader = core::utility::make_glowl_shader("tubeShader", shader_options,
            "protein_gl/protein/cartoon_cartoon.vert.glsl", "protein_gl/protein/cartoon_tube.geom.glsl",
            "protein_gl/protein/cartoon_cartoon.frag.glsl");
        glProgramParameteriEXT(tubeShader->getHandle(), GL_GEOMETRY_INPUT_TYPE_EXT, GL_TRIANGLES_ADJACENCY_EXT);
        glProgramParameteriEXT(tubeShader->getHandle(), GL_GEOMETRY_OUTPUT_TYPE_EXT, GL_TRIANGLE_STRIP);
        glProgramParameteriEXT(tubeShader->getHandle(), GL_GEOMETRY_VERTICES_OUT_EXT, 200);

        // TODO Note from shader factory migration: No idea if relink is ok or setting parameter must happen before initial link within glowl.
        glLinkProgram(tubeShader->getHandle());

        //////////////////////////////////////////////////
        // load the shader sources for the arrow shader //
        //////////////////////////////////////////////////

        this->arrowShader = core::utility::make_glowl_shader("arrowShader", shader_options,
            "protein_gl/protein/cartoon_cartoon.vert.glsl", "protein_gl/protein/cartoon_arrow.geom.glsl",
            "protein_gl/protein/cartoon_cartoon.frag.glsl");
        glProgramParameteriEXT(arrowShader->getHandle(), GL_GEOMETRY_INPUT_TYPE_EXT, GL_TRIANGLES_ADJACENCY_EXT);
        glProgramParameteriEXT(arrowShader->getHandle(), GL_GEOMETRY_OUTPUT_TYPE_EXT, GL_TRIANGLE_STRIP);
        glProgramParameteriEXT(arrowShader->getHandle(), GL_GEOMETRY_VERTICES_OUT_EXT, 200);

        // TODO Note from shader factory migration: No idea if relink is ok or setting parameter must happen before initial link within glowl.
        glLinkProgram(arrowShader->getHandle());

        /////////////////////////////////////////////////
        // load the shader sources for the helix shader //
        /////////////////////////////////////////////////

        this->helixShader = core::utility::make_glowl_shader("helixShader", shader_options,
            "protein_gl/protein/cartoon_cartoon.vert.glsl", "protein_gl/protein/cartoon_helix.geom.glsl",
            "protein_gl/protein/cartoon_cartoon.frag.glsl");
        glProgramParameteriEXT(helixShader->getHandle(), GL_GEOMETRY_INPUT_TYPE_EXT, GL_TRIANGLES_ADJACENCY_EXT);
        glProgramParameteriEXT(helixShader->getHandle(), GL_GEOMETRY_OUTPUT_TYPE_EXT, GL_TRIANGLE_STRIP);
        glProgramParameteriEXT(helixShader->getHandle(), GL_GEOMETRY_VERTICES_OUT_EXT, 200);

        // TODO Note from shader factory migration: No idea if relink is ok or setting parameter must happen before initial link within glowl.
        glLinkProgram(helixShader->getHandle());

        /////////////////////////////////////////////////
        // load the shader sources for the tube shader //
        /////////////////////////////////////////////////

        this->tubeSimpleShader = core::utility::make_glowl_shader("tubeSimpleShader", shader_options,
            "protein_gl/protein/cartoon_simple.vert.glsl", "protein_gl/protein/cartoon_simple_tube.geom.glsl",
            "protein_gl/protein/cartoon_simple.frag.glsl");
        glProgramParameteriEXT(tubeSimpleShader->getHandle(), GL_GEOMETRY_INPUT_TYPE_EXT, GL_TRIANGLES_ADJACENCY_EXT);
        glProgramParameteriEXT(tubeSimpleShader->getHandle(), GL_GEOMETRY_OUTPUT_TYPE_EXT, GL_TRIANGLE_STRIP);
        glProgramParameteriEXT(tubeSimpleShader->getHandle(), GL_GEOMETRY_VERTICES_OUT_EXT, 200);

        // TODO Note from shader factory migration: No idea if relink is ok or setting parameter must happen before initial link within glowl.
        glLinkProgram(tubeSimpleShader->getHandle());

        //////////////////////////////////////////////////
        // load the shader sources for the arrow shader //
        //////////////////////////////////////////////////

        this->arrowSimpleShader = core::utility::make_glowl_shader("arrowSimpleShader", shader_options,
            "protein_gl/protein/cartoon_simple.vert.glsl", "protein_gl/protein/cartoon_simple_arrow.geom.glsl",
            "protein_gl/protein/cartoon_simple.frag.glsl");
        glProgramParameteriEXT(arrowSimpleShader->getHandle(), GL_GEOMETRY_INPUT_TYPE_EXT, GL_TRIANGLES_ADJACENCY_EXT);
        glProgramParameteriEXT(arrowSimpleShader->getHandle(), GL_GEOMETRY_OUTPUT_TYPE_EXT, GL_TRIANGLE_STRIP);
        glProgramParameteriEXT(arrowSimpleShader->getHandle(), GL_GEOMETRY_VERTICES_OUT_EXT, 200);

        // TODO Note from shader factory migration: No idea if relink is ok or setting parameter must happen before initial link within glowl.
        glLinkProgram(arrowSimpleShader->getHandle());

        /////////////////////////////////////////////////
        // load the shader sources for the helix shader //
        /////////////////////////////////////////////////

        this->helixSimpleShader = core::utility::make_glowl_shader("helixSimpleShader", shader_options,
            "protein_gl/protein/cartoon_simple.vert.glsl", "protein_gl/protein/cartoon_simple_helix.geom.glsl",
            "protein_gl/protein/cartoon_simple.frag.glsl");
        glProgramParameteriEXT(helixSimpleShader->getHandle(), GL_GEOMETRY_INPUT_TYPE_EXT, GL_TRIANGLES_ADJACENCY_EXT);
        glProgramParameteriEXT(helixSimpleShader->getHandle(), GL_GEOMETRY_OUTPUT_TYPE_EXT, GL_TRIANGLE_STRIP);
        glProgramParameteriEXT(helixSimpleShader->getHandle(), GL_GEOMETRY_VERTICES_OUT_EXT, 200);

        // TODO Note from shader factory migration: No idea if relink is ok or setting parameter must happen before initial link within glowl.
        glLinkProgram(helixSimpleShader->getHandle());

        /////////////////////////////////////////////////////////
        // load the shader sources for the spline arrow shader //
        /////////////////////////////////////////////////////////

        this->arrowSplineShader = core::utility::make_glowl_shader("arrowSplineShader", shader_options,
            "protein_gl/protein/cartoon_spline.vert.glsl", "protein_gl/protein/cartoon_spline_arrow.geom.glsl",
            "protein_gl/protein/cartoon_spline.frag.glsl");
        glProgramParameteriEXT(arrowSplineShader->getHandle(), GL_GEOMETRY_INPUT_TYPE_EXT, GL_LINES_ADJACENCY_EXT);
        glProgramParameteriEXT(arrowSplineShader->getHandle(), GL_GEOMETRY_OUTPUT_TYPE_EXT, GL_TRIANGLE_STRIP);
        glProgramParameteriEXT(arrowSplineShader->getHandle(), GL_GEOMETRY_VERTICES_OUT_EXT, 1024);

        // TODO Note from shader factory migration: No idea if relink is ok or setting parameter must happen before initial link within glowl.
        glLinkProgram(arrowSplineShader->getHandle());

        ////////////////////////////////////////////////////////
        // load the shader sources for the spline tube shader //
        ////////////////////////////////////////////////////////

        this->tubeSplineShader = core::utility::make_glowl_shader("tubeSplineShader", shader_options,
            "protein_gl/protein/cartoon_spline.vert.glsl", "protein_gl/protein/cartoon_spline_tube.geom.glsl",
            "protein_gl/protein/cartoon_spline.frag.glsl");
        glProgramParameteriEXT(tubeSplineShader->getHandle(), GL_GEOMETRY_INPUT_TYPE_EXT, GL_LINES_ADJACENCY_EXT);
        glProgramParameteriEXT(tubeSplineShader->getHandle(), GL_GEOMETRY_OUTPUT_TYPE_EXT, GL_TRIANGLE_STRIP);
        glProgramParameteriEXT(tubeSplineShader->getHandle(), GL_GEOMETRY_VERTICES_OUT_EXT, 1024);

        // TODO Note from shader factory migration: No idea if relink is ok or setting parameter must happen before initial link within glowl.
        glLinkProgram(tubeSplineShader->getHandle());

        ////////////////////////////////////////////////////////
        // load the shader sources for the spline helix shader //
        ////////////////////////////////////////////////////////

        this->helixSplineShader = core::utility::make_glowl_shader("helixSplineShader", shader_options,
            "protein_gl/protein/cartoon_spline.vert.glsl", "protein_gl/protein/cartoon_spline_helix.geom.glsl",
            "protein_gl/protein/cartoon_spline.frag.glsl");
        glProgramParameteriEXT(helixSplineShader->getHandle(), GL_GEOMETRY_INPUT_TYPE_EXT, GL_LINES_ADJACENCY_EXT);
        glProgramParameteriEXT(helixSplineShader->getHandle(), GL_GEOMETRY_OUTPUT_TYPE_EXT, GL_TRIANGLE_STRIP);
        glProgramParameteriEXT(helixSplineShader->getHandle(), GL_GEOMETRY_VERTICES_OUT_EXT, 1024);

        // TODO Note from shader factory migration: No idea if relink is ok or setting parameter must happen before initial link within glowl.
        glLinkProgram(helixSplineShader->getHandle());

        //////////////////////////////////////////////////////
        // load the shader files for the per pixel lighting //
        //////////////////////////////////////////////////////

        this->lightShader = core::utility::make_glowl_shader("lightShader", shader_options,
            "protein_gl/protein/cartoon_perpixellight.vert.glsl", "protein_gl/protein/cartoon_perpixellight.frag.glsl");


        sphereShader_ = core::utility::make_shared_glowl_shader("sphere", shader_options,
            std::filesystem::path("protein_gl/simplemolecule/sm_sphere.vert.glsl"),
            std::filesystem::path("protein_gl/simplemolecule/sm_sphere.frag.glsl"));

        cylinderShader_ = core::utility::make_shared_glowl_shader("cylinder", shader_options,
            std::filesystem::path("protein_gl/simplemolecule/sm_cylinder.vert.glsl"),
            std::filesystem::path("protein_gl/simplemolecule/sm_cylinder.frag.glsl"));

    } catch (std::exception& e) {
        Log::DefaultLog.WriteError(("MoleculeCartoonRenderer: " + std::string(e.what())).c_str());
        return false;
    }

    // generate data buffers
    buffers_[static_cast<int>(Buffers::POSITION)] =
        std::make_unique<glowl::BufferObject>(GL_ARRAY_BUFFER, nullptr, 0, GL_DYNAMIC_DRAW);
    buffers_[static_cast<int>(Buffers::COLOR)] =
        std::make_unique<glowl::BufferObject>(GL_ARRAY_BUFFER, nullptr, 0, GL_DYNAMIC_DRAW);
    buffers_[static_cast<int>(Buffers::CYL_PARAMS)] =
        std::make_unique<glowl::BufferObject>(GL_ARRAY_BUFFER, nullptr, 0, GL_DYNAMIC_DRAW);
    buffers_[static_cast<int>(Buffers::CYL_QUAT)] =
        std::make_unique<glowl::BufferObject>(GL_ARRAY_BUFFER, nullptr, 0, GL_DYNAMIC_DRAW);
    buffers_[static_cast<int>(Buffers::CYL_COL1)] =
        std::make_unique<glowl::BufferObject>(GL_ARRAY_BUFFER, nullptr, 0, GL_DYNAMIC_DRAW);
    buffers_[static_cast<int>(Buffers::CYL_COL2)] =
        std::make_unique<glowl::BufferObject>(GL_ARRAY_BUFFER, nullptr, 0, GL_DYNAMIC_DRAW);
    buffers_[static_cast<int>(Buffers::FILTER)] =
        std::make_unique<glowl::BufferObject>(GL_ARRAY_BUFFER, nullptr, 0, GL_DYNAMIC_DRAW);
    buffers_[static_cast<int>(Buffers::LIGHT_POSITIONAL)] =
        std::make_unique<glowl::BufferObject>(GL_SHADER_STORAGE_BUFFER, nullptr, 0, GL_DYNAMIC_DRAW);
    buffers_[static_cast<int>(Buffers::LIGHT_DIRECTIONAL)] =
        std::make_unique<glowl::BufferObject>(GL_SHADER_STORAGE_BUFFER, nullptr, 0, GL_DYNAMIC_DRAW);

    glGenVertexArrays(1, &vertex_array_spheres_);
    glBindVertexArray(vertex_array_spheres_);

    buffers_[static_cast<int>(Buffers::POSITION)]->bind();
    glEnableVertexAttribArray(static_cast<int>(Buffers::POSITION));
    glVertexAttribPointer(static_cast<int>(Buffers::POSITION), 4, GL_FLOAT, GL_FALSE, 0, nullptr);

    buffers_[static_cast<int>(Buffers::COLOR)]->bind();
    glEnableVertexAttribArray(static_cast<int>(Buffers::COLOR));
    glVertexAttribPointer(static_cast<int>(Buffers::COLOR), 3, GL_FLOAT, GL_FALSE, 0, nullptr);

    buffers_[static_cast<int>(Buffers::CYL_PARAMS)]->bind();
    glEnableVertexAttribArray(static_cast<int>(Buffers::CYL_PARAMS));
    glVertexAttribPointer(static_cast<int>(Buffers::CYL_PARAMS), 2, GL_FLOAT, GL_FALSE, 0, nullptr);

    buffers_[static_cast<int>(Buffers::CYL_QUAT)]->bind();
    glEnableVertexAttribArray(static_cast<int>(Buffers::CYL_QUAT));
    glVertexAttribPointer(static_cast<int>(Buffers::CYL_QUAT), 4, GL_FLOAT, GL_FALSE, 0, nullptr);

    buffers_[static_cast<int>(Buffers::CYL_COL1)]->bind();
    glEnableVertexAttribArray(static_cast<int>(Buffers::CYL_COL1));
    glVertexAttribPointer(static_cast<int>(Buffers::CYL_COL1), 3, GL_FLOAT, GL_FALSE, 0, nullptr);

    buffers_[static_cast<int>(Buffers::CYL_COL2)]->bind();
    glEnableVertexAttribArray(static_cast<int>(Buffers::CYL_COL2));
    glVertexAttribPointer(static_cast<int>(Buffers::CYL_COL2), 3, GL_FLOAT, GL_FALSE, 0, nullptr);

    buffers_[static_cast<int>(Buffers::FILTER)]->bind();
    glEnableVertexAttribArray(static_cast<int>(Buffers::FILTER));
    glVertexAttribIPointer(static_cast<int>(Buffers::FILTER), 1, GL_INT, 0, nullptr);

    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glDisableVertexAttribArray(static_cast<int>(Buffers::POSITION));
    glDisableVertexAttribArray(static_cast<int>(Buffers::COLOR));
    glDisableVertexAttribArray(static_cast<int>(Buffers::CYL_PARAMS));
    glDisableVertexAttribArray(static_cast<int>(Buffers::CYL_QUAT));
    glDisableVertexAttribArray(static_cast<int>(Buffers::CYL_COL1));
    glDisableVertexAttribArray(static_cast<int>(Buffers::CYL_COL2));
    glDisableVertexAttribArray(static_cast<int>(Buffers::FILTER));

    return true;
}


/*
 * MoleculeCartoonRenderer::GetExtents
 */
bool MoleculeCartoonRenderer::GetExtents(mmstd_gl::CallRender3DGL& call) {
    core::view::AbstractCallRender* cr3d = dynamic_cast<core::view::AbstractCallRender*>(&call);
    if (cr3d == NULL)
        return false;

    MolecularDataCall* mol = this->molDataCallerSlot.CallAs<MolecularDataCall>();

    if (mol == NULL)
        return false;
    mol->SetFrameID(static_cast<unsigned int>(cr3d->Time()));
    if (!(*mol)(1))
        return false;

    float scale;
    if (!vislib::math::IsEqual(mol->AccessBoundingBoxes().ObjectSpaceBBox().LongestEdge(), 0.0f)) {
        scale = 2.0f / mol->AccessBoundingBoxes().ObjectSpaceBBox().LongestEdge();
    } else {
        scale = 1.0f;
    }

    cr3d->AccessBoundingBoxes() = mol->AccessBoundingBoxes();
    cr3d->SetTimeFramesCount(mol->FrameCount());

    // Get the pointer to CallRender3D (protein renderer) or CallRenderDeferred3D
    // if offscreen rendering is enabled

    mmstd_gl::CallRender3DGL* molrencr3d = this->molRendererCallerSlot.CallAs<mmstd_gl::CallRender3DGL>();
    if (molrencr3d) {
        (*molrencr3d)(core::view::AbstractCallRender::FnGetExtents);
    }

    return true;
}


/**********************************************************************
 * 'render'-functions
 **********************************************************************/

/*
 * MoleculeCartoonRenderer::Render
 */
bool MoleculeCartoonRenderer::Render(mmstd_gl::CallRender3DGL& call) {
    // cast the call to Render3D
    core::view::AbstractCallRender* cr3d = dynamic_cast<core::view::AbstractCallRender*>(&call);
    if (cr3d == NULL)
        return false;

    // get the pointer to AbstractCallRender (molecule renderer)
    core::view::AbstractCallRender* molrencr3d;
    molrencr3d = this->molRendererCallerSlot.CallAs<core::view::AbstractCallRender>();

    // get pointer to BindingSiteCall
    BindingSiteCall* bs = this->bsDataCallerSlot.CallAs<BindingSiteCall>();
    if (bs) {
        (*bs)(BindingSiteCall::CallForGetData);
    }

    // get camera information
    camera = call.GetCamera();
    view = camera.getViewMatrix();
    proj = camera.getProjectionMatrix();
    MVinv = glm::inverse(view);
    invProj = glm::inverse(proj);
    NormalM = glm::transpose(MVinv);
    MVP = proj * view;
    MVPinv = glm::inverse(MVP);
    MVPtransp = glm::transpose(MVP);

    planes = glm::vec2(0.0f, 0.0f);
    try {
        auto cam_intrinsics = camera.get<megamol::core::view::Camera::PerspectiveParameters>();
        planes.x = cam_intrinsics.near_plane;
        planes.y = cam_intrinsics.far_plane;
    } catch (...) {
        megamol::core::utility::log::Log::DefaultLog.WriteError(
            "MoleculeCartoonRenderer - Error when getting perspective camera intrinsics");
    }


    // get framebuffer information
    fbo = call.GetFramebuffer();

    // =============== Protein Rendering ===============
    if (molrencr3d) {
        // setup and call molecule renderer
        glPushMatrix();
        *molrencr3d = *cr3d;
        (*molrencr3d)();
        glPopMatrix();
    }

    float callTime = cr3d->Time();

    // get pointer to MolecularDataCall
    MolecularDataCall* mol = this->molDataCallerSlot.CallAs<MolecularDataCall>();


    if (mol == NULL)
        return false;

    mol->SetFrameID(static_cast<int>(callTime));
    if (!(*mol)(MolecularDataCall::CallForGetExtent))
        return false;
    if (!(*mol)(MolecularDataCall::CallForGetData))
        return false;

    // TODO store c-alpha of mol
    // TODO if mol2 == NULL : callForGetData with callTime + frameOffset of mol ; else : use mol2
    // TODO store c-alpha of mol2 (or new c-alpha of mol)


    // interpolate between frames ...
    int cnt;

    float* pos0 = new float[mol->AtomCount() * 3];
    memcpy(pos0, mol->AtomPositions(), mol->AtomCount() * 3 * sizeof(float));

    // set next frame ID and get positions of the second frame
    if (((static_cast<int>(callTime) + 1) < int(mol->FrameCount())) &&
        this->interpolParam.Param<param::BoolParam>()->Value())
        mol->SetFrameID(static_cast<int>(callTime) + 1);
    else
        mol->SetFrameID(static_cast<int>(callTime));
    if (!(*mol)(MolecularDataCall::CallForGetData)) {
        delete[] pos0;
        return false;
    }
    float* pos1 = new float[mol->AtomCount() * 3];
    memcpy(pos1, mol->AtomPositions(), mol->AtomCount() * 3 * sizeof(float));

    // interpolate atom positions between frames
    float* posInter = new float[mol->AtomCount() * 3];
    float inter = callTime - static_cast<float>(static_cast<int>(callTime));
    float threshold = vislib::math::Min(mol->AccessBoundingBoxes().ObjectSpaceBBox().Width(),
                          vislib::math::Min(mol->AccessBoundingBoxes().ObjectSpaceBBox().Height(),
                              mol->AccessBoundingBoxes().ObjectSpaceBBox().Depth())) *
                      0.75f;
#pragma omp parallel for
    for (cnt = 0; cnt < static_cast<int>(mol->AtomCount()); ++cnt) {
        if (std::sqrt(std::pow(pos0[3 * cnt + 0] - pos1[3 * cnt + 0], 2) +
                      std::pow(pos0[3 * cnt + 1] - pos1[3 * cnt + 1], 2) +
                      std::pow(pos0[3 * cnt + 2] - pos1[3 * cnt + 2], 2)) < threshold) {
            posInter[3 * cnt + 0] = (1.0f - inter) * pos0[3 * cnt + 0] + inter * pos1[3 * cnt + 0];
            posInter[3 * cnt + 1] = (1.0f - inter) * pos0[3 * cnt + 1] + inter * pos1[3 * cnt + 1];
            posInter[3 * cnt + 2] = (1.0f - inter) * pos0[3 * cnt + 2] + inter * pos1[3 * cnt + 2];
        } else if (inter < 0.5f) {
            posInter[3 * cnt + 0] = pos0[3 * cnt + 0];
            posInter[3 * cnt + 1] = pos0[3 * cnt + 1];
            posInter[3 * cnt + 2] = pos0[3 * cnt + 2];
        } else {
            posInter[3 * cnt + 0] = pos1[3 * cnt + 0];
            posInter[3 * cnt + 1] = pos1[3 * cnt + 1];
            posInter[3 * cnt + 2] = pos1[3 * cnt + 2];
        }
    }
    // ... interpolate between frames

    // check if the frame has changed
    if (this->currentFrameId != mol->FrameID()) {
        this->currentFrameId = mol->FrameID();
        this->RecomputeAll();
    }
    // check if the call time has changed
    if (this->oldCallTime != callTime) {
        this->oldCallTime = callTime;
        this->RecomputeAll();
    }

    // check last atom count with current atom count
    if (this->atomCount != mol->AtomCount()) {
        this->atomCount = mol->AtomCount();
        this->RecomputeAll();
    }

    if (this->recomputeAlwaysParam.Param<param::BoolParam>()->Value()) {
        this->atomCount = mol->AtomCount();
        this->RecomputeAll();
    }

    // force recomputation
    // this->RecomputeAll();

    this->colorLookupTable = {glm::make_vec3(this->minGradColorParam.Param<param::ColorParam>()->Value().data()),
        glm::make_vec3(this->midGradColorParam.Param<param::ColorParam>()->Value().data()),
        glm::make_vec3(this->maxGradColorParam.Param<param::ColorParam>()->Value().data())};

    // recompute colors
    // TODO Why is this done in every frame?
    ProteinColor::MakeWeightedColorTable(*mol, this->currentColoringMode0, this->currentColoringMode1,
        cmWeightParam.Param<param::FloatParam>()->Value(),        // weight for the first cm
        1.0f - cmWeightParam.Param<param::FloatParam>()->Value(), // weight for the second cm
        this->atomColorTable, this->colorLookupTable, this->fileLookupTable, this->rainbowColors, bs);

    // parameter refresh
    // Note this also recomputes the atom color table using the color call if necessary
    this->UpdateParameters(mol, mol->FrameID(), bs);


    // render...
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
    glDepthFunc(GL_LEQUAL);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    glEnable(GL_VERTEX_PROGRAM_TWO_SIDE);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glLoadMatrixf(glm::value_ptr(view));

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glLoadMatrixf(glm::value_ptr(proj));

    glDisable(GL_BLEND);
    glDisable(GL_CULL_FACE);

    if (this->currentRenderMode == CartoonRenderMode::CARTOON ||
        this->currentRenderMode == CartoonRenderMode::CARTOON_SIMPLE) {
        // ------------------------------------------------------------
        // --- CARTOON                                              ---
        // --- Hybrid Implementation using GLSL geometry shaders    ---
        // ------------------------------------------------------------
        this->RenderCartoonHybrid(mol, posInter);
    }

    if (this->currentRenderMode == CartoonRenderMode::CARTOON_CPU) {
        // ------------------------------------------------------------
        // --- CARTOON_CPU                                          ---
        // --- render the protein using OpenGL primitives           ---
        // ------------------------------------------------------------
        this->RenderCartoonCPU(mol, posInter);
    }

    if (this->currentRenderMode == CartoonRenderMode::CARTOON_LINE) {
        // ------------------------------------------------------------
        // --- CARTOON_LINE                                         ---
        // --- render the protein using OpenGL lines                ---
        // ------------------------------------------------------------
        this->RenderCartoonLineCPU(mol, posInter);
    }

    if (this->currentRenderMode == CartoonRenderMode::CARTOON_GPU) {
        // ------------------------------------------------------------
        // --- CARTOON_GPU                                          ---
        // --- render the protein using only GLSL geometry shaders  ---
        // ------------------------------------------------------------
        this->RenderCartoonGPU(mol, posInter);
    }

    if (this->currentRenderMode == CartoonRenderMode::CARTOON_TUBE_ONLY) {
        // ------------------------------------------------------------
        // --- CARTOON_TUBE_ONLY                                          ---
        // --- render the protein using only GLSL geometry shaders  ---
        // ------------------------------------------------------------
        this->RenderCartoonGPUTubeOnly(mol, posInter);
    }

    glDisable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS); // default depth function
    glDisable(GL_VERTEX_PROGRAM_POINT_SIZE);

    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();


    if (!this->proteinOnlyParam.Param<param::BoolParam>()->Value()) {
        // coloring mode for other molecules
        ProteinColor::MakeColorTable(*mol,
            static_cast<ProteinColor::ColoringMode>(
                int(this->stickColoringModeParam.Param<param::EnumParam>()->Value())),
            this->atomColorTable, this->colorLookupTable, this->fileLookupTable, this->rainbowColors, bs, nullptr,
            true);
        // render rest as stick
        this->RenderStick(mol, posInter, bs);
        // reset coloring mode
        this->currentColoringMode0 =
            static_cast<ProteinColor::ColoringMode>(int(this->coloringModeParam0.Param<param::EnumParam>()->Value()));
        this->currentColoringMode1 =
            static_cast<ProteinColor::ColoringMode>(int(this->coloringModeParam1.Param<param::EnumParam>()->Value()));
        ProteinColor::MakeWeightedColorTable(*mol, this->currentColoringMode0, this->currentColoringMode1,
            cmWeightParam.Param<param::FloatParam>()->Value(), 1.0f - cmWeightParam.Param<param::FloatParam>()->Value(),
            this->atomColorTable, this->colorLookupTable, this->fileLookupTable, this->rainbowColors, bs, nullptr,
            true);
    }

    delete[] pos0;
    delete[] pos1;
    delete[] posInter;

    // unlock the current frame
    mol->Unlock();

    return true;
}

/*
 * update parameters
 */
void MoleculeCartoonRenderer::UpdateParameters(
    MolecularDataCall* mol, unsigned int frameID, const BindingSiteCall* bs) {
    // color table param
    if (this->colorTableFileParam.IsDirty()) {
        ProteinColor::ReadColorTableFromFile(
            this->colorTableFileParam.Param<param::FilePathParam>()->Value(), this->fileLookupTable);
        this->colorTableFileParam.ResetDirty();
    }
    // parameter refresh
    if (this->renderingModeParam.IsDirty()) {
        this->currentRenderMode =
            static_cast<CartoonRenderMode>(int(this->renderingModeParam.Param<param::EnumParam>()->Value()));
        this->renderingModeParam.ResetDirty();
    }

    if (this->tubeRadiusParam.IsDirty()) {
        this->prepareCartoonHybrid = true;
        this->tubeRadiusParam.ResetDirty();
    }
    if (this->coloringModeParam0.IsDirty() || this->coloringModeParam1.IsDirty() || this->cmWeightParam.IsDirty() ||
        this->compareParam.IsDirty()) {

        this->currentColoringMode0 =
            static_cast<ProteinColor::ColoringMode>(int(this->coloringModeParam0.Param<param::EnumParam>()->Value()));
        this->currentColoringMode1 =
            static_cast<ProteinColor::ColoringMode>(int(this->coloringModeParam1.Param<param::EnumParam>()->Value()));
        this->compare = this->compareParam.Param<param::BoolParam>()->Value();
        RecomputeAll();

        if (!this->compare) {
            ProteinColor::MakeWeightedColorTable(*mol, this->currentColoringMode0, this->currentColoringMode1,
                cmWeightParam.Param<param::FloatParam>()->Value(),
                1.0f - cmWeightParam.Param<param::FloatParam>()->Value(), this->atomColorTable, this->colorLookupTable,
                this->fileLookupTable, this->rainbowColors, bs, nullptr, true);
            megamol::core::utility::log::Log::DefaultLog.WriteInfo(
                "Recomputing atom color table WITHOUT color module!");
        }

        this->coloringModeParam0.ResetDirty();
        this->coloringModeParam1.ResetDirty();
        this->cmWeightParam.ResetDirty();
        this->compareParam.ResetDirty();
    }
    if (this->smoothCartoonColoringParam.IsDirty()) {
        this->smoothCartoonColoringMode = this->smoothCartoonColoringParam.Param<param::BoolParam>()->Value();
        this->smoothCartoonColoringParam.ResetDirty();
        // if (this->currentRenderMode == CARTOON) {
        //    cartoonSplineCreated = false;
        //}
    }
}


/*
 * MoleculeCartoonRenderer::RenderCartoonHybrid
 */
void MoleculeCartoonRenderer::RenderCartoonHybrid(const MolecularDataCall* mol, float* atomPos) {
    // prepare hybrid cartoon representation, if necessary
    if (this->prepareCartoonHybrid) {
        unsigned int cntChain, cntS, cntAA, idx, firstSS, countSS, firstAA, countAA;
        // B-Spline
        BSpline bSpline;
        // control points for the first (center) b-spline
        std::vector<vislib::math::Vector<float, 3>> controlPoints;
        // control points for the second (direction) b-spline
        std::vector<vislib::math::Vector<float, 3>> controlPointsDir;
        // temporary vectors
        vislib::math::Vector<float, 3> vecCA, vecC, vecO, vecTmp, vecTmpOld;
        // temporary color
        // const float *color;
        // temporary color vector
        vislib::math::Vector<float, 3> colorVec;

        // coordinates of the first (center) b-spline (result of the spline computation)
        std::vector<std::vector<vislib::math::Vector<float, 3>>> bSplineCoords;
        // coordinates of the second (direction) b-spline (result of the spline computation)
        std::vector<std::vector<vislib::math::Vector<float, 3>>> bSplineCoordsDir;
        // secondary structure type for b-spline
        std::vector<std::vector<MolecularDataCall::SecStructure::ElementType>> bSplineSecStruct;
        // color of secondary structure b-spline
        std::vector<std::vector<vislib::math::Vector<float, 3>>> cartoonColor;

        // set the number of segments to create
        bSpline.setN(this->numberOfSplineSeg);

        // resize result vector for coordinates of first b-spline segments
        bSplineCoords.resize(mol->MoleculeCount());
        // resize result vector for coordinates of second b-spline segments
        bSplineCoordsDir.resize(mol->MoleculeCount());
        // resize vector for secondary structure
        bSplineSecStruct.resize(mol->MoleculeCount());
        // resize color vector
        cartoonColor.resize(mol->MoleculeCount());

        // --- compute the b-splines ---
        // loop over all chains
        MolecularDataCall::Molecule chain;
        MolecularDataCall::AminoAcid* aminoacid;
        for (cntChain = 0; cntChain < mol->MoleculeCount(); ++cntChain) {
            chain = mol->Molecules()[cntChain];
            controlPoints.clear();
            controlPointsDir.clear();
            // check if the first residue is an amino acid
            if (mol->Residues()[chain.FirstResidueIndex()]->Identifier() != MolecularDataCall::Residue::AMINOACID) {
                continue;
            }
            firstSS = chain.FirstSecStructIndex();
            countSS = firstSS + chain.SecStructCount();
            // loop over all secondary structure elements
            for (cntS = firstSS; cntS < countSS; ++cntS) {
                firstAA = mol->SecondaryStructures()[cntS].FirstAminoAcidIndex();
                countAA = firstAA + mol->SecondaryStructures()[cntS].AminoAcidCount();
                // loop over all amino acids in the current sec struct
                for (cntAA = firstAA; cntAA < countAA; ++cntAA) {
                    // add sec struct type
                    bSplineSecStruct[cntChain].push_back(mol->SecondaryStructures()[cntS].Type());
                    // get the index of the C-alpha atom
                    if (mol->Residues()[cntAA]->Identifier() == MolecularDataCall::Residue::AMINOACID)
                        aminoacid = (MolecularDataCall::AminoAcid*)(mol->Residues()[cntAA]);
                    else
                        continue;
                    idx = aminoacid->CAlphaIndex();
                    // get the coordinates of the C-alpha atom
                    vecCA.SetX(atomPos[idx * 3 + 0]);
                    vecCA.SetY(atomPos[idx * 3 + 1]);
                    vecCA.SetZ(atomPos[idx * 3 + 2]);
                    // add the C-alpha coords to the list of control points
                    controlPoints.push_back(vecCA);

                    // add the color of the C-alpha atom to the color vector
                    colorVec.SetX(this->atomColorTable[idx].x);
                    colorVec.SetY(this->atomColorTable[idx].y);
                    colorVec.SetZ(this->atomColorTable[idx].z);
                    cartoonColor[cntChain].push_back(colorVec);

                    // get the index of the C atom
                    idx = aminoacid->CCarbIndex();
                    // get the coordinates of the C-alpha atom
                    vecC.SetX(atomPos[idx * 3 + 0]);
                    vecC.SetY(atomPos[idx * 3 + 1]);
                    vecC.SetZ(atomPos[idx * 3 + 2]);

                    // get the index of the O atom
                    idx = aminoacid->OIndex();
                    // get the coordinates of the C-alpha atom
                    vecO.SetX(atomPos[idx * 3 + 0]);
                    vecO.SetY(atomPos[idx * 3 + 1]);
                    vecO.SetZ(atomPos[idx * 3 + 2]);

                    // compute control point of the second b-spline
                    vecTmp = vecO - vecC;
                    vecTmp.Normalise();
                    // check, if vector should be flipped
                    if (cntS > 0 && vecTmpOld.Dot(vecTmp) < 0.0f)
                        vecTmp = vecTmp * -1.0f;
                    vecTmpOld = vecTmp;
                    // add control point for the second b-spline to the list of control points
                    controlPointsDir.push_back(vecTmp + vecCA);
                }
            }
            // set the control points, compute the first spline and fetch the result
            bSpline.setBackbone(controlPoints);
            if (bSpline.computeSpline())
                bSpline.getResult(bSplineCoords[cntChain]);
            else
                continue; // --> return if spline could not be computed

            // set the control points, compute the second spline and fetch the result
            bSpline.setBackbone(controlPointsDir);
            if (bSpline.computeSpline())
                bSpline.getResult(bSplineCoordsDir[cntChain]);
            else
                continue; // --> return if spline could not be computed
        }

        // --- START store the vertices, colors and parameters ---
        this->totalCountTube = 0;
        this->totalCountArrow = 0;
        this->totalCountHelix = 0;
        for (unsigned int i = 0; i < bSplineCoords.size(); i++) {
            if (bSplineCoords[i].size() == 0)
                continue;
            for (unsigned int j = 2; j < bSplineCoords[i].size() - 1; j++) {
                if (bSplineSecStruct[i][j / this->numberOfSplineSeg] == MolecularDataCall::SecStructure::TYPE_SHEET)
                    this->totalCountArrow++;
                else if (bSplineSecStruct[i][j / this->numberOfSplineSeg] ==
                         MolecularDataCall::SecStructure::TYPE_HELIX)
                    this->totalCountHelix++;
                else
                    this->totalCountTube++;
            }
        }

        if (this->vertTube)
            delete[] this->vertTube;
        if (this->colorsParamsTube)
            delete[] this->colorsParamsTube;
        if (this->vertHelix)
            delete[] this->vertHelix;
        if (this->colorsParamsHelix)
            delete[] this->colorsParamsHelix;
        if (this->vertArrow)
            delete[] this->vertArrow;
        if (this->colorsParamsArrow)
            delete[] this->colorsParamsArrow;
        this->vertTube = new float[this->totalCountTube * 6 * 3];
        this->colorsParamsTube = new float[this->totalCountTube * 6 * 3];
        this->vertArrow = new float[this->totalCountArrow * 6 * 3];
        this->colorsParamsArrow = new float[this->totalCountArrow * 6 * 3];
        this->vertHelix = new float[this->totalCountHelix * 6 * 3];
        this->colorsParamsHelix = new float[this->totalCountHelix * 6 * 3];

        // auxiliary variables
        float start, end, f1, f2, type;
        unsigned int counterTube = 0;
        unsigned int counterArrow = 0;
        unsigned int counterHelix = 0;
        vislib::math::Vector<float, 3> col1, col2;
        // compute the inner b-spline (backbone)
        for (unsigned int i = 0; i < bSplineCoords.size(); i++) {
            if (bSplineCoords[i].size() == 0)
                continue;
            for (unsigned int j = 2; j < bSplineCoords[i].size() - 1; j++) {
                start = end = -1.0f;
                f1 = f2 = 1.0f;
                // set end caps --> if it is the first segment and the last sec struct was different
                if (j / this->numberOfSplineSeg > 0) {
                    if (bSplineSecStruct[i][j / this->numberOfSplineSeg] !=
                            bSplineSecStruct[i][j / this->numberOfSplineSeg - 1] &&
                        j % this->numberOfSplineSeg == 0)
                        end = 1.0f;
                } else if (j == 2)
                    end = 1.0f;
                // set start caps --> if its the last segment and the next sec struct is different
                if (j / this->numberOfSplineSeg < bSplineSecStruct[i].size() - 1) {
                    if (bSplineSecStruct[i][j / this->numberOfSplineSeg] !=
                            bSplineSecStruct[i][j / this->numberOfSplineSeg + 1] &&
                        j % this->numberOfSplineSeg == this->numberOfSplineSeg - 1)
                        start = 1.0f;
                } else if (j == bSplineCoords[i].size() - 2)
                    start = 1.0f;
                // set inParams --> set type and stretch factors of arrow head segments for the sheet
                if (bSplineSecStruct[i][j / this->numberOfSplineSeg] == MolecularDataCall::SecStructure::TYPE_SHEET) {
                    type = 1.0f;
                    if (bSplineSecStruct[i][j / this->numberOfSplineSeg + 1] !=
                        MolecularDataCall::SecStructure::TYPE_SHEET) {
                        if (j % this->numberOfSplineSeg == 0)
                            end = 1.0f;
                        f1 = 1.0f - float(j % this->numberOfSplineSeg) / float(this->numberOfSplineSeg - 1) +
                             1.0f / float(this->numberOfSplineSeg - 1) + 0.2f;
                        f2 = 1.0f - float(j % this->numberOfSplineSeg) / float(this->numberOfSplineSeg - 1) + 0.2f;
                    }
                } else if (bSplineSecStruct[i][j / this->numberOfSplineSeg] ==
                           MolecularDataCall::SecStructure::TYPE_HELIX)
                    type = 2.0f;
                else
                    type = 0.0f;
                // get the colors
                if (this->smoothCartoonColoringMode && j / this->numberOfSplineSeg > 0) {
                    col1 = cartoonColor[i][j / this->numberOfSplineSeg] * float(j % this->numberOfSplineSeg) /
                               float(this->numberOfSplineSeg - 1) +
                           cartoonColor[i][j / this->numberOfSplineSeg - 1] *
                               float((this->numberOfSplineSeg - 1) - j % this->numberOfSplineSeg) /
                               float(this->numberOfSplineSeg - 1);
                    int k = j + 1;
                    if (j % this->numberOfSplineSeg == this->numberOfSplineSeg - 1)
                        k = this->numberOfSplineSeg - 1;
                    col2 = cartoonColor[i][j / this->numberOfSplineSeg] * float(k % this->numberOfSplineSeg) /
                               float(this->numberOfSplineSeg - 1) +
                           cartoonColor[i][j / this->numberOfSplineSeg - 1] *
                               float((this->numberOfSplineSeg - 1) - k % this->numberOfSplineSeg) /
                               float(this->numberOfSplineSeg - 1);
                } else {
                    col1 = cartoonColor[i][j / this->numberOfSplineSeg];
                    col2 = cartoonColor[i][j / this->numberOfSplineSeg];
                }

                // store information in the apropriate arrays
                if (bSplineSecStruct[i][j / this->numberOfSplineSeg] == MolecularDataCall::SecStructure::TYPE_SHEET) {
                    this->colorsParamsArrow[counterArrow * 6 * 3 + 0] = col1.GetX();
                    this->colorsParamsArrow[counterArrow * 6 * 3 + 1] = col1.GetY();
                    this->colorsParamsArrow[counterArrow * 6 * 3 + 2] = col1.GetZ();
                    this->colorsParamsArrow[counterArrow * 6 * 3 + 3] = this->radiusCartoon;
                    this->colorsParamsArrow[counterArrow * 6 * 3 + 4] = f1;
                    this->colorsParamsArrow[counterArrow * 6 * 3 + 5] = f2;
                    this->colorsParamsArrow[counterArrow * 6 * 3 + 6] = type;
                    this->colorsParamsArrow[counterArrow * 6 * 3 + 7] = start;
                    this->colorsParamsArrow[counterArrow * 6 * 3 + 8] = end;
                    this->colorsParamsArrow[counterArrow * 6 * 3 + 9] = col2.GetX();
                    this->colorsParamsArrow[counterArrow * 6 * 3 + 10] = col2.GetY();
                    this->colorsParamsArrow[counterArrow * 6 * 3 + 11] = col2.GetZ();
                    this->colorsParamsArrow[counterArrow * 6 * 3 + 12] = 0.0f;
                    this->colorsParamsArrow[counterArrow * 6 * 3 + 13] = 0.0f;
                    this->colorsParamsArrow[counterArrow * 6 * 3 + 14] = 0.0f;
                    this->colorsParamsArrow[counterArrow * 6 * 3 + 15] = 0.0f;
                    this->colorsParamsArrow[counterArrow * 6 * 3 + 16] = 0.0f;
                    this->colorsParamsArrow[counterArrow * 6 * 3 + 17] = 0.0f;
                    this->vertArrow[counterArrow * 6 * 3 + 0] = bSplineCoords[i][j - 2].GetX();
                    this->vertArrow[counterArrow * 6 * 3 + 1] = bSplineCoords[i][j - 2].GetY();
                    this->vertArrow[counterArrow * 6 * 3 + 2] = bSplineCoords[i][j - 2].GetZ();
                    this->vertArrow[counterArrow * 6 * 3 + 3] = bSplineCoordsDir[i][j - 1].GetX();
                    this->vertArrow[counterArrow * 6 * 3 + 4] = bSplineCoordsDir[i][j - 1].GetY();
                    this->vertArrow[counterArrow * 6 * 3 + 5] = bSplineCoordsDir[i][j - 1].GetZ();
                    this->vertArrow[counterArrow * 6 * 3 + 6] = bSplineCoords[i][j - 1].GetX();
                    this->vertArrow[counterArrow * 6 * 3 + 7] = bSplineCoords[i][j - 1].GetY();
                    this->vertArrow[counterArrow * 6 * 3 + 8] = bSplineCoords[i][j - 1].GetZ();
                    this->vertArrow[counterArrow * 6 * 3 + 9] = bSplineCoords[i][j].GetX();
                    this->vertArrow[counterArrow * 6 * 3 + 10] = bSplineCoords[i][j].GetY();
                    this->vertArrow[counterArrow * 6 * 3 + 11] = bSplineCoords[i][j].GetZ();
                    this->vertArrow[counterArrow * 6 * 3 + 12] = bSplineCoordsDir[i][j].GetX();
                    this->vertArrow[counterArrow * 6 * 3 + 13] = bSplineCoordsDir[i][j].GetY();
                    this->vertArrow[counterArrow * 6 * 3 + 14] = bSplineCoordsDir[i][j].GetZ();
                    this->vertArrow[counterArrow * 6 * 3 + 15] = bSplineCoords[i][j + 1].GetX();
                    this->vertArrow[counterArrow * 6 * 3 + 16] = bSplineCoords[i][j + 1].GetY();
                    this->vertArrow[counterArrow * 6 * 3 + 17] = bSplineCoords[i][j + 1].GetZ();
                    counterArrow++;
                } else if (bSplineSecStruct[i][j / this->numberOfSplineSeg] ==
                           MolecularDataCall::SecStructure::TYPE_HELIX) {
                    this->colorsParamsHelix[counterHelix * 6 * 3 + 0] = col1.GetX();
                    this->colorsParamsHelix[counterHelix * 6 * 3 + 1] = col1.GetY();
                    this->colorsParamsHelix[counterHelix * 6 * 3 + 2] = col1.GetZ();
                    this->colorsParamsHelix[counterHelix * 6 * 3 + 3] = this->radiusCartoon;
                    this->colorsParamsHelix[counterHelix * 6 * 3 + 4] = f1;
                    this->colorsParamsHelix[counterHelix * 6 * 3 + 5] = f2;
                    this->colorsParamsHelix[counterHelix * 6 * 3 + 6] = type;
                    this->colorsParamsHelix[counterHelix * 6 * 3 + 7] = start;
                    this->colorsParamsHelix[counterHelix * 6 * 3 + 8] = end;
                    this->colorsParamsHelix[counterHelix * 6 * 3 + 9] = col2.GetX();
                    this->colorsParamsHelix[counterHelix * 6 * 3 + 10] = col2.GetY();
                    this->colorsParamsHelix[counterHelix * 6 * 3 + 11] = col2.GetZ();
                    this->colorsParamsHelix[counterHelix * 6 * 3 + 12] = 0.0f;
                    this->colorsParamsHelix[counterHelix * 6 * 3 + 13] = 0.0f;
                    this->colorsParamsHelix[counterHelix * 6 * 3 + 14] = 0.0f;
                    this->colorsParamsHelix[counterHelix * 6 * 3 + 15] = 0.0f;
                    this->colorsParamsHelix[counterHelix * 6 * 3 + 16] = 0.0f;
                    this->colorsParamsHelix[counterHelix * 6 * 3 + 17] = 0.0f;
                    this->vertHelix[counterHelix * 6 * 3 + 0] = bSplineCoords[i][j - 2].GetX();
                    this->vertHelix[counterHelix * 6 * 3 + 1] = bSplineCoords[i][j - 2].GetY();
                    this->vertHelix[counterHelix * 6 * 3 + 2] = bSplineCoords[i][j - 2].GetZ();
                    this->vertHelix[counterHelix * 6 * 3 + 3] = bSplineCoordsDir[i][j - 1].GetX();
                    this->vertHelix[counterHelix * 6 * 3 + 4] = bSplineCoordsDir[i][j - 1].GetY();
                    this->vertHelix[counterHelix * 6 * 3 + 5] = bSplineCoordsDir[i][j - 1].GetZ();
                    this->vertHelix[counterHelix * 6 * 3 + 6] = bSplineCoords[i][j - 1].GetX();
                    this->vertHelix[counterHelix * 6 * 3 + 7] = bSplineCoords[i][j - 1].GetY();
                    this->vertHelix[counterHelix * 6 * 3 + 8] = bSplineCoords[i][j - 1].GetZ();
                    this->vertHelix[counterHelix * 6 * 3 + 9] = bSplineCoords[i][j].GetX();
                    this->vertHelix[counterHelix * 6 * 3 + 10] = bSplineCoords[i][j].GetY();
                    this->vertHelix[counterHelix * 6 * 3 + 11] = bSplineCoords[i][j].GetZ();
                    this->vertHelix[counterHelix * 6 * 3 + 12] = bSplineCoordsDir[i][j].GetX();
                    this->vertHelix[counterHelix * 6 * 3 + 13] = bSplineCoordsDir[i][j].GetY();
                    this->vertHelix[counterHelix * 6 * 3 + 14] = bSplineCoordsDir[i][j].GetZ();
                    this->vertHelix[counterHelix * 6 * 3 + 15] = bSplineCoords[i][j + 1].GetX();
                    this->vertHelix[counterHelix * 6 * 3 + 16] = bSplineCoords[i][j + 1].GetY();
                    this->vertHelix[counterHelix * 6 * 3 + 17] = bSplineCoords[i][j + 1].GetZ();
                    counterHelix++;
                } else {
                    this->colorsParamsTube[counterTube * 6 * 3 + 0] = col1.GetX();
                    this->colorsParamsTube[counterTube * 6 * 3 + 1] = col1.GetY();
                    this->colorsParamsTube[counterTube * 6 * 3 + 2] = col1.GetZ();
                    this->colorsParamsTube[counterTube * 6 * 3 + 3] = this->radiusCartoon;
                    this->colorsParamsTube[counterTube * 6 * 3 + 4] = f1;
                    this->colorsParamsTube[counterTube * 6 * 3 + 5] = f2;
                    this->colorsParamsTube[counterTube * 6 * 3 + 6] = type;
                    this->colorsParamsTube[counterTube * 6 * 3 + 7] = start;
                    this->colorsParamsTube[counterTube * 6 * 3 + 8] = end;
                    this->colorsParamsTube[counterTube * 6 * 3 + 9] = col2.GetX();
                    this->colorsParamsTube[counterTube * 6 * 3 + 10] = col2.GetY();
                    this->colorsParamsTube[counterTube * 6 * 3 + 11] = col2.GetZ();
                    this->colorsParamsTube[counterTube * 6 * 3 + 12] = 0.0f;
                    this->colorsParamsTube[counterTube * 6 * 3 + 13] = 0.0f;
                    this->colorsParamsTube[counterTube * 6 * 3 + 14] = 0.0f;
                    this->colorsParamsTube[counterTube * 6 * 3 + 15] = 0.0f;
                    this->colorsParamsTube[counterTube * 6 * 3 + 16] = 0.0f;
                    this->colorsParamsTube[counterTube * 6 * 3 + 17] = 0.0f;
                    this->vertTube[counterTube * 6 * 3 + 0] = bSplineCoords[i][j - 2].GetX();
                    this->vertTube[counterTube * 6 * 3 + 1] = bSplineCoords[i][j - 2].GetY();
                    this->vertTube[counterTube * 6 * 3 + 2] = bSplineCoords[i][j - 2].GetZ();
                    this->vertTube[counterTube * 6 * 3 + 3] = bSplineCoordsDir[i][j - 1].GetX();
                    this->vertTube[counterTube * 6 * 3 + 4] = bSplineCoordsDir[i][j - 1].GetY();
                    this->vertTube[counterTube * 6 * 3 + 5] = bSplineCoordsDir[i][j - 1].GetZ();
                    this->vertTube[counterTube * 6 * 3 + 6] = bSplineCoords[i][j - 1].GetX();
                    this->vertTube[counterTube * 6 * 3 + 7] = bSplineCoords[i][j - 1].GetY();
                    this->vertTube[counterTube * 6 * 3 + 8] = bSplineCoords[i][j - 1].GetZ();
                    this->vertTube[counterTube * 6 * 3 + 9] = bSplineCoords[i][j].GetX();
                    this->vertTube[counterTube * 6 * 3 + 10] = bSplineCoords[i][j].GetY();
                    this->vertTube[counterTube * 6 * 3 + 11] = bSplineCoords[i][j].GetZ();
                    this->vertTube[counterTube * 6 * 3 + 12] = bSplineCoordsDir[i][j].GetX();
                    this->vertTube[counterTube * 6 * 3 + 13] = bSplineCoordsDir[i][j].GetY();
                    this->vertTube[counterTube * 6 * 3 + 14] = bSplineCoordsDir[i][j].GetZ();
                    this->vertTube[counterTube * 6 * 3 + 15] = bSplineCoords[i][j + 1].GetX();
                    this->vertTube[counterTube * 6 * 3 + 16] = bSplineCoords[i][j + 1].GetY();
                    this->vertTube[counterTube * 6 * 3 + 17] = bSplineCoords[i][j + 1].GetZ();
                    counterTube++;
                }
            }
        }

        // --- END store vertex/color/inparams ---

        // set cartoon as created
        this->prepareCartoonHybrid = false;
    }

    float spec[4] = {1.0f, 1.0f, 1.0f, 1.0f};
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, spec);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 50.0f);
    glDisable(GL_COLOR_MATERIAL);

    // Get current window size
    float curVP[4];
    glGetFloatv(GL_VIEWPORT, curVP);

    // enable tube shader
    if (this->currentRenderMode == CartoonRenderMode::CARTOON)
        this->tubeShader->use();
    else
        this->tubeSimpleShader->use();

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
    glVertexPointer(3, GL_FLOAT, 0, this->vertTube);
    glColorPointer(3, GL_FLOAT, 0, this->colorsParamsTube);
    glDrawArrays(GL_TRIANGLES_ADJACENCY_EXT, 0, this->totalCountTube * 6);
    // disable tube shader
    glUseProgram(0);

    // enable arrow shader
    if (this->currentRenderMode == CartoonRenderMode::CARTOON)
        this->arrowShader->use();
    else
        this->arrowSimpleShader->use();
    glVertexPointer(3, GL_FLOAT, 0, this->vertArrow);
    glColorPointer(3, GL_FLOAT, 0, this->colorsParamsArrow);
    glDrawArrays(GL_TRIANGLES_ADJACENCY_EXT, 0, this->totalCountArrow * 6);
    // disable arrow shader
    glUseProgram(0);

    // enable helix shader
    if (this->currentRenderMode == CartoonRenderMode::CARTOON)
        this->helixShader->use();
    else
        this->helixSimpleShader->use();
    glVertexPointer(3, GL_FLOAT, 0, this->vertHelix);
    glColorPointer(3, GL_FLOAT, 0, this->colorsParamsHelix);
    glDrawArrays(GL_TRIANGLES_ADJACENCY_EXT, 0, this->totalCountHelix * 6);
    // disable helix shader
    glUseProgram(0);

    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
    glDisable(GL_COLOR_MATERIAL);
}


/*
 * MoleculeCartoonRenderer::RenderCartoonCPU
 */
void MoleculeCartoonRenderer::RenderCartoonCPU(const MolecularDataCall* mol, float* atomPos) {

    // prepare CPU cartoon representation, if necessary
    if (this->prepareCartoonCPU) {
        unsigned int cntChain, cntS, cntAA, idx, firstSS, countSS, firstAA, countAA;
        // B-Spline
        BSpline bSpline;
        // control points for the first (center) b-spline
        std::vector<vislib::math::Vector<float, 3>> controlPoints;
        // control points for the second (direction) b-spline
        std::vector<vislib::math::Vector<float, 3>> controlPointsDir;
        // temporary vectors
        vislib::math::Vector<float, 3> vecCA, vecC, vecO, vecTmp, vecTmpOld(0.0f, 0.0f, 0.0f);
        // temporary color vector
        vislib::math::Vector<float, 3> colorVec;

        // secondary structure type for b-spline
        std::vector<std::vector<MolecularDataCall::SecStructure::ElementType>> bSplineSecStruct;


        // set the number of segments to create
        bSpline.setN(this->numberOfSplineSeg);

        // resize result vector for coordinates of first b-spline segments
        bSplineCoordsCPU.resize(mol->MoleculeCount());
        // resize result vector for coordinates of first b-spline segments
        bSplineCoordsDirCPU.resize(mol->MoleculeCount());
        // resize vector for secondary structure
        bSplineSecStruct.resize(mol->MoleculeCount());
        // resize color vector
        this->cartoonColorCPU.resize(mol->MoleculeCount());

        // --- compute the b-splines ---
        // loop over all chains
        MolecularDataCall::Molecule chain;
        MolecularDataCall::AminoAcid* aminoacid;
        for (cntChain = 0; cntChain < mol->MoleculeCount(); ++cntChain) {
            chain = mol->Molecules()[cntChain];
            controlPoints.clear();
            controlPointsDir.clear();
            this->cartoonColorCPU[cntChain].clear();
            // check if the first residue is an amino acid
            if (mol->Residues()[chain.FirstResidueIndex()]->Identifier() != MolecularDataCall::Residue::AMINOACID) {
                continue;
            }
            firstSS = chain.FirstSecStructIndex();
            countSS = firstSS + chain.SecStructCount();
            // loop over all secondary structure elements
            for (cntS = firstSS; cntS < countSS; ++cntS) {
                firstAA = mol->SecondaryStructures()[cntS].FirstAminoAcidIndex();
                countAA = firstAA + mol->SecondaryStructures()[cntS].AminoAcidCount();
                // loop over all amino acids in the current sec struct
                for (cntAA = firstAA; cntAA < countAA; ++cntAA) {
                    // add sec struct type
                    bSplineSecStruct[cntChain].push_back(mol->SecondaryStructures()[cntS].Type());

                    // get the index of the C-alpha atom
                    if (mol->Residues()[cntAA]->Identifier() == MolecularDataCall::Residue::AMINOACID)
                        aminoacid = (MolecularDataCall::AminoAcid*)(mol->Residues()[cntAA]);
                    else
                        continue;

                    // get the index of the C-alpha atom
                    idx = aminoacid->CAlphaIndex();
                    // get the coordinates of the C-alpha atom
                    vecCA.SetX(atomPos[idx * 3 + 0]);
                    vecCA.SetY(atomPos[idx * 3 + 1]);
                    vecCA.SetZ(atomPos[idx * 3 + 2]);
                    // add the C-alpha coords to the list of control points
                    controlPoints.push_back(vecCA);

                    // add the color of the C-alpha atom to the color vector
                    colorVec.SetX(this->atomColorTable[idx].x);
                    colorVec.SetY(this->atomColorTable[idx].y);
                    colorVec.SetZ(this->atomColorTable[idx].z);
                    cartoonColorCPU[cntChain].push_back(colorVec);

                    // get the index of the C atom
                    idx = aminoacid->CCarbIndex();
                    // get the coordinates of the C-alpha atom
                    vecC.SetX(atomPos[idx * 3 + 0]);
                    vecC.SetY(atomPos[idx * 3 + 1]);
                    vecC.SetZ(atomPos[idx * 3 + 2]);

                    // get the index of the O atom
                    idx = aminoacid->OIndex();
                    // get the coordinates of the C-alpha atom
                    vecO.SetX(atomPos[idx * 3 + 0]);
                    vecO.SetY(atomPos[idx * 3 + 1]);
                    vecO.SetZ(atomPos[idx * 3 + 2]);

                    // compute control point of the second b-spline
                    vecTmp = vecO - vecC;
                    vecTmp.Normalise();
                    // check, if vector should be flipped
                    if (cntS > 0 && vecTmpOld.Dot(vecTmp) < 0.0f)
                        vecTmp = vecTmp * -1.0f;
                    vecTmpOld = vecTmp;
                    // add control point for the second b-spline to the list of control points
                    controlPointsDir.push_back(vecTmp + vecCA);
                }
            }
            // set the control points, compute the first spline and fetch the result
            bSpline.setBackbone(controlPoints);
            if (bSpline.computeSpline())
                bSpline.getResult(bSplineCoordsCPU[cntChain]);
            else
                continue; // --> return if spline could not be computed

            // set the control points, compute the second spline and fetch the result
            bSpline.setBackbone(controlPointsDir);
            if (bSpline.computeSpline())
                bSpline.getResult(bSplineCoordsDirCPU[cntChain]);
            else
                continue; // --> return if spline could not be computed
        }

        // --- START store the vertices, colors and parameters ---
        this->totalCountTube = 0;
        this->totalCountArrow = 0;
        this->totalCountHelix = 0;
        for (unsigned int i = 0; i < bSplineCoordsCPU.size(); i++) {
            if (bSplineCoordsCPU[i].size() == 0)
                continue;
            for (unsigned int j = 2; j < bSplineCoordsCPU[i].size() - 1; j++) {
                if (bSplineSecStruct[i][j / this->numberOfSplineSeg] == MolecularDataCall::SecStructure::TYPE_SHEET)
                    this->totalCountArrow++;
                else if (bSplineSecStruct[i][j / this->numberOfSplineSeg] ==
                         MolecularDataCall::SecStructure::TYPE_HELIX)
                    this->totalCountHelix++;
                else
                    this->totalCountTube++;
            }
        }

        if (this->vertTube)
            delete[] this->vertTube;
        if (this->colorsParamsTube)
            delete[] this->colorsParamsTube;
        if (this->vertHelix)
            delete[] this->vertHelix;
        if (this->colorsParamsHelix)
            delete[] this->colorsParamsHelix;
        if (this->vertArrow)
            delete[] this->vertArrow;
        if (this->colorsParamsArrow)
            delete[] this->colorsParamsArrow;
        if (this->normalTube)
            delete[] this->normalTube;
        if (this->normalHelix)
            delete[] this->normalHelix;
        if (this->normalArrow)
            delete[] this->normalArrow;
        this->vertTube = new float[this->totalCountTube * 3 * 4 * this->numberOfTubeSeg];
        this->colorsParamsTube = new float[this->totalCountTube * 3 * 4 * this->numberOfTubeSeg];
        // 4 3D-Punkte pro Quad, 4 Quads pro Segment, d.h. 16 3D-Punkte pro Segment
        this->vertArrow = new float[this->totalCountArrow * 3 * 16];
        this->colorsParamsArrow = new float[this->totalCountArrow * 3 * 16];
        this->vertHelix = new float[this->totalCountHelix * 3 * 16];
        this->colorsParamsHelix = new float[this->totalCountHelix * 3 * 16];
        this->normalTube = new float[this->totalCountTube * 3 * 4 * this->numberOfTubeSeg];
        this->normalHelix = new float[this->totalCountHelix * 3 * 16];
        this->normalArrow = new float[this->totalCountArrow * 3 * 16];

        // auxiliary variables
        float start, end, f1, f2, type;
        unsigned int counterTube = 0;
        unsigned int counterArrow = 0;
        unsigned int counterHelix = 0;
        vislib::math::Vector<float, 3> col1, col2;

        vislib::math::Vector<float, 3> v0;
        vislib::math::Vector<float, 3> v1;
        vislib::math::Vector<float, 3> v2;
        vislib::math::Vector<float, 3> v3;
        vislib::math::Vector<float, 3> v4;
        vislib::math::Vector<float, 3> v5;
        vislib::math::Vector<float, 3> dir20;
        vislib::math::Vector<float, 3> dir12;
        vislib::math::Vector<float, 3> dir32;
        vislib::math::Vector<float, 3> dir43;
        vislib::math::Vector<float, 3> dir53;
        vislib::math::Vector<float, 3> res1;
        vislib::math::Vector<float, 3> res2;
        float scale;
        float stretch1;
        float stretch2;
        vislib::math::Vector<float, 3> ortho1;
        vislib::math::Vector<float, 3> ortho2;
        vislib::math::Vector<float, 3> dir1;
        vislib::math::Vector<float, 3> dir2;
        vislib::math::Vector<float, 3> norm1;
        vislib::math::Vector<float, 3> norm2;
        vislib::math::Quaternion<float> q1;
        vislib::math::Quaternion<float> q2;
        // angle for the rotation
        float alpha;

        // compute the geometry
        for (unsigned int i = 0; i < bSplineCoordsCPU.size(); i++) {
            if (bSplineCoordsCPU[i].size() == 0)
                continue;
            for (unsigned int j = 2; j < bSplineCoordsCPU[i].size() - 1; j++) {
                start = end = -1.0f;
                f1 = f2 = 1.0f;
                // set end caps --> if it is the first segment and the last sec struct was different
                if (j / this->numberOfSplineSeg > 0) {
                    if (bSplineSecStruct[i][j / this->numberOfSplineSeg] !=
                            bSplineSecStruct[i][j / this->numberOfSplineSeg - 1] &&
                        j % this->numberOfSplineSeg == 0)
                        end = 1.0f;
                } else if (j == 2)
                    end = 1.0f;
                // set start caps --> if its the last segment and the next sec struct is different
                if (j / this->numberOfSplineSeg < bSplineSecStruct[i].size() - 1) {
                    if (bSplineSecStruct[i][j / this->numberOfSplineSeg] !=
                            bSplineSecStruct[i][j / this->numberOfSplineSeg + 1] &&
                        j % this->numberOfSplineSeg == this->numberOfSplineSeg - 1)
                        start = 1.0f;
                } else if (j == bSplineCoordsCPU[i].size() - 2)
                    start = 1.0f;
                // set inParams --> set type and stretch factors of arrow head segments for the sheet
                if (bSplineSecStruct[i][j / this->numberOfSplineSeg] == MolecularDataCall::SecStructure::TYPE_SHEET) {
                    type = 1.0f;
                    if (bSplineSecStruct[i][j / this->numberOfSplineSeg + 1] !=
                        MolecularDataCall::SecStructure::TYPE_SHEET) {
                        if (j % this->numberOfSplineSeg == 0)
                            end = 1.0f;
                        f1 = 1.0f - float(j % this->numberOfSplineSeg) / float(this->numberOfSplineSeg - 1) +
                             1.0f / float(this->numberOfSplineSeg - 1) + 0.2f;
                        f2 = 1.0f - float(j % this->numberOfSplineSeg) / float(this->numberOfSplineSeg - 1) + 0.2f;
                    }
                } else if (bSplineSecStruct[i][j / this->numberOfSplineSeg] ==
                           MolecularDataCall::SecStructure::TYPE_HELIX)
                    type = 2.0f;
                else
                    type = 0.0f;
                // get the colors
                if (this->smoothCartoonColoringMode && j / this->numberOfSplineSeg > 0) {
                    col1 = cartoonColorCPU[i][j / this->numberOfSplineSeg] * float(j % this->numberOfSplineSeg) /
                               float(this->numberOfSplineSeg - 1) +
                           cartoonColorCPU[i][j / this->numberOfSplineSeg - 1] *
                               float((this->numberOfSplineSeg - 1) - j % this->numberOfSplineSeg) /
                               float(this->numberOfSplineSeg - 1);
                    int k = j + 1;
                    if (j % this->numberOfSplineSeg == this->numberOfSplineSeg - 1)
                        k = this->numberOfSplineSeg - 1;
                    col2 = cartoonColorCPU[i][j / this->numberOfSplineSeg] * float(k % this->numberOfSplineSeg) /
                               float(this->numberOfSplineSeg - 1) +
                           cartoonColorCPU[i][j / this->numberOfSplineSeg - 1] *
                               float((this->numberOfSplineSeg - 1) - k % this->numberOfSplineSeg) /
                               float(this->numberOfSplineSeg - 1);
                } else {
                    col1 = cartoonColorCPU[i][j / this->numberOfSplineSeg];
                    col2 = cartoonColorCPU[i][j / this->numberOfSplineSeg];
                }

                // -------------------------------------
                // --- START computation from shader ---
                // -------------------------------------

                // get all vertex positions
                v0 = bSplineCoordsCPU[i][j - 2];
                v1 = bSplineCoordsDirCPU[i][j - 1];
                v2 = bSplineCoordsCPU[i][j - 1];
                v3 = bSplineCoordsCPU[i][j];
                v4 = bSplineCoordsDirCPU[i][j];
                v5 = bSplineCoordsCPU[i][j + 1];
                // compute all needed directions
                dir20 = v2 - v0;
                dir12 = v1 - v2;
                dir32 = v3 - v2;
                dir43 = v4 - v3;
                dir53 = v5 - v3;
                // scale factor for the width of the tube

                //      this->colorsParamsTube[counterTube*6*3+7] = start;
                //      this->colorsParamsTube[counterTube*6*3+8] = end;
                scale = this->radiusCartoon;
                stretch1 = f1;
                stretch2 = f2;
                ortho1 = (dir20 + dir32);
                ortho1.Normalise();
                ortho2 = (dir32 + dir53);
                ortho2.Normalise();
                dir1 = (dir12.Cross(ortho1));
                dir1.Normalise();
                dir2 = (dir43.Cross(ortho2));
                dir2.Normalise();

                dir1 = (dir1.Cross(ortho1));
                dir1.Normalise();
                dir1 = dir1 * stretch1;
                dir2 = (dir2.Cross(ortho2));
                dir2.Normalise();
                dir2 = dir2 * stretch2;
                norm1 = (dir1.Cross(ortho1));
                norm1.Normalise();
                norm2 = (dir2.Cross(ortho2));
                norm2.Normalise();

                // -----------------------------------
                // --- END computation from shader ---
                // -----------------------------------

                // store information in the apropriate arrays
                if (bSplineSecStruct[i][j / this->numberOfSplineSeg] == MolecularDataCall::SecStructure::TYPE_SHEET) {
                    this->vertArrow[counterArrow * 3 * 16 + 0] = (v2 - dir1 + norm1 * scale).GetX();
                    this->vertArrow[counterArrow * 3 * 16 + 1] = (v2 - dir1 + norm1 * scale).GetY();
                    this->vertArrow[counterArrow * 3 * 16 + 2] = (v2 - dir1 + norm1 * scale).GetZ();
                    this->vertArrow[counterArrow * 3 * 16 + 3] = (v2 + dir1 + norm1 * scale).GetX();
                    this->vertArrow[counterArrow * 3 * 16 + 4] = (v2 + dir1 + norm1 * scale).GetY();
                    this->vertArrow[counterArrow * 3 * 16 + 5] = (v2 + dir1 + norm1 * scale).GetZ();
                    this->vertArrow[counterArrow * 3 * 16 + 6] = (v3 + dir2 + norm2 * scale).GetX();
                    this->vertArrow[counterArrow * 3 * 16 + 7] = (v3 + dir2 + norm2 * scale).GetY();
                    this->vertArrow[counterArrow * 3 * 16 + 8] = (v3 + dir2 + norm2 * scale).GetZ();
                    this->vertArrow[counterArrow * 3 * 16 + 9] = (v3 - dir2 + norm2 * scale).GetX();
                    this->vertArrow[counterArrow * 3 * 16 + 10] = (v3 - dir2 + norm2 * scale).GetY();
                    this->vertArrow[counterArrow * 3 * 16 + 11] = (v3 - dir2 + norm2 * scale).GetZ();
                    this->vertArrow[counterArrow * 3 * 16 + 12] = (v2 - dir1 - norm1 * scale).GetX();
                    this->vertArrow[counterArrow * 3 * 16 + 13] = (v2 - dir1 - norm1 * scale).GetY();
                    this->vertArrow[counterArrow * 3 * 16 + 14] = (v2 - dir1 - norm1 * scale).GetZ();
                    this->vertArrow[counterArrow * 3 * 16 + 15] = (v2 + dir1 - norm1 * scale).GetX();
                    this->vertArrow[counterArrow * 3 * 16 + 16] = (v2 + dir1 - norm1 * scale).GetY();
                    this->vertArrow[counterArrow * 3 * 16 + 17] = (v2 + dir1 - norm1 * scale).GetZ();
                    this->vertArrow[counterArrow * 3 * 16 + 18] = (v3 + dir2 - norm2 * scale).GetX();
                    this->vertArrow[counterArrow * 3 * 16 + 19] = (v3 + dir2 - norm2 * scale).GetY();
                    this->vertArrow[counterArrow * 3 * 16 + 20] = (v3 + dir2 - norm2 * scale).GetZ();
                    this->vertArrow[counterArrow * 3 * 16 + 21] = (v3 - dir2 - norm2 * scale).GetX();
                    this->vertArrow[counterArrow * 3 * 16 + 22] = (v3 - dir2 - norm2 * scale).GetY();
                    this->vertArrow[counterArrow * 3 * 16 + 23] = (v3 - dir2 - norm2 * scale).GetZ();
                    this->vertArrow[counterArrow * 3 * 16 + 24] = (v2 + dir1 + norm1 * scale).GetX();
                    this->vertArrow[counterArrow * 3 * 16 + 25] = (v2 + dir1 + norm1 * scale).GetY();
                    this->vertArrow[counterArrow * 3 * 16 + 26] = (v2 + dir1 + norm1 * scale).GetZ();
                    this->vertArrow[counterArrow * 3 * 16 + 27] = (v2 + dir1 - norm1 * scale).GetX();
                    this->vertArrow[counterArrow * 3 * 16 + 28] = (v2 + dir1 - norm1 * scale).GetY();
                    this->vertArrow[counterArrow * 3 * 16 + 29] = (v2 + dir1 - norm1 * scale).GetZ();
                    this->vertArrow[counterArrow * 3 * 16 + 30] = (v3 + dir2 - norm2 * scale).GetX();
                    this->vertArrow[counterArrow * 3 * 16 + 31] = (v3 + dir2 - norm2 * scale).GetY();
                    this->vertArrow[counterArrow * 3 * 16 + 32] = (v3 + dir2 - norm2 * scale).GetZ();
                    this->vertArrow[counterArrow * 3 * 16 + 33] = (v3 + dir2 + norm2 * scale).GetX();
                    this->vertArrow[counterArrow * 3 * 16 + 34] = (v3 + dir2 + norm2 * scale).GetY();
                    this->vertArrow[counterArrow * 3 * 16 + 35] = (v3 + dir2 + norm2 * scale).GetZ();
                    this->vertArrow[counterArrow * 3 * 16 + 36] = (v2 - dir1 + norm1 * scale).GetX();
                    this->vertArrow[counterArrow * 3 * 16 + 37] = (v2 - dir1 + norm1 * scale).GetY();
                    this->vertArrow[counterArrow * 3 * 16 + 38] = (v2 - dir1 + norm1 * scale).GetZ();
                    this->vertArrow[counterArrow * 3 * 16 + 39] = (v2 - dir1 - norm1 * scale).GetX();
                    this->vertArrow[counterArrow * 3 * 16 + 40] = (v2 - dir1 - norm1 * scale).GetY();
                    this->vertArrow[counterArrow * 3 * 16 + 41] = (v2 - dir1 - norm1 * scale).GetZ();
                    this->vertArrow[counterArrow * 3 * 16 + 42] = (v3 - dir2 - norm2 * scale).GetX();
                    this->vertArrow[counterArrow * 3 * 16 + 43] = (v3 - dir2 - norm2 * scale).GetY();
                    this->vertArrow[counterArrow * 3 * 16 + 44] = (v3 - dir2 - norm2 * scale).GetZ();
                    this->vertArrow[counterArrow * 3 * 16 + 45] = (v3 - dir2 + norm2 * scale).GetX();
                    this->vertArrow[counterArrow * 3 * 16 + 46] = (v3 - dir2 + norm2 * scale).GetY();
                    this->vertArrow[counterArrow * 3 * 16 + 47] = (v3 - dir2 + norm2 * scale).GetZ();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 0] = col1.GetX();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 1] = col1.GetY();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 2] = col1.GetZ();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 3] = col1.GetX();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 4] = col1.GetY();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 5] = col1.GetZ();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 6] = col2.GetX();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 7] = col2.GetY();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 8] = col2.GetZ();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 9] = col2.GetX();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 10] = col2.GetY();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 11] = col2.GetZ();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 12] = col1.GetX();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 13] = col1.GetY();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 14] = col1.GetZ();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 15] = col1.GetX();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 16] = col1.GetY();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 17] = col1.GetZ();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 18] = col2.GetX();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 19] = col2.GetY();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 20] = col2.GetZ();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 21] = col2.GetX();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 22] = col2.GetY();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 23] = col2.GetZ();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 24] = col1.GetX();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 25] = col1.GetY();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 26] = col1.GetZ();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 27] = col1.GetX();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 28] = col1.GetY();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 29] = col1.GetZ();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 30] = col2.GetX();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 31] = col2.GetY();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 32] = col2.GetZ();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 33] = col2.GetX();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 34] = col2.GetY();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 35] = col2.GetZ();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 36] = col1.GetX();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 37] = col1.GetY();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 38] = col1.GetZ();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 39] = col1.GetX();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 40] = col1.GetY();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 41] = col1.GetZ();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 42] = col2.GetX();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 43] = col2.GetY();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 44] = col2.GetZ();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 45] = col2.GetX();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 46] = col2.GetY();
                    this->colorsParamsArrow[counterArrow * 3 * 16 + 47] = col2.GetZ();
                    norm1.Normalise();
                    norm2.Normalise();
                    dir1.Normalise();
                    dir2.Normalise();
                    this->normalArrow[counterArrow * 3 * 16 + 0] = (norm1).GetX();
                    this->normalArrow[counterArrow * 3 * 16 + 1] = (norm1).GetY();
                    this->normalArrow[counterArrow * 3 * 16 + 2] = (norm1).GetZ();
                    this->normalArrow[counterArrow * 3 * 16 + 3] = (norm1).GetX();
                    this->normalArrow[counterArrow * 3 * 16 + 4] = (norm1).GetY();
                    this->normalArrow[counterArrow * 3 * 16 + 5] = (norm1).GetZ();
                    this->normalArrow[counterArrow * 3 * 16 + 6] = (norm2).GetX();
                    this->normalArrow[counterArrow * 3 * 16 + 7] = (norm2).GetY();
                    this->normalArrow[counterArrow * 3 * 16 + 8] = (norm2).GetZ();
                    this->normalArrow[counterArrow * 3 * 16 + 9] = (norm2).GetX();
                    this->normalArrow[counterArrow * 3 * 16 + 10] = (norm2).GetY();
                    this->normalArrow[counterArrow * 3 * 16 + 11] = (norm2).GetZ();
                    this->normalArrow[counterArrow * 3 * 16 + 12] = (-norm1).GetX();
                    this->normalArrow[counterArrow * 3 * 16 + 13] = (-norm1).GetY();
                    this->normalArrow[counterArrow * 3 * 16 + 14] = (-norm1).GetZ();
                    this->normalArrow[counterArrow * 3 * 16 + 15] = (-norm1).GetX();
                    this->normalArrow[counterArrow * 3 * 16 + 16] = (-norm1).GetY();
                    this->normalArrow[counterArrow * 3 * 16 + 17] = (-norm1).GetZ();
                    this->normalArrow[counterArrow * 3 * 16 + 18] = (-norm2).GetX();
                    this->normalArrow[counterArrow * 3 * 16 + 19] = (-norm2).GetY();
                    this->normalArrow[counterArrow * 3 * 16 + 20] = (-norm2).GetZ();
                    this->normalArrow[counterArrow * 3 * 16 + 21] = (-norm2).GetX();
                    this->normalArrow[counterArrow * 3 * 16 + 22] = (-norm2).GetY();
                    this->normalArrow[counterArrow * 3 * 16 + 23] = (-norm2).GetZ();
                    this->normalArrow[counterArrow * 3 * 16 + 24] = (dir1).GetX();
                    this->normalArrow[counterArrow * 3 * 16 + 25] = (dir1).GetY();
                    this->normalArrow[counterArrow * 3 * 16 + 26] = (dir1).GetZ();
                    this->normalArrow[counterArrow * 3 * 16 + 27] = (dir1).GetX();
                    this->normalArrow[counterArrow * 3 * 16 + 28] = (dir1).GetY();
                    this->normalArrow[counterArrow * 3 * 16 + 29] = (dir1).GetZ();
                    this->normalArrow[counterArrow * 3 * 16 + 30] = (dir2).GetX();
                    this->normalArrow[counterArrow * 3 * 16 + 31] = (dir2).GetY();
                    this->normalArrow[counterArrow * 3 * 16 + 32] = (dir2).GetZ();
                    this->normalArrow[counterArrow * 3 * 16 + 33] = (dir2).GetX();
                    this->normalArrow[counterArrow * 3 * 16 + 34] = (dir2).GetY();
                    this->normalArrow[counterArrow * 3 * 16 + 35] = (dir2).GetZ();
                    this->normalArrow[counterArrow * 3 * 16 + 36] = (-dir1).GetX();
                    this->normalArrow[counterArrow * 3 * 16 + 37] = (-dir1).GetY();
                    this->normalArrow[counterArrow * 3 * 16 + 38] = (-dir1).GetZ();
                    this->normalArrow[counterArrow * 3 * 16 + 39] = (-dir1).GetX();
                    this->normalArrow[counterArrow * 3 * 16 + 40] = (-dir1).GetY();
                    this->normalArrow[counterArrow * 3 * 16 + 41] = (-dir1).GetZ();
                    this->normalArrow[counterArrow * 3 * 16 + 42] = (-dir2).GetX();
                    this->normalArrow[counterArrow * 3 * 16 + 43] = (-dir2).GetY();
                    this->normalArrow[counterArrow * 3 * 16 + 44] = (-dir2).GetZ();
                    this->normalArrow[counterArrow * 3 * 16 + 45] = (-dir2).GetX();
                    this->normalArrow[counterArrow * 3 * 16 + 46] = (-dir2).GetY();
                    this->normalArrow[counterArrow * 3 * 16 + 47] = (-dir2).GetZ();
                    counterArrow++;
                } else if (bSplineSecStruct[i][j / this->numberOfSplineSeg] ==
                           MolecularDataCall::SecStructure::TYPE_HELIX) {
                    this->vertHelix[counterHelix * 3 * 16 + 0] = (v2 - dir1 + norm1 * scale).GetX();
                    this->vertHelix[counterHelix * 3 * 16 + 1] = (v2 - dir1 + norm1 * scale).GetY();
                    this->vertHelix[counterHelix * 3 * 16 + 2] = (v2 - dir1 + norm1 * scale).GetZ();
                    this->vertHelix[counterHelix * 3 * 16 + 3] = (v2 + dir1 + norm1 * scale).GetX();
                    this->vertHelix[counterHelix * 3 * 16 + 4] = (v2 + dir1 + norm1 * scale).GetY();
                    this->vertHelix[counterHelix * 3 * 16 + 5] = (v2 + dir1 + norm1 * scale).GetZ();
                    this->vertHelix[counterHelix * 3 * 16 + 6] = (v3 + dir2 + norm2 * scale).GetX();
                    this->vertHelix[counterHelix * 3 * 16 + 7] = (v3 + dir2 + norm2 * scale).GetY();
                    this->vertHelix[counterHelix * 3 * 16 + 8] = (v3 + dir2 + norm2 * scale).GetZ();
                    this->vertHelix[counterHelix * 3 * 16 + 9] = (v3 - dir2 + norm2 * scale).GetX();
                    this->vertHelix[counterHelix * 3 * 16 + 10] = (v3 - dir2 + norm2 * scale).GetY();
                    this->vertHelix[counterHelix * 3 * 16 + 11] = (v3 - dir2 + norm2 * scale).GetZ();
                    this->vertHelix[counterHelix * 3 * 16 + 12] = (v2 - dir1 - norm1 * scale).GetX();
                    this->vertHelix[counterHelix * 3 * 16 + 13] = (v2 - dir1 - norm1 * scale).GetY();
                    this->vertHelix[counterHelix * 3 * 16 + 14] = (v2 - dir1 - norm1 * scale).GetZ();
                    this->vertHelix[counterHelix * 3 * 16 + 15] = (v2 + dir1 - norm1 * scale).GetX();
                    this->vertHelix[counterHelix * 3 * 16 + 16] = (v2 + dir1 - norm1 * scale).GetY();
                    this->vertHelix[counterHelix * 3 * 16 + 17] = (v2 + dir1 - norm1 * scale).GetZ();
                    this->vertHelix[counterHelix * 3 * 16 + 18] = (v3 + dir2 - norm2 * scale).GetX();
                    this->vertHelix[counterHelix * 3 * 16 + 19] = (v3 + dir2 - norm2 * scale).GetY();
                    this->vertHelix[counterHelix * 3 * 16 + 20] = (v3 + dir2 - norm2 * scale).GetZ();
                    this->vertHelix[counterHelix * 3 * 16 + 21] = (v3 - dir2 - norm2 * scale).GetX();
                    this->vertHelix[counterHelix * 3 * 16 + 22] = (v3 - dir2 - norm2 * scale).GetY();
                    this->vertHelix[counterHelix * 3 * 16 + 23] = (v3 - dir2 - norm2 * scale).GetZ();
                    this->vertHelix[counterHelix * 3 * 16 + 24] = (v2 + dir1 + norm1 * scale).GetX();
                    this->vertHelix[counterHelix * 3 * 16 + 25] = (v2 + dir1 + norm1 * scale).GetY();
                    this->vertHelix[counterHelix * 3 * 16 + 26] = (v2 + dir1 + norm1 * scale).GetZ();
                    this->vertHelix[counterHelix * 3 * 16 + 27] = (v2 + dir1 - norm1 * scale).GetX();
                    this->vertHelix[counterHelix * 3 * 16 + 28] = (v2 + dir1 - norm1 * scale).GetY();
                    this->vertHelix[counterHelix * 3 * 16 + 29] = (v2 + dir1 - norm1 * scale).GetZ();
                    this->vertHelix[counterHelix * 3 * 16 + 30] = (v3 + dir2 - norm2 * scale).GetX();
                    this->vertHelix[counterHelix * 3 * 16 + 31] = (v3 + dir2 - norm2 * scale).GetY();
                    this->vertHelix[counterHelix * 3 * 16 + 32] = (v3 + dir2 - norm2 * scale).GetZ();
                    this->vertHelix[counterHelix * 3 * 16 + 33] = (v3 + dir2 + norm2 * scale).GetX();
                    this->vertHelix[counterHelix * 3 * 16 + 34] = (v3 + dir2 + norm2 * scale).GetY();
                    this->vertHelix[counterHelix * 3 * 16 + 35] = (v3 + dir2 + norm2 * scale).GetZ();
                    this->vertHelix[counterHelix * 3 * 16 + 36] = (v2 - dir1 + norm1 * scale).GetX();
                    this->vertHelix[counterHelix * 3 * 16 + 37] = (v2 - dir1 + norm1 * scale).GetY();
                    this->vertHelix[counterHelix * 3 * 16 + 38] = (v2 - dir1 + norm1 * scale).GetZ();
                    this->vertHelix[counterHelix * 3 * 16 + 39] = (v2 - dir1 - norm1 * scale).GetX();
                    this->vertHelix[counterHelix * 3 * 16 + 40] = (v2 - dir1 - norm1 * scale).GetY();
                    this->vertHelix[counterHelix * 3 * 16 + 41] = (v2 - dir1 - norm1 * scale).GetZ();
                    this->vertHelix[counterHelix * 3 * 16 + 42] = (v3 - dir2 - norm2 * scale).GetX();
                    this->vertHelix[counterHelix * 3 * 16 + 43] = (v3 - dir2 - norm2 * scale).GetY();
                    this->vertHelix[counterHelix * 3 * 16 + 44] = (v3 - dir2 - norm2 * scale).GetZ();
                    this->vertHelix[counterHelix * 3 * 16 + 45] = (v3 - dir2 + norm2 * scale).GetX();
                    this->vertHelix[counterHelix * 3 * 16 + 46] = (v3 - dir2 + norm2 * scale).GetY();
                    this->vertHelix[counterHelix * 3 * 16 + 47] = (v3 - dir2 + norm2 * scale).GetZ();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 0] = col1.GetX();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 1] = col1.GetY();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 2] = col1.GetZ();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 3] = col1.GetX();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 4] = col1.GetY();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 5] = col1.GetZ();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 6] = col2.GetX();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 7] = col2.GetY();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 8] = col2.GetZ();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 9] = col2.GetX();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 10] = col2.GetY();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 11] = col2.GetZ();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 12] = col1.GetX();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 13] = col1.GetY();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 14] = col1.GetZ();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 15] = col1.GetX();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 16] = col1.GetY();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 17] = col1.GetZ();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 18] = col2.GetX();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 19] = col2.GetY();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 20] = col2.GetZ();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 21] = col2.GetX();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 22] = col2.GetY();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 23] = col2.GetZ();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 24] = col1.GetX();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 25] = col1.GetY();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 26] = col1.GetZ();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 27] = col1.GetX();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 28] = col1.GetY();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 29] = col1.GetZ();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 30] = col2.GetX();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 31] = col2.GetY();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 32] = col2.GetZ();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 33] = col2.GetX();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 34] = col2.GetY();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 35] = col2.GetZ();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 36] = col1.GetX();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 37] = col1.GetY();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 38] = col1.GetZ();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 39] = col1.GetX();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 40] = col1.GetY();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 41] = col1.GetZ();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 42] = col2.GetX();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 43] = col2.GetY();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 44] = col2.GetZ();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 45] = col2.GetX();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 46] = col2.GetY();
                    this->colorsParamsHelix[counterHelix * 3 * 16 + 47] = col2.GetZ();
                    norm1.Normalise();
                    norm2.Normalise();
                    dir1.Normalise();
                    dir2.Normalise();
                    this->normalHelix[counterHelix * 3 * 16 + 0] = (norm1).GetX();
                    this->normalHelix[counterHelix * 3 * 16 + 1] = (norm1).GetY();
                    this->normalHelix[counterHelix * 3 * 16 + 2] = (norm1).GetZ();
                    this->normalHelix[counterHelix * 3 * 16 + 3] = (norm1).GetX();
                    this->normalHelix[counterHelix * 3 * 16 + 4] = (norm1).GetY();
                    this->normalHelix[counterHelix * 3 * 16 + 5] = (norm1).GetZ();
                    this->normalHelix[counterHelix * 3 * 16 + 6] = (norm2).GetX();
                    this->normalHelix[counterHelix * 3 * 16 + 7] = (norm2).GetY();
                    this->normalHelix[counterHelix * 3 * 16 + 8] = (norm2).GetZ();
                    this->normalHelix[counterHelix * 3 * 16 + 9] = (norm2).GetX();
                    this->normalHelix[counterHelix * 3 * 16 + 10] = (norm2).GetY();
                    this->normalHelix[counterHelix * 3 * 16 + 11] = (norm2).GetZ();
                    this->normalHelix[counterHelix * 3 * 16 + 12] = (-norm1).GetX();
                    this->normalHelix[counterHelix * 3 * 16 + 13] = (-norm1).GetY();
                    this->normalHelix[counterHelix * 3 * 16 + 14] = (-norm1).GetZ();
                    this->normalHelix[counterHelix * 3 * 16 + 15] = (-norm1).GetX();
                    this->normalHelix[counterHelix * 3 * 16 + 16] = (-norm1).GetY();
                    this->normalHelix[counterHelix * 3 * 16 + 17] = (-norm1).GetZ();
                    this->normalHelix[counterHelix * 3 * 16 + 18] = (-norm2).GetX();
                    this->normalHelix[counterHelix * 3 * 16 + 19] = (-norm2).GetY();
                    this->normalHelix[counterHelix * 3 * 16 + 20] = (-norm2).GetZ();
                    this->normalHelix[counterHelix * 3 * 16 + 21] = (-norm2).GetX();
                    this->normalHelix[counterHelix * 3 * 16 + 22] = (-norm2).GetY();
                    this->normalHelix[counterHelix * 3 * 16 + 23] = (-norm2).GetZ();
                    this->normalHelix[counterHelix * 3 * 16 + 24] = (dir1).GetX();
                    this->normalHelix[counterHelix * 3 * 16 + 25] = (dir1).GetY();
                    this->normalHelix[counterHelix * 3 * 16 + 26] = (dir1).GetZ();
                    this->normalHelix[counterHelix * 3 * 16 + 27] = (dir1).GetX();
                    this->normalHelix[counterHelix * 3 * 16 + 28] = (dir1).GetY();
                    this->normalHelix[counterHelix * 3 * 16 + 29] = (dir1).GetZ();
                    this->normalHelix[counterHelix * 3 * 16 + 30] = (dir2).GetX();
                    this->normalHelix[counterHelix * 3 * 16 + 31] = (dir2).GetY();
                    this->normalHelix[counterHelix * 3 * 16 + 32] = (dir2).GetZ();
                    this->normalHelix[counterHelix * 3 * 16 + 33] = (dir2).GetX();
                    this->normalHelix[counterHelix * 3 * 16 + 34] = (dir2).GetY();
                    this->normalHelix[counterHelix * 3 * 16 + 35] = (dir2).GetZ();
                    this->normalHelix[counterHelix * 3 * 16 + 36] = (-dir1).GetX();
                    this->normalHelix[counterHelix * 3 * 16 + 37] = (-dir1).GetY();
                    this->normalHelix[counterHelix * 3 * 16 + 38] = (-dir1).GetZ();
                    this->normalHelix[counterHelix * 3 * 16 + 39] = (-dir1).GetX();
                    this->normalHelix[counterHelix * 3 * 16 + 40] = (-dir1).GetY();
                    this->normalHelix[counterHelix * 3 * 16 + 41] = (-dir1).GetZ();
                    this->normalHelix[counterHelix * 3 * 16 + 42] = (-dir2).GetX();
                    this->normalHelix[counterHelix * 3 * 16 + 43] = (-dir2).GetY();
                    this->normalHelix[counterHelix * 3 * 16 + 44] = (-dir2).GetZ();
                    this->normalHelix[counterHelix * 3 * 16 + 45] = (-dir2).GetX();
                    this->normalHelix[counterHelix * 3 * 16 + 46] = (-dir2).GetY();
                    this->normalHelix[counterHelix * 3 * 16 + 47] = (-dir2).GetZ();
                    counterHelix++;
                } else {
                    dir1 = dir1 * scale;
                    dir2 = dir2 * scale;

                    for (unsigned int k = 0; k < this->numberOfTubeSeg; k++) {
                        alpha = (float(2.0 * M_PI) / float(this->numberOfTubeSeg)) * float(k);
                        q1.Set(alpha, ortho1);
                        q2.Set(alpha, ortho2);
                        res1 = q1 * dir1;
                        res2 = q2 * dir2;

                        // v1
                        this->vertTube[counterTube] = (v2 + res1).GetX();
                        this->colorsParamsTube[counterTube] = col1.GetX();
                        counterTube++;
                        this->vertTube[counterTube] = (v2 + res1).GetY();
                        this->colorsParamsTube[counterTube] = col1.GetY();
                        counterTube++;
                        this->vertTube[counterTube] = (v2 + res1).GetZ();
                        this->colorsParamsTube[counterTube] = col1.GetZ();
                        counterTube++;
                        res1.Normalise();
                        this->normalTube[counterTube - 3] = res1.GetX();
                        this->normalTube[counterTube - 2] = res1.GetY();
                        this->normalTube[counterTube - 1] = res1.GetZ();
                        // v2
                        this->vertTube[counterTube] = (v3 + res2).GetX();
                        this->colorsParamsTube[counterTube] = col2.GetX();
                        counterTube++;
                        this->vertTube[counterTube] = (v3 + res2).GetY();
                        this->colorsParamsTube[counterTube] = col2.GetY();
                        counterTube++;
                        this->vertTube[counterTube] = (v3 + res2).GetZ();
                        this->colorsParamsTube[counterTube] = col2.GetZ();
                        counterTube++;
                        res2.Normalise();
                        this->normalTube[counterTube - 3] = res2.GetX();
                        this->normalTube[counterTube - 2] = res2.GetY();
                        this->normalTube[counterTube - 1] = res2.GetZ();

                        alpha = (float(2.0f * M_PI) / float(this->numberOfTubeSeg)) * float(k + 1);
                        q1.Set(alpha, ortho1);
                        q2.Set(alpha, ortho2);
                        res1 = q1 * dir1;
                        res2 = q2 * dir2;

                        // v3
                        this->vertTube[counterTube] = (v3 + res2).GetX();
                        this->colorsParamsTube[counterTube] = col2.GetX();
                        counterTube++;
                        this->vertTube[counterTube] = (v3 + res2).GetY();
                        this->colorsParamsTube[counterTube] = col2.GetY();
                        counterTube++;
                        this->vertTube[counterTube] = (v3 + res2).GetZ();
                        this->colorsParamsTube[counterTube] = col2.GetZ();
                        counterTube++;
                        res2.Normalise();
                        this->normalTube[counterTube - 3] = res2.GetX();
                        this->normalTube[counterTube - 2] = res2.GetY();
                        this->normalTube[counterTube - 1] = res2.GetZ();
                        // v4
                        this->vertTube[counterTube] = (v2 + res1).GetX();
                        this->colorsParamsTube[counterTube] = col1.GetX();
                        counterTube++;
                        this->vertTube[counterTube] = (v2 + res1).GetY();
                        this->colorsParamsTube[counterTube] = col1.GetY();
                        counterTube++;
                        this->vertTube[counterTube] = (v2 + res1).GetZ();
                        this->colorsParamsTube[counterTube] = col1.GetZ();
                        counterTube++;
                        res1.Normalise();
                        this->normalTube[counterTube - 3] = res1.GetX();
                        this->normalTube[counterTube - 2] = res1.GetY();
                        this->normalTube[counterTube - 1] = res1.GetZ();
                    }
                }
            }
        }
        // --- END store vertex/color/inparams ---

        // set cartoon CPU as created
        this->prepareCartoonCPU = false;
    }

    float spec[4] = {1.0f, 1.0f, 1.0f, 1.0f};
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, spec);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 50.0f);
    glEnable(GL_COLOR_MATERIAL);
    glDisable(GL_LIGHTING);

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);

    this->lightShader->use();
    // glDisable(GL_LIGHTING);

    // tube
    glVertexPointer(3, GL_FLOAT, 0, this->vertTube);
    glNormalPointer(GL_FLOAT, 0, this->normalTube);
    glColorPointer(3, GL_FLOAT, 0, this->colorsParamsTube);
    glDrawArrays(GL_QUADS, 0, this->totalCountTube * 4 * this->numberOfTubeSeg);

    // arrow
    glVertexPointer(3, GL_FLOAT, 0, this->vertArrow);
    glNormalPointer(GL_FLOAT, 0, this->normalArrow);
    glColorPointer(3, GL_FLOAT, 0, this->colorsParamsArrow);
    glDrawArrays(GL_QUADS, 0, this->totalCountArrow * 16);

    // helix
    glVertexPointer(3, GL_FLOAT, 0, this->vertHelix);
    glColorPointer(3, GL_FLOAT, 0, this->colorsParamsHelix);
    glNormalPointer(GL_FLOAT, 0, this->normalHelix);
    glDrawArrays(GL_QUADS, 0, this->totalCountHelix * 16);

    glUseProgram(0);

    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisable(GL_COLOR_MATERIAL);
}


/*
 * MoleculeCartoonRenderer::RenderCartoonLineCPU
 */
void MoleculeCartoonRenderer::RenderCartoonLineCPU(const MolecularDataCall* mol, float* atomPos) {

    // prepare cartoon line representation, if necessary
    if (this->prepareCartoonLine) {
        unsigned int cntChain, cntS, cntAA, idx, firstSS, countSS, firstAA, countAA;
        // B-Spline
        BSpline bSpline;
        // control points for the first (center) b-spline
        std::vector<vislib::math::Vector<float, 3>> controlPoints;
        // temporary vectors
        vislib::math::Vector<float, 3> vecCA, vecC, vecO, vecTmp, vecTmpOld;
        // temporary color
        // const float *color;
        // temporary color vector
        vislib::math::Vector<float, 3> colorVec;

        // secondary structure type for b-spline
        std::vector<std::vector<MolecularDataCall::SecStructure::ElementType>> bSplineSecStruct;


        // set the number of segments to create
        bSpline.setN(this->numberOfSplineSeg);

        // resize result vector for coordinates of first b-spline segments
        bSplineCoordsCPU.resize(mol->MoleculeCount());
        // resize vector for secondary structure
        bSplineSecStruct.resize(mol->MoleculeCount());
        // resize color vector
        this->cartoonColorCPU.resize(mol->MoleculeCount());

        // --- compute the b-splines ---
        // loop over all chains
        MolecularDataCall::Molecule chain;
        MolecularDataCall::AminoAcid* aminoacid;
        for (cntChain = 0; cntChain < mol->MoleculeCount(); ++cntChain) {
            chain = mol->Molecules()[cntChain];
            controlPoints.clear();
            this->cartoonColorCPU[cntChain].clear();
            // check if the first residue is an amino acid
            if (mol->Residues()[chain.FirstResidueIndex()]->Identifier() != MolecularDataCall::Residue::AMINOACID) {
                continue;
            }
            firstSS = chain.FirstSecStructIndex();
            countSS = firstSS + chain.SecStructCount();
            // loop over all secondary structure elements
            for (cntS = firstSS; cntS < countSS; ++cntS) {
                firstAA = mol->SecondaryStructures()[cntS].FirstAminoAcidIndex();
                countAA = firstAA + mol->SecondaryStructures()[cntS].AminoAcidCount();
                // loop over all amino acids in the current sec struct
                for (cntAA = firstAA; cntAA < countAA; ++cntAA) {
                    // add sec struct type
                    bSplineSecStruct[cntChain].push_back(mol->SecondaryStructures()[cntS].Type());
                    // get the index of the C-alpha atom
                    if (mol->Residues()[cntAA]->Identifier() == MolecularDataCall::Residue::AMINOACID)
                        aminoacid = (MolecularDataCall::AminoAcid*)(mol->Residues()[cntAA]);
                    else
                        continue;
                    idx = aminoacid->CAlphaIndex();
                    // get the coordinates of the C-alpha atom
                    vecCA.SetX(atomPos[idx * 3 + 0]);
                    vecCA.SetY(atomPos[idx * 3 + 1]);
                    vecCA.SetZ(atomPos[idx * 3 + 2]);
                    // add the C-alpha coords to the list of control points
                    controlPoints.push_back(vecCA);

                    // add the color of the C-alpha atom to the color vector
                    colorVec.SetX(this->atomColorTable[idx].x);
                    colorVec.SetY(this->atomColorTable[idx].y);
                    colorVec.SetZ(this->atomColorTable[idx].z);
                    cartoonColorCPU[cntChain].push_back(colorVec);
                }
            }
            // set the control points, compute the first spline and fetch the result
            bSpline.setBackbone(controlPoints);
            if (bSpline.computeSpline())
                bSpline.getResult(bSplineCoordsCPU[cntChain]);
            else
                continue; // --> return if spline could not be computed
        }

        this->prepareCartoonLine = false;
    }


    if (!this->prepareCartoonLine) {
        glDisable(GL_LIGHTING);
        glColor3f(1.0f, 1.0f, 1.0f);
        for (unsigned int j = 0; j < bSplineCoordsCPU.size(); j++) {
            glBegin(GL_LINE_STRIP);
            for (unsigned int i = 0; i < bSplineCoordsCPU[j].size(); i++) {
                glColor3fv(cartoonColorCPU[j][i / this->numberOfSplineSeg].PeekComponents());
                glVertex3fv(bSplineCoordsCPU[j][i].PeekComponents());
            }
            glEnd(); // GL_LINE_STRIP
        }
    }
}


/*
 * Render protein in geometry shader CARTOON_GPU mode
 */
void MoleculeCartoonRenderer::RenderCartoonGPU(const MolecularDataCall* mol, float* atomPos) {
    // ------------------------------------------------------------
    // --- CARTOON SPLINE                                       ---
    // --- GPU Implementation                                   ---
    // --- use geometry shader for whole computation            ---
    // ------------------------------------------------------------

    unsigned int cntChain, cntS, cntAA, idx, firstSS, countSS, firstAA, countAA;

    float spec[4] = {1.0f, 1.0f, 1.0f, 1.0f};
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, spec);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 50.0f);
    glEnable(GL_COLOR_MATERIAL);

    vislib::math::Vector<float, 3> v0, v1, v2, v3, v4, v5;
    vislib::math::Vector<float, 3> n1, n2, n3, n4;
    vislib::math::Vector<float, 3> color;
    float flip = 1.0f;
    float factor = 0.0f;

    // loop over all molecules
    MolecularDataCall::Molecule chain;
    MolecularDataCall::AminoAcid* aminoacid;
    for (cntChain = 0; cntChain < mol->MoleculeCount(); ++cntChain) {
        chain = mol->Molecules()[cntChain];
        // check if the first residue is an amino acid
        if (mol->Residues()[chain.FirstResidueIndex()]->Identifier() != MolecularDataCall::Residue::AMINOACID) {
            continue;
        }
        firstSS = chain.FirstSecStructIndex();
        countSS = firstSS + chain.SecStructCount();
        // loop over all secondary structure elements
        for (cntS = firstSS; cntS < countSS; ++cntS) {
            firstAA = mol->SecondaryStructures()[cntS].FirstAminoAcidIndex();
            countAA = firstAA + mol->SecondaryStructures()[cntS].AminoAcidCount();
            // loop over all amino acids in the current sec struct
            for (cntAA = firstAA; cntAA < countAA; ++cntAA) {
                idx = chain.FirstResidueIndex() + chain.ResidueCount();
                if (cntAA + 3 >= idx)
                    continue;
                // do nothing if the current residue is no amino acid
                if (mol->Residues()[cntAA]->Identifier() != MolecularDataCall::Residue::AMINOACID)
                    continue;

                idx = cntS;
                while (idx + 1 < mol->SecondaryStructureCount() &&
                       cntAA + 2 >= mol->SecondaryStructures()[idx].FirstAminoAcidIndex() +
                                        mol->SecondaryStructures()[idx].AminoAcidCount()) {
                    idx++;
                }

                if (mol->SecondaryStructures()[idx].Type() == MolecularDataCall::SecStructure::TYPE_HELIX) {
                    this->helixSplineShader->use();
                } else if (mol->SecondaryStructures()[idx].Type() == MolecularDataCall::SecStructure::TYPE_SHEET) {
                    this->arrowSplineShader->use();
                    if ((cntAA + 3) == countAA) {
                        factor = 1.0f;
                    } else {
                        factor = 0.0f;
                    }
                } else {
                    this->tubeSplineShader->use();
                }

                glBegin(GL_LINES_ADJACENCY_EXT);

                // vertex 1
                aminoacid = (MolecularDataCall::AminoAcid*)(mol->Residues()[cntAA]);
                idx = aminoacid->CAlphaIndex();
                v1.SetX(atomPos[idx * 3 + 0]);
                v1.SetY(atomPos[idx * 3 + 1]);
                v1.SetZ(atomPos[idx * 3 + 2]);
                idx = aminoacid->CCarbIndex();
                v5.SetX(atomPos[idx * 3 + 0]);
                v5.SetY(atomPos[idx * 3 + 1]);
                v5.SetZ(atomPos[idx * 3 + 2]);
                idx = aminoacid->OIndex();
                n1.SetX(atomPos[idx * 3 + 0]);
                n1.SetY(atomPos[idx * 3 + 1]);
                n1.SetZ(atomPos[idx * 3 + 2]);
                n1 = n1 - v5;
                n1.Normalise();
                if (cntAA > firstAA && n3.Dot(n1) < 0.0f)
                    flip = -1.0;
                else
                    flip = 1.0;
                n1 *= flip;
                aminoacid = (MolecularDataCall::AminoAcid*)(mol->Residues()[cntAA + 2]);
                idx = aminoacid->CAlphaIndex();
                glSecondaryColor3f(
                    this->atomColorTable[idx].x, this->atomColorTable[idx].y, this->atomColorTable[idx].z);
                glColor3fv(n1.PeekComponents());
                glVertex3fv(v1.PeekComponents());

                // vertex 2
                aminoacid = (MolecularDataCall::AminoAcid*)(mol->Residues()[cntAA + 1]);
                idx = aminoacid->CAlphaIndex();
                v2.SetX(atomPos[idx * 3 + 0]);
                v2.SetY(atomPos[idx * 3 + 1]);
                v2.SetZ(atomPos[idx * 3 + 2]);
                idx = aminoacid->CCarbIndex();
                v5.SetX(atomPos[idx * 3 + 0]);
                v5.SetY(atomPos[idx * 3 + 1]);
                v5.SetZ(atomPos[idx * 3 + 2]);
                idx = aminoacid->OIndex();
                n2.SetX(atomPos[idx * 3 + 0]);
                n2.SetY(atomPos[idx * 3 + 1]);
                n2.SetZ(atomPos[idx * 3 + 2]);
                n2 = n2 - v5;
                n2.Normalise();
                if (n1.Dot(n2) < 0.0f)
                    flip = -1.0;
                else
                    flip = 1.0;
                n2 *= flip;
                glSecondaryColor3f(0.2f, 1.0f, factor);
                glColor3fv(n2.PeekComponents());
                glVertex3fv(v2.PeekComponents());

                // vertex 3
                aminoacid = (MolecularDataCall::AminoAcid*)(mol->Residues()[cntAA + 2]);
                idx = aminoacid->CAlphaIndex();
                v3.SetX(atomPos[idx * 3 + 0]);
                v3.SetY(atomPos[idx * 3 + 1]);
                v3.SetZ(atomPos[idx * 3 + 2]);
                idx = aminoacid->CCarbIndex();
                v5.SetX(atomPos[idx * 3 + 0]);
                v5.SetY(atomPos[idx * 3 + 1]);
                v5.SetZ(atomPos[idx * 3 + 2]);
                idx = aminoacid->OIndex();
                n3.SetX(atomPos[idx * 3 + 0]);
                n3.SetY(atomPos[idx * 3 + 1]);
                n3.SetZ(atomPos[idx * 3 + 2]);
                n3 = n3 - v5;
                n3.Normalise();
                if (n2.Dot(n3) < 0.0f)
                    flip = -1.0;
                else
                    flip = 1.0;
                n3 *= flip;
                glColor3fv(n3.PeekComponents());
                glVertex3fv(v3.PeekComponents());

                // vertex 4
                aminoacid = (MolecularDataCall::AminoAcid*)(mol->Residues()[cntAA + 3]);
                idx = aminoacid->CAlphaIndex();
                v4.SetX(atomPos[idx * 3 + 0]);
                v4.SetY(atomPos[idx * 3 + 1]);
                v4.SetZ(atomPos[idx * 3 + 2]);
                idx = aminoacid->CCarbIndex();
                v5.SetX(atomPos[idx * 3 + 0]);
                v5.SetY(atomPos[idx * 3 + 1]);
                v5.SetZ(atomPos[idx * 3 + 2]);
                idx = aminoacid->OIndex();
                n4.SetX(atomPos[idx * 3 + 0]);
                n4.SetY(atomPos[idx * 3 + 1]);
                n4.SetZ(atomPos[idx * 3 + 2]);
                n4 = n4 - v5;
                n4.Normalise();
                if (n3.Dot(n4) < 0.0f)
                    flip = -1.0;
                else
                    flip = 1.0;
                n4 *= flip;
                glColor3fv(n4.PeekComponents());
                glVertex3fv(v4.PeekComponents());

                // store last vertex for comparison (flip)
                n3 = n1;

                glEnd(); // GL_LINES_ADJACENCY_EXT

                glUseProgram(0);
            }
        }
    }
    glDisable(GL_COLOR_MATERIAL);
}


/*
 * Render protein in geometry shader CARTOON_GPU mode
 */
void MoleculeCartoonRenderer::RenderCartoonGPUTubeOnly(const MolecularDataCall* mol, float* atomPos) {
    // prepare hybrid cartoon representation, if necessary
    if (this->prepareCartoonHybrid) {
        unsigned int cntChain, cntS, cntAA, idx, firstSS, countSS, firstAA, countAA;
        // B-Spline
        BSpline bSpline;
        // control points for the first (center) b-spline
        std::vector<vislib::math::Vector<float, 3>> controlPoints;
        // control points for the second (direction) b-spline
        std::vector<vislib::math::Vector<float, 3>> controlPointsDir;
        // temporary vectors
        vislib::math::Vector<float, 3> vecCA, vecC, vecO, vecTmp, vecTmpOld;
        // temporary color
        // const float *color;
        // temporary color vector
        vislib::math::Vector<float, 3> colorVec;

        // coordinates of the first (center) b-spline (result of the spline computation)
        std::vector<std::vector<vislib::math::Vector<float, 3>>> bSplineCoords;
        // coordinates of the second (direction) b-spline (result of the spline computation)
        std::vector<std::vector<vislib::math::Vector<float, 3>>> bSplineCoordsDir;
        // secondary structure type for b-spline
        std::vector<std::vector<MolecularDataCall::SecStructure::ElementType>> bSplineSecStruct;
        // color of secondary structure b-spline
        std::vector<std::vector<vislib::math::Vector<float, 3>>> cartoonColor;

        // set the number of segments to create
        bSpline.setN(this->numberOfSplineSeg);

        // resize result vector for coordinates of first b-spline segments
        bSplineCoords.resize(mol->MoleculeCount());
        // resize result vector for coordinates of second b-spline segments
        bSplineCoordsDir.resize(mol->MoleculeCount());
        // resize vector for secondary structure
        bSplineSecStruct.resize(mol->MoleculeCount());
        // resize color vector
        cartoonColor.resize(mol->MoleculeCount());

        // --- compute the b-splines ---
        // loop over all chains
        MolecularDataCall::Molecule chain;
        MolecularDataCall::AminoAcid* aminoacid;
        for (cntChain = 0; cntChain < mol->MoleculeCount(); ++cntChain) {
            chain = mol->Molecules()[cntChain];
            controlPoints.clear();
            controlPointsDir.clear();
            // check if the first residue is an amino acid
            if (mol->Residues()[chain.FirstResidueIndex()]->Identifier() != MolecularDataCall::Residue::AMINOACID) {
                continue;
            }
            firstSS = chain.FirstSecStructIndex();
            countSS = firstSS + chain.SecStructCount();
            // loop over all secondary structure elements
            for (cntS = firstSS; cntS < countSS; ++cntS) {
                firstAA = mol->SecondaryStructures()[cntS].FirstAminoAcidIndex();
                countAA = firstAA + mol->SecondaryStructures()[cntS].AminoAcidCount();
                // loop over all amino acids in the current sec struct
                for (cntAA = firstAA; cntAA < countAA; ++cntAA) {
                    // add sec struct type
                    bSplineSecStruct[cntChain].push_back(mol->SecondaryStructures()[cntS].Type());
                    // get the index of the C-alpha atom
                    if (mol->Residues()[cntAA]->Identifier() == MolecularDataCall::Residue::AMINOACID)
                        aminoacid = (MolecularDataCall::AminoAcid*)(mol->Residues()[cntAA]);
                    else
                        continue;
                    idx = aminoacid->CAlphaIndex();
                    // get the coordinates of the C-alpha atom
                    vecCA.SetX(atomPos[idx * 3 + 0]);
                    vecCA.SetY(atomPos[idx * 3 + 1]);
                    vecCA.SetZ(atomPos[idx * 3 + 2]);
                    // add the C-alpha coords to the list of control points
                    controlPoints.push_back(vecCA);

                    // add the color of the C-alpha atom to the color vector
                    colorVec.SetX(this->atomColorTable[idx].x);
                    colorVec.SetY(this->atomColorTable[idx].y);
                    colorVec.SetZ(this->atomColorTable[idx].z);
                    cartoonColor[cntChain].push_back(colorVec);

                    // get the index of the C atom
                    idx = aminoacid->CCarbIndex();
                    // get the coordinates of the C-alpha atom
                    vecC.SetX(atomPos[idx * 3 + 0]);
                    vecC.SetY(atomPos[idx * 3 + 1]);
                    vecC.SetZ(atomPos[idx * 3 + 2]);

                    // get the index of the O atom
                    idx = aminoacid->OIndex();
                    // get the coordinates of the C-alpha atom
                    vecO.SetX(atomPos[idx * 3 + 0]);
                    vecO.SetY(atomPos[idx * 3 + 1]);
                    vecO.SetZ(atomPos[idx * 3 + 2]);

                    // compute control point of the second b-spline
                    vecTmp = vecO - vecC;
                    vecTmp.Normalise();
                    // check, if vector should be flipped
                    if (cntS > 0 && vecTmpOld.Dot(vecTmp) < 0.0f)
                        vecTmp = vecTmp * -1.0f;
                    vecTmpOld = vecTmp;
                    // add control point for the second b-spline to the list of control points
                    controlPointsDir.push_back(vecTmp + vecCA);
                }
            }
            // set the control points, compute the first spline and fetch the result
            bSpline.setBackbone(controlPoints);
            if (bSpline.computeSpline())
                bSpline.getResult(bSplineCoords[cntChain]);
            else
                continue; // --> return if spline could not be computed

            // set the control points, compute the second spline and fetch the result
            bSpline.setBackbone(controlPointsDir);
            if (bSpline.computeSpline())
                bSpline.getResult(bSplineCoordsDir[cntChain]);
            else
                continue; // --> return if spline could not be computed
        }

        // --- START store the vertices, colors and parameters ---
        this->totalCountTube = 0;
        for (unsigned int i = 0; i < bSplineCoords.size(); i++) {
            this->totalCountTube += (unsigned int)bSplineCoords[i].size();
        }

        if (this->vertTube)
            delete[] this->vertTube;
        if (this->colorsParamsTube)
            delete[] this->colorsParamsTube;

        this->vertTube = new float[this->totalCountTube * 6 * 3];
        this->colorsParamsTube = new float[this->totalCountTube * 6 * 3];

        // auxiliary variables
        float start, end, f1, f2, type = 0.0f;
        unsigned int counterTube = 0;
        vislib::math::Vector<float, 3> col1, col2;
        // compute the inner b-spline (backbone)
        for (unsigned int i = 0; i < bSplineCoords.size(); i++) {
            if (bSplineCoords[i].size() == 0)
                continue;
            for (unsigned int j = 2; j < bSplineCoords[i].size() - 1; j++) {
                start = end = -1.0f;
                f1 = f2 = 1.0f;

                // get the colors
                if (this->smoothCartoonColoringMode && j / this->numberOfSplineSeg > 0) {
                    col1 = cartoonColor[i][j / this->numberOfSplineSeg] * float(j % this->numberOfSplineSeg) /
                               float(this->numberOfSplineSeg - 1) +
                           cartoonColor[i][j / this->numberOfSplineSeg - 1] *
                               float((this->numberOfSplineSeg - 1) - j % this->numberOfSplineSeg) /
                               float(this->numberOfSplineSeg - 1);
                    int k = j + 1;
                    if (j % this->numberOfSplineSeg == this->numberOfSplineSeg - 1)
                        k = this->numberOfSplineSeg - 1;
                    col2 = cartoonColor[i][j / this->numberOfSplineSeg] * float(k % this->numberOfSplineSeg) /
                               float(this->numberOfSplineSeg - 1) +
                           cartoonColor[i][j / this->numberOfSplineSeg - 1] *
                               float((this->numberOfSplineSeg - 1) - k % this->numberOfSplineSeg) /
                               float(this->numberOfSplineSeg - 1);
                } else {
                    col1 = cartoonColor[i][j / this->numberOfSplineSeg];
                    col2 = cartoonColor[i][j / this->numberOfSplineSeg];
                }

                this->colorsParamsTube[counterTube * 6 * 3 + 0] = col1.GetX();
                this->colorsParamsTube[counterTube * 6 * 3 + 1] = col1.GetY();
                this->colorsParamsTube[counterTube * 6 * 3 + 2] = col1.GetZ();
                this->colorsParamsTube[counterTube * 6 * 3 + 3] =
                    this->tubeRadiusParam.Param<param::FloatParam>()->Value();
                this->colorsParamsTube[counterTube * 6 * 3 + 4] = f1;
                this->colorsParamsTube[counterTube * 6 * 3 + 5] = f2;
                this->colorsParamsTube[counterTube * 6 * 3 + 6] = type;
                this->colorsParamsTube[counterTube * 6 * 3 + 7] = start;
                this->colorsParamsTube[counterTube * 6 * 3 + 8] = end;
                this->colorsParamsTube[counterTube * 6 * 3 + 9] = col2.GetX();
                this->colorsParamsTube[counterTube * 6 * 3 + 10] = col2.GetY();
                this->colorsParamsTube[counterTube * 6 * 3 + 11] = col2.GetZ();
                this->colorsParamsTube[counterTube * 6 * 3 + 12] = 0.0f;
                this->colorsParamsTube[counterTube * 6 * 3 + 13] = 0.0f;
                this->colorsParamsTube[counterTube * 6 * 3 + 14] = 0.0f;
                this->colorsParamsTube[counterTube * 6 * 3 + 15] = 0.0f;
                this->colorsParamsTube[counterTube * 6 * 3 + 16] = 0.0f;
                this->colorsParamsTube[counterTube * 6 * 3 + 17] = 0.0f;
                this->vertTube[counterTube * 6 * 3 + 0] = bSplineCoords[i][j - 2].GetX();
                this->vertTube[counterTube * 6 * 3 + 1] = bSplineCoords[i][j - 2].GetY();
                this->vertTube[counterTube * 6 * 3 + 2] = bSplineCoords[i][j - 2].GetZ();
                this->vertTube[counterTube * 6 * 3 + 3] = bSplineCoordsDir[i][j - 1].GetX();
                this->vertTube[counterTube * 6 * 3 + 4] = bSplineCoordsDir[i][j - 1].GetY();
                this->vertTube[counterTube * 6 * 3 + 5] = bSplineCoordsDir[i][j - 1].GetZ();
                this->vertTube[counterTube * 6 * 3 + 6] = bSplineCoords[i][j - 1].GetX();
                this->vertTube[counterTube * 6 * 3 + 7] = bSplineCoords[i][j - 1].GetY();
                this->vertTube[counterTube * 6 * 3 + 8] = bSplineCoords[i][j - 1].GetZ();
                this->vertTube[counterTube * 6 * 3 + 9] = bSplineCoords[i][j].GetX();
                this->vertTube[counterTube * 6 * 3 + 10] = bSplineCoords[i][j].GetY();
                this->vertTube[counterTube * 6 * 3 + 11] = bSplineCoords[i][j].GetZ();
                this->vertTube[counterTube * 6 * 3 + 12] = bSplineCoordsDir[i][j].GetX();
                this->vertTube[counterTube * 6 * 3 + 13] = bSplineCoordsDir[i][j].GetY();
                this->vertTube[counterTube * 6 * 3 + 14] = bSplineCoordsDir[i][j].GetZ();
                this->vertTube[counterTube * 6 * 3 + 15] = bSplineCoords[i][j + 1].GetX();
                this->vertTube[counterTube * 6 * 3 + 16] = bSplineCoords[i][j + 1].GetY();
                this->vertTube[counterTube * 6 * 3 + 17] = bSplineCoords[i][j + 1].GetZ();
                counterTube++;
                //                }
            }
        }

        // --- END store vertex/color/inparams ---

        // set cartoon as created
        this->prepareCartoonHybrid = false;
    }

    float spec[4] = {1.0f, 1.0f, 1.0f, 1.0f};
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, spec);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 50.0f);
    glDisable(GL_COLOR_MATERIAL);

    // Get current window size
    float curVP[4];
    glGetFloatv(GL_VIEWPORT, curVP);

    // enable tube shader
    if (this->currentRenderMode == CartoonRenderMode::CARTOON)
        this->tubeShader->use();
    else
        this->tubeSimpleShader->use();
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
    glVertexPointer(3, GL_FLOAT, 0, this->vertTube);
    glColorPointer(3, GL_FLOAT, 0, this->colorsParamsTube);
    glDrawArrays(GL_TRIANGLES_ADJACENCY_EXT, 0, this->totalCountTube * 6);
    // disable tube shader
    glUseProgram(0);

    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
    glDisable(GL_COLOR_MATERIAL);
}


/*
 * MoleculeCartoonRenderer::RecomputeAll
 */
void MoleculeCartoonRenderer::RecomputeAll() {
    this->prepareCartoonHybrid = true;
    this->prepareCartoonCPU = true;
    this->prepareCartoonLine = true;

    this->atomColorTable.clear();
}


/*
 * Render the molecular data in stick mode.
 */
void MoleculeCartoonRenderer::RenderStick(
    const MolecularDataCall* mol, const float* atomPos, const BindingSiteCall* bs) {
    // ----- prepare stick raycasting -----

    /** vertex array for spheres */
    vislib::Array<float> vertSpheres;
    /** color array for spheres */
    vislib::Array<float> colorSpheres;
    /** vertex array for cylinders */
    vislib::Array<float> vertCylinders;
    /** attribute array for quaterinons of the cylinders */
    vislib::Array<float> quatCylinders;
    /** attribute array for inParam of the cylinders (radius and length) */
    vislib::Array<float> inParaCylinders;
    /** first color array for cylinder */
    vislib::Array<float> color1Cylinders;
    /** second color array for cylinder */
    vislib::Array<float> color2Cylinders;

    unsigned int totalAtomCnt = 0;
    unsigned int totalCylinderCnt = 0;

    vertSpheres.SetCount(mol->AtomCount() * 4);
    colorSpheres.SetCount(mol->AtomCount() * 3);
    vertCylinders.SetCount(mol->ConnectionCount() * 4);
    quatCylinders.SetCount(mol->ConnectionCount() * 4);
    inParaCylinders.SetCount(mol->ConnectionCount() * 2);
    color1Cylinders.SetCount(mol->ConnectionCount() * 3);
    color2Cylinders.SetCount(mol->ConnectionCount() * 3);
    this->atomVisible.SetCount(mol->AtomCount());

    unsigned int cnt, idx, molCnt, resCnt, atomCnt, conCnt, idxAtom, cntAtom;

    // reset atom visiblity
    for (cnt = 0; cnt < mol->AtomCount(); cnt++) {
        this->atomVisible[cnt] = false;
    }

    for (molCnt = 0; molCnt < mol->MoleculeCount(); ++molCnt) {
        idx = mol->Molecules()[molCnt].FirstResidueIndex();
        cnt = idx + mol->Molecules()[molCnt].ResidueCount();
        // do nothing if molecule is a protein
        if (mol->Residues()[idx]->Identifier() == MolecularDataCall::Residue::AMINOACID)
            continue;
        for (resCnt = idx; resCnt < cnt; ++resCnt) {
            idxAtom = mol->Residues()[resCnt]->FirstAtomIndex();
            cntAtom = idxAtom + mol->Residues()[resCnt]->AtomCount();
            for (atomCnt = idxAtom; atomCnt < cntAtom; ++atomCnt) {
                vertSpheres[4 * totalAtomCnt + 0] = atomPos[3 * atomCnt + 0];
                vertSpheres[4 * totalAtomCnt + 1] = atomPos[3 * atomCnt + 1];
                vertSpheres[4 * totalAtomCnt + 2] = atomPos[3 * atomCnt + 2];
                vertSpheres[4 * totalAtomCnt + 3] = this->stickRadiusParam.Param<param::FloatParam>()->Value();
                colorSpheres[3 * totalAtomCnt + 0] = this->atomColorTable[atomCnt].x;
                colorSpheres[3 * totalAtomCnt + 1] = this->atomColorTable[atomCnt].y;
                colorSpheres[3 * totalAtomCnt + 2] = this->atomColorTable[atomCnt].z;
                this->atomVisible[atomCnt] = true;
                totalAtomCnt++;
            }
        }
    }


    // search for binding sites if BindingSiteCall is available
    if (bs) {
        // temporary variables
        unsigned int firstMol;
        unsigned int firstRes;
        unsigned int firstAtom;
        unsigned int atomIdx;
        for (unsigned int cCnt = 0; cCnt < mol->ChainCount(); cCnt++) {
            firstMol = mol->Chains()[cCnt].FirstMoleculeIndex();
            for (unsigned int mCnt = firstMol; mCnt < firstMol + mol->Chains()[cCnt].MoleculeCount(); mCnt++) {
                firstRes = mol->Molecules()[mCnt].FirstResidueIndex();
                for (unsigned int rCnt = 0; rCnt < mol->Molecules()[mCnt].ResidueCount(); rCnt++) {
                    // try to match binding sites
                    std::pair<char, unsigned int> bsRes;
                    // loop over all binding sites
                    for (unsigned int bsCnt = 0; bsCnt < bs->GetBindingSiteCount(); bsCnt++) {
                        for (unsigned int bsResCnt = 0; bsResCnt < bs->GetBindingSite(bsCnt)->size(); bsResCnt++) {
                            bsRes = bs->GetBindingSite(bsCnt)->operator[](bsResCnt);
                            if (mol->Chains()[cCnt].Name() == bsRes.first &&
                                mol->Residues()[firstRes + rCnt]->OriginalResIndex() == bsRes.second &&
                                mol->ResidueTypeNames()[mol->Residues()[firstRes + rCnt]->Type()].PeekBuffer() ==
                                    bs->GetBindingSiteResNames(bsCnt)->operator[](bsResCnt)) {
                                // TODO loop over all atoms and add the color
                                firstAtom = mol->Residues()[firstRes + rCnt]->FirstAtomIndex();
                                for (unsigned int aCnt = 0; aCnt < mol->Residues()[firstRes + rCnt]->AtomCount();
                                     aCnt++) {
                                    atomIdx = firstAtom + aCnt;
                                    // position
                                    vertSpheres[4 * totalAtomCnt + 0] = atomPos[3 * atomIdx];
                                    vertSpheres[4 * totalAtomCnt + 1] = atomPos[3 * atomIdx + 1];
                                    vertSpheres[4 * totalAtomCnt + 2] = atomPos[3 * atomIdx + 2];
                                    vertSpheres[4 * totalAtomCnt + 3] =
                                        this->stickRadiusParam.Param<param::FloatParam>()->Value();
                                    // color
                                    colorSpheres[3 * totalAtomCnt + 0] = this->atomColorTable[atomIdx].x;
                                    colorSpheres[3 * totalAtomCnt + 1] = this->atomColorTable[atomIdx].y;
                                    colorSpheres[3 * totalAtomCnt + 2] = this->atomColorTable[atomIdx].z;
                                    // set visiblity
                                    this->atomVisible[atomIdx] = true;
                                    totalAtomCnt++;
                                }
                            }
                        }
                    }
                } // residues
            }     // molecules
        }         // chains
    }             // BindingSiteCall available

    unsigned int idx0, idx1;
    vislib::math::Vector<float, 3> firstAtomPos, secondAtomPos;
    vislib::math::Quaternion<float> quatC(0, 0, 0, 1);
    vislib::math::Vector<float, 3> tmpVec, ortho, dir, position;
    float angle;
    // loop over all connections and compute cylinder parameters
    for (molCnt = 0; molCnt < mol->MoleculeCount(); ++molCnt) {
        idx = mol->Molecules()[molCnt].FirstConnectionIndex();
        cnt = mol->Molecules()[molCnt].ConnectionCount();
        // do nothing if molecule is a protein and no binding sites are set
        if (!bs && mol->Residues()[mol->Molecules()[molCnt].FirstResidueIndex()]->Identifier() ==
                       MolecularDataCall::Residue::AMINOACID)
            continue;
        for (conCnt = 0; conCnt < cnt; ++conCnt) {
            idx0 = mol->Connection()[idx + 2 * conCnt];
            idx1 = mol->Connection()[idx + 2 * conCnt + 1];

            // do noting for invisible atoms
            if (!this->atomVisible[idx0] || !this->atomVisible[idx1])
                continue;

            firstAtomPos.SetX(atomPos[3 * idx0 + 0]);
            firstAtomPos.SetY(atomPos[3 * idx0 + 1]);
            firstAtomPos.SetZ(atomPos[3 * idx0 + 2]);

            secondAtomPos.SetX(atomPos[3 * idx1 + 0]);
            secondAtomPos.SetY(atomPos[3 * idx1 + 1]);
            secondAtomPos.SetZ(atomPos[3 * idx1 + 2]);

            // compute the quaternion for the rotation of the cylinder
            dir = secondAtomPos - firstAtomPos;
            tmpVec.Set(1.0f, 0.0f, 0.0f);
            angle = -tmpVec.Angle(dir);
            ortho = tmpVec.Cross(dir);
            ortho.Normalise();
            quatC.Set(angle, ortho);
            // compute the absolute position 'position' of the cylinder (center point)
            position = firstAtomPos + (dir / 2.0f);

            inParaCylinders[2 * totalCylinderCnt] = this->stickRadiusParam.Param<param::FloatParam>()->Value();
            inParaCylinders[2 * totalCylinderCnt + 1] = (firstAtomPos - secondAtomPos).Length();

            quatCylinders[4 * totalCylinderCnt + 0] = quatC.GetX();
            quatCylinders[4 * totalCylinderCnt + 1] = quatC.GetY();
            quatCylinders[4 * totalCylinderCnt + 2] = quatC.GetZ();
            quatCylinders[4 * totalCylinderCnt + 3] = quatC.GetW();

            color1Cylinders[3 * totalCylinderCnt + 0] = this->atomColorTable[idx0].x;
            color1Cylinders[3 * totalCylinderCnt + 1] = this->atomColorTable[idx0].y;
            color1Cylinders[3 * totalCylinderCnt + 2] = this->atomColorTable[idx0].z;

            color2Cylinders[3 * totalCylinderCnt + 0] = this->atomColorTable[idx1].x;
            color2Cylinders[3 * totalCylinderCnt + 1] = this->atomColorTable[idx1].y;
            color2Cylinders[3 * totalCylinderCnt + 2] = this->atomColorTable[idx1].z;

            vertCylinders[4 * totalCylinderCnt + 0] = position.X();
            vertCylinders[4 * totalCylinderCnt + 1] = position.Y();
            vertCylinders[4 * totalCylinderCnt + 2] = position.Z();
            vertCylinders[4 * totalCylinderCnt + 3] = 0.0f;

            totalCylinderCnt++;
        }
    }

    // ---------- actual rendering ----------
    // get viewpoint parameters for raycasting
    glm::vec4 viewportStuff = glm::vec4(0, 0, fbo->getWidth(), fbo->getHeight());
    if (viewportStuff[2] < 1.0f)
        viewportStuff[2] = 1.0f;
    if (viewportStuff[3] < 1.0f)
        viewportStuff[3] = 1.0f;
    viewportStuff[2] = 2.0f / viewportStuff[2];
    viewportStuff[3] = 2.0f / viewportStuff[3];

    auto camera_pose = camera.get<core::view::Camera::Pose>();

    buffers_[static_cast<int>(Buffers::POSITION)]->rebuffer(
        vertSpheres.PeekElements(), vertSpheres.Count() * sizeof(float));
    buffers_[static_cast<int>(Buffers::COLOR)]->rebuffer(
        colorSpheres.PeekElements(), colorSpheres.Count() * sizeof(float));

    glBindVertexArray(vertex_array_spheres_);

    // enable sphere shader
    sphereShader_->use();
    // set shader variables
    sphereShader_->setUniform("viewAttr", viewportStuff);
    sphereShader_->setUniform("camIn", camera_pose.direction);
    sphereShader_->setUniform("camRight", camera_pose.right);
    sphereShader_->setUniform("camUp", camera_pose.up);
    sphereShader_->setUniform("MVP", MVP);
    sphereShader_->setUniform("MVinv", MVinv);
    sphereShader_->setUniform("MVPinv", MVPinv);
    sphereShader_->setUniform("MVPtransp", MVPtransp);
    sphereShader_->setUniform("NormalM", NormalM);
    sphereShader_->setUniform("planes", planes);
    // draw them
    glDrawArrays(GL_POINTS, 0, totalAtomCnt);
    // disable sphere shader
    glUseProgram(0);
    glBindVertexArray(0);

    buffers_[static_cast<int>(Buffers::POSITION)]->rebuffer(
        vertCylinders.PeekElements(), vertCylinders.Count() * sizeof(float));
    buffers_[static_cast<int>(Buffers::COLOR)]->rebuffer(
        colorSpheres.PeekElements(), colorSpheres.Count() * sizeof(float));
    buffers_[static_cast<int>(Buffers::CYL_PARAMS)]->rebuffer(
        inParaCylinders.PeekElements(), inParaCylinders.Count() * sizeof(float));
    buffers_[static_cast<int>(Buffers::CYL_QUAT)]->rebuffer(
        quatCylinders.PeekElements(), quatCylinders.Count() * sizeof(float));
    buffers_[static_cast<int>(Buffers::CYL_COL1)]->rebuffer(
        color1Cylinders.PeekElements(), color1Cylinders.Count() * sizeof(float));
    buffers_[static_cast<int>(Buffers::CYL_COL2)]->rebuffer(
        color2Cylinders.PeekElements(), color2Cylinders.Count() * sizeof(float));

    glBindVertexArray(vertex_array_spheres_);

    // enable cylinder shader
    cylinderShader_->use();
    // set shader variables
    cylinderShader_->setUniform("viewAttr", viewportStuff);
    cylinderShader_->setUniform("camIn", camera_pose.direction);
    cylinderShader_->setUniform("camRight", camera_pose.right);
    cylinderShader_->setUniform("camUp", camera_pose.up);
    cylinderShader_->setUniform("MVP", MVP);
    cylinderShader_->setUniform("MVinv", MVinv);
    cylinderShader_->setUniform("MVPinv", MVPinv);
    cylinderShader_->setUniform("MVPtransp", MVPtransp);
    cylinderShader_->setUniform("NormalM", NormalM);
    cylinderShader_->setUniform("planes", planes);

    // draw them
    glDrawArrays(GL_POINTS, 0, totalCylinderCnt);

    // disable cylinder shader
    glUseProgram(0);
}
