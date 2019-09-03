#include "stdafx.h"
#include "LocalLighting.h"

#include "mmcore/CoreInstance.h"

#include "vislib/graphics/gl/ShaderSource.h"

#include "compositing/CompositingCalls.h"

megamol::compositing::LocalLighting::LocalLighting() 
    : Renderer3DModule_2()
    , m_output_texture(nullptr)
    , m_lights_buffer(nullptr)
    , m_output_tex_slot("OutputTexture", "Gives access to resulting output texture")
    , m_albedo_tex_slot("AlbedoTexture", "Connect to the albedo render target texture")
    , m_normal_tex_slot("NormalTexture", "Connects to the normals render target texture")
    , m_depth_tex_slot("DepthTexture", "Connects to the depth render target texture")
    , m_roughness_metalness_tex_slot("RoughMetalTexture","Connects to the roughness/metalness render target texture")
{
    this->m_output_tex_slot.SetCallback(CallTexture2D::ClassName(), "GetData", &LocalLighting::getDataCallback);
    this->m_output_tex_slot.SetCallback(CallTexture2D::ClassName(), "GetMetaData", &LocalLighting::getMetaDataCallback);
    this->MakeSlotAvailable(&this->m_output_tex_slot);

    this->m_albedo_tex_slot.SetCompatibleCall<CallTexture2DDescription>();
    this->MakeSlotAvailable(&this->m_albedo_tex_slot);

    this->m_normal_tex_slot.SetCompatibleCall<CallTexture2DDescription>();
    this->MakeSlotAvailable(&this->m_normal_tex_slot);

    this->m_depth_tex_slot.SetCompatibleCall<CallTexture2DDescription>();
    this->MakeSlotAvailable(&this->m_depth_tex_slot);

    this->m_roughness_metalness_tex_slot.SetCompatibleCall<CallTexture2DDescription>();
    this->MakeSlotAvailable(&this->m_roughness_metalness_tex_slot);
}

megamol::compositing::LocalLighting::~LocalLighting() {}

bool megamol::compositing::LocalLighting::create() {

    try {
        // create shader program
        m_lighting_prgm = std::make_unique<GLSLComputeShader>();

        vislib::graphics::gl::ShaderSource compute_src;

        if (!instance()->ShaderSourceFactory().MakeShaderSource("Compositing::lambert", compute_src))
            return false;
        if (!m_lighting_prgm->Compile(compute_src.Code(), compute_src.Count())) return false;
        if (!m_lighting_prgm->Link()) return false;

    } catch (vislib::graphics::gl::AbstractOpenGLShader::CompileException ce) {
        vislib::sys::Log::DefaultLog.WriteMsg(vislib::sys::Log::LEVEL_ERROR, "Unable to compile shader (@%s): %s\n",
            vislib::graphics::gl::AbstractOpenGLShader::CompileException::CompileActionName(ce.FailedAction()),
            ce.GetMsgA());
        return false;
    } catch (vislib::Exception e) {
        vislib::sys::Log::DefaultLog.WriteMsg(
            vislib::sys::Log::LEVEL_ERROR, "Unable to compile shader: %s\n", e.GetMsgA());
        return false;
    } catch (...) {
        vislib::sys::Log::DefaultLog.WriteMsg(
            vislib::sys::Log::LEVEL_ERROR, "Unable to compile shader: Unknown exception\n");
        return false;
    }

    glowl::TextureLayout tx_layout(GL_RGBA16F, 1, 1, 1, GL_RGBA, GL_HALF_FLOAT, 1);
    m_output_texture = std::make_shared<glowl::Texture2D>("lighting_output", tx_layout, nullptr);

    m_lights_buffer = std::make_unique<glowl::BufferObject>(GL_SHADER_STORAGE_BUFFER, nullptr, 0, GL_DYNAMIC_DRAW);

    return true;
}

void megamol::compositing::LocalLighting::release() {}

bool megamol::compositing::LocalLighting::GetExtents(core::view::CallRender3D_2& call) { return true; }

bool megamol::compositing::LocalLighting::Render(core::view::CallRender3D_2& call) {

    megamol::core::view::CallRender3D_2* cr = &call;
    auto call_albedo = m_albedo_tex_slot.CallAs<CallTexture2D>();
    auto call_normal = m_normal_tex_slot.CallAs<CallTexture2D>();
    auto call_depth = m_normal_tex_slot.CallAs<CallTexture2D>();

    if (cr == NULL) return false;
    if (call_albedo == NULL) return false;
    if (call_normal == NULL) return false;
    if (call_depth == NULL) return false;

    if (!(*call_albedo)(0)) return false;
    if (!(*call_normal)(0)) return false;
    if (!(*call_depth)(0)) return false;

    // set output texture size to primary input texture
    auto albedo_tx2D = call_albedo->getData();
    auto normal_tx2D = call_normal->getData();
    auto depth_tx2D = call_normal->getData();
    std::array<float, 2> texture_res = {
        static_cast<float>(albedo_tx2D->getWidth()), static_cast<float>(albedo_tx2D->getHeight())};

    if (m_output_texture->getWidth() != std::get<0>(texture_res) ||
        m_output_texture->getHeight() != std::get<1>(texture_res)) {
        glowl::TextureLayout tx_layout(
            GL_RGBA16F, std::get<0>(texture_res), std::get<1>(texture_res), 1, GL_RGBA, GL_HALF_FLOAT, 1);
        m_output_texture->reload(tx_layout, nullptr);
    }

    // obtain camera information
    core::view::Camera_2 cam(cr->GetCamera());
    cam_type::snapshot_type snapshot;
    cam_type::matrix_type view_tmp, proj_tmp;
    cam.calc_matrices(snapshot, view_tmp, proj_tmp, core::thecam::snapshot_content::all);
    glm::mat4 view_mx = view_tmp;
    glm::mat4 proj_mx = proj_tmp;

    auto light_update = this->GetLights();
    if (light_update) {
        struct LightParams {
            float x, y, z, intensity;
        };

        auto light_cnt = lightMap.size();

        std::vector<LightParams> lights;
        lights.reserve(light_cnt);

        for (const auto element : this->lightMap) {
            auto light = element.second;
            lights.push_back({light.pl_position[0], light.pl_position[1], light.pl_position[2], light.lightIntensity});
        }

        m_lights_buffer->rebuffer(lights);
    }

    if (m_lighting_prgm != nullptr && m_lights_buffer != nullptr) {
        m_lighting_prgm->Enable();

        m_lights_buffer->bind(1);
        glUniform1i(m_lighting_prgm->ParameterLocation("light_cnt"), 1);

        glActiveTexture(GL_TEXTURE0);
        albedo_tx2D->bindTexture();
        glUniform1i(m_lighting_prgm->ParameterLocation("albedo_tx2D"), 0);

        glActiveTexture(GL_TEXTURE1);
        normal_tx2D->bindTexture();
        glUniform1i(m_lighting_prgm->ParameterLocation("normal_tx2D"), 1);

        glActiveTexture(GL_TEXTURE2);
        depth_tx2D->bindTexture();
        glUniform1i(m_lighting_prgm->ParameterLocation("depth_tx2D"), 2);

        auto inv_view_mx = glm::inverse(view_mx);
        auto inv_proj_mx = glm::inverse(proj_mx);
        glUniformMatrix4fv(m_lighting_prgm->ParameterLocation("inv_view_mx"), 1, GL_FALSE, glm::value_ptr(inv_view_mx));
        glUniformMatrix4fv(m_lighting_prgm->ParameterLocation("inv_proj_mx"), 1, GL_FALSE, glm::value_ptr(inv_proj_mx));

        m_output_texture->bindImage(0, GL_WRITE_ONLY);

         m_lighting_prgm->Dispatch(static_cast<int>(std::ceil(std::get<0>(texture_res) / 8.0f)),
            static_cast<int>(std::ceil(std::get<1>(texture_res) / 8.0f)), 1);

        m_lighting_prgm->Disable();
    }

    return true;
}

void megamol::compositing::LocalLighting::PreRender(core::view::CallRender3D_2& call) {}

bool megamol::compositing::LocalLighting::getDataCallback(core::Call& caller) {
    auto lhs_tc = dynamic_cast<CallTexture2D*>(&caller);

    if (lhs_tc->getData() == nullptr) {
        lhs_tc->setData(m_output_texture);
    }

    return true; 
}

bool megamol::compositing::LocalLighting::getMetaDataCallback(core::Call& caller) { return true; }
