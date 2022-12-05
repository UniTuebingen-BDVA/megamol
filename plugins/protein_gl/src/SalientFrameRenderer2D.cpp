/**
 * @author Nils Lichtenberg
 */

#include "SalientFrameRenderer2D.h"
#include "protein_calls/SalientFrameCall.h"

using namespace megamol;
using namespace megamol::core;
using namespace megamol::protein_gl;
using megamol::core::utility::log::Log;

SalientFrameRenderer2D::SalientFrameRenderer2D(void)
        : mmstd_gl::Renderer2DModuleGL()
        , dataCallerSlot("getData", "Connects the rendering with data storage.") {

    this->dataCallerSlot.SetCompatibleCall<protein_calls::SalientFrameCallDescription>();
    this->dataCallerSlot.SetNecessity(core::AbstractCallSlotPresentation::SLOT_REQUIRED);
    this->MakeSlotAvailable(&this->dataCallerSlot);
}

bool SalientFrameRenderer2D::create(void) {
    return true;
}

void SalientFrameRenderer2D::release(void) {}

bool SalientFrameRenderer2D::MouseEvent(float x, float y, megamol::core::view::MouseFlags flags) {
    return RendererModule::MouseEvent(x, y, flags);
}

bool SalientFrameRenderer2D::GetExtents(mmstd_gl::CallRender2DGL& call) {
    return true;
}

bool SalientFrameRenderer2D::Render(mmstd_gl::CallRender2DGL& call) {
    protein_calls::SalientFrameCall* pCall = this->dataCallerSlot.CallAs<protein_calls::SalientFrameCall>();
    if (pCall == NULL)
        return false;
    // execute the call
    if (!(*pCall)(protein_calls::SalientFrameCall::CallForGetData))
        return false;
    return true;
}

SalientFrameRenderer2D::~SalientFrameRenderer2D(void) {}
