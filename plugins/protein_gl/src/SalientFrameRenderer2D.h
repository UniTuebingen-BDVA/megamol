/**
 * @author Nils Lichtenberg
 */

#include "mmstd_gl/renderer/Renderer2DModuleGL.h"

#ifndef _SalientFrameRenderer2D_H_
#define _SalientFrameRenderer2D_H_

namespace megamol {
namespace protein_gl {
class SalientFrameRenderer2D : public megamol::mmstd_gl::Renderer2DModuleGL {
public:
    /**
     * Answer the name of this module.
     *
     * @return The name of this module.
     */
    static const char* ClassName(void) {
        return "SalientFrameRenderer";
    }

    /**
     * Answer a human readable description of this module.
     *
     * @return A human readable description of this module.
     */
    static const char* Description(void) {
        return "Renderer for MD frame saliency.";
    }

    /**
     * Answers whether this module is available on the current system.
     *
     * @return 'true' if the module is available, 'false' otherwise.
     */
    static bool IsAvailable(void) {
        return true;
    }

    /** ctor */
    SalientFrameRenderer2D(void);

    /** dtor */
    ~SalientFrameRenderer2D(void);

protected:
    /**
     * Implementation of 'Create'.
     *
     * @return 'true' on success, 'false' otherwise.
     */
    virtual bool create(void);

    /**
     * Implementation of 'Release'.
     */
    virtual void release(void);

    /**
     * Callback for mouse events (move, press, and release)
     *
     * @param x The x coordinate of the mouse in world space
     * @param y The y coordinate of the mouse in world space
     * @param flags The mouse flags
     */
    virtual bool MouseEvent(float x, float y, megamol::core::view::MouseFlags flags);

private:
    /**
     * The get extents callback. The module should set the members of
     * 'call' to tell the caller the extents of its data (bounding boxes
     * and times).
     *
     * @param call The calling call.
     *
     * @return The return value of the function.
     */
    virtual bool GetExtents(mmstd_gl::CallRender2DGL& call);

    /**
     * The Open GL Render callback.
     *
     * @param call The calling call.
     * @return The return value of the function.
     */
    virtual bool Render(mmstd_gl::CallRender2DGL& call);

    /** caller slot */
    core::CallerSlot dataCallerSlot;
};

} // namespace protein_gl
} /* end namespace megamol */
#endif //_SalientFrameRenderer2D_H_
