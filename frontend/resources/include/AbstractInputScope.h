/*
 * AbstractInputScope.h
 *
 * Copyright (C) 2018 by VISUS (Universitaet Stuttgart).
 * Alle Rechte vorbehalten.
 */

#ifndef MEGAMOLCORE_ABSTRACTINPUTSCOPE_H_INCLUDED
#define MEGAMOLCORE_ABSTRACTINPUTSCOPE_H_INCLUDED
#if (defined(_MSC_VER) && (_MSC_VER > 1000))
#pragma once
#endif /* (defined(_MSC_VER) && (_MSC_VER > 1000)) */

#include "KeyboardMouseInput.h"
// TODO: do this include correctly via CMake.
// the structs used here are now located in the frontend_resources CMake module.

namespace megamol {
namespace frontend_resources {

using megamol::frontend_resources::Key;
using megamol::frontend_resources::KeyAction;
using megamol::frontend_resources::Modifiers;
using megamol::frontend_resources::MouseButton;
using megamol::frontend_resources::MouseButtonAction;

class AbstractInputScope {
public:
    /**
     * This event handler can be reimplemented to receive key code events.
     *
     * @return true to stop propagation.
     */
    virtual bool OnKey(Key key, KeyAction action, Modifiers mods) {
        return false;
    }

    /**
     * This event handler can be reimplemented to receive unicode events.
     *
     * @return Returns true if the event was accepted (stopping propagation), otherwise false.
     */
    virtual bool OnChar(unsigned int codePoint) {
        return false;
    }

    /**
     * This event handler can be reimplemented to receive mouse button events.
     *
     * @return Returns true if the event was accepted (stopping propagation), otherwise false.
     */
    virtual bool OnMouseButton(MouseButton button, MouseButtonAction action, Modifiers mods) {
        return false;
    }

    /**
     * This event handler can be reimplemented to receive mouse move events.
     *
     * @return Returns true if the event was accepted (stopping propagation), otherwise false.
     */
    virtual bool OnMouseMove(double x, double y) {
        return false;
    }

    /**
     * This event handler can be reimplemented to receive mouse scroll events.
     *
     * @return Returns true if the event was accepted (stopping propagation), otherwise false.
     */
    virtual bool OnMouseScroll(double dx, double dy) {
        return false;
    }

protected:
    AbstractInputScope() = default;
    virtual ~AbstractInputScope() = default;
};

} /* end namespace frontend_resources */
} /* end namespace megamol */

#endif /* MEGAMOLCORE_ABSTRACTINPUTSCOPE_H_INCLUDED */
