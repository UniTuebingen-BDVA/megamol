/*
* OSPRayThinGlassMaterial.h
* Copyright (C) 2009-2017 by MegaMol Team
* Alle Rechte vorbehalten.
*/
#pragma once

#include "AbstractOSPRayMaterial.h"

namespace megamol {
namespace ospray {

class OSPRayThinGlassMaterial : public AbstractOSPRayMaterial {

    /**
    * Answer the name of this module.
    *
    * @return The name of this module.
    */
    static const char *ClassName(void) {
        return "OSPRayThinGlassMaterial";
    }

    /**
    * Answer a human readable description of this module.
    *
    * @return A human readable description of this module.
    */
    static const char *Description(void) {
        return "Configuration module for an OSPRay thin glass material";
    }

    /**
    * Answers whether this module is available on the current system.
    *
    * @return 'true' if the module is available, 'false' otherwise.
    */
    static bool IsAvailable(void) {
        return true;
    }

    /** Ctor. */
    OSPRayThinGlassMaterial(void);

    /** Dtor. */
    virtual ~OSPRayThinGlassMaterial(void);

private:

    //THINGLASS
    core::param::ParamSlot thinglassTransmission;
    core::param::ParamSlot thinglassEta;
    core::param::ParamSlot thinglassThickness;


    virtual bool InterfaceIsDirty();
    virtual void readParams();

};


} // namespace ospray
} // namespace megamol