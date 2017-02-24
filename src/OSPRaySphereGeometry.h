/*
* OSPRaySphereRenderer.h
* Copyright (C) 2009-2017 by MegaMol Team
* Alle Rechte vorbehalten.
*/
#pragma once

#include "mmcore/Module.h"
#include "mmcore/param/ParamSlot.h"
#include "mmcore/CallerSlot.h"
#include "mmcore/CalleeSlot.h"
#include "CallOSPRayStructure.h"

namespace megamol {
namespace ospray {

class OSPRaySphereGeometry : core::Module {

public:

    /**
    * Answer the name of this module.
    *
    * @return The name of this module.
    */
    static const char *ClassName(void) {
        return "OSPRaySphereGeometry";
    }

    /**
    * Answer a human readable description of this module.
    *
    * @return A human readable description of this module.
    */
    static const char *Description(void) {
        return "Creator for OSPRay sphere geometries.";
    }

    /**
    * Answers whether this module is available on the current system.
    *
    * @return 'true' if the module is available, 'false' otherwise.
    */
    static bool IsAvailable(void) {
        return true;
    }

    /** Dtor. */
    virtual ~OSPRaySphereGeometry(void);

    /** Ctor. */
    OSPRaySphereGeometry(void);

protected:
    /**
    * color transfer helper
    * @param array with gray scales
    * @param transferfunction table/texture
    * @param transferfunction table/texture size
    * @param target array (rgba)
    */
    void colorTransferGray(std::vector<float> &grayArray, float const* transferTable, unsigned int tableSize, std::vector<float> &rgbaArray);

    virtual bool create();
    virtual void release();

    bool readData(core::Call &call);

    bool getStructureCallback(core::Call& call);
    bool InterfaceIsDirty();
    void OSPRaySphereGeometry::getClipData(float *clipDat, float *clipCol);

    core::param::ParamSlot particleList;

    /** The call for data */
    core::CallerSlot getDataSlot;

    /** The call for clipping plane */
    core::CallerSlot getClipPlaneSlot;

    /** The call for Transfer function */
    core::CallerSlot getTFSlot;

    /** The callee for Structure */
    core::CalleeSlot deployStructureSlot;

    /** The call for additional structures */
    core::CallerSlot getStructureSlot;

    /** The call for materials */
    core::CallerSlot getMaterialSlot;

private:


    //data objects
    std::vector<float> cd_rgba;
    std::vector<float> vd;


    // color transfer data
    unsigned int tex_size;
    SIZE_T datahash;

    OSPRayStructureContainer structureContainer;

    OSPRayMaterialContainer* materialContainer;

};

} // namespace ospray
} // namespace megamol