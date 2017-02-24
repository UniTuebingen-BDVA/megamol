/*
* CallOSPRayStructure.h
*
* Copyright (C) 2017 by Universitaet Stuttgart (VISUS).
* Alle Rechte vorbehalten.
*/

#pragma once
#include "mmcore/factories/CallAutoDescription.h"
#include "mmcore/Call.h"
#include <map>
#include <vector>
#include "CallOSPRayMaterial.h"


namespace megamol {
namespace ospray {

enum structureTypeEnum {
    GEOMETRY,
    VOLUME
};

enum geometryTypeEnum {
    SPHERES,
    TRIANGLES,
    STREAMLINES,
    CYLINDERS
};

enum volumeTypeEnum {
    STRUCTUREDVOLUME,
    BLOCKBRICKEDVOLUME,
    GHOSTBLOCKBRICKEDVOLUME,
    ISOSURFACE,
    SLICE
};


class OSPRayStructureContainer {
public:
    structureTypeEnum type;
    OSPRayMaterialContainer* materialContainer;
    geometryTypeEnum geometryType;
    volumeTypeEnum volumeType;

    std::shared_ptr<std::vector<float>> vertexData;
    std::shared_ptr<std::vector<float>> colorData;

    bool dataChanged;
    bool isValid;

    OSPRayStructureContainer();
    ~OSPRayStructureContainer();

};


class CallOSPRayStructure;
typedef std::map<CallOSPRayStructure*, OSPRayStructureContainer> OSPRayStrcutrureMap;


class CallOSPRayStructure : public core::Call {
public:

    /**
    * Answer the name of the objects of this description.
    *
    * @return The name of the objects of this description.
    */
    static const char *ClassName(void) {
        return "CallOSPRayStructure";
    }

    /**
    * Gets a human readable description of the module.
    *
    * @return A human readable description of the module.
    */
    static const char *Description(void) {
        return "Call for an OSPRay structure";
    }

    /**
    * Answer the number of functions used for this call.
    *
    * @return The number of functions used for this call.
    */
    static unsigned int FunctionCount(void) {
        return 1;
    }

    /**
    * Answer the name of the function used for this call.
    *
    * @param idx The index of the function to return it's name.
    *
    * @return The name of the requested function.
    */
    static const char * FunctionName(unsigned int idx) {
        switch (idx) {
        case 0: return "GetDataCall";
        default: return NULL;
        }
    }

    /** Ctor. */
    CallOSPRayStructure();

    /** Dtor. */
    virtual ~CallOSPRayStructure(void);

    /**
    * Assignment operator
    *
    * @param rhs The right hand side operand
    *
    * @return A reference to this
    */
    CallOSPRayStructure& operator=(const CallOSPRayStructure& rhs);

    void setStructureMap(OSPRayStrcutrureMap*sm);
    void addStructure(OSPRayStructureContainer &sc);
    void fillStructureMap();

    OSPRayStrcutrureMap *structureMap;

    void setTime(float time);
    float getTime();


private:

    float time;

};
typedef core::factories::CallAutoDescription<CallOSPRayStructure> CallOSPRayStructureDescription;
} // namespace ospray
} // namespace megamol