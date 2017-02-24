/*
* OSPRayVelvetMaterial.cpp
* Copyright (C) 2009-2017 by MegaMol Team
* Alle Rechte vorbehalten.
*/

#include "stdafx.h"
#include "OSPRayVelvetMaterial.h"
#include "mmcore/param/FloatParam.h"
#include "mmcore/param/Vector3fParam.h"

using namespace megamol::ospray;


OSPRayVelvetMaterial::OSPRayVelvetMaterial(void) :
    AbstractOSPRayMaterial(),
    // VELVET
    velvetReflectance("Material::Velvet::", "Reflectance"),
    velvetBackScattering("Material::Velvet::", "BackScattering"),
    velvetHorizonScatteringColor("Material::Velvet::", "Scattering color"),
    velvetHorizonScatteringFallOff("Material::Velvet::", "Scattering fall off") {

    this->velvetReflectance << new core::param::Vector3fParam(vislib::math::Vector<float, 3>(0.4f, 0.0f, 0.0f));
    this->velvetHorizonScatteringColor << new core::param::Vector3fParam(vislib::math::Vector<float, 3>(0.75f, 0.1f, 0.1f));
    this->velvetBackScattering << new core::param::FloatParam(0.5f);
    this->velvetHorizonScatteringFallOff << new core::param::FloatParam(10.0f);
    this->MakeSlotAvailable(&this->velvetBackScattering);
    this->MakeSlotAvailable(&this->velvetHorizonScatteringColor);
    this->MakeSlotAvailable(&this->velvetHorizonScatteringFallOff);
    this->MakeSlotAvailable(&this->velvetReflectance);
}

OSPRayVelvetMaterial::~OSPRayVelvetMaterial(void) {
    // empty
}

void OSPRayVelvetMaterial::readParams() {
    materialContainer.materialType = materialTypeEnum::VELVET;

    auto reflect = this->velvetReflectance.Param<core::param::Vector3fParam>()->Value().PeekComponents();
    materialContainer.velvetReflectance.assign(reflect, reflect + 3);

    auto color = this->velvetHorizonScatteringColor.Param<core::param::Vector3fParam>()->Value().PeekComponents();
    materialContainer.velvetHorizonScatteringColor.assign(color, color + 3);

    materialContainer.velvetBackScattering = this->velvetBackScattering.Param<core::param::FloatParam>()->Value();

    materialContainer.velvetHorizonScatteringFallOff = this->velvetHorizonScatteringFallOff.Param<core::param::FloatParam>()->Value();
}

bool OSPRayVelvetMaterial::InterfaceIsDirty() {
    if (
        this->velvetBackScattering.IsDirty() ||
        this->velvetHorizonScatteringColor.IsDirty() ||
        this->velvetHorizonScatteringFallOff.IsDirty() ||
        this->velvetReflectance.IsDirty() 
        ) {
        this->velvetBackScattering.ResetDirty();
        this->velvetHorizonScatteringColor.ResetDirty();
        this->velvetHorizonScatteringFallOff.ResetDirty();
        this->velvetReflectance.ResetDirty();
        return true;
    } else {
        return false;
    }
}