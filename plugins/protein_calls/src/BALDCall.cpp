/**
 * BALDCall.cpp
 *
 * Copyright (C) 2010 by VISUS (University of Stuttgart)
 * Alle Rechte vorbehalten.
 *
 * @author Nils Lichtenberg
 */
#include "protein_calls/BALDCall.h"

using namespace megamol;
using namespace megamol::protein_calls;

BALDCall::BALDCall(void) : core::AbstractGetData3DCall() {
    // intentionally empty
}

BALDCall::~BALDCall(void) {
    this->Unlock(); // just for paranoia reasons
}
