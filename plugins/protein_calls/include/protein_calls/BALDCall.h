/**
 * BALDCall.h
 *
 * Copyright (C) 2010 by VISUS (University of Stuttgart)
 * Alle Rechte vorbehalten.
 *
 * @author Nils Lichtenberg
 */
#ifndef MEGAMOL_PROTEIN_BALDCALL_H_INCLUDED
#define MEGAMOL_PROTEIN_BALDCALL_H_INCLUDED
#if (defined(_MSC_VER) && (_MSC_VER > 1000))
#pragma once
#endif /* (defined(_MSC_VER) && (_MSC_VER > 1000)) */

#include "mmstd/data/AbstractGetData3DCall.h"


namespace megamol::protein_calls {

/**
 * Get data call for BALD input.
 */
class BALDCall : public megamol::core::AbstractGetData3DCall {
public:
    /**
     * Answer the name of the objects of this description.
     *
     * @return The name of the objects of this description.
     */
    static const char* ClassName() {
        return "BALDCall";
    }

    /**
     * Gets a human readable description of the module.
     *
     * @return A human readable description of the module.
     */
    static const char* Description() {
        return "Call to get BALD information from BALDLoader";
    }

    /**
     * Answer the number of functions used for this call.
     *
     * @return The number of functions used for this call.
     */
    static unsigned int FunctionCount() {
        return megamol::core::AbstractGetData3DCall::FunctionCount();
    }

    /**
     * Answer the name of the function used for this call.
     *
     * @param idx The index of the function to return it's name.
     *
     * @return The name of the requested function.
     */
    static const char* FunctionName(unsigned int idx) {
        return megamol::core::AbstractGetData3DCall::FunctionName(idx);
    }

    /** Ctor */
    BALDCall();

    /** Dtor. */
    ~BALDCall() override;

    /** Index of the 'GetData' function */
    static const unsigned int CallForGetData = 0;

    /** Index of the 'GetExtent' function */
    static const unsigned int CallForGetExtent = 1;
};

/** Description class typedef */
typedef megamol::core::factories::CallAutoDescription<BALDCall> BALDCallDescription;

} // namespace megamol::protein_calls

#endif /*  MEGAMOL_PROTEIN_BALDCALL_H_INCLUDED */
